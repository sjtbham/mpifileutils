#include <dirent.h>
#include <limits.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h> /* asctime / localtime */

#include <stdarg.h> /* variable length args */

#include <pwd.h> /* for getpwent */
#include <grp.h> /* for getgrent */
#include <errno.h>

#include "libcircle.h"
#include "dtcmp.h"
#include "bayer.h"

// getpwent getgrent to read user and group entries

/* TODO: change globals to struct */
static int verbose   = 1;
static int walk_stat = 1;

static char* outname = NULL;
static FILE* outfile = NULL;

/* we'll set this to the youngest age still acceptable */
static uint64_t cutoff;

/* keep stats during walk */
uint64_t total_dirs    = 0;
uint64_t total_files   = 0;
uint64_t total_links   = 0;
uint64_t total_unknown = 0;
uint64_t total_bytes   = 0;

/****************************************
 * Global counter and callbacks for LIBCIRCLE reductions
 ***************************************/

uint64_t reduce_items;

static void reduce_init(void)
{
    CIRCLE_reduce(&reduce_items, sizeof(uint64_t));
}

static void reduce_exec(const void* buf1, size_t size1, const void* buf2, size_t size2)
{
    const uint64_t* a = (const uint64_t*) buf1;
    const uint64_t* b = (const uint64_t*) buf2;
    uint64_t val = a[0] + b[0];
    CIRCLE_reduce(&val, sizeof(uint64_t));
}

static void reduce_fini(const void* buf, size_t size)
{
    /* get current time */
    time_t walk_start_t = time(NULL);
    if (walk_start_t == (time_t) - 1) {
        /* TODO: ERROR! */
    }

    /* format timestamp string */
    char walk_s[30];
    size_t rc = strftime(walk_s, sizeof(walk_s) - 1, "%FT%T", localtime(&walk_start_t));
    if (rc == 0) {
        walk_s[0] = '\0';
    }

    /* get result of reduction */
    const uint64_t* a = (const uint64_t*) buf;
    unsigned long long val = (unsigned long long) a[0];

    /* print status to stdout */
    printf("%s: Items walked %llu ...\n", walk_s, val);
    fflush(stdout);
}

static char CURRENT_DIR[PATH_MAX];

/* insert a file given its mode and optional stat data */
static void list_insert_stat(const char* fpath, mode_t mode, const struct stat* sb)
{
    /* set file type */
    if (S_ISDIR(mode)) {
        total_dirs++;
    }
    else if (S_ISREG(mode)) {
        total_files++;
    }
    else if (S_ISLNK(mode)) {
        total_links++;
    }
    else {
        /* unknown file type */
        total_unknown++;
    }

    /* copy stat info */
    if (sb != NULL) {
#if 0
        elem->mode  = (uint64_t) sb->st_mode;
        elem->uid   = (uint64_t) sb->st_uid;
        elem->gid   = (uint64_t) sb->st_gid;
        elem->atime = (uint64_t) sb->st_atime;
        elem->mtime = (uint64_t) sb->st_mtime;
        elem->ctime = (uint64_t) sb->st_ctime;
        elem->size  = (uint64_t) sb->st_size;
#endif

        total_bytes += (uint64_t) sb->st_size;

        if (S_ISREG(mode) || S_ISLNK(mode)) {
            uint64_t atime  = (uint64_t) sb->st_atime;
            uint64_t mtime  = (uint64_t) sb->st_mtime;
            uint64_t crtime = (uint64_t) sb->st_ctime;
            if (atime < cutoff && mtime < cutoff && crtime < cutoff) {
                /* write path to file */
                fprintf(outfile, "%s\n", fpath);
            }
        }
    } else {
        /* write path to file */
        fprintf(outfile, "%s\n", fpath);
    }

    return;
}

/****************************************
 * Walk directory tree using stat at top level and readdir
 ***************************************/

static void walk_readdir_process_dir(char* dir, CIRCLE_handle* handle)
{
    /* TODO: may need to try these functions multiple times */
    DIR* dirp = bayer_opendir(dir);

    if (! dirp) {
        /* TODO: print error */
    }
    else {
        /* Read all directory entries */
        while (1) {
            /* read next directory entry */
            struct dirent* entry = bayer_readdir(dirp);
            if (entry == NULL) {
                break;
            }

            /* process component, unless it's "." or ".." */
            char* name = entry->d_name;
            if ((strncmp(name, ".", 2)) && (strncmp(name, "..", 3))) {
                /* <dir> + '/' + <name> + '/0' */
                char newpath[CIRCLE_MAX_STRING_LEN];
                size_t len = strlen(dir) + 1 + strlen(name) + 1;
                if (len < sizeof(newpath)) {
                    /* build full path to item */
                    strcpy(newpath, dir);
                    strcat(newpath, "/");
                    strcat(newpath, name);

#ifdef _DIRENT_HAVE_D_TYPE
                    /* record info for item */
                    mode_t mode;
                    int have_mode = 0;
                    if (entry->d_type != DT_UNKNOWN) {
                        /* we can read object type from directory entry */
                        have_mode = 1;
                        mode = DTTOIF(entry->d_type);
                        list_insert_stat(newpath, mode, NULL);
                    }
                    else {
                        /* type is unknown, we need to stat it */
                        struct stat st;
                        int status = bayer_lstat(newpath, &st);
                        if (status == 0) {
                            have_mode = 1;
                            mode = st.st_mode;
                            list_insert_stat(newpath, mode, &st);
                        }
                        else {
                            /* error */
                        }
                    }

                    /* increment our item count */
                    reduce_items++;

                    /* recurse into directories */
                    if (have_mode && S_ISDIR(mode)) {
                        handle->enqueue(newpath);
                    }
#endif
                }
                else {
                    /* TODO: print error in correct format */
                    /* name is too long */
                    printf("Path name is too long: %lu chars exceeds limit %lu\n", len, sizeof(newpath));
                    fflush(stdout);
                }
            }
        }
    }

    bayer_closedir(dirp);

    return;
}

/** Call back given to initialize the dataset. */
static void walk_readdir_create(CIRCLE_handle* handle)
{
    char* path = CURRENT_DIR;

    /* stat top level item */
    struct stat st;
    int status = bayer_lstat(path, &st);
    if (status != 0) {
        /* TODO: print error */
        return;
    }

    /* increment our item count */
    reduce_items++;

    /* record item info */
    list_insert_stat(path, st.st_mode, &st);

    /* recurse into directory */
    if (S_ISDIR(st.st_mode)) {
        walk_readdir_process_dir(path, handle);
    }

    return;
}

/** Callback given to process the dataset. */
static void walk_readdir_process(CIRCLE_handle* handle)
{
    /* in this case, only items on queue are directories */
    char path[CIRCLE_MAX_STRING_LEN];
    handle->dequeue(path);
    walk_readdir_process_dir(path, handle);
    return;
}

/****************************************
 * Walk directory tree using stat on every object
 ***************************************/

static void walk_stat_process_dir(char* dir, CIRCLE_handle* handle)
{
    /* TODO: may need to try these functions multiple times */
    DIR* dirp = bayer_opendir(dir);

    if (! dirp) {
        /* TODO: print error */
    }
    else {
        while (1) {
            /* read next directory entry */
            struct dirent* entry = bayer_readdir(dirp);
            if (entry == NULL) {
                break;
            }

            /* We don't care about . or .. */
            char* name = entry->d_name;
            if ((strncmp(name, ".", 2)) && (strncmp(name, "..", 3))) {
                /* <dir> + '/' + <name> + '/0' */
                char newpath[CIRCLE_MAX_STRING_LEN];
                size_t len = strlen(dir) + 1 + strlen(name) + 1;
                if (len < sizeof(newpath)) {
                    /* build full path to item */
                    strcpy(newpath, dir);
                    strcat(newpath, "/");
                    strcat(newpath, name);

                    /* add item to queue */
                    handle->enqueue(newpath);
                }
                else {
                    /* TODO: print error in correct format */
                    /* name is too long */
                    printf("Path name is too long: %lu chars exceeds limit %lu\n", len, sizeof(newpath));
                    fflush(stdout);
                }
            }
        }
    }

    bayer_closedir(dirp);

    return;
}

/** Call back given to initialize the dataset. */
static void walk_stat_create(CIRCLE_handle* handle)
{
    /* we'll call stat on every item */
    handle->enqueue(CURRENT_DIR);
}

/** Callback given to process the dataset. */
static void walk_stat_process(CIRCLE_handle* handle)
{
    /* get path from queue */
    char path[CIRCLE_MAX_STRING_LEN];
    handle->dequeue(path);

    /* stat item */
    struct stat st;
    int status = bayer_lstat(path, &st);
    if (status != 0) {
        /* print error */
        return;
    }

    /* increment our item count */
    reduce_items++;

    /* TODO: filter items by stat info */

    /* record info for item in list */
    list_insert_stat(path, st.st_mode, &st);

    /* recurse into directory */
    if (S_ISDIR(st.st_mode)) {
        /* TODO: check that we can recurse into directory */
        walk_stat_process_dir(path, handle);
    }

    return;
}

/* Set up and execute directory walk */
static void walk_path(const char* dirpath, int use_stat)
{
    /* initialize libcircle */
    CIRCLE_init(0, NULL, CIRCLE_SPLIT_EQUAL);

    /* set libcircle verbosity level */
    enum CIRCLE_loglevel loglevel = CIRCLE_LOG_WARN;
    CIRCLE_enable_logging(loglevel);

    /* set some global variables to do the file walk */
    strncpy(CURRENT_DIR, dirpath, PATH_MAX);

    /* we lookup users and groups first in case we can use
     * them to filter the walk */
#if 0
    flist->detail = 0;
    if (use_stat) {
        flist->detail = 1;
        if (flist->have_users == 0) {
            get_users(&flist->users);
            create_map(&flist->users, flist->user_id2name);
            flist->have_users = 1;
        }
        if (flist->have_groups == 0) {
            get_groups(&flist->groups);
            create_map(&flist->groups, flist->group_id2name);
            flist->have_groups = 1;
        }
    }
#endif

    /* register callbacks */
    if (use_stat) {
        /* walk directories by calling stat on every item */
        CIRCLE_cb_create(&walk_stat_create);
        CIRCLE_cb_process(&walk_stat_process);
//        CIRCLE_cb_create(&walk_lustrestat_create);
//        CIRCLE_cb_process(&walk_lustrestat_process);
    }
    else {
        /* walk directories using file types in readdir */
        CIRCLE_cb_create(&walk_readdir_create);
        CIRCLE_cb_process(&walk_readdir_process);
//        CIRCLE_cb_create(&walk_getdents_create);
//        CIRCLE_cb_process(&walk_getdents_process);
    }

    /* prepare callbacks and initialize variables for reductions */
    reduce_items = 0;
    CIRCLE_cb_reduce_init(&reduce_init);
    CIRCLE_cb_reduce_op(&reduce_exec);
    CIRCLE_cb_reduce_fini(&reduce_fini);

    /* run the libcircle job */
    CIRCLE_begin();
    CIRCLE_finalize();

    return;
}

static void print_summary(int use_stat)
{
    /* get our rank and the size of comm_world */
    int rank, ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);

    /* get total directories, files, links, and bytes */
    uint64_t all_dirs, all_files, all_links, all_unknown, all_bytes, all_count;
    MPI_Allreduce(&total_dirs,    &all_dirs,     1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&total_files,   &all_files,    1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&total_links,   &all_links,    1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&total_unknown, &all_unknown,  1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&total_bytes,   &all_bytes,    1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&reduce_items,  &all_count,    1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

    /* convert total size to units */
    if (verbose && rank == 0) {
        printf("Items: %llu\n", (unsigned long long) all_count);
        printf("  Directories: %llu\n", (unsigned long long) all_dirs);
        printf("  Files: %llu\n", (unsigned long long) all_files);
        printf("  Links: %llu\n", (unsigned long long) all_links);
        /* printf("  Unknown: %lu\n", (unsigned long long) all_unknown); */

        if (use_stat) {
            double agg_size_tmp;
            const char* agg_size_units;
            bayer_format_bytes(all_bytes, &agg_size_tmp, &agg_size_units);

            uint64_t size_per_file = 0.0;
            if (all_files > 0) {
                size_per_file = (uint64_t)((double)all_bytes / (double)all_files);
            }
            double size_per_file_tmp;
            const char* size_per_file_units;
            bayer_format_bytes(size_per_file, &size_per_file_tmp, &size_per_file_units);

            printf("Data: %.3lf %s (%.3lf %s per file)\n", agg_size_tmp, agg_size_units, size_per_file_tmp, size_per_file_units);
        }
    }

    return;
}

static uint32_t gettime(void)
{
    uint32_t t = 0;
    time_t now = time(NULL);
    if (now != (time_t)-1) {
        t = (uint32_t) now;
    }
    return t;
}

static void print_usage(void)
{
    printf("\n");
    printf("Usage: purgelist [options] -o file <path> ...\n");
    printf("\n");
    printf("Options:\n");
    printf("  -d, --days <n>      - set purge threshold at n days ([cma]time)\n");
    printf("  -o, --output <file> - write processed list to file\n");
    printf("  -l, --lite          - walk file system without stat\n");
//    printf("  -v, --verbose       - verbose output\n");
    printf("  -h, --help          - print usage\n");
    printf("\n");
    fflush(stdout);
    return;
}

int main(int argc, char** argv)
{
    int i;

    /* initialize MPI */
    MPI_Init(&argc, &argv);
    bayer_init();

    /* get our rank and the size of comm_world */
    int rank, ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);

    char* outputname = NULL;
    int walk = 0;
    int days = 60;

    int option_index = 0;
    static struct option long_options[] = {
        {"days",     1, 0, 'd'},
        {"output",   1, 0, 'o'},
        {"lite",     0, 0, 'l'},
        {"help",     0, 0, 'h'},
        {"verbose",  0, 0, 'v'},
        {0, 0, 0, 0}
    };

    int usage = 0;
    while (1) {
        int c = getopt_long(
                    argc, argv, "d:o:lhv",
                    long_options, &option_index
                );

        if (c == -1) {
            break;
        }

        switch (c) {
            case 'd':
                days = atoi(optarg);
                break;
            case 'o':
                outputname = BAYER_STRDUP(optarg);
                break;
            case 'l':
                walk_stat = 0;
                break;
            case 'h':
                usage = 1;
                break;
            case 'v':
                verbose = 1;
                break;
            case '?':
                usage = 1;
                break;
            default:
                if (rank == 0) {
                    printf("?? getopt returned character code 0%o ??\n", c);
                }
        }
    }

    /* paths to walk come after the options */
    int numpaths = 0;
    bayer_param_path* paths = NULL;
    if (optind < argc) {
        /* got a path to walk */
        walk = 1;

        /* determine number of paths specified by user */
        numpaths = argc - optind;

        /* allocate space for each path */
        paths = (bayer_param_path*) BAYER_MALLOC((size_t)numpaths * sizeof(bayer_param_path));

        /* process each path */
        for (i = 0; i < numpaths; i++) {
            const char* path = argv[optind];
            bayer_param_path_set(path, &paths[i]);
            optind++;
        }
    }

    /* check that days is valid */
    if (days < 0) {
        BAYER_LOG(BAYER_LOG_ERR, "--days option must be non-negative `%d'\n", days);
        usage = 1;
    }

    /* require user to give us an output file name */
    if (outputname == NULL) {
        usage = 1;
    }

    if (usage) {
        if (rank == 0) {
            print_usage();
        }
        MPI_Finalize();
        return 0;
    }

    /* open out output file */
    outname = BAYER_STRDUPF("%s.%d", outputname, rank);
    outfile = fopen(outname, "w");
    if (outfile == NULL) {
        BAYER_ABORT("Error opening file `%s': errno=%d %s\n", outname, errno, strerror(errno));
    }

    /* determine cutoff point we're willing to keep */
    uint64_t limit = days * 24 * 3600; /* 60 days */
    uint64_t now = (uint64_t) gettime();
    cutoff = now - limit;

    uint64_t walk_start, walk_end;

    time_t walk_start_t = time(NULL);
    if (walk_start_t == (time_t)-1) {
        /* TODO: ERROR! */
    }
    walk_start = (uint64_t) walk_start_t;

    /* report walk count, time, and rate */
    double start_walk = MPI_Wtime();
    for (i = 0; i < numpaths; i++) {
        /* get path for this step */
        const char* target = paths[i].path;

        /* print message to user that we're starting */
        if (verbose && rank == 0) {
            char walk_s[30];
            size_t rc = strftime(walk_s, sizeof(walk_s) - 1, "%FT%T", localtime(&walk_start_t));
            if (rc == 0) {
                walk_s[0] = '\0';
            }
            printf("%s: Walking %s\n", walk_s, target);
            fflush(stdout);
        }

        /* walk file tree and record stat data for each file */
        walk_path(target, walk_stat);
    }
    double end_walk = MPI_Wtime();

    time_t walk_end_t = time(NULL);
    if (walk_end_t == (time_t)-1) {
        /* TODO: ERROR! */
    }
    walk_end = (uint64_t) walk_end_t;

    /* compute total number of items walked */
    uint64_t all_count = 0;
    MPI_Allreduce(&reduce_items, &all_count, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

    /* report walk count, time, and rate */
    if (verbose && rank == 0) {
        double secs = end_walk - start_walk;
        double rate = 0.0;
        if (secs > 0.0) {
            rate = ((double)all_count) / secs;
        }
        printf("Walked %lu files in %f seconds (%f files/sec)\n",
               all_count, secs, rate
              );
    }

    print_summary(walk_stat);

    /* close our output file */
    int close_rc = fclose(outfile);
    if (close_rc == EOF) {
        BAYER_LOG(BAYER_LOG_ERR, "Error closing file `%s': errno=%d %s\n", outname, errno, strerror(errno));
    }
    bayer_free(&outname);

    /* free memory allocated for options */
    bayer_free(&outputname);

    /* free the path parameters */
    for (i = 0; i < numpaths; i++) {
        bayer_param_path_free(&paths[i]);
    }
    bayer_free(&paths);

    /* shut down MPI */
    bayer_finalize();
    MPI_Finalize();

    return 0;
}
