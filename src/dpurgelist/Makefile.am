include $(top_srcdir)/common.mk

bin_PROGRAMS = dpurgelist

dpurgelist_SOURCES = dpurgelist.c
dpurgelist_CPPFLAGS = -I../common/ $(MPI_CFLAGS) $(libcircle_CFLAGS)
dpurgelist_LDFLAGS = $(MPI_CLDFLAGS)
dpurgelist_LDADD = ../common/libfileutils_common.a $(MPI_CLDFLAGS) $(libcircle_LIBS)
