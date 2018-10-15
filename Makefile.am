ACLOCAL_AMFLAGS = -I m4

lib_LTLIBRARIES = libsharp.la

src_sharp = \
  libsharp/sharp.c \
  libsharp/sharp_almhelpers.c \
  libsharp/sharp_announce.c \
  libsharp/sharp_core.c \
  libsharp/sharp_geomhelpers.c \
  libsharp/sharp_legendre_roots.c \
  libsharp/sharp_ylmgen_c.c \
  libsharp/sharp_announce.h \
  libsharp/sharp_complex_hacks.h \
  libsharp/sharp_core.h \
  libsharp/sharp_internal.h \
  libsharp/sharp_legendre_roots.h \
  libsharp/sharp_lowlevel.h \
  libsharp/sharp_vecsupport.h \
  libsharp/sharp_vecutil.h \
  libsharp/sharp_ylmgen_c.h

include_HEADERS = \
  libsharp/sharp.h \
  libsharp/sharp_lowlevel.h \
  libsharp/sharp_geomhelpers.h \
  libsharp/sharp_almhelpers.h \
  libsharp/sharp_cxx.h

EXTRA_DIST = \
  libsharp/sharp_core_inc.c \
  libsharp/sharp_core_inc2.c \
  libsharp/sharp_core_inchelper.c

libsharp_la_SOURCES = $(src_sharp)

#check_PROGRAMS = ffttest
#ffttest_SOURCES = ffttest.c
#ffttest_LDADD = libpocketfft.la -lm

#TESTS = ffttest

AM_CFLAGS = -I$(top_srcdir)

pkgconfigdir = $(libdir)/pkgconfig
nodist_pkgconfig_DATA = @PACKAGE_NAME@.pc

DISTCLEANFILES=@PACKAGE_NAME@.pc @PACKAGE_NAME@.pc.in  @PACKAGE_NAME@-uninstalled.pc  @PACKAGE_NAME@-uninstalled.sh
