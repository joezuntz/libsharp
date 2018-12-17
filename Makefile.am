ACLOCAL_AMFLAGS = -I m4

lib_LTLIBRARIES = libsharp.la

src_sharp = \
  c_utils/c_utils.c \
  c_utils/c_utils.h \
  pocketfft/pocketfft.c \
  pocketfft/pocketfft.h \
  libsharp/sharp.c \
  libsharp/sharp_almhelpers.c \
  libsharp/sharp_announce.c \
  libsharp/sharp_core.c \
  libsharp/sharp_geomhelpers.c \
  libsharp/sharp_legendre_roots.c \
  libsharp/sharp_ylmgen_c.c \
  libsharp/sharp_announce.h \
  libsharp/sharp_complex_hacks.h \
  libsharp/sharp_internal.h \
  libsharp/sharp_legendre_roots.h \
  libsharp/sharp_vecsupport.h \
  libsharp/sharp_vecutil.h \
  libsharp/sharp_ylmgen_c.h

include_HEADERS = \
  libsharp/sharp.h \
  libsharp/sharp_geomhelpers.h \
  libsharp/sharp_almhelpers.h \
  libsharp/sharp_cxx.h

libsharp_la_SOURCES = $(src_sharp)

check_PROGRAMS = sharp_testsuite
sharp_testsuite_SOURCES = libsharp/sharp_testsuite.c c_utils/memusage.c c_utils/memusage.h c_utils/walltime_c.c c_utils/walltime_c.h
sharp_testsuite_LDADD = libsharp.la

#TESTS = ffttest

AM_CFLAGS = -I$(top_srcdir)/c_utils -I$(top_srcdir)/libsharp @AM_CFLAGS@

pkgconfigdir = $(libdir)/pkgconfig
nodist_pkgconfig_DATA = @PACKAGE_NAME@.pc

DISTCLEANFILES=@PACKAGE_NAME@.pc @PACKAGE_NAME@.pc.in  @PACKAGE_NAME@-uninstalled.pc  @PACKAGE_NAME@-uninstalled.sh
