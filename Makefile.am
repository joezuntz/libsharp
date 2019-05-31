ACLOCAL_AMFLAGS = -I m4

lib_LTLIBRARIES = libsharp.la

libsharp_la_SOURCES = \
  pocketfft/pocketfft.c \
  pocketfft/pocketfft.h \
  libsharp/sharp_utils.c \
  libsharp/sharp_utils.h \
  libsharp/sharp.c \
  libsharp/sharp_almhelpers.c \
  libsharp/sharp_core.c \
  libsharp/sharp_geomhelpers.c \
  libsharp/sharp_legendre_roots.c \
  libsharp/sharp_ylmgen_c.c \
  libsharp/sharp_internal.h \
  libsharp/sharp_legendre_roots.h \
  libsharp/sharp_vecsupport.h \
  libsharp/sharp_ylmgen_c.h

# format is "current:revision:age"
# any change: increase revision
# any interface change: increase current, revision=0
# any backward-compatible change: increase age
# any backward-incompatible change: age=0
# ==> age <= current
libsharp_la_LDFLAGS = -version-info 0:0:0

AM_CFLAGS = @AM_CFLAGS@

if HAVE_MULTIARCH

libavx_la_SOURCES = libsharp/sharp_core_inc.c
libavx2_la_SOURCES = libsharp/sharp_core_inc.c
libfma_la_SOURCES = libsharp/sharp_core_inc.c
libfma4_la_SOURCES = libsharp/sharp_core_inc.c
libavx512f_la_SOURCES = libsharp/sharp_core_inc.c

noinst_LTLIBRARIES = libavx.la libavx2.la libfma.la libfma4.la libavx512f.la

libavx_la_CFLAGS = ${AM_CFLAGS} -mavx -DARCH=avx
libavx2_la_CFLAGS = ${AM_CFLAGS} -mavx2 -DARCH=avx2
libfma_la_CFLAGS = ${AM_CFLAGS} -mfma -DARCH=fma
libfma4_la_CFLAGS = ${AM_CFLAGS} -mfma4 -DARCH=fma4
libavx512f_la_CFLAGS = ${AM_CFLAGS} -mavx512f -DARCH=avx512f

libsharp_la_LIBADD = libavx.la libavx2.la libfma.la libfma4.la libavx512f.la

endif

nobase_include_HEADERS = \
  libsharp/sharp.h \
  libsharp/sharp_mpi.h \
  libsharp/sharp_geomhelpers.h \
  libsharp/sharp_almhelpers.h \
  libsharp/sharp_cxx.h

EXTRA_DIST = \
  runtest.sh fortran/sharp.f90 fortran/test_sharp.f90 libsharp/sharp_mpi.c

check_PROGRAMS = sharp_testsuite
sharp_testsuite_SOURCES = test/sharp_testsuite.c test/memusage.c test/memusage.h
sharp_testsuite_LDADD = libsharp.la -lm

TESTS = runtest.sh

pkgconfigdir = $(libdir)/pkgconfig
nodist_pkgconfig_DATA = @PACKAGE_NAME@.pc

DISTCLEANFILES=@PACKAGE_NAME@.pc @PACKAGE_NAME@.pc.in  @PACKAGE_NAME@-uninstalled.pc  @PACKAGE_NAME@-uninstalled.sh
