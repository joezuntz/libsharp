ACLOCAL_AMFLAGS = -I m4

lib_LTLIBRARIES = libsharp2.la

libsharp2_la_SOURCES = \
  libsharp2/pocketfft.c \
  libsharp2/pocketfft.h \
  libsharp2/sharp_utils.c \
  libsharp2/sharp_utils.h \
  libsharp2/sharp.c \
  libsharp2/sharp_almhelpers.c \
  libsharp2/sharp_core.c \
  libsharp2/sharp_geomhelpers.c \
  libsharp2/sharp_legendre_roots.c \
  libsharp2/sharp_ylmgen_c.c \
  libsharp2/sharp_internal.h \
  libsharp2/sharp_legendre_roots.h \
  libsharp2/sharp_vecsupport.h \
  libsharp2/sharp_ylmgen_c.h

# format is "current:revision:age"
# any change: increase revision
# any interface change: increase current, revision=0
# any backward-compatible change: increase age
# any backward-incompatible change: age=0
# ==> age <= current
libsharp2_la_LDFLAGS = -version-info 0:0:0

AM_CFLAGS = @AM_CFLAGS@

if HAVE_MULTIARCH

libavx_la_SOURCES = libsharp2/sharp_core_inc.c
libavx2_la_SOURCES = libsharp2/sharp_core_inc.c
libfma_la_SOURCES = libsharp2/sharp_core_inc.c
libfma4_la_SOURCES = libsharp2/sharp_core_inc.c
libavx512f_la_SOURCES = libsharp2/sharp_core_inc.c

noinst_LTLIBRARIES = libavx.la libavx2.la libfma.la libfma4.la libavx512f.la

libavx_la_CFLAGS = ${AM_CFLAGS} -mavx -DARCH=avx
libavx2_la_CFLAGS = ${AM_CFLAGS} -mavx2 -DARCH=avx2
libfma_la_CFLAGS = ${AM_CFLAGS} -mfma -DARCH=fma
libfma4_la_CFLAGS = ${AM_CFLAGS} -mfma4 -DARCH=fma4
libavx512f_la_CFLAGS = ${AM_CFLAGS} -mavx512f -DARCH=avx512f

libsharp2_la_LIBADD = libavx.la libavx2.la libfma.la libfma4.la libavx512f.la

endif

nobase_include_HEADERS = \
  libsharp2/sharp.h \
  libsharp2/sharp_mpi.h \
  libsharp2/sharp_geomhelpers.h \
  libsharp2/sharp_almhelpers.h \
  libsharp2/sharp_cxx.h

EXTRA_DIST = \
  runtest.sh fortran/sharp.f90 fortran/test_sharp.f90 libsharp2/sharp_mpi.c

check_PROGRAMS = sharp2_testsuite
sharp2_testsuite_SOURCES = test/sharp2_testsuite.c test/memusage.c test/memusage.h
sharp2_testsuite_LDADD = libsharp2.la -lm

TESTS = runtest.sh

pkgconfigdir = $(libdir)/pkgconfig
nodist_pkgconfig_DATA = @PACKAGE_NAME@.pc

DISTCLEANFILES=@PACKAGE_NAME@.pc @PACKAGE_NAME@.pc.in  @PACKAGE_NAME@-uninstalled.pc  @PACKAGE_NAME@-uninstalled.sh
