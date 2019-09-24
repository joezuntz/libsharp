/*
 *  This file is part of libsharp.
 *
 *  libsharp is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libsharp is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libsharp; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/* libsharp is being developed at the Max-Planck-Institut fuer Astrophysik */

/*
 *  Convenience functions
 *
 *  Copyright (C) 2008-2019 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include <stdio.h>
#include "libsharp/sharp_utils.h"

void sharp_fail_ (const char *file, int line, const char *func, const char *msg)
  {
  fprintf(stderr,"%s, %i (%s):\n%s\n",file,line,func,msg);
  exit(1);
  }

/* This function tries to avoid allocations with a total size close to a high
   power of two (called the "critical stride" here), by adding a few more bytes
   if necessary. This lowers the probability that two arrays differ by a multiple
   of the critical stride in their starting address, which in turn lowers the
   risk of cache line contention. */
static size_t manipsize(size_t sz)
  {
  const size_t critical_stride=4096, cacheline=64, overhead=32;
  if (sz < (critical_stride/2)) return sz;
  if (((sz+overhead)%critical_stride)>(2*cacheline)) return sz;
  return sz+2*cacheline;
  }

#ifdef __SSE__
#include <xmmintrin.h>
void *sharp_malloc_ (size_t sz)
  {
  void *res;
  if (sz==0) return NULL;
  res = _mm_malloc(manipsize(sz),32);
  UTIL_ASSERT(res,"_mm_malloc() failed");
  return res;
  }
void sharp_free_ (void *ptr)
  { if ((ptr)!=NULL) _mm_free(ptr); }
#else
void *sharp_malloc_ (size_t sz)
  {
  void *res;
  if (sz==0) return NULL;
  res = malloc(manipsize(sz));
  UTIL_ASSERT(res,"malloc() failed");
  return res;
  }
void sharp_free_ (void *ptr)
  { if ((ptr)!=NULL) free(ptr); }
#endif

#if defined (_OPENMP)
#include <omp.h>
#elif defined (USE_MPI)
#include "mpi.h"
#elif defined (_WIN32)
#include <Windows.h>
#else
#include <sys/time.h>
#include <stdlib.h>
#endif

double sharp_wallTime(void)
  {
#if defined (_OPENMP)
  return omp_get_wtime();
#elif defined (USE_MPI)
  return MPI_Wtime();
#elif defined (_WIN32)
  static double inv_freq = -1.;
  if (inv_freq<0)
    {
    LARGE_INTEGER freq;
    QueryPerformanceFrequency(&freq);
    inv_freq = 1. / double(freq.QuadPart);
    }
  LARGE_INTEGER count;
  QueryPerformanceCounter(&count);
  return count.QuadPart*inv_freq;
#else
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + 1e-6*t.tv_usec;
#endif
  }
