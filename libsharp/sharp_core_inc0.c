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

/*
 *  libsharp is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*! \file sharp_core_inc0.c
 *  Computational core
 *
 *  Copyright (C) 2012-2018 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include <complex.h>
#include <math.h>
#include <string.h>
#include "sharp_vecsupport.h"
#include "sharp_complex_hacks.h"
#include "sharp.h"
#include "sharp_core.h"
#include "c_utils.h"

typedef complex double dcmplx;

// must be in the range [0;6]
#define MAXJOB_SPECIAL 2

#define XCONCATX(a,b) a##b
#define CONCATX(a,b) XCONCATX(a,b)
#define XCONCAT2(a,b) a##_##b
#define CONCAT2(a,b) XCONCAT2(a,b)
#define XCONCAT3(a,b,c) a##_##b##_##c
#define CONCAT3(a,b,c) XCONCAT3(a,b,c)

#define nvec 1
#include "sharp_core_inchelper.c"
#undef nvec

#define nvec 2
#include "sharp_core_inchelper.c"
#undef nvec

#define nvec 3
#include "sharp_core_inchelper.c"
#undef nvec

#define nvec 4
#include "sharp_core_inchelper.c"
#undef nvec

#define nvec 5
#include "sharp_core_inchelper.c"
#undef nvec

#define nvec 6
#include "sharp_core_inchelper.c"
#undef nvec

void CONCATX(inner_loop,ARCH) (sharp_job *job, const int *ispair,const double *cth,
  const double *sth, int llim, int ulim, sharp_Ylmgen_C *gen, int mi,
  const int *mlim)
  {
  int nv=job->flags&SHARP_NVMAX;
  switch (nv)
    {
    case 0x1:
      CONCAT2(inner_loop,1) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
      return;
    case 0x2:
      CONCAT2(inner_loop,2) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
      return;
    case 0x3:
      CONCAT2(inner_loop,3) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
      return;
    case 0x4:
      CONCAT2(inner_loop,4) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
      return;
    case 0x5:
      CONCAT2(inner_loop,5) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
      return;
    case 0x6:
      CONCAT2(inner_loop,6) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
      return;
    }
  UTIL_FAIL("Incorrect vector parameters");
  }
