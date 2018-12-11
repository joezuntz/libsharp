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

#define XCONCATX(a,b) a##b
#define CONCATX(a,b) XCONCATX(a,b)
#define XCONCAT2(a,b) a##_##b
#define CONCAT2(a,b) XCONCAT2(a,b)

#define nvec 6
#include "sharp_core_inchelper.c"
#undef nvec

void CONCATX(inner_loop,ARCH) (sharp_job *job, const int *ispair,const double *cth,
  const double *sth, int llim, int ulim, sharp_Ylmgen_C *gen, int mi,
  const int *mlim)
  {
  CONCAT2(inner_loop,6) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
  }
