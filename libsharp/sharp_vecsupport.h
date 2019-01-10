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

/*  \file sharp_vecsupport.h
 *  Convenience functions for vector arithmetics
 *
 *  Copyright (C) 2012-2019 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#ifndef SHARP_VECSUPPORT_H
#define SHARP_VECSUPPORT_H

#include <math.h>

#ifndef VLEN

#if (defined(__AVX512F__))
#define VLEN 8
#elif (defined (__AVX__))
#define VLEN 4
#elif (defined (__SSE2__))
#define VLEN 2
#else
#define VLEN 1
#endif

#endif

typedef double Ts;

#if (VLEN==1)

typedef double Tv;
typedef int Tm;

#define vload(a) (a)
#define vzero 0.
#define vone 1.

#define vaddeq_mask(mask,a,b) if (mask) (a)+=(b);
#define vsubeq_mask(mask,a,b) if (mask) (a)-=(b);
#define vmuleq_mask(mask,a,b) if (mask) (a)*=(b);
#define vneg(a) (-(a))
#define vabs(a) fabs(a)
#define vsqrt(a) sqrt(a)
#define vlt(a,b) ((a)<(b))
#define vgt(a,b) ((a)>(b))
#define vge(a,b) ((a)>=(b))
#define vne(a,b) ((a)!=(b))
#define vand_mask(a,b) ((a)&&(b))
#define vor_mask(a,b) ((a)||(b))
static inline Tv vmin (Tv a, Tv b) { return (a<b) ? a : b; }
static inline Tv vmax (Tv a, Tv b) { return (a>b) ? a : b; }
#define vanyTrue(a) (a)
#define vallTrue(a) (a)

#endif

#if (VLEN==2)

#include <emmintrin.h>

#if defined (__SSE3__)
#include <pmmintrin.h>
#endif
#if defined (__SSE4_1__)
#include <smmintrin.h>
#endif

typedef __m128d Tv;
typedef __m128d Tm;

#if defined(__SSE4_1__)
#define vblend__(m,a,b) _mm_blendv_pd(b,a,m)
#else
static inline Tv vblend__(Tv m, Tv a, Tv b)
  { return _mm_or_pd(_mm_and_pd(a,m),_mm_andnot_pd(m,b)); }
#endif
#define vload(a) _mm_set1_pd(a)
#define vzero _mm_setzero_pd()
#define vone vload(1.)

#define vaddeq_mask(mask,a,b) a+=vblend__(mask,b,vzero)
#define vsubeq_mask(mask,a,b) a-=vblend__(mask,b,vzero)
#define vmuleq_mask(mask,a,b) a*=vblend__(mask,b,vone)
#define vneg(a) _mm_xor_pd(vload(-0.),a)
#define vabs(a) _mm_andnot_pd(vload(-0.),a)
#define vsqrt(a) _mm_sqrt_pd(a)
#define vlt(a,b) _mm_cmplt_pd(a,b)
#define vgt(a,b) _mm_cmpgt_pd(a,b)
#define vge(a,b) _mm_cmpge_pd(a,b)
#define vne(a,b) _mm_cmpneq_pd(a,b)
#define vand_mask(a,b) _mm_and_pd(a,b)
#define vor_mask(a,b) _mm_or_pd(a,b)
#define vmin(a,b) _mm_min_pd(a,b)
#define vmax(a,b) _mm_max_pd(a,b);
#define vanyTrue(a) (_mm_movemask_pd(a)!=0)
#define vallTrue(a) (_mm_movemask_pd(a)==3)

#endif

#if (VLEN==4)

#include <immintrin.h>

typedef __m256d Tv;
typedef __m256d Tm;

#define vblend__(m,a,b) _mm256_blendv_pd(b,a,m)
#define vload(a) _mm256_set1_pd(a)
#define vzero _mm256_setzero_pd()
#define vone vload(1.)

#define vaddeq_mask(mask,a,b) a+=vblend__(mask,b,vzero)
#define vsubeq_mask(mask,a,b) a-=vblend__(mask,b,vzero)
#define vmuleq_mask(mask,a,b) a*=vblend__(mask,b,vone)
#define vneg(a) _mm256_xor_pd(vload(-0.),a)
#define vabs(a) _mm256_andnot_pd(vload(-0.),a)
#define vsqrt(a) _mm256_sqrt_pd(a)
#define vlt(a,b) _mm256_cmp_pd(a,b,_CMP_LT_OQ)
#define vgt(a,b) _mm256_cmp_pd(a,b,_CMP_GT_OQ)
#define vge(a,b) _mm256_cmp_pd(a,b,_CMP_GE_OQ)
#define vne(a,b) _mm256_cmp_pd(a,b,_CMP_NEQ_OQ)
#define vand_mask(a,b) _mm256_and_pd(a,b)
#define vor_mask(a,b) _mm256_or_pd(a,b)
#define vmin(a,b) _mm256_min_pd(a,b)
#define vmax(a,b) _mm256_max_pd(a,b)
#define vanyTrue(a) (_mm256_movemask_pd(a)!=0)
#define vallTrue(a) (_mm256_movemask_pd(a)==15)

#endif

#if (VLEN==8)

#include <immintrin.h>

typedef __m512d Tv;
typedef __mmask8 Tm;

#define vload(a) _mm512_set1_pd(a)
#define vzero _mm512_setzero_pd()
#define vone vload(1.)

#define vaddeq_mask(mask,a,b) a=_mm512_mask_add_pd(a,mask,a,b);
#define vsubeq_mask(mask,a,b) a=_mm512_mask_sub_pd(a,mask,a,b);
#define vmuleq_mask(mask,a,b) a=_mm512_mask_mul_pd(a,mask,a,b);
#define vneg(a) _mm512_mul_pd(a,vload(-1.))
#define vabs(a) (__m512d)_mm512_andnot_epi64((__m512i)vload(-0.),(__m512i)a)
#define vsqrt(a) _mm512_sqrt_pd(a)
#define vlt(a,b) _mm512_cmp_pd_mask(a,b,_CMP_LT_OQ)
#define vgt(a,b) _mm512_cmp_pd_mask(a,b,_CMP_GT_OQ)
#define vge(a,b) _mm512_cmp_pd_mask(a,b,_CMP_GE_OQ)
#define vne(a,b) _mm512_cmp_pd_mask(a,b,_CMP_NEQ_OQ)
#define vand_mask(a,b) ((a)&(b))
#define vor_mask(a,b) ((a)|(b))
#define vmin(a,b) _mm512_min_pd(a,b)
#define vmax(a,b) _mm512_max_pd(a,b)
#define vanyTrue(a) (a!=0)
#define vallTrue(a) (a==255)

#endif

#endif
