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

/*! \file sharp_core.c
 *  Computational core
 *
 *  Copyright (C) 2012-2013 Max-Planck-Society
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

#define nvec (256/VLEN)

typedef union
  { Tv v; double s[VLEN]; } Tvu;

typedef struct
  { Tv v[nvec]; } Tb;

typedef union
  { Tb b; double s[VLEN*nvec]; } Tbu;

typedef struct
  { Tb r, i; } Tbri;

typedef struct
  { Tb qr, qi, ur, ui; } Tbqu;

typedef struct
  { double r[VLEN*nvec], i[VLEN*nvec]; } Tsri;

typedef struct
  { double qr[VLEN*nvec],qi[VLEN*nvec],ur[VLEN*nvec],ui[VLEN*nvec]; } Tsqu;

typedef union
  { Tbri b; Tsri s; } Tburi;

typedef union
  { Tbqu b; Tsqu s; } Tbuqu;

static inline Tb Tbconst(double val)
  {
  Tv v=vload(val);
  Tb res;
  for (int i=0; i<nvec; ++i) res.v[i]=v;
  return res;
  }

static inline void Tbmuleq1(Tb * restrict a, double b)
  { Tv v=vload(b); for (int i=0; i<nvec; ++i) vmuleq(a->v[i],v); }

static inline Tb Tbprod(Tb a, Tb b)
  { Tb r; for (int i=0; i<nvec; ++i) r.v[i]=vmul(a.v[i],b.v[i]); return r; }

static inline void Tbmuleq(Tb * restrict a, Tb b)
  { for (int i=0; i<nvec; ++i) vmuleq(a->v[i],b.v[i]); }

static void Tbnormalize (Tv * restrict val, Tv * restrict scale,
  double maxval)
  {
  const Tv vfmin=vload(sharp_fsmall*maxval), vfmax=vload(maxval);
  const Tv vfsmall=vload(sharp_fsmall), vfbig=vload(sharp_fbig);
  Tm mask = vgt(vabs(*val),vfmax);
  while (vanyTrue(mask))
    {
    vmuleq_mask(mask,*val,vfsmall);
    vaddeq_mask(mask,*scale,vone);
    mask = vgt(vabs(*val),vfmax);
    }
  mask = vand_mask(vlt(vabs(*val),vfmin),vne(*val,vzero));
  while (vanyTrue(mask))
    {
    vmuleq_mask(mask,*val,vfbig);
    vsubeq_mask(mask,*scale,vone);
    mask = vand_mask(vlt(vabs(*val),vfmin),vne(*val,vzero));
    }
  }

NOINLINE static void mypow (Tb val, int npow, const double * restrict powlimit,
  Tb * restrict resd, Tb * restrict ress)
  {
  Tv vminv=vload(powlimit[npow]);
  int npsave=npow;
  for (int i=0;i<nvec; ++i)
    {
    npow=npsave;
    Tv res=vone;
    Tm mask = vlt(vabs(val.v[i]),vminv);
    if (!vanyTrue(mask)) // no underflows possible, use quick algoritm
      {
      Tv res=vone;
      do
        {
        if (npow&1)
          vmuleq(res,val.v[i]);
        vmuleq(val.v[i],val.v[i]);
        }
      while(npow>>=1);
      resd->v[i]=res;
      ress->v[i]=vzero;
      }
    else
      {
      Tv scale=vzero, scaleint=vzero, res=vone;
      Tbnormalize(&val.v[i],&scaleint,sharp_fbighalf);
      do
        {
        if (npow&1)
          {
          vmuleq(res,val.v[i]);
          vaddeq(scale,scaleint);
          Tbnormalize(&res,&scale,sharp_fbighalf);
          }
        vmuleq(val.v[i],val.v[i]);
        vaddeq(scaleint,scaleint);
        Tbnormalize(&val.v[i],&scaleint,sharp_fbighalf);
        }
      while(npow>>=1);
      resd->v[i]=res;
      ress->v[i]=scale;
      }
    }
  }

static inline int TballLt(Tb a,double b)
  {
  Tv vb=vload(b);
  Tm res=vlt(a.v[0],vb);
  for (int i=1; i<nvec; ++i)
    res=vand_mask(res,vlt(a.v[i],vb));
  return vallTrue(res);
  }
static inline int TballGt(Tb a,double b)
  {
  Tv vb=vload(b);
  Tm res=vgt(a.v[0],vb);
  for (int i=1; i<nvec; ++i)
    res=vand_mask(res,vgt(a.v[i],vb));
  return vallTrue(res);
  }
static inline int TballGe(Tb a,double b)
  {
  Tv vb=vload(b);
  Tm res=vge(a.v[0],vb);
  for (int i=1; i<nvec; ++i)
    res=vand_mask(res,vge(a.v[i],vb));
  return vallTrue(res);
  }

static void getCorfac(Tb scale, Tb * restrict corfac,
  const double * restrict cf)
  {
  Tbu sc, corf;
  sc.b=scale;
  for (int i=0; i<VLEN*nvec; ++i)
    corf.s[i] = (sc.s[i]<sharp_minscale) ?
      0. : cf[(int)(sc.s[i])-sharp_minscale];
  *corfac=corf.b;
  }

NOINLINE static void iter_to_ieee (const Tb sth, Tb cth, int *l_,
  Tb * restrict lam_1, Tb * restrict lam_2, Tb * restrict scale,
  const sharp_Ylmgen_C * restrict gen)
  {
  int l=gen->m;
  for (int i=0; i<nvec; ++i) lam_1->v[i]=vzero;
  mypow(sth,l,gen->powlimit,lam_2,scale);
  Tbmuleq1(lam_2,(gen->m&1) ? -gen->mfac[gen->m]:gen->mfac[gen->m]);
  for (int i=0; i<nvec; ++i)
    Tbnormalize(&lam_2->v[i],&scale->v[i],sharp_ftol);
  Tv fsmall=vload(sharp_fsmall), limscale=vload(sharp_limscale);

  int below_limit = TballLt(*scale,sharp_limscale);
  while (below_limit)
    {
    if (l+2>gen->lmax) {*l_=gen->lmax+1;return;}
    below_limit=1;
    Tv r10=vload(gen->rf[l  ].f[0]), r11=vload(gen->rf[l  ].f[1]),
       r20=vload(gen->rf[l+1].f[0]), r21=vload(gen->rf[l+1].f[1]);
    for (int i=0; i<nvec; ++i)
      {
      lam_1->v[i] = r10*cth.v[i]*lam_2->v[i] - r11*lam_1->v[i];
      lam_2->v[i] = r20*cth.v[i]*lam_1->v[i] - r21*lam_2->v[i];
      Tm mask = vgt(vabs(lam_2->v[i]),vload(sharp_ftol));
      if (vanyTrue(mask))
        {
        vmuleq_mask(mask,lam_1->v[i],fsmall);
        vmuleq_mask(mask,lam_2->v[i],fsmall);
        vaddeq_mask(mask,scale->v[i],vone);
        below_limit &= vallTrue(vlt(scale->v[i],limscale));
        }
      }
    l+=2;
    }
  *l_=l;
  }


NOINLINE static void alm2map_kernel(const Tb cth, Tbri * restrict p1,
  Tbri * restrict p2, Tb lam_1, Tb lam_2,
  const sharp_ylmgen_dbl2 * restrict rf, const dcmplx * restrict alm,
  int l, int lmax)
  {
  while (l<=lmax)
    {
    Tv ar1=vload(creal(alm[l  ])), ai1=vload(cimag(alm[l  ]));
    Tv ar2=vload(creal(alm[l+1])), ai2=vload(cimag(alm[l+1]));
    Tv f10=vload(rf[l  ].f[0]), f11=vload(rf[l  ].f[1]),
       f20=vload(rf[l+1].f[0]), f21=vload(rf[l+1].f[1]);
    for (int i=0; i<nvec; ++i)
      {
      lam_1.v[i] = f10*cth.v[i]*lam_2.v[i] - f11*lam_1.v[i];
      vfmaeq(p1->r.v[i],lam_2.v[i],ar1);
      vfmaeq(p1->i.v[i],lam_2.v[i],ai1);
      lam_2.v[i] = f20*cth.v[i]*lam_1.v[i] - f21*lam_2.v[i];
      vfmaeq(p2->r.v[i],lam_1.v[i],ar2);
      vfmaeq(p2->i.v[i],lam_1.v[i],ai2);
      }
    l+=2;
    }
  }

NOINLINE static void map2alm_kernel (const Tb cth,
  const Tbri * restrict p1, const Tbri * restrict p2, Tb lam_1, Tb lam_2,
  const sharp_ylmgen_dbl2 * restrict rf, int l, int lmax, Tv *restrict atmp)
  {
  while (l<=lmax)
    {
    Tv f10=vload(rf[l  ].f[0]), f11=vload(rf[l  ].f[1]),
       f20=vload(rf[l+1].f[0]), f21=vload(rf[l+1].f[1]);
    for (int i=0; i<nvec; ++i)
      {
      lam_1.v[i] = f10*cth.v[i]*lam_2.v[i] - f11*lam_1.v[i];
      vfmaeq(atmp[2*l  ],lam_2.v[i],p1->r.v[i]);
      vfmaeq(atmp[2*l+1],lam_2.v[i],p1->i.v[i]);
      lam_2.v[i] = f20*cth.v[i]*lam_1.v[i] - f21*lam_2.v[i];
      vfmaeq(atmp[2*(l+1)  ],lam_1.v[i],p2->r.v[i]);
      vfmaeq(atmp[2*(l+1)+1],lam_1.v[i],p2->i.v[i]);
      }
    l+=2;
    }
  }

NOINLINE static void calc_alm2map (const Tb cth, const Tb sth,
  const sharp_Ylmgen_C *gen, sharp_job *job, Tbri * restrict p1,
  Tbri * restrict p2)
  {
  int l,lmax=gen->lmax;
  Tb lam_1=Tbconst(0.),lam_2=Tbconst(0.),scale;
  iter_to_ieee(sth,cth,&l,&lam_1,&lam_2,&scale,gen);
  job->opcnt += (l-gen->m) * 4*VLEN*nvec;
  if (l>lmax) return;
  job->opcnt += (lmax+1-l) * 8*VLEN*nvec;

  Tb corfac;
  getCorfac(scale,&corfac,gen->cf);
  const sharp_ylmgen_dbl2 * restrict rf = gen->rf;
  const dcmplx * restrict alm=job->almtmp;
  int full_ieee = TballGe(scale,sharp_minscale);
  while((!full_ieee) && (l<=lmax))
    {
    Tv ar1=vload(creal(alm[l  ])), ai1=vload(cimag(alm[l  ]));
    Tv ar2=vload(creal(alm[l+1])), ai2=vload(cimag(alm[l+1]));
    Tv f10=vload(rf[l  ].f[0]), f11=vload(rf[l  ].f[1]),
       f20=vload(rf[l+1].f[0]), f21=vload(rf[l+1].f[1]);
    full_ieee=1;
    for (int i=0; i<nvec; ++i)
      {
      lam_1.v[i] = f10*cth.v[i]*lam_2.v[i] - f11*lam_1.v[i];
      vfmaeq(p1->r.v[i],lam_2.v[i]*corfac.v[i],ar1);
      vfmaeq(p1->i.v[i],lam_2.v[i]*corfac.v[i],ai1);
      lam_2.v[i] = f20*cth.v[i]*lam_1.v[i] - f21*lam_2.v[i];
      Tm mask = vgt(vabs(lam_2.v[i]),vload(sharp_ftol));
      if (vanyTrue(mask))
        {
        vmuleq_mask(mask,lam_1.v[i],vload(sharp_fsmall));
        vmuleq_mask(mask,lam_2.v[i],vload(sharp_fsmall));
        vaddeq_mask(mask,scale.v[i],vone);
        Tvu sc, corf;
        sc.v=scale.v[i];
        for (int j=0; j<VLEN; ++j)
          corf.s[j] = (sc.s[j]<sharp_minscale) ?
            0. : gen->cf[(int)(sc.s[j])-sharp_minscale];
        corfac.v[i]=corf.v;
        full_ieee &= vallTrue(vge(scale.v[i],vload(sharp_minscale)));
        }
      vfmaeq(p2->r.v[i],lam_1.v[i]*corfac.v[i],ar2);
      vfmaeq(p2->i.v[i],lam_1.v[i]*corfac.v[i],ai2);
      }
    l+=2;
    }
  if (l>lmax) return;

  Tbmuleq(&lam_1,corfac); Tbmuleq(&lam_2,corfac);
  alm2map_kernel(cth, p1, p2, lam_1, lam_2, rf, alm, l, lmax);
  }

NOINLINE static void calc_map2alm(const Tb cth, const Tb sth,
  const sharp_Ylmgen_C *gen, sharp_job *job, const Tbri * restrict p1,
  const Tbri * restrict p2, Tv *restrict atmp)
  {
  int lmax=gen->lmax;
  Tb lam_1=Tbconst(0.),lam_2=Tbconst(0.),scale;
  int l=gen->m;
  iter_to_ieee(sth,cth,&l,&lam_1,&lam_2,&scale,gen);
  job->opcnt += (l-gen->m) * 4*VLEN*nvec;
  if (l>lmax) return;
  job->opcnt += (lmax+1-l) * 8*VLEN*nvec;

  const sharp_ylmgen_dbl2 * restrict rf = gen->rf;
  Tb corfac;
  getCorfac(scale,&corfac,gen->cf);
  int full_ieee = TballGe(scale,sharp_minscale);
  while ((!full_ieee) && (l<=lmax))
    {
    full_ieee=1;
    Tv f10=vload(rf[l  ].f[0]), f11=vload(rf[l  ].f[1]),
       f20=vload(rf[l+1].f[0]), f21=vload(rf[l+1].f[1]);
    for (int i=0; i<nvec; ++i)
      {
      lam_1.v[i] = f10*cth.v[i]*lam_2.v[i] - f11*lam_1.v[i];
      vfmaeq(atmp[2*l  ],lam_2.v[i]*corfac.v[i],p1->r.v[i]);
      vfmaeq(atmp[2*l+1],lam_2.v[i]*corfac.v[i],p1->i.v[i]);
      lam_2.v[i] = f20*cth.v[i]*lam_1.v[i] - f21*lam_2.v[i];
      Tm mask = vgt(vabs(lam_2.v[i]),vload(sharp_ftol));
      if (vanyTrue(mask))
        {
        vmuleq_mask(mask,lam_1.v[i],vload(sharp_fsmall));
        vmuleq_mask(mask,lam_2.v[i],vload(sharp_fsmall));
        vaddeq_mask(mask,scale.v[i],vone);
        Tvu sc, corf;
        sc.v=scale.v[i];
        for (int j=0; j<VLEN; ++j)
          corf.s[j] = (sc.s[j]<sharp_minscale) ?
            0. : gen->cf[(int)(sc.s[j])-sharp_minscale];
        corfac.v[i]=corf.v;
        full_ieee &= vallTrue(vge(scale.v[i],vload(sharp_minscale)));
        }
      vfmaeq(atmp[2*(l+1)  ],lam_1.v[i]*corfac.v[i],p2->r.v[i]);
      vfmaeq(atmp[2*(l+1)+1],lam_1.v[i]*corfac.v[i],p2->i.v[i]);
      }
    l+=2;
    }

  Tbmuleq(&lam_1,corfac); Tbmuleq(&lam_2,corfac);
  map2alm_kernel(cth, p1, p2, lam_1, lam_2, rf, l, lmax, atmp);
  }


#define VZERO(var) do { memset(&(var),0,sizeof(var)); } while(0)

NOINLINE static void inner_loop_a2m(sharp_job *job, const int *ispair,
  const double *cth_, const double *sth_, int llim, int ulim,
  sharp_Ylmgen_C *gen, int mi, const int *mlim)
  {
  const int nval=nvec*VLEN;
  const int m = job->ainfo->mval[mi];
  sharp_Ylmgen_prepare (gen, m);

  switch (job->type)
    {
    case SHARP_ALM2MAP:
    case SHARP_ALM2MAP_DERIV1:
      {
      if (job->spin==0)
        {
        for (int ith=0; ith<ulim-llim; ith+=nval)
          {
          Tburi p1,p2; VZERO(p1); VZERO(p2);
          Tbu cth, sth;

          int skip=1;
          for (int i=0; i<nval; ++i)
            {
            int itot=i+ith;
            if (itot>=ulim-llim) itot=ulim-llim-1;
            if (mlim[itot]>=m) skip=0;
            cth.s[i]=cth_[itot]; sth.s[i]=sth_[itot];
            }
          if (!skip)
            calc_alm2map (cth.b,sth.b,gen,job,&p1.b,&p2.b);

          for (int i=0; i<nval; ++i)
            {
            int itot=i+ith;
            if (itot<ulim-llim)
              {
              int phas_idx = itot*job->s_th + mi*job->s_m;
              complex double r1 = p1.s.r[i] + p1.s.i[i]*_Complex_I,
                             r2 = p2.s.r[i] + p2.s.i[i]*_Complex_I;
              job->phase[phas_idx] = r1+r2;
              if (ispair[itot])
                job->phase[phas_idx+1] = r1-r2;
              }
            }
          }
        }
      else
        {
        UTIL_FAIL("only spin==0 allowed at the moment");
        }
      break;
      }
    default:
      {
      UTIL_FAIL("must not happen");
      break;
      }
    }
  }

NOINLINE static void inner_loop_m2a(sharp_job *job, const int *ispair,
  const double *cth_, const double *sth_, int llim, int ulim,
  sharp_Ylmgen_C *gen, int mi, const int *mlim)
  {
  const int nval=nvec*VLEN;
  const int m = job->ainfo->mval[mi];
  sharp_Ylmgen_prepare (gen, m);

  switch (job->type)
    {
    case SHARP_MAP2ALM:
      {
      if (job->spin==0)
        {
        Tv atmp[2*(gen->lmax+2)];
        memset (&atmp[2*m],0,2*(gen->lmax+2-m)*sizeof(Tv));
        for (int ith=0; ith<ulim-llim; ith+=nval)
          {
          Tburi p1, p2; VZERO(p1); VZERO(p2);
          Tbu cth, sth;
          int skip=1;

          for (int i=0; i<nval; ++i)
            {
            int itot=i+ith;
            if (itot>=ulim-llim) itot=ulim-llim-1;
            if (mlim[itot]>=m) skip=0;
            cth.s[i]=cth_[itot]; sth.s[i]=sth_[itot];
            if ((i+ith<ulim-llim)&&(mlim[itot]>=m))
              {
              int phas_idx = itot*job->s_th + mi*job->s_m;
              dcmplx ph1=job->phase[phas_idx];
              dcmplx ph2=ispair[itot] ? job->phase[phas_idx+1] : 0.;
              p1.s.r[i]=creal(ph1+ph2); p1.s.i[i]=cimag(ph1+ph2);
              p2.s.r[i]=creal(ph1-ph2); p2.s.i[i]=cimag(ph1-ph2);
              }
            }
          if (!skip)
            calc_map2alm(cth.b,sth.b,gen,job,&p1.b,&p2.b, atmp);
          }
        {
        int istart=m, istop=gen->lmax+1;
        for(; istart<istop-2; istart+=2)
          vhsum_cmplx_special(atmp[2*istart],atmp[2*istart+1],atmp[2*istart+2],atmp[2*istart+3],&(job->almtmp[istart]));
        for(; istart<istop; istart++)
          job->almtmp[istart]+=vhsum_cmplx(atmp[2*istart],atmp[2*istart+1]);
        }
        }
      else
        {
        UTIL_FAIL("only spin==0 allowed at the moment");
        }
      break;
      }
    default:
      {
      UTIL_FAIL("must not happen");
      break;
      }
    }
  }

void inner_loop (sharp_job *job, const int *ispair,
  const double *cth_, const double *sth_, int llim, int ulim,
  sharp_Ylmgen_C *gen, int mi, const int *mlim)
  {
  (job->type==SHARP_MAP2ALM) ?
    inner_loop_m2a(job,ispair,cth_,sth_,llim,ulim,gen,mi,mlim) :
    inner_loop_a2m(job,ispair,cth_,sth_,llim,ulim,gen,mi,mlim);
  }

#undef VZERO

int sharp_veclen(void)
  {
  return VLEN;
  }

int sharp_max_nvec(void)
  {
  return nvec;
  }
