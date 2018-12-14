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
  { double r[VLEN*nvec], i[VLEN*nvec]; } Tsri;

typedef union
  { Tbri b; Tsri s; } Tburi;

static inline void Tvnormalize (Tv * restrict val, Tv * restrict scale,
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

static void mypow(Tv val, int npow, const double * restrict powlimit,
  Tv * restrict resd, Tv * restrict ress)
  {
  Tv vminv=vload(powlimit[npow]);
  Tm mask = vlt(vabs(val),vminv);
  if (!vanyTrue(mask)) // no underflows possible, use quick algoritm
    {
    Tv res=vone;
    do
      {
      if (npow&1)
        vmuleq(res,val);
      vmuleq(val,val);
      }
    while(npow>>=1);
    *resd=res;
    *ress=vzero;
    }
  else
    {
    Tv scale=vzero, scaleint=vzero, res=vone;
    Tvnormalize(&val,&scaleint,sharp_fbighalf);
    do
      {
      if (npow&1)
        {
        vmuleq(res,val);
        vaddeq(scale,scaleint);
        Tvnormalize(&res,&scale,sharp_fbighalf);
        }
      vmuleq(val,val);
      vaddeq(scaleint,scaleint);
      Tvnormalize(&val,&scaleint,sharp_fbighalf);
      }
    while(npow>>=1);
    *resd=res;
    *ress=scale;
    }
  }

static inline void getCorfac(Tv scale, Tv * restrict corfac,
  const double * restrict cf)
  {
  Tvu sc, corf;
  sc.v=scale;
  for (int i=0; i<VLEN; ++i)
    corf.s[i] = (sc.s[i]<sharp_minscale) ?
      0. : cf[(int)(sc.s[i])-sharp_minscale];
  *corfac=corf.v;
  }

static inline int rescale(Tv * restrict v1, Tv * restrict v2, Tv * restrict s, Tv eps)
  {
  Tm mask = vgt(vabs(*v2),eps);
  if (vanyTrue(mask))
    {
    vmuleq_mask(mask,*v1,vload(sharp_fsmall));
    vmuleq_mask(mask,*v2,vload(sharp_fsmall));
    vaddeq_mask(mask,*s,vone);
    return 1;
    }
  return 0;
  }

NOINLINE static void iter_to_ieee (const Tb sth, Tb cth, int *l_,
  Tb * restrict lam_1, Tb * restrict lam_2, Tb * restrict scale,
  const sharp_Ylmgen_C * restrict gen, int nv2)
  {
  int l=gen->m;
  Tv mfac = vload((gen->m&1) ? -gen->mfac[gen->m]:gen->mfac[gen->m]);
  Tv limscale=vload(sharp_limscale);
  int below_limit = 1;
  for (int i=0; i<nv2; ++i)
    {
    lam_1->v[i]=vzero;
    mypow(sth.v[i],l,gen->powlimit,&lam_2->v[i],&scale->v[i]);
    lam_2->v[i] *= mfac;
    Tvnormalize(&lam_2->v[i],&scale->v[i],sharp_ftol);
    below_limit &= vallTrue(vlt(scale->v[i],limscale));
    }

  while (below_limit)
    {
    if (l+2>gen->lmax) {*l_=gen->lmax+1;return;}
    below_limit=1;
    Tv r10=vload(gen->rf[l  ].f[0]), r11=vload(gen->rf[l  ].f[1]),
       r20=vload(gen->rf[l+1].f[0]), r21=vload(gen->rf[l+1].f[1]);
    for (int i=0; i<nv2; ++i)
      {
      lam_1->v[i] = r10*cth.v[i]*lam_2->v[i] - r11*lam_1->v[i];
      lam_2->v[i] = r20*cth.v[i]*lam_1->v[i] - r21*lam_2->v[i];
      if (rescale(&lam_1->v[i], &lam_2->v[i], &scale->v[i], vload(sharp_ftol)))
        below_limit &= vallTrue(vlt(scale->v[i],limscale));
      }
    l+=2;
    }
  *l_=l;
  }

#if 1
static inline void rec_step (Tv * restrict rxp, Tv * restrict rxm,
  Tv * restrict ryp, Tv * restrict rym, const Tv cth,
  const sharp_ylmgen_dbl3 fx)
  {
  Tv fx0=vload(fx.f[0]),fx1=vload(fx.f[1]),fx2=vload(fx.f[2]);
  *rxp = (cth-fx1)*fx0* *ryp - fx2* *rxp;
  *rxm = (cth+fx1)*fx0* *rym - fx2* *rxm;
  }

NOINLINE static void iter_to_ieee_spin (const Tb cth, const Tb sth, int *l_,
  Tb * rec1p, Tb * rec1m, Tb * rec2p, Tb * rec2m,
  Tb * scalep, Tb * scalem, const sharp_Ylmgen_C * restrict gen, int nv2)
  {
  const sharp_ylmgen_dbl3 * restrict fx = gen->fx;
  Tb ccp, ccps, ssp, ssps, csp, csps, scp, scps;
  Tv prefac=vload(gen->prefac[gen->m]),
     prescale=vload(gen->fscale[gen->m]);
  Tv limscale=vload(sharp_limscale);
  int below_limit=1;
  for (int i=0; i<nv2; ++i)
    {
    Tv cth2, sth2;
    cth2=vsqrt(vmul(vadd(vone,cth.v[i]),vload(0.5)));
    cth2=vmax(cth2,vload(1e-15));
    sth2=vsqrt(vmul(vsub(vone,cth.v[i]),vload(0.5)));
    sth2=vmax(sth2,vload(1e-15));
    Tm mask=vlt(sth.v[i],vzero);
    Tm cmask=vand_mask(mask,vlt(cth.v[i],vzero));
    vmuleq_mask(cmask,cth2,vload(-1.));
    Tm smask=vand_mask(mask,vgt(cth.v[i],vzero));
    vmuleq_mask(smask,sth2,vload(-1.));

    mypow(cth2,gen->cosPow,gen->powlimit,&ccp.v[i],&ccps.v[i]);
    mypow(sth2,gen->sinPow,gen->powlimit,&ssp.v[i],&ssps.v[i]);
    mypow(cth2,gen->sinPow,gen->powlimit,&csp.v[i],&csps.v[i]);
    mypow(sth2,gen->cosPow,gen->powlimit,&scp.v[i],&scps.v[i]);

    rec1p->v[i] = vzero;
    rec1m->v[i] = vzero;
    rec2p->v[i]=vmul(prefac,ccp.v[i]);
    scalep->v[i]=vadd(prescale,ccps.v[i]);
    rec2m->v[i]=vmul(prefac,csp.v[i]);
    scalem->v[i]=vadd(prescale,csps.v[i]);
    Tvnormalize(&rec2m->v[i],&scalem->v[i],sharp_fbighalf);
    Tvnormalize(&rec2p->v[i],&scalep->v[i],sharp_fbighalf);

    rec2p->v[i]=vmul(rec2p->v[i],ssp.v[i]);
    scalep->v[i]=vadd(scalep->v[i],ssps.v[i]);
    rec2m.v[i]=vmul(rec2m.v[i],scp.v[i]);
    scalem.v[i]=vadd(scalem.v[i],scps.v[i]);
    if (gen->preMinus_p)
      rec2p.v[i]=vneg(rec2p.v[i]);
    if (gen->preMinus_m)
      rec2m.v[i]=vneg(rec2m.v[i]);
    if (gen->s&1)
      rec2p.v[i]=vneg(rec2p.v[i]);

    Tvnormalize(&rec2m.v[i],&scalem.v[i],sharp_ftol);
    Tvnormalize(&rec2p.v[i],&scalep.v[i],sharp_ftol);

    below_limit &= vallTrue(vand_mask(vlt(scalem.v[i],limscale),vlt(scalep.v[i],limscale)));
    }

  int l=gen->mhi;

  while (below_limit)
    {
    if (l+2>gen->lmax) {*l_=gen->lmax+1;return;}
    for (int i=0; i<nv2; ++i)
      {
      rec_step(&rec1p.v[i],&rec1m.v[i],&rec2p.v[i],&rec2m.v[i],cth.v[i],fx[l+1]);
      rec_step(&rec2p.v[i],&rec2m.v[i],&rec1p.v[i],&rec1m.v[i],cth.v[i],fx[l+2]);
      if (rescale(&rec1p.v[i],&rec2p.v[i],&scalep.v[i],vload(sharp_ftol)) ||
          rescale(&rec1m.v[i],&rec2m.v[i],&scalem.v[i],vload(sharp_ftol)))
      below_limit &= vallTrue(vlt(scalep.v[i],limscale)) &&
                     vallTrue(vlt(scalem.v[i],limscale));
      }
    l+=2;
    }

  *l_=l;
  *rec1p_=rec1p; *rec2p_=rec2p; *scalep_=scalep;
  *rec1m_=rec1m; *rec2m_=rec2m; *scalem_=scalem;
  }

NOINLINE static void alm2map_spin_kernel(Tb cth, Tbqu * restrict p1,
  Tbqu * restrict p2, Tb rec1p, Tb rec1m, Tb rec2p, Tb rec2m,
  const sharp_ylmgen_dbl3 * restrict fx, const dcmplx * restrict alm, int l,
  int lmax, int nv2)
  {
  while (l<=lmax)
    {
    Tv fx10=vload(fx[l+1].f[0]),fx1=v1load(fx[l+1].f[1]),
       fx12=vload(fx[l+1].f[2]);
    for (int i=0; i<nvec; ++i)
      {
      rec1p.v[i] = (cth.v[i]-fx1)*fx0*rec2p.v[i] - fx2*rec1p.v[i];
      rec1m.v[i] = (cth.v[i]+fx1)*fx0*rec2m.v[i] - fx2*rec1m.v[i];
      }
    Z(saddstepb)(p1,p2,rec1p,rec1m,rec2p,rec2m,&alm[2*njobs*l],
      &alm[2*njobs*(l+1)] NJ2);
    fx0=vload(fx[l+2].f[0]);fx1=vload(fx[l+2].f[1]);
    fx2=vload(fx[l+2].f[2]);
    for (int i=0; i<nvec; ++i)
      {
      rec2p.v[i] = (cth.v[i]-fx1)*fx0*rec1p.v[i] - fx2*rec2p.v[i];
      rec2m.v[i] = (cth.v[i]+fx1)*fx0*rec1m.v[i] - fx2*rec2m.v[i];
      }
    l+=2;
    }
  if (l==lmax)
    Z(saddstep)(p1, p2, rec2p, rec2m, &alm[2*njobs*l] NJ2);
  }
#endif

NOINLINE static void alm2map_kernel(const Tb cth, Tbri * restrict p1,
  Tbri * restrict p2, Tb lam_1, Tb lam_2,
  const sharp_ylmgen_dbl2 * restrict rf, const dcmplx * restrict alm,
  int l, int lmax, int nv2)
  {
  while (l<=lmax)
    {
    Tv ar1=vload(creal(alm[l  ])), ai1=vload(cimag(alm[l  ]));
    Tv ar2=vload(creal(alm[l+1])), ai2=vload(cimag(alm[l+1]));
    Tv f10=vload(rf[l  ].f[0]), f11=vload(rf[l  ].f[1]),
       f20=vload(rf[l+1].f[0]), f21=vload(rf[l+1].f[1]);
    for (int i=0; i<nv2; ++i)
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
  const sharp_ylmgen_dbl2 * restrict rf, dcmplx * restrict alm, int l,
  int lmax, int nv2)
  {
  while (l<=lmax)
    {
    Tv f10=vload(rf[l  ].f[0]), f11=vload(rf[l  ].f[1]),
       f20=vload(rf[l+1].f[0]), f21=vload(rf[l+1].f[1]);
    Tv atmp[4] = {vzero, vzero, vzero, vzero};
    for (int i=0; i<nv2; ++i)
      {
      lam_1.v[i] = f10*cth.v[i]*lam_2.v[i] - f11*lam_1.v[i];
      vfmaeq(atmp[0],lam_2.v[i],p1->r.v[i]);
      vfmaeq(atmp[1],lam_2.v[i],p1->i.v[i]);
      lam_2.v[i] = f20*cth.v[i]*lam_1.v[i] - f21*lam_2.v[i];
      vfmaeq(atmp[2],lam_1.v[i],p2->r.v[i]);
      vfmaeq(atmp[3],lam_1.v[i],p2->i.v[i]);
      }
    vhsum_cmplx_special (atmp[0], atmp[1], atmp[2], atmp[3], &alm[l]);
    l+=2;
    }
  }

NOINLINE static void calc_alm2map (const Tb cth, const Tb sth,
  const sharp_Ylmgen_C *gen, sharp_job *job, Tbri * restrict p1,
  Tbri * restrict p2, int nth)
  {
  int l,lmax=gen->lmax;
  Tb lam_1,lam_2,scale;
  int nv2 = (nth+VLEN-1)/VLEN;
  iter_to_ieee(sth,cth,&l,&lam_1,&lam_2,&scale,gen,nv2);
  job->opcnt += (l-gen->m) * 4*nth;
  if (l>lmax) return;
  job->opcnt += (lmax+1-l) * 8*nth;

  const sharp_ylmgen_dbl2 * restrict rf = gen->rf;
  const dcmplx * restrict alm=job->almtmp;
  int full_ieee=1;
  Tb corfac;
  for (int i=0; i<nv2; ++i)
    {
    getCorfac(scale.v[i], &corfac.v[i], gen->cf);
    full_ieee &= vallTrue(vge(scale.v[i],vload(sharp_minscale)));
    }

  while((!full_ieee) && (l<=lmax))
    {
    Tv ar1=vload(creal(alm[l  ])), ai1=vload(cimag(alm[l  ]));
    Tv ar2=vload(creal(alm[l+1])), ai2=vload(cimag(alm[l+1]));
    Tv f10=vload(rf[l  ].f[0]), f11=vload(rf[l  ].f[1]),
       f20=vload(rf[l+1].f[0]), f21=vload(rf[l+1].f[1]);
    full_ieee=1;
    for (int i=0; i<nv2; ++i)
      {
      lam_1.v[i] = f10*cth.v[i]*lam_2.v[i] - f11*lam_1.v[i];
      vfmaeq(p1->r.v[i],lam_2.v[i]*corfac.v[i],ar1);
      vfmaeq(p1->i.v[i],lam_2.v[i]*corfac.v[i],ai1);
      lam_2.v[i] = f20*cth.v[i]*lam_1.v[i] - f21*lam_2.v[i];
      if (rescale(&lam_1.v[i], &lam_2.v[i], &scale.v[i], vload(sharp_ftol)))
        {
        getCorfac(scale.v[i], &corfac.v[i], gen->cf);
        full_ieee &= vallTrue(vge(scale.v[i],vload(sharp_minscale)));
        }
      vfmaeq(p2->r.v[i],lam_1.v[i]*corfac.v[i],ar2);
      vfmaeq(p2->i.v[i],lam_1.v[i]*corfac.v[i],ai2);
      }
    l+=2;
    }
  if (l>lmax) return;

  for (int i=0; i<nv2; ++i)
    {
    lam_1.v[i] *= corfac.v[i];
    lam_2.v[i] *= corfac.v[i];
    }
  alm2map_kernel(cth, p1, p2, lam_1, lam_2, rf, alm, l, lmax, nv2);
  }

NOINLINE static void calc_map2alm(const Tb cth, const Tb sth,
  const sharp_Ylmgen_C *gen, sharp_job *job, const Tbri * restrict p1,
  const Tbri * restrict p2, int nth)
  {
  int lmax=gen->lmax;
  Tb lam_1,lam_2,scale;
  int l=gen->m;
  int nv2 = (nth+VLEN-1)/VLEN;
  iter_to_ieee(sth,cth,&l,&lam_1,&lam_2,&scale,gen,nv2);
  job->opcnt += (l-gen->m) * 4*nth;
  if (l>lmax) return;
  job->opcnt += (lmax+1-l) * 8*nth;

  const sharp_ylmgen_dbl2 * restrict rf = gen->rf;
  dcmplx * restrict alm=job->almtmp;
  int full_ieee=1;
  Tb corfac;
  for (int i=0; i<nv2; ++i)
    {
    getCorfac(scale.v[i], &corfac.v[i], gen->cf);
    full_ieee &= vallTrue(vge(scale.v[i],vload(sharp_minscale)));
    }

  while ((!full_ieee) && (l<=lmax))
    {
    full_ieee=1;
    Tv f10=vload(rf[l  ].f[0]), f11=vload(rf[l  ].f[1]),
       f20=vload(rf[l+1].f[0]), f21=vload(rf[l+1].f[1]);
    Tv atmp[4] = {vzero, vzero, vzero, vzero};
    for (int i=0; i<nv2; ++i)
      {
      lam_1.v[i] = f10*cth.v[i]*lam_2.v[i] - f11*lam_1.v[i];
      vfmaeq(atmp[0],lam_2.v[i]*corfac.v[i],p1->r.v[i]);
      vfmaeq(atmp[1],lam_2.v[i]*corfac.v[i],p1->i.v[i]);
      lam_2.v[i] = f20*cth.v[i]*lam_1.v[i] - f21*lam_2.v[i];
      if (rescale(&lam_1.v[i], &lam_2.v[i], &scale.v[i], vload(sharp_ftol)))
        {
        getCorfac(scale.v[i], &corfac.v[i], gen->cf);
        full_ieee &= vallTrue(vge(scale.v[i],vload(sharp_minscale)));
        }
      vfmaeq(atmp[2],lam_1.v[i]*corfac.v[i],p2->r.v[i]);
      vfmaeq(atmp[3],lam_1.v[i]*corfac.v[i],p2->i.v[i]);
      }
    vhsum_cmplx_special (atmp[0], atmp[1], atmp[2], atmp[3], &alm[l]);
    l+=2;
    }

  for (int i=0; i<nv2; ++i)
    {
    lam_1.v[i] *= corfac.v[i];
    lam_2.v[i] *= corfac.v[i];
    }
  map2alm_kernel(cth, p1, p2, lam_1, lam_2, rf, alm, l, lmax, nv2);
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
        int ith=0;
        int itgt[nvec*VLEN];
        while (ith<ulim-llim)
          {
          Tburi p1,p2; VZERO(p1); VZERO(p2);
          Tbu cth, sth;
          int nth=0;
          while ((nth<nval)&&(ith<ulim-llim))
            {
            if (mlim[ith]>=m)
              {
              itgt[nth] = ith;
              cth.s[nth]=cth_[ith]; sth.s[nth]=sth_[ith];
              ++nth;
              }
            else
              {
              int phas_idx = ith*job->s_th + mi*job->s_m;
              job->phase[phas_idx] = job->phase[phas_idx+1] = 0;
              }
            ++ith;
            }
          if (nth>0)
            {
            int i2=((nth+VLEN-1)/VLEN)*VLEN;
            for (int i=nth; i<i2; ++i)
              {
              cth.s[i]=cth.s[nth-1];
              sth.s[i]=sth.s[nth-1];
              }
            calc_alm2map (cth.b,sth.b,gen,job,&p1.b,&p2.b,nth);
            for (int i=0; i<nth; ++i)
              {
              int tgt=itgt[i];
              int phas_idx = tgt*job->s_th + mi*job->s_m;
              complex double r1 = p1.s.r[i] + p1.s.i[i]*_Complex_I,
                             r2 = p2.s.r[i] + p2.s.i[i]*_Complex_I;
              job->phase[phas_idx] = r1+r2;
              if (ispair[tgt])
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
        int ith=0;
        while (ith<ulim-llim)
          {
          Tburi p1,p2; VZERO(p1); VZERO(p2);
          Tbu cth, sth;
          int nth=0;
          while ((nth<nval)&&(ith<ulim-llim))
            {
            if (mlim[ith]>=m)
              {
              cth.s[nth]=cth_[ith]; sth.s[nth]=sth_[ith];
              int phas_idx = ith*job->s_th + mi*job->s_m;
              dcmplx ph1=job->phase[phas_idx];
              dcmplx ph2=ispair[ith] ? job->phase[phas_idx+1] : 0.;
              p1.s.r[nth]=creal(ph1+ph2); p1.s.i[nth]=cimag(ph1+ph2);
              p2.s.r[nth]=creal(ph1-ph2); p2.s.i[nth]=cimag(ph1-ph2);
              ++nth;
              }
            ++ith;
            }
          if (nth>0)
            {
            int i2=((nth+VLEN-1)/VLEN)*VLEN;
            for (int i=nth; i<i2; ++i)
              {
              cth.s[i]=cth.s[nth-1];
              sth.s[i]=sth.s[nth-1];
              p1.s.r[i]=p1.s.i[i]=p2.s.r[i]=p2.s.i[i]=0.;
              }
            calc_map2alm(cth.b,sth.b,gen,job,&p1.b,&p2.b, nth);
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
