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
#include "sharp_internal.h"
#include "c_utils.h"

typedef complex double dcmplx;

#define nv0 (128/VLEN)
#define nvx (64/VLEN)

typedef Tv Tbv0[nv0];
typedef double Tbs0[nv0*VLEN];

typedef struct
  {
  Tbv0 sth, corfac, scale, lam1, lam2, cth, p1r, p1i, p2r, p2i;
  } s0data_v;

typedef struct
  {
  Tbs0 sth, corfac, scale, lam1, lam2, cth, p1r, p1i, p2r, p2i;
  } s0data_s;

typedef union
  {
  s0data_v v;
  s0data_s s;
  } s0data_u;

typedef Tv Tbvx[nvx];
typedef double Tbsx[nvx*VLEN];

typedef struct
  {
  Tbvx sth, cfp, cfm, scp, scm, l1p, l2p, l1m, l2m, cth,
       p1pr, p1pi, p2pr, p2pi, p1mr, p1mi, p2mr, p2mi;
  } sxdata_v;

typedef struct
  {
  Tbsx sth, cfp, cfm, scp, scm, l1p, l2p, l1m, l2m, cth,
       p1pr, p1pi, p2pr, p2pi, p1mr, p1mi, p2mr, p2mi;
  } sxdata_s;

typedef union
  {
  sxdata_v v;
  sxdata_s s;
  } sxdata_u;

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
        res*=val;
      val*=val;
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
        res*=val;
        scale+=scaleint;
        Tvnormalize(&res,&scale,sharp_fbighalf);
        }
      val*=val;
      scaleint+=scaleint;
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
  typedef union
    { Tv v; double s[VLEN]; } Tvu;

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

NOINLINE static void iter_to_ieee(const sharp_Ylmgen_C * restrict gen,
  s0data_v * restrict d, int * restrict l_, int nv2)
  {
  int l=gen->m;
  Tv mfac = vload((gen->m&1) ? -gen->mfac[gen->m]:gen->mfac[gen->m]);
  Tv limscale=vload(sharp_limscale);
  int below_limit = 1;
  for (int i=0; i<nv2; ++i)
    {
    d->lam1[i]=vzero;
    mypow(d->sth[i],l,gen->powlimit,&d->lam2[i],&d->scale[i]);
    d->lam2[i] *= mfac;
    Tvnormalize(&d->lam2[i],&d->scale[i],sharp_ftol);
    below_limit &= vallTrue(vlt(d->scale[i],limscale));
    }

  while (below_limit)
    {
    if (l+2>gen->lmax) {*l_=gen->lmax+1;return;}
    below_limit=1;
    Tv r10=vload(gen->rf[l  ].f[0]), r11=vload(gen->rf[l  ].f[1]),
       r20=vload(gen->rf[l+1].f[0]), r21=vload(gen->rf[l+1].f[1]);
    for (int i=0; i<nv2; ++i)
      {
      d->lam1[i] = r10*d->cth[i]*d->lam2[i] - r11*d->lam1[i];
      d->lam2[i] = r20*d->cth[i]*d->lam1[i] - r21*d->lam2[i];
      if (rescale(&d->lam1[i], &d->lam2[i], &d->scale[i], vload(sharp_ftol)))
        below_limit &= vallTrue(vlt(d->scale[i],limscale));
      }
    l+=2;
    }
  *l_=l;
  }

NOINLINE static void alm2map_kernel(s0data_v * restrict d,
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
      d->lam1[i] = f10*d->cth[i]*d->lam2[i] - f11*d->lam1[i];
      d->p1r[i] += d->lam2[i]*ar1;
      d->p1i[i] += d->lam2[i]*ai1;
      d->lam2[i] = f20*d->cth[i]*d->lam1[i] - f21*d->lam2[i];
      d->p2r[i] += d->lam1[i]*ar2;
      d->p2i[i] += d->lam1[i]*ai2;
      }
    l+=2;
    }
  }

NOINLINE static void calc_alm2map (sharp_job * restrict job,
  const sharp_Ylmgen_C * restrict gen, s0data_v * restrict d, int nth)
  {
  int l,lmax=gen->lmax;
  int nv2 = (nth+VLEN-1)/VLEN;
  iter_to_ieee(gen, d, &l, nv2);
  job->opcnt += (l-gen->m) * 4*nth;
  if (l>lmax) return;
  job->opcnt += (lmax+1-l) * 8*nth;

  const sharp_ylmgen_dbl2 * restrict rf = gen->rf;
  const dcmplx * restrict alm=job->almtmp;
  int full_ieee=1;
  for (int i=0; i<nv2; ++i)
    {
    getCorfac(d->scale[i], &d->corfac[i], gen->cf);
    full_ieee &= vallTrue(vge(d->scale[i],vload(sharp_minscale)));
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
      d->lam1[i] = f10*d->cth[i]*d->lam2[i] - f11*d->lam1[i];
      d->p1r[i] += d->lam2[i]*d->corfac[i]*ar1;
      d->p1i[i] += d->lam2[i]*d->corfac[i]*ai1;
      d->lam2[i] = f20*d->cth[i]*d->lam1[i] - f21*d->lam2[i];
      if (rescale(&d->lam1[i], &d->lam2[i], &d->scale[i], vload(sharp_ftol)))
        getCorfac(d->scale[i], &d->corfac[i], gen->cf);
      full_ieee &= vallTrue(vge(d->scale[i],vload(sharp_minscale)));
      d->p2r[i] += d->lam1[i]*d->corfac[i]*ar2;
      d->p2i[i] += d->lam1[i]*d->corfac[i]*ai2;
      }
    l+=2;
    }
  if (l>lmax) return;

  for (int i=0; i<nv2; ++i)
    {
    d->lam1[i] *= d->corfac[i];
    d->lam2[i] *= d->corfac[i];
    }
  alm2map_kernel(d, rf, alm, l, lmax, nv2);
  }

NOINLINE static void map2alm_kernel(s0data_v * restrict d,
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
      d->lam1[i] = f10*d->cth[i]*d->lam2[i] - f11*d->lam1[i];
      atmp[0] += d->lam2[i]*d->p1r[i];
      atmp[1] += d->lam2[i]*d->p1i[i];
      d->lam2[i] = f20*d->cth[i]*d->lam1[i] - f21*d->lam2[i];
      atmp[2] += d->lam1[i]*d->p2r[i];
      atmp[3] += d->lam1[i]*d->p2i[i];
      }
    vhsum_cmplx_special (atmp[0], atmp[1], atmp[2], atmp[3], &alm[l]);
    l+=2;
    }
  }

NOINLINE static void calc_map2alm(sharp_job * restrict job,
  const sharp_Ylmgen_C *gen, s0data_v * restrict d, int nth)
  {
  int lmax=gen->lmax;
  int l=gen->m;
  int nv2 = (nth+VLEN-1)/VLEN;
  iter_to_ieee(gen, d, &l, nv2);
  job->opcnt += (l-gen->m) * 4*nth;
  if (l>lmax) return;
  job->opcnt += (lmax+1-l) * 8*nth;

  const sharp_ylmgen_dbl2 * restrict rf = gen->rf;
  dcmplx * restrict alm=job->almtmp;
  int full_ieee=1;
  for (int i=0; i<nv2; ++i)
    {
    getCorfac(d->scale[i], &d->corfac[i], gen->cf);
    full_ieee &= vallTrue(vge(d->scale[i],vload(sharp_minscale)));
    }

  while ((!full_ieee) && (l<=lmax))
    {
    full_ieee=1;
    Tv f10=vload(rf[l  ].f[0]), f11=vload(rf[l  ].f[1]),
       f20=vload(rf[l+1].f[0]), f21=vload(rf[l+1].f[1]);
    Tv atmp[4] = {vzero, vzero, vzero, vzero};
    for (int i=0; i<nv2; ++i)
      {
      d->lam1[i] = f10*d->cth[i]*d->lam2[i] - f11*d->lam1[i];
      atmp[0] += d->lam2[i]*d->corfac[i]*d->p1r[i];
      atmp[1] += d->lam2[i]*d->corfac[i]*d->p1i[i];
      d->lam2[i] = f20*d->cth[i]*d->lam1[i] - f21*d->lam2[i];
      if (rescale(&d->lam1[i], &d->lam2[i], &d->scale[i], vload(sharp_ftol)))
        getCorfac(d->scale[i], &d->corfac[i], gen->cf);
      full_ieee &= vallTrue(vge(d->scale[i],vload(sharp_minscale)));
      atmp[2] += d->lam1[i]*d->corfac[i]*d->p2r[i];
      atmp[3] += d->lam1[i]*d->corfac[i]*d->p2i[i];
      }
    vhsum_cmplx_special (atmp[0], atmp[1], atmp[2], atmp[3], &alm[l]);
    l+=2;
    }
  if (l>lmax) return;

  for (int i=0; i<nv2; ++i)
    {
    d->lam1[i] *= d->corfac[i];
    d->lam2[i] *= d->corfac[i];
    }
  map2alm_kernel(d, rf, alm, l, lmax, nv2);
  }

NOINLINE static void iter_to_ieee_spin (const sharp_Ylmgen_C * restrict gen,
  sxdata_v * restrict d, int * restrict l_, int nv2)
  {
  const sharp_ylmgen_dbl3 * restrict fx = gen->fx;
  Tv prefac=vload(gen->prefac[gen->m]),
     prescale=vload(gen->fscale[gen->m]);
  Tv limscale=vload(sharp_limscale);
  int below_limit=1;
  for (int i=0; i<nv2; ++i)
    {
    Tv cth2=vmax(vload(1e-15),vsqrt((vone+d->cth[i])*vload(0.5)));
    Tv sth2=vmax(vload(1e-15),vsqrt((vone-d->cth[i])*vload(0.5)));
    Tm mask=vlt(d->sth[i],vzero);
    vmuleq_mask(vand_mask(mask,vlt(d->cth[i],vzero)),cth2,vload(-1.));
    vmuleq_mask(vand_mask(mask,vgt(d->cth[i],vzero)),sth2,vload(-1.));

    Tv ccp, ccps, ssp, ssps, csp, csps, scp, scps;
    mypow(cth2,gen->cosPow,gen->powlimit,&ccp,&ccps);
    mypow(sth2,gen->sinPow,gen->powlimit,&ssp,&ssps);
    mypow(cth2,gen->sinPow,gen->powlimit,&csp,&csps);
    mypow(sth2,gen->cosPow,gen->powlimit,&scp,&scps);

    d->l1p[i] = vzero;
    d->l1m[i] = vzero;
    d->l2p[i] = prefac*ccp;
    d->scp[i] = prescale+ccps;
    d->l2m[i] = prefac*csp;
    d->scm[i] = prescale+csps;
    Tvnormalize(&d->l2m[i],&d->scm[i],sharp_fbighalf);
    Tvnormalize(&d->l2p[i],&d->scp[i],sharp_fbighalf);
    d->l2p[i] *= ssp;
    d->scp[i] += ssps;
    d->l2m[i] *= scp;
    d->scm[i] += scps;
    if (gen->preMinus_p)
      d->l2p[i] = vneg(d->l2p[i]);
    if (gen->preMinus_m)
      d->l2m[i] = vneg(d->l2m[i]);
    if (gen->s&1)
      d->l2p[i] = vneg(d->l2p[i]);

    Tvnormalize(&d->l2m[i],&d->scm[i],sharp_ftol);
    Tvnormalize(&d->l2p[i],&d->scp[i],sharp_ftol);

    below_limit &= vallTrue(vlt(d->scm[i],limscale)) &&
                   vallTrue(vlt(d->scp[i],limscale));
    }

  int l=gen->mhi;

  while (below_limit)
    {
    if (l+2>gen->lmax) {*l_=gen->lmax+1;return;}
    below_limit=1;
    Tv fx10=vload(fx[l+1].f[0]),fx11=vload(fx[l+1].f[1]),
       fx12=vload(fx[l+1].f[2]);
    Tv fx20=vload(fx[l+2].f[0]),fx21=vload(fx[l+2].f[1]),
       fx22=vload(fx[l+2].f[2]);
    for (int i=0; i<nv2; ++i)
      {
      d->l1p[i] = (d->cth[i]-fx11)*fx10*d->l2p[i] - fx12*d->l1p[i];
      d->l1m[i] = (d->cth[i]+fx11)*fx10*d->l2m[i] - fx12*d->l1m[i];
      d->l2p[i] = (d->cth[i]-fx21)*fx20*d->l1p[i] - fx22*d->l2p[i];
      d->l2m[i] = (d->cth[i]+fx21)*fx20*d->l1m[i] - fx22*d->l2m[i];
      if (rescale(&d->l1p[i],&d->l2p[i],&d->scp[i],vload(sharp_ftol)) ||
          rescale(&d->l1m[i],&d->l2m[i],&d->scm[i],vload(sharp_ftol)))
        below_limit &= vallTrue(vlt(d->scp[i],limscale)) &&
                       vallTrue(vlt(d->scm[i],limscale));
      }
    l+=2;
    }

  *l_=l;
  }

NOINLINE static void alm2map_spin_kernel(sxdata_v * restrict d,
  const sharp_ylmgen_dbl3 * restrict fx, const dcmplx * restrict alm,
  int l, int lmax, int nv2)
  {
  while (l<=lmax)
    {
    Tv fx10=vload(fx[l+1].f[0]),fx11=vload(fx[l+1].f[1]),
       fx12=vload(fx[l+1].f[2]);
    Tv fx20=vload(fx[l+2].f[0]),fx21=vload(fx[l+2].f[1]),
       fx22=vload(fx[l+2].f[2]);
    Tv agr1=vload(creal(alm[2*l  ])), agi1=vload(cimag(alm[2*l  ])),
       acr1=vload(creal(alm[2*l+1])), aci1=vload(cimag(alm[2*l+1]));
    Tv agr2=vload(creal(alm[2*l+2])), agi2=vload(cimag(alm[2*l+2])),
       acr2=vload(creal(alm[2*l+3])), aci2=vload(cimag(alm[2*l+3]));
    for (int i=0; i<nv2; ++i)
      {
      d->l1p[i] = (d->cth[i]-fx11)*fx10*d->l2p[i] - fx12*d->l1p[i];
      d->l1m[i] = (d->cth[i]+fx11)*fx10*d->l2m[i] - fx12*d->l1m[i];
      Tv lw1=d->l2p[i]+d->l2m[i];
      Tv lx2=d->l1m[i]-d->l1p[i];
      d->p1pr[i] += agr1*lw1 - aci2*lx2;
      d->p1pi[i] += agi1*lw1 + acr2*lx2;
      d->p1mr[i] += acr1*lw1 + agi2*lx2;
      d->p1mi[i] += aci1*lw1 - agr2*lx2;
      Tv lx1=d->l2m[i]-d->l2p[i];
      Tv lw2=d->l1p[i]+d->l1m[i];
      d->p2pr[i] += agr2*lw2 - aci1*lx1;
      d->p2pi[i] += agi2*lw2 + acr1*lx1;
      d->p2mr[i] += acr2*lw2 + agi1*lx1;
      d->p2mi[i] += aci2*lw2 - agr1*lx1;
      d->l2p[i] = (d->cth[i]-fx21)*fx20*d->l1p[i] - fx22*d->l2p[i];
      d->l2m[i] = (d->cth[i]+fx21)*fx20*d->l1m[i] - fx22*d->l2m[i];
      }
    l+=2;
    }
  }

NOINLINE static void calc_alm2map_spin (sharp_job * restrict job,
  const sharp_Ylmgen_C * restrict gen, sxdata_v * restrict d, int nth)
  {
  int l,lmax=gen->lmax;
  int nv2 = (nth+VLEN-1)/VLEN;
  iter_to_ieee_spin(gen, d, &l, nv2);
  job->opcnt += (l-gen->mhi) * 10*nth;
  if (l>lmax) return;
  job->opcnt += (lmax+1-l) * 28*nth;

  const sharp_ylmgen_dbl3 * restrict fx = gen->fx;
  const dcmplx * restrict alm=job->almtmp;
  int full_ieee=1;
  for (int i=0; i<nv2; ++i)
    {
    getCorfac(d->scp[i], &d->cfp[i], gen->cf);
    getCorfac(d->scm[i], &d->cfm[i], gen->cf);
    full_ieee &= vallTrue(vge(d->scp[i],vload(sharp_minscale))) &&
                 vallTrue(vge(d->scm[i],vload(sharp_minscale)));
    }

  while((!full_ieee) && (l<=lmax))
    {
    Tv fx10=vload(fx[l+1].f[0]),fx11=vload(fx[l+1].f[1]),
       fx12=vload(fx[l+1].f[2]);
    Tv fx20=vload(fx[l+2].f[0]),fx21=vload(fx[l+2].f[1]),
       fx22=vload(fx[l+2].f[2]);
    Tv agr1=vload(creal(alm[2*l  ])), agi1=vload(cimag(alm[2*l  ])),
       acr1=vload(creal(alm[2*l+1])), aci1=vload(cimag(alm[2*l+1]));
    Tv agr2=vload(creal(alm[2*l+2])), agi2=vload(cimag(alm[2*l+2])),
       acr2=vload(creal(alm[2*l+3])), aci2=vload(cimag(alm[2*l+3]));
    full_ieee=1;
    for (int i=0; i<nv2; ++i)
      {
      d->l1p[i] = (d->cth[i]-fx11)*fx10*d->l2p[i] - fx12*d->l1p[i];
      d->l1m[i] = (d->cth[i]+fx11)*fx10*d->l2m[i] - fx12*d->l1m[i];
      Tv lw1=d->l2p[i]*d->cfp[i] + d->l2m[i]*d->cfm[i];
      Tv lx2=d->l1m[i]*d->cfm[i] - d->l1p[i]*d->cfp[i];
      d->p1pr[i] += agr1*lw1 - aci2*lx2;
      d->p1pi[i] += agi1*lw1 + acr2*lx2;
      d->p1mr[i] += acr1*lw1 + agi2*lx2;
      d->p1mi[i] += aci1*lw1 - agr2*lx2;
      Tv lx1=d->l2m[i]*d->cfm[i] - d->l2p[i]*d->cfp[i];
      Tv lw2=d->l1p[i]*d->cfp[i] + d->l1m[i]*d->cfm[i];
      d->p2pr[i] += agr2*lw2 - aci1*lx1;
      d->p2pi[i] += agi2*lw2 + acr1*lx1;
      d->p2mr[i] += acr2*lw2 + agi1*lx1;
      d->p2mi[i] += aci2*lw2 - agr1*lx1;
      d->l2p[i] = (d->cth[i]-fx21)*fx20*d->l1p[i] - fx22*d->l2p[i];
      d->l2m[i] = (d->cth[i]+fx21)*fx20*d->l1m[i] - fx22*d->l2m[i];
      if (rescale(&d->l1p[i], &d->l2p[i], &d->scp[i], vload(sharp_ftol)))
        getCorfac(d->scp[i], &d->cfp[i], gen->cf);
      full_ieee &= vallTrue(vge(d->scp[i],vload(sharp_minscale)));
      if (rescale(&d->l1m[i], &d->l2m[i], &d->scm[i], vload(sharp_ftol)))
        getCorfac(d->scm[i], &d->cfm[i], gen->cf);
      full_ieee &= vallTrue(vge(d->scm[i],vload(sharp_minscale)));
      }
    l+=2;
    }
  if (l>lmax) return;

  for (int i=0; i<nv2; ++i)
    {
    d->l1p[i] *= d->cfp[i];
    d->l2p[i] *= d->cfp[i];
    d->l1m[i] *= d->cfm[i];
    d->l2m[i] *= d->cfm[i];
    }
  alm2map_spin_kernel(d, fx, alm, l, lmax, nv2);
  }

NOINLINE static void map2alm_spin_kernel(sxdata_v * restrict d,
  const sharp_ylmgen_dbl3 * restrict fx, dcmplx * restrict alm,
  int l, int lmax, int nv2)
  {
  while (l<=lmax)
    {
    Tv fx10=vload(fx[l+1].f[0]),fx11=vload(fx[l+1].f[1]),
       fx12=vload(fx[l+1].f[2]);
    Tv fx20=vload(fx[l+2].f[0]),fx21=vload(fx[l+2].f[1]),
       fx22=vload(fx[l+2].f[2]);
    Tv agr1=vzero, agi1=vzero, acr1=vzero, aci1=vzero;
    Tv agr2=vzero, agi2=vzero, acr2=vzero, aci2=vzero;
    for (int i=0; i<nv2; ++i)
      {
      d->l1p[i] = (d->cth[i]-fx11)*fx10*d->l2p[i] - fx12*d->l1p[i];
      d->l1m[i] = (d->cth[i]+fx11)*fx10*d->l2m[i] - fx12*d->l1m[i];
      Tv lw = d->l2p[i] + d->l2m[i];
      Tv lx = d->l2m[i] - d->l2p[i];
      agr1 += d->p1pr[i]*lw - d->p2mi[i]*lx;;
      agi1 += d->p1pi[i]*lw + d->p2mr[i]*lx;
      acr1 += d->p1mr[i]*lw + d->p2pi[i]*lx;
      aci1 += d->p1mi[i]*lw - d->p2pr[i]*lx;
      d->l2p[i] = (d->cth[i]-fx21)*fx20*d->l1p[i] - fx22*d->l2p[i];
      d->l2m[i] = (d->cth[i]+fx21)*fx20*d->l1m[i] - fx22*d->l2m[i];
      lw = d->l1p[i] + d->l1m[i];
      lx = d->l1m[i] - d->l1p[i];
      agr2 += d->p2pr[i]*lw - d->p1mi[i]*lx;
      agi2 += d->p2pi[i]*lw + d->p1mr[i]*lx;
      acr2 += d->p2mr[i]*lw + d->p1pi[i]*lx;
      aci2 += d->p2mi[i]*lw - d->p1pr[i]*lx;
      }
    vhsum_cmplx_special (agr1,agi1,acr1,aci1,&alm[2*l]);
    vhsum_cmplx_special (agr2,agi2,acr2,aci2,&alm[2*l+2]);
    l+=2;
    }
  }

NOINLINE static void calc_map2alm_spin (sharp_job * restrict job,
  const sharp_Ylmgen_C * restrict gen, sxdata_v * restrict d, int nth)
  {
  int l,lmax=gen->lmax;
  int nv2 = (nth+VLEN-1)/VLEN;
  iter_to_ieee_spin(gen, d, &l, nv2);
  job->opcnt += (l-gen->mhi) * 10*nth;
  if (l>lmax) return;
  job->opcnt += (lmax+1-l) * 28*nth;

  const sharp_ylmgen_dbl3 * restrict fx = gen->fx;
  dcmplx * restrict alm=job->almtmp;
  int full_ieee=1;
  for (int i=0; i<nv2; ++i)
    {
    getCorfac(d->scp[i], &d->cfp[i], gen->cf);
    getCorfac(d->scm[i], &d->cfm[i], gen->cf);
    full_ieee &= vallTrue(vge(d->scp[i],vload(sharp_minscale))) &&
                 vallTrue(vge(d->scm[i],vload(sharp_minscale)));
    }

  while((!full_ieee) && (l<=lmax))
    {
    Tv fx10=vload(fx[l+1].f[0]),fx11=vload(fx[l+1].f[1]),
       fx12=vload(fx[l+1].f[2]);
    Tv fx20=vload(fx[l+2].f[0]),fx21=vload(fx[l+2].f[1]),
       fx22=vload(fx[l+2].f[2]);
    Tv agr1=vzero, agi1=vzero, acr1=vzero, aci1=vzero;
    Tv agr2=vzero, agi2=vzero, acr2=vzero, aci2=vzero;
    full_ieee=1;
    for (int i=0; i<nv2; ++i)
      {
      d->l1p[i] = (d->cth[i]-fx11)*fx10*d->l2p[i] - fx12*d->l1p[i];
      d->l1m[i] = (d->cth[i]+fx11)*fx10*d->l2m[i] - fx12*d->l1m[i];
      Tv lw = d->l2p[i]*d->cfp[i] + d->l2m[i]*d->cfm[i];
      Tv lx = d->l2m[i]*d->cfm[i] - d->l2p[i]*d->cfp[i];
      agr1 += d->p1pr[i]*lw - d->p2mi[i]*lx;
      agi1 += d->p1pi[i]*lw + d->p2mr[i]*lx;
      acr1 += d->p1mr[i]*lw + d->p2pi[i]*lx;
      aci1 += d->p1mi[i]*lw - d->p2pr[i]*lx;
      d->l2p[i] = (d->cth[i]-fx21)*fx20*d->l1p[i] - fx22*d->l2p[i];
      d->l2m[i] = (d->cth[i]+fx21)*fx20*d->l1m[i] - fx22*d->l2m[i];
      lw = d->l1p[i]*d->cfp[i] + d->l1m[i]*d->cfm[i];
      lx = d->l1m[i]*d->cfm[i] - d->l1p[i]*d->cfp[i];
      agr2 += d->p2pr[i]*lw - d->p1mi[i]*lx;
      agi2 += d->p2pi[i]*lw + d->p1mr[i]*lx;
      acr2 += d->p2mr[i]*lw + d->p1pi[i]*lx;
      aci2 += d->p2mi[i]*lw - d->p1pr[i]*lx;
      if (rescale(&d->l1p[i], &d->l2p[i], &d->scp[i], vload(sharp_ftol)))
        getCorfac(d->scp[i], &d->cfp[i], gen->cf);
      full_ieee &= vallTrue(vge(d->scp[i],vload(sharp_minscale)));
      if (rescale(&d->l1m[i], &d->l2m[i], &d->scm[i], vload(sharp_ftol)))
        getCorfac(d->scm[i], &d->cfm[i], gen->cf);
      full_ieee &= vallTrue(vge(d->scm[i],vload(sharp_minscale)));
      }
    vhsum_cmplx_special (agr1,agi1,acr1,aci1,&alm[2*l]);
    vhsum_cmplx_special (agr2,agi2,acr2,aci2,&alm[2*l+2]);
    l+=2;
    }
  if (l>lmax) return;

  for (int i=0; i<nv2; ++i)
    {
    d->l1p[i] *= d->cfp[i];
    d->l2p[i] *= d->cfp[i];
    d->l1m[i] *= d->cfm[i];
    d->l2m[i] *= d->cfm[i];
    }
  map2alm_spin_kernel(d, fx, alm, l, lmax, nv2);
  }


NOINLINE static void alm2map_deriv1_kernel(sxdata_v * restrict d,
  const sharp_ylmgen_dbl3 * restrict fx, const dcmplx * restrict alm,
  int l, int lmax, int nv2)
  {
  while (l<=lmax)
    {
    Tv fx10=vload(fx[l+1].f[0]),fx11=vload(fx[l+1].f[1]),
       fx12=vload(fx[l+1].f[2]);
    Tv fx20=vload(fx[l+2].f[0]),fx21=vload(fx[l+2].f[1]),
       fx22=vload(fx[l+2].f[2]);
    Tv ar1=vload(creal(alm[l  ])), ai1=vload(cimag(alm[l  ])),
       ar2=vload(creal(alm[l+1])), ai2=vload(cimag(alm[l+1]));
    for (int i=0; i<nv2; ++i)
      {
      d->l1p[i] = (d->cth[i]-fx11)*fx10*d->l2p[i] - fx12*d->l1p[i];
      d->l1m[i] = (d->cth[i]+fx11)*fx10*d->l2m[i] - fx12*d->l1m[i];
      Tv lw=d->l2p[i]+d->l2m[i];
      d->p1pr[i] += ar1*lw;
      d->p1pi[i] += ai1*lw;
      Tv lx=d->l2m[i]-d->l2p[i];
      d->p2mr[i] += ai1*lx;
      d->p2mi[i] -= ar1*lx;
      lw=d->l1p[i]+d->l1m[i];
      d->p2pr[i] += ar2*lw;
      d->p2pi[i] += ai2*lw;
      lx=d->l1m[i]-d->l1p[i];
      d->p1mr[i] += ai2*lx;
      d->p1mi[i] -= ar2*lx;
      d->l2p[i] = (d->cth[i]-fx21)*fx20*d->l1p[i] - fx22*d->l2p[i];
      d->l2m[i] = (d->cth[i]+fx21)*fx20*d->l1m[i] - fx22*d->l2m[i];
      }
    l+=2;
    }
  }

NOINLINE static void calc_alm2map_deriv1(sharp_job * restrict job,
  const sharp_Ylmgen_C * restrict gen, sxdata_v * restrict d, int nth)
  {
  int l,lmax=gen->lmax;
  int nv2 = (nth+VLEN-1)/VLEN;
  iter_to_ieee_spin(gen, d, &l, nv2);
  job->opcnt += (l-gen->mhi) * 10*nth;
  if (l>lmax) return;
  job->opcnt += (lmax+1-l) * 20*nth;

  const sharp_ylmgen_dbl3 * restrict fx = gen->fx;
  const dcmplx * restrict alm=job->almtmp;
  int full_ieee=1;
  for (int i=0; i<nv2; ++i)
    {
    getCorfac(d->scp[i], &d->cfp[i], gen->cf);
    getCorfac(d->scm[i], &d->cfm[i], gen->cf);
    full_ieee &= vallTrue(vge(d->scp[i],vload(sharp_minscale))) &&
                 vallTrue(vge(d->scm[i],vload(sharp_minscale)));
    }

  while((!full_ieee) && (l<=lmax))
    {
    Tv fx10=vload(fx[l+1].f[0]),fx11=vload(fx[l+1].f[1]),
       fx12=vload(fx[l+1].f[2]);
    Tv fx20=vload(fx[l+2].f[0]),fx21=vload(fx[l+2].f[1]),
       fx22=vload(fx[l+2].f[2]);
    Tv ar1=vload(creal(alm[l  ])), ai1=vload(cimag(alm[l  ])),
       ar2=vload(creal(alm[l+1])), ai2=vload(cimag(alm[l+1]));
    full_ieee=1;
    for (int i=0; i<nv2; ++i)
      {
      d->l1p[i] = (d->cth[i]-fx11)*fx10*d->l2p[i] - fx12*d->l1p[i];
      d->l1m[i] = (d->cth[i]+fx11)*fx10*d->l2m[i] - fx12*d->l1m[i];
      Tv lw=d->l2p[i]*d->cfp[i]+d->l2m[i]*d->cfm[i];
      d->p1pr[i] += ar1*lw;
      d->p1pi[i] += ai1*lw;
      Tv lx=d->l2m[i]*d->cfm[i]-d->l2p[i]*d->cfp[i];
      d->p2mr[i] += ai1*lx;
      d->p2mi[i] -= ar1*lx;
      lw=d->l1p[i]*d->cfp[i]+d->l1m[i]*d->cfm[i];
      d->p2pr[i] += ar2*lw;
      d->p2pi[i] += ai2*lw;
      lx=d->l1m[i]*d->cfm[i]-d->l1p[i]*d->cfp[i];
      d->p1mr[i] += ai2*lx;
      d->p1mi[i] -= ar2*lx;
      d->l2p[i] = (d->cth[i]-fx21)*fx20*d->l1p[i] - fx22*d->l2p[i];
      d->l2m[i] = (d->cth[i]+fx21)*fx20*d->l1m[i] - fx22*d->l2m[i];
      if (rescale(&d->l1p[i], &d->l2p[i], &d->scp[i], vload(sharp_ftol)))
        {
        getCorfac(d->scp[i], &d->cfp[i], gen->cf);
        full_ieee &= vallTrue(vge(d->scp[i],vload(sharp_minscale)));
        }
      if (rescale(&d->l1m[i], &d->l2m[i], &d->scm[i], vload(sharp_ftol)))
        {
        getCorfac(d->scm[i], &d->cfm[i], gen->cf);
        full_ieee &= vallTrue(vge(d->scm[i],vload(sharp_minscale)));
        }
      }
    l+=2;
    }
  if (l>lmax) return;

  for (int i=0; i<nv2; ++i)
    {
    d->l1p[i] *= d->cfp[i];
    d->l2p[i] *= d->cfp[i];
    d->l1m[i] *= d->cfm[i];
    d->l2m[i] *= d->cfm[i];
    }
  alm2map_deriv1_kernel(d, fx, alm, l, lmax, nv2);
  }


#define VZERO(var) do { memset(&(var),0,sizeof(var)); } while(0)

NOINLINE static void inner_loop_a2m(sharp_job *job, const int *ispair,
  const double *cth_, const double *sth_, int llim, int ulim,
  sharp_Ylmgen_C *gen, int mi, const int *mlim)
  {
  const int m = job->ainfo->mval[mi];
  sharp_Ylmgen_prepare (gen, m);

  switch (job->type)
    {
    case SHARP_ALM2MAP:
    case SHARP_ALM2MAP_DERIV1:
      {
      if (job->spin==0)
        {
        const int nval=nv0*VLEN;
        int ith=0;
        int itgt[nval];
        while (ith<ulim-llim)
          {
          s0data_u d;
          VZERO(d.s.p1r); VZERO(d.s.p1i); VZERO(d.s.p2r); VZERO(d.s.p2i);
          int nth=0;
          while ((nth<nval)&&(ith<ulim-llim))
            {
            if (mlim[ith]>=m)
              {
              itgt[nth] = ith;
              d.s.cth[nth]=cth_[ith]; d.s.sth[nth]=sth_[ith];
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
              d.s.cth[i]=d.s.cth[nth-1];
              d.s.sth[i]=d.s.sth[nth-1];
              d.s.p1r[i]=d.s.p1i[i]=d.s.p2r[i]=d.s.p2i[i]=0.;
              }
            calc_alm2map (job, gen, &d.v, nth);
            for (int i=0; i<nth; ++i)
              {
              int tgt=itgt[i];
              int phas_idx = tgt*job->s_th + mi*job->s_m;
              complex double r1 = d.s.p1r[i] + d.s.p1i[i]*_Complex_I,
                             r2 = d.s.p2r[i] + d.s.p2i[i]*_Complex_I;
              job->phase[phas_idx] = r1+r2;
              if (ispair[tgt])
                job->phase[phas_idx+1] = r1-r2;
              }
            }
          }
        }
      else
        {
        const int nval=nvx*VLEN;
        int ith=0;
        int itgt[nval];
        while (ith<ulim-llim)
          {
          sxdata_u d;
          VZERO(d.s.p1pr); VZERO(d.s.p1pi); VZERO(d.s.p2pr); VZERO(d.s.p2pi);
          VZERO(d.s.p1mr); VZERO(d.s.p1mi); VZERO(d.s.p2mr); VZERO(d.s.p2mi);
          int nth=0;
          while ((nth<nval)&&(ith<ulim-llim))
            {
            if (mlim[ith]>=m)
              {
              itgt[nth] = ith;
              d.s.cth[nth]=cth_[ith]; d.s.sth[nth]=sth_[ith];
              ++nth;
              }
            else
              {
              int phas_idx = ith*job->s_th + mi*job->s_m;
              job->phase[phas_idx  ] = job->phase[phas_idx+1] = 0;
              job->phase[phas_idx+2] = job->phase[phas_idx+3] = 0;
              }
            ++ith;
            }
          if (nth>0)
            {
            int i2=((nth+VLEN-1)/VLEN)*VLEN;
            for (int i=nth; i<i2; ++i)
              {
              d.s.cth[i]=d.s.cth[nth-1];
              d.s.sth[i]=d.s.sth[nth-1];
              d.s.p1pr[i]=d.s.p1pi[i]=d.s.p2pr[i]=d.s.p2pi[i]=0.;
              d.s.p1mr[i]=d.s.p1mi[i]=d.s.p2mr[i]=d.s.p2mi[i]=0.;
              }
            (job->type==SHARP_ALM2MAP) ?
              calc_alm2map_spin  (job, gen, &d.v, nth) :
              calc_alm2map_deriv1(job, gen, &d.v, nth);
            for (int i=0; i<nth; ++i)
              {
              int tgt=itgt[i];
              int phas_idx = tgt*job->s_th + mi*job->s_m;
              complex double q1 = d.s.p1pr[i] + d.s.p1pi[i]*_Complex_I,
                             q2 = d.s.p2pr[i] + d.s.p2pi[i]*_Complex_I,
                             u1 = d.s.p1mr[i] + d.s.p1mi[i]*_Complex_I,
                             u2 = d.s.p2mr[i] + d.s.p2mi[i]*_Complex_I;
              job->phase[phas_idx  ] = q1+q2;
              job->phase[phas_idx+2] = u1+u2;
              if (ispair[tgt])
                {
                dcmplx *phQ = &(job->phase[phas_idx+1]),
                       *phU = &(job->phase[phas_idx+3]);
                *phQ = q1-q2;
                *phU = u1-u2;
                if ((gen->mhi-gen->m+gen->s)&1)
                  { *phQ=-(*phQ); *phU=-(*phU); }
                }
              }
            }
          }
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
  const int m = job->ainfo->mval[mi];
  sharp_Ylmgen_prepare (gen, m);

  switch (job->type)
    {
    case SHARP_MAP2ALM:
      {
      if (job->spin==0)
        {
        const int nval=nv0*VLEN;
        int ith=0;
        while (ith<ulim-llim)
          {
          s0data_u d;
          int nth=0;
          while ((nth<nval)&&(ith<ulim-llim))
            {
            if (mlim[ith]>=m)
              {
              d.s.cth[nth]=cth_[ith]; d.s.sth[nth]=sth_[ith];
              int phas_idx = ith*job->s_th + mi*job->s_m;
              dcmplx ph1=job->phase[phas_idx];
              dcmplx ph2=ispair[ith] ? job->phase[phas_idx+1] : 0.;
              d.s.p1r[nth]=creal(ph1+ph2); d.s.p1i[nth]=cimag(ph1+ph2);
              d.s.p2r[nth]=creal(ph1-ph2); d.s.p2i[nth]=cimag(ph1-ph2);
              ++nth;
              }
            ++ith;
            }
          if (nth>0)
            {
            int i2=((nth+VLEN-1)/VLEN)*VLEN;
            for (int i=nth; i<i2; ++i)
              {
              d.s.cth[i]=d.s.cth[nth-1];
              d.s.sth[i]=d.s.sth[nth-1];
              d.s.p1r[i]=d.s.p1i[i]=d.s.p2r[i]=d.s.p2i[i]=0.;
              }
            calc_map2alm (job, gen, &d.v, nth);
            }
          }
        }
      else
        {
        const int nval=nvx*VLEN;
        int ith=0;
        while (ith<ulim-llim)
          {
          sxdata_u d;
          int nth=0;
          while ((nth<nval)&&(ith<ulim-llim))
            {
            if (mlim[ith]>=m)
              {
              d.s.cth[nth]=cth_[ith]; d.s.sth[nth]=sth_[ith];
              int phas_idx = ith*job->s_th + mi*job->s_m;
              dcmplx p1Q=job->phase[phas_idx],
                     p1U=job->phase[phas_idx+2],
                     p2Q=ispair[ith] ? job->phase[phas_idx+1]:0.,
                     p2U=ispair[ith] ? job->phase[phas_idx+3]:0.;
              if ((gen->mhi-gen->m+gen->s)&1)
                { p2Q=-p2Q; p2U=-p2U; }
              d.s.p1pr[nth]=creal(p1Q+p2Q); d.s.p1pi[nth]=cimag(p1Q+p2Q);
              d.s.p1mr[nth]=creal(p1U+p2U); d.s.p1mi[nth]=cimag(p1U+p2U);
              d.s.p2pr[nth]=creal(p1Q-p2Q); d.s.p2pi[nth]=cimag(p1Q-p2Q);
              d.s.p2mr[nth]=creal(p1U-p2U); d.s.p2mi[nth]=cimag(p1U-p2U);
              ++nth;
              }
            ++ith;
            }
          if (nth>0)
            {
            int i2=((nth+VLEN-1)/VLEN)*VLEN;
            for (int i=nth; i<i2; ++i)
              {
              d.s.cth[i]=d.s.cth[nth-1];
              d.s.sth[i]=d.s.sth[nth-1];
              d.s.p1pr[i]=d.s.p1pi[i]=d.s.p2pr[i]=d.s.p2pi[i]=0.;
              d.s.p1mr[i]=d.s.p1mi[i]=d.s.p2mr[i]=d.s.p2mi[i]=0.;
              }
            calc_map2alm_spin(job, gen, &d.v, nth);
            }
          }
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
  return nv0;
  }
