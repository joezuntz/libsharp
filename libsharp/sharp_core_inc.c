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

/*! \file sharp_core_inc.c
 *  Type-dependent code for the computational core
 *
 *  Copyright (C) 2012-2017 Max-Planck-Society
 *  \author Martin Reinecke
 */

typedef struct
  { Tv v[nvec]; } Tb;

typedef union
  { Tb b; double s[VLEN*nvec]; } Y(Tbu);

typedef struct
  { Tb r, i; } Y(Tbri);

typedef struct
  { Tb qr, qi, ur, ui; } Y(Tbqu);

typedef struct
  { double r[VLEN*nvec], i[VLEN*nvec]; } Y(Tsri);

typedef struct
  { double qr[VLEN*nvec],qi[VLEN*nvec],ur[VLEN*nvec],ui[VLEN*nvec]; } Y(Tsqu);

typedef union
  { Y(Tbri) b; Y(Tsri)s; } Y(Tburi);

typedef union
  { Y(Tbqu) b; Y(Tsqu)s; } Y(Tbuqu);

static inline Tb Y(Tbconst)(double val)
  {
  Tv v=vload(val);
  Tb res;
  for (int i=0; i<nvec; ++i) res.v[i]=v;
  return res;
  }

static inline void Y(Tbmuleq1)(Tb * restrict a, double b)
  { Tv v=vload(b); for (int i=0; i<nvec; ++i) vmuleq(a->v[i],v); }

static inline Tb Y(Tbprod)(Tb a, Tb b)
  { Tb r; for (int i=0; i<nvec; ++i) r.v[i]=vmul(a.v[i],b.v[i]); return r; }

static inline void Y(Tbmuleq)(Tb * restrict a, Tb b)
  { for (int i=0; i<nvec; ++i) vmuleq(a->v[i],b.v[i]); }

static void Y(Tbnormalize) (Tb * restrict val, Tb * restrict scale,
  double maxval)
  {
  const Tv vfmin=vload(sharp_fsmall*maxval), vfmax=vload(maxval);
  const Tv vfsmall=vload(sharp_fsmall), vfbig=vload(sharp_fbig);
  for (int i=0;i<nvec; ++i)
    {
    Tm mask = vgt(vabs(val->v[i]),vfmax);
    while (vanyTrue(mask))
      {
      vmuleq_mask(mask,val->v[i],vfsmall);
      vaddeq_mask(mask,scale->v[i],vone);
      mask = vgt(vabs(val->v[i]),vfmax);
      }
    mask = vand_mask(vlt(vabs(val->v[i]),vfmin),vne(val->v[i],vzero));
    while (vanyTrue(mask))
      {
      vmuleq_mask(mask,val->v[i],vfbig);
      vsubeq_mask(mask,scale->v[i],vone);
      mask = vand_mask(vlt(vabs(val->v[i]),vfmin),vne(val->v[i],vzero));
      }
    }
  }

NOINLINE static void Y(mypow) (Tb val, int npow, const double * restrict powlimit,
  Tb * restrict resd, Tb * restrict ress)
  {
  Tv vminv=vload(powlimit[npow]);
  Tm mask = vlt(vabs(val.v[0]),vminv);
  for (int i=1;i<nvec; ++i)
    mask=vor_mask(mask,vlt(vabs(val.v[i]),vminv));
  if (!vanyTrue(mask)) // no underflows possible, use quick algoritm
    {
    Tb res=Y(Tbconst)(1.);
    do
      {
      if (npow&1)
        for (int i=0; i<nvec; ++i)
          {
          vmuleq(res.v[i],val.v[i]);
          vmuleq(val.v[i],val.v[i]);
          }
      else
        for (int i=0; i<nvec; ++i)
          vmuleq(val.v[i],val.v[i]);
      }
    while(npow>>=1);
    *resd=res;
    *ress=Y(Tbconst)(0.);
    }
  else
    {
    Tb scale=Y(Tbconst)(0.), scaleint=Y(Tbconst)(0.), res=Y(Tbconst)(1.);
    Y(Tbnormalize)(&val,&scaleint,sharp_fbighalf);
    do
      {
      if (npow&1)
        {
        for (int i=0; i<nvec; ++i)
          {
          vmuleq(res.v[i],val.v[i]);
          vaddeq(scale.v[i],scaleint.v[i]);
          }
        Y(Tbnormalize)(&res,&scale,sharp_fbighalf);
        }
      for (int i=0; i<nvec; ++i)
        {
        vmuleq(val.v[i],val.v[i]);
        vaddeq(scaleint.v[i],scaleint.v[i]);
        }
      Y(Tbnormalize)(&val,&scaleint,sharp_fbighalf);
      }
    while(npow>>=1);
    *resd=res;
    *ress=scale;
    }
  }

static inline int Y(rescale) (Tb * restrict lam1, Tb * restrict lam2,
  Tb * restrict scale)
  {
  int did_scale=0;
  for (int i=0;i<nvec; ++i)
    {
    Tm mask = vgt(vabs(lam2->v[i]),vload(sharp_ftol));
    if (vanyTrue(mask))
      {
      did_scale=1;
      vmuleq_mask(mask,lam1->v[i],vload(sharp_fsmall));
      vmuleq_mask(mask,lam2->v[i],vload(sharp_fsmall));
      vaddeq_mask(mask,scale->v[i],vone);
      }
    }
  return did_scale;
  }

static inline int Y(TballLt)(Tb a,double b)
  {
  Tv vb=vload(b);
  Tm res=vlt(a.v[0],vb);
  for (int i=1; i<nvec; ++i)
    res=vand_mask(res,vlt(a.v[i],vb));
  return vallTrue(res);
  }
static inline int Y(TballGt)(Tb a,double b)
  {
  Tv vb=vload(b);
  Tm res=vgt(a.v[0],vb);
  for (int i=1; i<nvec; ++i)
    res=vand_mask(res,vgt(a.v[i],vb));
  return vallTrue(res);
  }
static inline int Y(TballGe)(Tb a,double b)
  {
  Tv vb=vload(b);
  Tm res=vge(a.v[0],vb);
  for (int i=1; i<nvec; ++i)
    res=vand_mask(res,vge(a.v[i],vb));
  return vallTrue(res);
  }

static void Y(getCorfac)(Tb scale, Tb * restrict corfac,
  const double * restrict cf)
  {
  Y(Tbu) sc, corf;
  sc.b=scale;
  for (int i=0; i<VLEN*nvec; ++i)
    corf.s[i] = (sc.s[i]<sharp_minscale) ?
      0. : cf[(int)(sc.s[i])-sharp_minscale];
  *corfac=corf.b;
  }

NOINLINE static void Y(iter_to_ieee) (const Tb sth, Tb cth, int *l_,
  Tb * restrict lam_1_, Tb * restrict lam_2_, Tb * restrict scale_,
  const sharp_Ylmgen_C * restrict gen)
  {
  int l=gen->m;
  Tb lam_1=Y(Tbconst)(0.), lam_2, scale;
  Y(mypow) (sth,l,gen->powlimit,&lam_2,&scale);
  Y(Tbmuleq1) (&lam_2,(gen->m&1) ? -gen->mfac[gen->m]:gen->mfac[gen->m]);
  Y(Tbnormalize)(&lam_2,&scale,sharp_ftol);

  int below_limit = Y(TballLt)(scale,sharp_limscale);
  while (below_limit)
    {
    if (l+2>gen->lmax) {*l_=gen->lmax+1;return;}
    for (int i=0; i<nvec; ++i)
      {
      lam_1.v[i] = vload(gen->rf[l].f[0])*(cth.v[i]*lam_2.v[i])
                 - vload(gen->rf[l].f[1])*lam_1.v[i];
      lam_2.v[i] = vload(gen->rf[l+1].f[0])*(cth.v[i]*lam_1.v[i])
                 - vload(gen->rf[l+1].f[1])*lam_2.v[i];
      }
    if (Y(rescale)(&lam_1,&lam_2,&scale))
      below_limit = Y(TballLt)(scale,sharp_limscale);
    l+=2;
    }
  *l_=l; *lam_1_=lam_1; *lam_2_=lam_2; *scale_=scale;
  }

static inline void Y(rec_step) (Tb * restrict rxp, Tb * restrict rxm,
  Tb * restrict ryp, Tb * restrict rym, const Tb cth,
  const sharp_ylmgen_dbl3 fx)
  {
  Tv fx0=vload(fx.f[0]),fx1=vload(fx.f[1]),fx2=vload(fx.f[2]);
  for (int i=0; i<nvec; ++i)
    {
    rxp->v[i] = (cth.v[i]-fx1)*fx0*ryp->v[i] - fx2*rxp->v[i];
    rxm->v[i] = (cth.v[i]+fx1)*fx0*rym->v[i] - fx2*rxm->v[i];
    }
  }

static void Y(iter_to_ieee_spin) (const Tb cth, const Tb sth, int *l_,
  Tb * rec1p_, Tb * rec1m_, Tb * rec2p_, Tb * rec2m_,
  Tb * scalep_, Tb * scalem_, const sharp_Ylmgen_C * restrict gen)
  {
  const sharp_ylmgen_dbl3 * restrict fx = gen->fx;
  Tb cth2, sth2;
  for (int i=0; i<nvec; ++i)
    {
    cth2.v[i]=vsqrt(vmul(vadd(vone,cth.v[i]),vload(0.5)));
    cth2.v[i]=vmax(cth2.v[i],vload(1e-15));
    sth2.v[i]=vsqrt(vmul(vsub(vone,cth.v[i]),vload(0.5)));
    sth2.v[i]=vmax(sth2.v[i],vload(1e-15));
    Tm mask=vlt(sth.v[i],vzero);
    Tm cmask=vand_mask(mask,vlt(cth.v[i],vzero));
    vmuleq_mask(cmask,cth2.v[i],vload(-1.));
    Tm smask=vand_mask(mask,vgt(cth.v[i],vzero));
    vmuleq_mask(smask,sth2.v[i],vload(-1.));
    }

  Tb ccp, ccps, ssp, ssps, csp, csps, scp, scps;
  Y(mypow)(cth2,gen->cosPow,gen->powlimit,&ccp,&ccps);
  Y(mypow)(sth2,gen->sinPow,gen->powlimit,&ssp,&ssps);
  Y(mypow)(cth2,gen->sinPow,gen->powlimit,&csp,&csps);
  Y(mypow)(sth2,gen->cosPow,gen->powlimit,&scp,&scps);

  Tb rec2p, rec2m, scalep, scalem;
  Tb rec1p=Y(Tbconst)(0.), rec1m=Y(Tbconst)(0.);
  Tv prefac=vload(gen->prefac[gen->m]),
     prescale=vload(gen->fscale[gen->m]);
  for (int i=0; i<nvec; ++i)
    {
    rec2p.v[i]=vmul(prefac,ccp.v[i]);
    scalep.v[i]=vadd(prescale,ccps.v[i]);
    rec2m.v[i]=vmul(prefac,csp.v[i]);
    scalem.v[i]=vadd(prescale,csps.v[i]);
    }
  Y(Tbnormalize)(&rec2m,&scalem,sharp_fbighalf);
  Y(Tbnormalize)(&rec2p,&scalep,sharp_fbighalf);
  for (int i=0; i<nvec; ++i)
    {
    rec2p.v[i]=vmul(rec2p.v[i],ssp.v[i]);
    scalep.v[i]=vadd(scalep.v[i],ssps.v[i]);
    rec2m.v[i]=vmul(rec2m.v[i],scp.v[i]);
    scalem.v[i]=vadd(scalem.v[i],scps.v[i]);
    if (gen->preMinus_p)
      rec2p.v[i]=vneg(rec2p.v[i]);
    if (gen->preMinus_m)
      rec2m.v[i]=vneg(rec2m.v[i]);
    if (gen->s&1)
      rec2p.v[i]=vneg(rec2p.v[i]);
    }
  Y(Tbnormalize)(&rec2m,&scalem,sharp_ftol);
  Y(Tbnormalize)(&rec2p,&scalep,sharp_ftol);

  int l=gen->mhi;

  int below_limit = Y(TballLt)(scalep,sharp_limscale)
                 && Y(TballLt)(scalem,sharp_limscale);
  while (below_limit)
    {
    if (l+2>gen->lmax) {*l_=gen->lmax+1;return;}
    Y(rec_step)(&rec1p,&rec1m,&rec2p,&rec2m,cth,fx[l+1]);
    Y(rec_step)(&rec2p,&rec2m,&rec1p,&rec1m,cth,fx[l+2]);
    if (Y(rescale)(&rec1p,&rec2p,&scalep) | Y(rescale)(&rec1m,&rec2m,&scalem))
      below_limit = Y(TballLt)(scalep,sharp_limscale)
                 && Y(TballLt)(scalem,sharp_limscale);
    l+=2;
    }

  *l_=l;
  *rec1p_=rec1p; *rec2p_=rec2p; *scalep_=scalep;
  *rec1m_=rec1m; *rec2m_=rec2m; *scalem_=scalem;
  }


NOINLINE static void Y(alm2map_kernel) (const Tb cth, Y(Tbri) * restrict p1,
  Y(Tbri) * restrict p2, Tb lam_1, Tb lam_2,
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
      lam_1.v[i] = f10*(cth.v[i]*lam_2.v[i]) - f11*lam_1.v[i];
      p1->r.v[i] += lam_2.v[i]*ar1;
      p1->i.v[i] += lam_2.v[i]*ai1;
      lam_2.v[i] = f20*(cth.v[i]*lam_1.v[i]) - f21*lam_2.v[i];
      p2->r.v[i] += lam_1.v[i]*ar2;
      p2->i.v[i] += lam_1.v[i]*ai2;
      }
    l+=2;
    }
  }

NOINLINE static void Y(map2alm_kernel) (const Tb cth,
  const Y(Tbri) * restrict p1, const Y(Tbri) * restrict p2, Tb lam_1, Tb lam_2,
  const sharp_ylmgen_dbl2 * restrict rf, int l, int lmax, Tv *restrict atmp)
  {
  while (l<=lmax)
    {
    Tv f10=vload(rf[l  ].f[0]), f11=vload(rf[l  ].f[1]),
       f20=vload(rf[l+1].f[0]), f21=vload(rf[l+1].f[1]);
    for (int i=0; i<nvec; ++i)
      {
      lam_1.v[i] = f10*(cth.v[i]*lam_2.v[i]) - f11*lam_1.v[i];
      vfmaeq(atmp[2*l  ],lam_2.v[i],p1->r.v[i]);
      vfmaeq(atmp[2*l+1],lam_2.v[i],p1->i.v[i]);
      lam_2.v[i] = f20*(cth.v[i]*lam_1.v[i]) - f21*lam_2.v[i];
      vfmaeq(atmp[2*(l+1)  ],lam_1.v[i],p2->r.v[i]);
      vfmaeq(atmp[2*(l+1)+1],lam_1.v[i],p2->i.v[i]);
      }
    l+=2;
    }
  }

NOINLINE static void Y(calc_alm2map) (const Tb cth, const Tb sth,
  const sharp_Ylmgen_C *gen, sharp_job *job, Y(Tbri) * restrict p1,
  Y(Tbri) * restrict p2)
  {
  int l,lmax=gen->lmax;
  Tb lam_1=Y(Tbconst)(0.),lam_2=Y(Tbconst)(0.),scale;
  Y(iter_to_ieee) (sth,cth,&l,&lam_1,&lam_2,&scale,gen);
  job->opcnt += (l-gen->m) * 4*VLEN*nvec;
  if (l>lmax) return;
  job->opcnt += (lmax+1-l) * 8*VLEN*nvec;

  Tb corfac;
  Y(getCorfac)(scale,&corfac,gen->cf);
  const sharp_ylmgen_dbl2 * restrict rf = gen->rf;
  const dcmplx * restrict alm=job->almtmp;
  int full_ieee = Y(TballGe)(scale,sharp_minscale);
  while (!full_ieee)
    {
    {
    Tv ar=vload(creal(alm[l])),ai=vload(cimag(alm[l]));
    for (int i=0; i<nvec; ++i)
      {
      Tv tmp=vmul(lam_2.v[i],corfac.v[i]);
      vfmaeq(p1->r.v[i],tmp,ar);
      vfmaeq(p1->i.v[i],tmp,ai);
      }
    }
    if (++l>lmax) break;
    Tv r0=vload(rf[l-1].f[0]),r1=vload(rf[l-1].f[1]);
    for (int i=0; i<nvec; ++i)
      lam_1.v[i] = vsub(vmul(vmul(cth.v[i],lam_2.v[i]),r0),vmul(lam_1.v[i],r1));
    {
    Tv ar=vload(creal(alm[l])),ai=vload(cimag(alm[l]));
    for (int i=0; i<nvec; ++i)
      {
      Tv tmp=vmul(lam_1.v[i],corfac.v[i]);
      vfmaeq(p2->r.v[i],tmp,ar);
      vfmaeq(p2->i.v[i],tmp,ai);
      }
    }
    if (++l>lmax) break;
    r0=vload(rf[l-1].f[0]); r1=vload(rf[l-1].f[1]);
    for (int i=0; i<nvec; ++i)
      lam_2.v[i] = vsub(vmul(vmul(cth.v[i],lam_1.v[i]),r0),vmul(lam_2.v[i],r1));
    if (Y(rescale)(&lam_1,&lam_2,&scale))
      {
      Y(getCorfac)(scale,&corfac,gen->cf);
      full_ieee = Y(TballGe)(scale,sharp_minscale);
      }
    }
  if (l>lmax) return;

  Y(Tbmuleq)(&lam_1,corfac); Y(Tbmuleq)(&lam_2,corfac);
  Y(alm2map_kernel) (cth, p1, p2, lam_1, lam_2, rf, alm, l, lmax);
  }

NOINLINE static void Y(calc_map2alm) (const Tb cth, const Tb sth,
  const sharp_Ylmgen_C *gen, sharp_job *job, const Y(Tbri) * restrict p1,
  const Y(Tbri) * restrict p2, Tv *restrict atmp)
  {
  int lmax=gen->lmax;
  Tb lam_1=Y(Tbconst)(0.),lam_2=Y(Tbconst)(0.),scale;
  int l=gen->m;
  Y(iter_to_ieee) (sth,cth,&l,&lam_1,&lam_2,&scale,gen);
  job->opcnt += (l-gen->m) * 4*VLEN*nvec;
  if (l>lmax) return;
  job->opcnt += (lmax+1-l) * 8*VLEN*nvec;

  const sharp_ylmgen_dbl2 * restrict rf = gen->rf;
  Tb corfac;
  Y(getCorfac)(scale,&corfac,gen->cf);
  int full_ieee = Y(TballGe)(scale,sharp_minscale);
  while (!full_ieee)
    {
    for (int i=0; i<nvec; ++i)
      {
      Tv tmp=lam_2.v[i]*corfac.v[i];
      atmp[2*l  ]+=tmp*p1->r.v[i];
      atmp[2*l+1]+=tmp*p1->i.v[i];
      }
    if (++l>lmax) return;
    for (int i=0; i<nvec; ++i)
      {
      lam_1.v[i] = vload(rf[l-1].f[0])*(cth.v[i]*lam_2.v[i])
                 - vload(rf[l-1].f[1])*lam_1.v[i];
      Tv tmp=lam_1.v[i]*corfac.v[i];
      atmp[2*l  ]+=tmp*p2->r.v[i];
      atmp[2*l+1]+=tmp*p2->i.v[i];
      }
    if (++l>lmax) return;
    for (int i=0; i<nvec; ++i)
      lam_2.v[i] = vload(rf[l-1].f[0])*(cth.v[i]*lam_1.v[i])
                 - vload(rf[l-1].f[1])*lam_2.v[i];
    if (Y(rescale)(&lam_1,&lam_2,&scale))
      {
      Y(getCorfac)(scale,&corfac,gen->cf);
      full_ieee = Y(TballGe)(scale,sharp_minscale);
      }
    }

  Y(Tbmuleq)(&lam_1,corfac); Y(Tbmuleq)(&lam_2,corfac);
  Y(map2alm_kernel) (cth, p1, p2, lam_1, lam_2, rf, l, lmax, atmp);
  }

static inline void Y(saddstep) (Y(Tbqu) * restrict px, Y(Tbqu) * restrict py,
  const Tb rxp, const Tb rxm, const dcmplx * restrict alm)
  {
  Tv agr=vload(creal(alm[0])), agi=vload(cimag(alm[0])),
     acr=vload(creal(alm[1])), aci=vload(cimag(alm[1]));
  for (int i=0; i<nvec; ++i)
    {
    Tv lw=vadd(rxp.v[i],rxm.v[i]);
    vfmaeq(px->qr.v[i],agr,lw);
    vfmaeq(px->qi.v[i],agi,lw);
    vfmaeq(px->ur.v[i],acr,lw);
    vfmaeq(px->ui.v[i],aci,lw);
    Tv lx=vsub(rxm.v[i],rxp.v[i]);
    vfmseq(py->qr.v[i],aci,lx);
    vfmaeq(py->qi.v[i],acr,lx);
    vfmaeq(py->ur.v[i],agi,lx);
    vfmseq(py->ui.v[i],agr,lx);
    }
  }

static inline void Y(saddstepb) (Y(Tbqu) * restrict p1, Y(Tbqu) * restrict p2,
  const Tb r1p, const Tb r1m, const Tb r2p, const Tb r2m,
  const dcmplx * restrict alm1, const dcmplx * restrict alm2)
  {
  Tv agr1=vload(creal(alm1[0])), agi1=vload(cimag(alm1[0])),
     acr1=vload(creal(alm1[1])), aci1=vload(cimag(alm1[1]));
  Tv agr2=vload(creal(alm2[0])), agi2=vload(cimag(alm2[0])),
     acr2=vload(creal(alm2[1])), aci2=vload(cimag(alm2[1]));
  for (int i=0; i<nvec; ++i)
    {
    Tv lw1=r2p.v[i]+r2m.v[i];
    Tv lx2=r1m.v[i]-r1p.v[i];
    vfmaseq(p1->qr.v[i],agr1,lw1,aci2,lx2);
    vfmaaeq(p1->qi.v[i],agi1,lw1,acr2,lx2);
    vfmaaeq(p1->ur.v[i],acr1,lw1,agi2,lx2);
    vfmaseq(p1->ui.v[i],aci1,lw1,agr2,lx2);
    Tv lx1=r2m.v[i]-r2p.v[i];
    Tv lw2=r1p.v[i]+r1m.v[i];
    vfmaseq(p2->qr.v[i],agr2,lw2,aci1,lx1);
    vfmaaeq(p2->qi.v[i],agi2,lw2,acr1,lx1);
    vfmaaeq(p2->ur.v[i],acr2,lw2,agi1,lx1);
    vfmaseq(p2->ui.v[i],aci2,lw2,agr1,lx1);
    }
  }

static inline void Y(saddstep2) (const Y(Tbqu) * restrict px,
  const Y(Tbqu) * restrict py, const Tb * restrict rxp,
  const Tb * restrict rxm, dcmplx * restrict alm)
  {
  Tv agr=vzero, agi=vzero, acr=vzero, aci=vzero;
  for (int i=0; i<nvec; ++i)
    {
    Tv lw=vadd(rxp->v[i],rxm->v[i]);
    vfmaeq(agr,px->qr.v[i],lw);
    vfmaeq(agi,px->qi.v[i],lw);
    vfmaeq(acr,px->ur.v[i],lw);
    vfmaeq(aci,px->ui.v[i],lw);
    Tv lx=vsub(rxm->v[i],rxp->v[i]);
    vfmseq(agr,py->ui.v[i],lx);
    vfmaeq(agi,py->ur.v[i],lx);
    vfmaeq(acr,py->qi.v[i],lx);
    vfmseq(aci,py->qr.v[i],lx);
    }
  vhsum_cmplx_special(agr,agi,acr,aci,alm);
  }

NOINLINE static void Y(alm2map_spin_kernel) (Tb cth, Y(Tbqu) * restrict p1,
  Y(Tbqu) * restrict p2, Tb rec1p, Tb rec1m, Tb rec2p, Tb rec2m,
  const sharp_ylmgen_dbl3 * restrict fx, const dcmplx * restrict alm, int l,
  int lmax)
  {
  while (l<lmax)
    {
    Tv fx0=vload(fx[l+1].f[0]),fx1=vload(fx[l+1].f[1]),
       fx2=vload(fx[l+1].f[2]);
    for (int i=0; i<nvec; ++i)
      {
      rec1p.v[i] = (cth.v[i]-fx1)*fx0*rec2p.v[i] - fx2*rec1p.v[i];
      rec1m.v[i] = (cth.v[i]+fx1)*fx0*rec2m.v[i] - fx2*rec1m.v[i];
      }
    Y(saddstepb)(p1,p2,rec1p,rec1m,rec2p,rec2m,&alm[2*l],
      &alm[2*(l+1)]);
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
    Y(saddstep)(p1, p2, rec2p, rec2m, &alm[2*l]);
  }

NOINLINE static void Y(map2alm_spin_kernel) (Tb cth, const Y(Tbqu) * restrict p1,
  const Y(Tbqu) * restrict p2, Tb rec1p, Tb rec1m, Tb rec2p, Tb rec2m,
  const sharp_ylmgen_dbl3 * restrict fx, dcmplx * restrict alm, int l, int lmax)
  {
  while (l<lmax)
    {
    Tv fx0=vload(fx[l+1].f[0]),fx1=vload(fx[l+1].f[1]),
       fx2=vload(fx[l+1].f[2]);
    for (int i=0; i<nvec; ++i)
      {
      rec1p.v[i] = vsub(vmul(vsub(cth.v[i],fx1),vmul(fx0,rec2p.v[i])),
                        vmul(fx2,rec1p.v[i]));
      rec1m.v[i] = vsub(vmul(vadd(cth.v[i],fx1),vmul(fx0,rec2m.v[i])),
                        vmul(fx2,rec1m.v[i]));
      }
    Y(saddstep2)(p1, p2, &rec2p, &rec2m, &alm[2*l]);
    Y(saddstep2)(p2, p1, &rec1p, &rec1m, &alm[2*(l+1)]);
    fx0=vload(fx[l+2].f[0]);fx1=vload(fx[l+2].f[1]);
    fx2=vload(fx[l+2].f[2]);
    for (int i=0; i<nvec; ++i)
      {
      rec2p.v[i] = vsub(vmul(vsub(cth.v[i],fx1),vmul(fx0,rec1p.v[i])),
                        vmul(fx2,rec2p.v[i]));
      rec2m.v[i] = vsub(vmul(vadd(cth.v[i],fx1),vmul(fx0,rec1m.v[i])),
                        vmul(fx2,rec2m.v[i]));
      }
    l+=2;
    }
  if (l==lmax)
    Y(saddstep2)(p1, p2, &rec2p, &rec2m, &alm[2*l]);
  }

NOINLINE static void Y(calc_alm2map_spin) (const Tb cth, const Tb sth,
  const sharp_Ylmgen_C *gen, sharp_job *job, Y(Tbqu) * restrict p1,
  Y(Tbqu) * restrict p2)
  {
  int l, lmax=gen->lmax;
  Tb rec1p, rec1m, rec2p, rec2m, scalem, scalep;
  Y(iter_to_ieee_spin)
    (cth,sth,&l,&rec1p,&rec1m,&rec2p,&rec2m,&scalep,&scalem,gen);
  job->opcnt += (l-gen->m) * 10*VLEN*nvec;
  if (l>lmax) return;
  job->opcnt += (lmax+1-l) * 28*VLEN*nvec;

  const sharp_ylmgen_dbl3 * restrict fx = gen->fx;
  Tb corfacp,corfacm;
  Y(getCorfac)(scalep,&corfacp,gen->cf);
  Y(getCorfac)(scalem,&corfacm,gen->cf);
  const dcmplx * restrict alm=job->almtmp;
  int full_ieee = Y(TballGe)(scalep,sharp_minscale)
               && Y(TballGe)(scalem,sharp_minscale);
  while (!full_ieee)
    {
    Y(saddstep)(p1, p2, Y(Tbprod)(rec2p,corfacp), Y(Tbprod)(rec2m,corfacm),
      &alm[2*l]);
    if (++l>lmax) break;
    Y(rec_step)(&rec1p,&rec1m,&rec2p,&rec2m,cth,fx[l]);
    Y(saddstep)(p2, p1, Y(Tbprod)(rec1p,corfacp), Y(Tbprod)(rec1m,corfacm),
      &alm[2*l]);
    if (++l>lmax) break;
    Y(rec_step)(&rec2p,&rec2m,&rec1p,&rec1m,cth,fx[l]);
    if (Y(rescale)(&rec1p,&rec2p,&scalep) | Y(rescale)(&rec1m,&rec2m,&scalem))
      {
      Y(getCorfac)(scalep,&corfacp,gen->cf);
      Y(getCorfac)(scalem,&corfacm,gen->cf);
      full_ieee = Y(TballGe)(scalep,sharp_minscale)
               && Y(TballGe)(scalem,sharp_minscale);
      }
    }

  if (l>lmax) return;

  Y(Tbmuleq)(&rec1p,corfacp); Y(Tbmuleq)(&rec2p,corfacp);
  Y(Tbmuleq)(&rec1m,corfacm); Y(Tbmuleq)(&rec2m,corfacm);
  Y(alm2map_spin_kernel) (cth, p1, p2, rec1p, rec1m, rec2p, rec2m, fx, alm, l,
    lmax);
  }

NOINLINE static void Y(calc_map2alm_spin) (Tb cth, Tb sth,
  const sharp_Ylmgen_C * restrict gen, sharp_job *job,
  const Y(Tbqu) * restrict p1, const Y(Tbqu) * restrict p2)
  {
  int l, lmax=gen->lmax;
  Tb rec1p, rec1m, rec2p, rec2m, scalem, scalep;
  Y(iter_to_ieee_spin)
    (cth,sth,&l,&rec1p,&rec1m,&rec2p,&rec2m,&scalep,&scalem,gen);
  job->opcnt += (l-gen->m) * 10*VLEN*nvec;
  if (l>lmax) return;
  job->opcnt += (lmax+1-l) * 28*VLEN*nvec;

  const sharp_ylmgen_dbl3 * restrict fx = gen->fx;
  Tb corfacp,corfacm;
  Y(getCorfac)(scalep,&corfacp,gen->cf);
  Y(getCorfac)(scalem,&corfacm,gen->cf);
  dcmplx * restrict alm=job->almtmp;
  int full_ieee = Y(TballGe)(scalep,sharp_minscale)
               && Y(TballGe)(scalem,sharp_minscale);
  while (!full_ieee)
    {
    Tb t1=Y(Tbprod)(rec2p,corfacp), t2=Y(Tbprod)(rec2m,corfacm);
    Y(saddstep2)(p1, p2, &t1, &t2, &alm[2*l]);
    if (++l>lmax) return;
    Y(rec_step)(&rec1p,&rec1m,&rec2p,&rec2m,cth,fx[l]);
    t1=Y(Tbprod)(rec1p,corfacp); t2=Y(Tbprod)(rec1m,corfacm);
    Y(saddstep2)(p2, p1, &t1, &t2, &alm[2*l]);
    if (++l>lmax) return;
    Y(rec_step)(&rec2p,&rec2m,&rec1p,&rec1m,cth,fx[l]);
    if (Y(rescale)(&rec1p,&rec2p,&scalep) | Y(rescale)(&rec1m,&rec2m,&scalem))
      {
      Y(getCorfac)(scalep,&corfacp,gen->cf);
      Y(getCorfac)(scalem,&corfacm,gen->cf);
      full_ieee = Y(TballGe)(scalep,sharp_minscale)
               && Y(TballGe)(scalem,sharp_minscale);
      }
    }

  Y(Tbmuleq)(&rec1p,corfacp); Y(Tbmuleq)(&rec2p,corfacp);
  Y(Tbmuleq)(&rec1m,corfacm); Y(Tbmuleq)(&rec2m,corfacm);
  Y(map2alm_spin_kernel)(cth,p1,p2,rec1p,rec1m,rec2p,rec2m,fx,alm,l,lmax);
  }

static inline void Y(saddstep_d) (Y(Tbqu) * restrict px, Y(Tbqu) * restrict py,
  const Tb rxp, const Tb rxm, const dcmplx * restrict alm)
  {
  Tv ar=vload(creal(alm[0])), ai=vload(cimag(alm[0]));
  for (int i=0; i<nvec; ++i)
    {
    Tv lw=vadd(rxp.v[i],rxm.v[i]);
    vfmaeq(px->qr.v[i],ar,lw);
    vfmaeq(px->qi.v[i],ai,lw);
    Tv lx=vsub(rxm.v[i],rxp.v[i]);
    vfmaeq(py->ur.v[i],ai,lx);
    vfmseq(py->ui.v[i],ar,lx);
    }
  }

NOINLINE static void Y(alm2map_deriv1_kernel) (Tb cth, Y(Tbqu) * restrict p1,
  Y(Tbqu) * restrict p2, Tb rec1p, Tb rec1m, Tb rec2p, Tb rec2m,
  const sharp_ylmgen_dbl3 * restrict fx, const dcmplx * restrict alm, int l,
  int lmax)
  {
  while (l<lmax)
    {
    Tv fx0=vload(fx[l+1].f[0]),fx1=vload(fx[l+1].f[1]),
       fx2=vload(fx[l+1].f[2]);
    for (int i=0; i<nvec; ++i)
      {
      rec1p.v[i] = vsub(vmul(vsub(cth.v[i],fx1),vmul(fx0,rec2p.v[i])),
                        vmul(fx2,rec1p.v[i]));
      rec1m.v[i] = vsub(vmul(vadd(cth.v[i],fx1),vmul(fx0,rec2m.v[i])),
                        vmul(fx2,rec1m.v[i]));
      }
    Y(saddstep_d)(p1,p2,rec2p,rec2m,&alm[l]);
    Y(saddstep_d)(p2,p1,rec1p,rec1m,&alm[l+1]);
    fx0=vload(fx[l+2].f[0]);fx1=vload(fx[l+2].f[1]);
    fx2=vload(fx[l+2].f[2]);
    for (int i=0; i<nvec; ++i)
      {
      rec2p.v[i] = vsub(vmul(vsub(cth.v[i],fx1),vmul(fx0,rec1p.v[i])),
                        vmul(fx2,rec2p.v[i]));
      rec2m.v[i] = vsub(vmul(vadd(cth.v[i],fx1),vmul(fx0,rec1m.v[i])),
                        vmul(fx2,rec2m.v[i]));
      }
    l+=2;
    }
  if (l==lmax)
    Y(saddstep_d)(p1, p2, rec2p, rec2m, &alm[l]);
  }

NOINLINE static void Y(calc_alm2map_deriv1) (const Tb cth, const Tb sth,
  const sharp_Ylmgen_C *gen, sharp_job *job, Y(Tbqu) * restrict p1,
  Y(Tbqu) * restrict p2)
  {
  int l, lmax=gen->lmax;
  Tb rec1p, rec1m, rec2p, rec2m, scalem, scalep;
  Y(iter_to_ieee_spin)
    (cth,sth,&l,&rec1p,&rec1m,&rec2p,&rec2m,&scalep,&scalem,gen);
  job->opcnt += (l-gen->m) * 10*VLEN*nvec;
  if (l>lmax) return;
  job->opcnt += (lmax+1-l) * 20*VLEN*nvec;

  const sharp_ylmgen_dbl3 * restrict fx = gen->fx;
  Tb corfacp,corfacm;
  Y(getCorfac)(scalep,&corfacp,gen->cf);
  Y(getCorfac)(scalem,&corfacm,gen->cf);
  const dcmplx * restrict alm=job->almtmp;
  int full_ieee = Y(TballGe)(scalep,sharp_minscale)
               && Y(TballGe)(scalem,sharp_minscale);
  while (!full_ieee)
    {
    Y(saddstep_d)(p1, p2, Y(Tbprod)(rec2p,corfacp), Y(Tbprod)(rec2m,corfacm),
      &alm[l]);
    if (++l>lmax) break;
    Y(rec_step)(&rec1p,&rec1m,&rec2p,&rec2m,cth,fx[l]);
    Y(saddstep_d)(p2, p1, Y(Tbprod)(rec1p,corfacp), Y(Tbprod)(rec1m,corfacm),
      &alm[l]);
    if (++l>lmax) break;
    Y(rec_step)(&rec2p,&rec2m,&rec1p,&rec1m,cth,fx[l]);
    if (Y(rescale)(&rec1p,&rec2p,&scalep) | Y(rescale)(&rec1m,&rec2m,&scalem))
      {
      Y(getCorfac)(scalep,&corfacp,gen->cf);
      Y(getCorfac)(scalem,&corfacm,gen->cf);
      full_ieee = Y(TballGe)(scalep,sharp_minscale)
               && Y(TballGe)(scalem,sharp_minscale);
      }
    }

  if (l>lmax) return;

  Y(Tbmuleq)(&rec1p,corfacp); Y(Tbmuleq)(&rec2p,corfacp);
  Y(Tbmuleq)(&rec1m,corfacm); Y(Tbmuleq)(&rec2m,corfacm);
  Y(alm2map_deriv1_kernel) (cth, p1, p2, rec1p, rec1m, rec2p, rec2m, fx, alm, l,
    lmax);
  }


#define VZERO(var) do { memset(&(var),0,sizeof(var)); } while(0)

NOINLINE static void Y(inner_loop_a2m) (sharp_job *job, const int *ispair,
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
          Y(Tburi) p1,p2; VZERO(p1); VZERO(p2);
          Y(Tbu) cth, sth;

          int skip=1;
          for (int i=0; i<nval; ++i)
            {
            int itot=i+ith;
            if (itot>=ulim-llim) itot=ulim-llim-1;
            if (mlim[itot]>=m) skip=0;
            cth.s[i]=cth_[itot]; sth.s[i]=sth_[itot];
            }
          if (!skip)
            Y(calc_alm2map) (cth.b,sth.b,gen,job,&p1.b,&p2.b);

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
        for (int ith=0; ith<ulim-llim; ith+=nval)
          {
          Y(Tbuqu) p1,p2; VZERO(p1); VZERO(p2);
          Y(Tbu) cth, sth;
          int skip=1;

          for (int i=0; i<nval; ++i)
            {
            int itot=i+ith;
            if (itot>=ulim-llim) itot=ulim-llim-1;
            if (mlim[itot]>=m) skip=0;
            cth.s[i]=cth_[itot]; sth.s[i]=sth_[itot];
            }
          if (!skip)
            (job->type==SHARP_ALM2MAP) ?
              Y(calc_alm2map_spin  )
                (cth.b,sth.b,gen,job,&p1.b,&p2.b) :
              Y(calc_alm2map_deriv1)
                (cth.b,sth.b,gen,job,&p1.b,&p2.b);

          for (int i=0; i<nval; ++i)
            {
            int itot=i+ith;
            if (itot<ulim-llim)
              {
              int phas_idx = itot*job->s_th + mi*job->s_m;
              complex double q1 = p1.s.qr[i] + p1.s.qi[i]*_Complex_I,
                             q2 = p2.s.qr[i] + p2.s.qi[i]*_Complex_I,
                             u1 = p1.s.ur[i] + p1.s.ui[i]*_Complex_I,
                             u2 = p2.s.ur[i] + p2.s.ui[i]*_Complex_I;
              job->phase[phas_idx] = q1+q2;
              job->phase[phas_idx+2] = u1+u2;
              if (ispair[itot])
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

NOINLINE static void Y(inner_loop_m2a) (sharp_job *job, const int *ispair,
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
          Y(Tburi) p1, p2; VZERO(p1); VZERO(p2);
          Y(Tbu) cth, sth;
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
            Y(calc_map2alm)(cth.b,sth.b,gen,job,&p1.b,&p2.b, atmp);
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
        for (int ith=0; ith<ulim-llim; ith+=nval)
          {
          Y(Tbuqu) p1, p2; VZERO(p1); VZERO(p2);
          Y(Tbu) cth, sth;
          int skip=1;

          for (int i=0; i<nval; ++i)
            {
            int itot=i+ith;
            if (itot>=ulim-llim) itot=ulim-llim-1;
            if (mlim[itot]>=m) skip=0;
            cth.s[i]=cth_[itot]; sth.s[i]=sth_[itot];
            if (i+ith<ulim-llim)
              {
              int phas_idx = itot*job->s_th + mi*job->s_m;
              dcmplx p1Q=job->phase[phas_idx],
                     p1U=job->phase[phas_idx+2],
                     p2Q=ispair[itot] ? job->phase[phas_idx+1]:0.,
                     p2U=ispair[itot] ? job->phase[phas_idx+3]:0.;
              if ((gen->mhi-gen->m+gen->s)&1)
                { p2Q=-p2Q; p2U=-p2U; }
              p1.s.qr[i]=creal(p1Q+p2Q); p1.s.qi[i]=cimag(p1Q+p2Q);
              p1.s.ur[i]=creal(p1U+p2U); p1.s.ui[i]=cimag(p1U+p2U);
              p2.s.qr[i]=creal(p1Q-p2Q); p2.s.qi[i]=cimag(p1Q-p2Q);
              p2.s.ur[i]=creal(p1U-p2U); p2.s.ui[i]=cimag(p1U-p2U);
              }
            }
          if (!skip)
            Y(calc_map2alm_spin) (cth.b,sth.b,gen,job,&p1.b,&p2.b);
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

static void Y(inner_loop) (sharp_job *job, const int *ispair,
  const double *cth_, const double *sth_, int llim, int ulim,
  sharp_Ylmgen_C *gen, int mi, const int *mlim)
  {
  (job->type==SHARP_MAP2ALM) ?
    Y(inner_loop_m2a)(job,ispair,cth_,sth_,llim,ulim,gen,mi,mlim) :
    Y(inner_loop_a2m)(job,ispair,cth_,sth_,llim,ulim,gen,mi,mlim);
  }

#undef VZERO
