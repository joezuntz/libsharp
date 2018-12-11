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

static void Tbnormalize (Tb * restrict val, Tb * restrict scale,
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

NOINLINE static void mypow (Tb val, int npow, const double * restrict powlimit,
  Tb * restrict resd, Tb * restrict ress)
  {
  Tv vminv=vload(powlimit[npow]);
  Tm mask = vlt(vabs(val.v[0]),vminv);
  for (int i=1;i<nvec; ++i)
    mask=vor_mask(mask,vlt(vabs(val.v[i]),vminv));
  if (!vanyTrue(mask)) // no underflows possible, use quick algoritm
    {
    Tb res=Tbconst(1.);
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
    *ress=Tbconst(0.);
    }
  else
    {
    Tb scale=Tbconst(0.), scaleint=Tbconst(0.), res=Tbconst(1.);
    Tbnormalize(&val,&scaleint,sharp_fbighalf);
    do
      {
      if (npow&1)
        {
        for (int i=0; i<nvec; ++i)
          {
          vmuleq(res.v[i],val.v[i]);
          vaddeq(scale.v[i],scaleint.v[i]);
          }
        Tbnormalize(&res,&scale,sharp_fbighalf);
        }
      for (int i=0; i<nvec; ++i)
        {
        vmuleq(val.v[i],val.v[i]);
        vaddeq(scaleint.v[i],scaleint.v[i]);
        }
      Tbnormalize(&val,&scaleint,sharp_fbighalf);
      }
    while(npow>>=1);
    *resd=res;
    *ress=scale;
    }
  }

static inline int rescale(Tb * restrict lam1, Tb * restrict lam2,
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
  Tb * restrict lam_1_, Tb * restrict lam_2_, Tb * restrict scale_,
  const sharp_Ylmgen_C * restrict gen)
  {
  int l=gen->m;
  Tb lam_1=Tbconst(0.), lam_2, scale;
  mypow(sth,l,gen->powlimit,&lam_2,&scale);
  Tbmuleq1(&lam_2,(gen->m&1) ? -gen->mfac[gen->m]:gen->mfac[gen->m]);
  Tbnormalize(&lam_2,&scale,sharp_ftol);

  int below_limit = TballLt(scale,sharp_limscale);
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
    if (rescale(&lam_1,&lam_2,&scale))
      below_limit = TballLt(scale,sharp_limscale);
    l+=2;
    }
  *l_=l; *lam_1_=lam_1; *lam_2_=lam_2; *scale_=scale;
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
    if (rescale(&lam_1,&lam_2,&scale))
      {
      getCorfac(scale,&corfac,gen->cf);
      full_ieee = TballGe(scale,sharp_minscale);
      }
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
    if (rescale(&lam_1,&lam_2,&scale))
      {
      getCorfac(scale,&corfac,gen->cf);
      full_ieee = TballGe(scale,sharp_minscale);
      }
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

static void inner_loop_ (sharp_job *job, const int *ispair,
  const double *cth_, const double *sth_, int llim, int ulim,
  sharp_Ylmgen_C *gen, int mi, const int *mlim)
  {
  (job->type==SHARP_MAP2ALM) ?
    inner_loop_m2a(job,ispair,cth_,sth_,llim,ulim,gen,mi,mlim) :
    inner_loop_a2m(job,ispair,cth_,sth_,llim,ulim,gen,mi,mlim);
  }

#undef VZERO
