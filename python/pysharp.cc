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
 *
 *  For more information about HEALPix, see http://healpix.sourceforge.net
 */

/*
 *  libsharp is being developed at the Max-Planck-Institut fuer Astrophysik
 */

/*
 *  Copyright (C) 2017 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <iostream>
#include <vector>
#include <complex>
#include <string>

#include "libsharp/sharp_cxx.h"
#include "libsharp/sharp_legendre_roots.h"

using namespace std;

namespace py = pybind11;

namespace {

template<typename T> using pyarr = py::array_t<T>;
template<typename T> using pyarr_c
  = py::array_t<T, py::array::c_style | py::array::forcecast>;

using a_d = py::array_t<double>;
using a_c = py::array_t<complex<double>>;
using a_i = py::array_t<int64_t>;
using a_d_c = py::array_t<double, py::array::c_style | py::array::forcecast>;
using a_c_c = py::array_t<complex<double>,
  py::array::c_style | py::array::forcecast>;

void myassert(bool cond, const char *msg)
  { if (!cond) throw runtime_error(msg); }

template<typename T> class py_sharpjob
  {
  private:
    sharp_cxxjob<T> job;
    int64_t lmax_, mmax_, npix_;

  public:
    py_sharpjob () : lmax_(0), mmax_(0), npix_(0) {}

    void set_Gauss_geometry(int64_t nrings, int64_t nphi)
      {
      myassert((nrings>0)&&(nphi>0),"bad grid dimensions");
      npix_=nrings*nphi;
      job.set_Gauss_geometry(nrings, nphi);
      }
    void set_Healpix_geometry(int64_t nside)
      {
      myassert(nside>0,"bad Nside value");
      npix_=12*nside*nside;
      job.set_Healpix_geometry(nside);
      }
    void set_triangular_alm_info (int64_t lmax, int64_t mmax)
      {
      myassert(mmax>=0,"negative mmax");
      myassert(mmax<=lmax,"mmax must not be larger than lmax");
      lmax_=lmax; mmax_=mmax;
      job.set_triangular_alm_info (lmax,mmax);
      }

    int64_t n_alm() const
      { return ((mmax_+1)*(mmax_+2))/2 + (mmax_+1)*(lmax_-mmax_); }

    a_d_c alm2map (const a_c_c &alm) const
      {
      myassert(npix_>0,"no map geometry specified");
      myassert (alm.size()==n_alm(),
        "incorrect size of a_lm array");
      a_d_c map(npix_);
      auto mr=map.mutable_unchecked<1>();
      auto ar=alm.unchecked<1>();
      job.alm2map(&ar[0],&mr[0],false);
      return map;
      }
    a_c_c alm2map_adjoint (const a_d_c &map) const
      {
      myassert(npix_>0,"no map geometry specified");
      myassert (map.size()==npix_,"incorrect size of map array");
      a_c_c alm(n_alm());
      auto mr=map.unchecked<1>();
      auto ar=alm.mutable_unchecked<1>();
      job.alm2map_adjoint(&mr[0],&ar[0],false);
      return alm;
      }
    a_c_c map2alm (const a_d_c &map) const
      {
      myassert(npix_>0,"no map geometry specified");
      myassert (map.size()==npix_,"incorrect size of map array");
      a_c_c alm(n_alm());
      auto mr=map.unchecked<1>();
      auto ar=alm.mutable_unchecked<1>();
      job.map2alm(&mr[0],&ar[0],false);
      return alm;
      }
    a_d_c alm2map_spin (const a_c_c &alm, int64_t spin) const
      {
      myassert(npix_>0,"no map geometry specified");
      auto ar=alm.unchecked<2>();
      myassert((ar.shape(0)==2)&&(ar.shape(1)==n_alm()),
        "incorrect size of a_lm array");
      a_d_c map(vector<size_t>{2,size_t(npix_)});
      auto mr=map.mutable_unchecked<2>();
      job.alm2map_spin(&ar(0,0),&ar(1,0),&mr(0,0),&mr(1,0),spin,false);
      return map;
      }
    a_c_c map2alm_spin (const a_d_c &map, int64_t spin) const
      {
      myassert(npix_>0,"no map geometry specified");
      auto mr=map.unchecked<2>();
      myassert ((mr.shape(0)==2)&&(mr.shape(1)==npix_),
        "incorrect size of map array");
      a_c_c alm(vector<size_t>{2,size_t(n_alm())});
      auto ar=alm.mutable_unchecked<2>();
      job.map2alm_spin(&mr(0,0),&mr(1,0),&ar(0,0),&ar(1,0),spin,false);
      return alm;
      }
  };

a_d_c GL_weights(int64_t nlat, int64_t nlon)
  {
  constexpr double twopi=6.283185307179586476925286766559005768394;
  a_d_c res(nlat);
  auto rr=res.mutable_unchecked<1>();
  vector<double> dummy_roots(nlat);
  sharp_legendre_roots(nlat, dummy_roots.data(), &rr[0]);
  for (size_t i=0; i<size_t(rr.shape(0)); ++i)
    rr[i]*=twopi/nlon;
  return res;
  }

a_d_c GL_thetas(int64_t nlat)
  {
  a_d_c res(nlat);
  auto rr=res.mutable_unchecked<1>();
  vector<double> dummy_weights(nlat);
  sharp_legendre_roots(nlat, &rr[0], dummy_weights.data());
  for (size_t i=0; i<size_t(rr.shape(0)); ++i)
    rr[i]=acos(-rr[i]);
  return res;
  }


const char *pysharp_DS = R"""(
Python interface for libsharp

All angles are interpreted as radians.
The theta coordinate is measured as co-latitude, ranging from 0 (North Pole)
to pi (South Pole).

Error conditions are reported by raising exceptions.
)""";


} // unnamed namespace

PYBIND11_MODULE(pysharp, m)
  {
  using namespace pybind11::literals;

  m.doc() = pysharp_DS;

  py::class_<py_sharpjob<double>> (m, "sharpjob_d")
    .def(py::init<>())
    .def("set_Gauss_geometry", &py_sharpjob<double>::set_Gauss_geometry,
      "nrings"_a,"nphi"_a)
    .def("set_Healpix_geometry", &py_sharpjob<double>::set_Healpix_geometry,
      "nside"_a)
    .def("set_triangular_alm_info",
      &py_sharpjob<double>::set_triangular_alm_info, "lmax"_a, "mmax"_a)
    .def("n_alm", &py_sharpjob<double>::n_alm)
    .def("alm2map", &py_sharpjob<double>::alm2map,"alm"_a)
    .def("alm2map_adjoint", &py_sharpjob<double>::alm2map_adjoint,"map"_a)
    .def("map2alm", &py_sharpjob<double>::map2alm,"map"_a)
    .def("alm2map_spin", &py_sharpjob<double>::alm2map_spin,"alm"_a,"spin"_a)
    .def("map2alm_spin", &py_sharpjob<double>::map2alm_spin,"map"_a,"spin"_a)
    ;

  m.def("GL_weights",&GL_weights, "nlat"_a, "nlon"_a);
  m.def("GL_thetas",&GL_thetas, "nlat"_a);
  }
