#ifndef CCTBX_BOOST_PYTHON_FLEX_FWD_H
#define CCTBX_BOOST_PYTHON_FLEX_FWD_H

#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <cctbx/maptbx/accessors/c_grid_p1.h>
#include <cctbx/maptbx/accessors/c_grid_padded_p1.h>
#include <cctbx/miller.h>
#include <cctbx/hendrickson_lattman.h>
#include <cctbx/xray/scatterer.h>

#if defined(__sgi) && !defined(__GNUC__)

namespace cctbx { namespace boost_python {

  template <typename T>
  struct flex_fwd
  {
    friend void f(scitbx::af::versa<T, maptbx::c_grid_p1<3> > const&);
    friend void f(scitbx::af::versa<T, maptbx::c_grid_padded_p1<3> > const&);
    friend void f(scitbx::af::ref<T, maptbx::c_grid_p1<3> > const&);
    friend void f(scitbx::af::ref<T, maptbx::c_grid_padded_p1<3> > const&);
    friend void f(scitbx::af::const_ref<T, maptbx::c_grid_p1<3> > const&);
    friend void f(scitbx::af::const_ref<T, maptbx::c_grid_padded_p1<3> > const&);
  };

  inline void
  flex_fwd_types()
  {
    flex_fwd<bool>();
    flex_fwd<int>();
    flex_fwd<long>();
    flex_fwd<std::size_t>();
    flex_fwd<float>();
    flex_fwd<double>();
    flex_fwd<std::complex<double> >();
    flex_fwd<std::string>();

    scitbx::af::boost_python::flex_fwd<scitbx::af::tiny<std::size_t, 2> >();
    scitbx::af::boost_python::flex_fwd<scitbx::sym_mat3<double> >();
    scitbx::af::boost_python::flex_fwd<miller::index<> >();
    scitbx::af::boost_python::flex_fwd<hendrickson_lattman<> >();
    scitbx::af::boost_python::flex_fwd<xray::scatterer<> >();
  }

}} // namespace cctbx::boost_python

#endif // defined(__sgi) && !defined(__GNUC__)

#endif // CCTBX_BOOST_PYTHON_FLEX_FWD_H
