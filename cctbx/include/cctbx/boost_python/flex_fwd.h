/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Nov: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_BOOST_PYTHON_FLEX_FWD_H
#define CCTBX_BOOST_PYTHON_FLEX_FWD_H

#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <cctbx/maptbx/accessors/c_grid_p1.h>
#include <cctbx/maptbx/accessors/c_grid_padded_p1.h>
#include <cctbx/miller.h>

namespace cctbx { namespace boost_python {

  template <typename T>
  struct flex_fwd
  {
    friend void f(scitbx::af::versa<T, maptbx::c_grid_p1<3> > const&);
    friend void f(scitbx::af::versa<T, maptbx::c_grid_padded_p1<3> > const&);
  };

  inline void
  flex_fwd_types()
  {
    flex_fwd<int>();
    flex_fwd<long>();
    flex_fwd<float>();
    flex_fwd<double>();
    flex_fwd<std::complex<double> >();

    scitbx::af::boost_python::flex_fwd<scitbx::vec3<double> >();
    scitbx::af::boost_python::flex_fwd<miller::index<> >();
    scitbx::af::boost_python::flex_fwd<scitbx::af::tiny<std::size_t, 2> >();
  }

}} // namespace cctbx::boost_python

#endif // CCTBX_BOOST_PYTHON_FLEX_FWD_H
