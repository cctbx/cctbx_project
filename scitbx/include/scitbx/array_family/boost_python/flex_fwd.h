/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Nov: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_FLEX_FWD_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_FLEX_FWD_H

#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/array_family/accessors/c_grid_padded.h>
#include <complex>

namespace scitbx { namespace af { namespace boost_python {

  template <typename T>
  struct flex_fwd
  {
    friend void f(shared_plain<T> const&);
    friend void f(shared<T> const&);
    friend void f(versa<T> const&);
    friend void f(versa<T, flex_grid<> > const&);
    friend void f(versa<T, c_grid<2> > const&);
    friend void f(versa<T, c_grid<3> > const&);
    friend void f(versa<T, c_grid_padded<2> > const&);
    friend void f(versa<T, c_grid_padded<3> > const&);
    friend void f(ref<T> const&);
    friend void f(ref<T, flex_grid<> > const&);
    friend void f(ref<T, c_grid<2> > const&);
    friend void f(ref<T, c_grid<3> > const&);
    friend void f(ref<T, c_grid_padded<2> > const&);
    friend void f(ref<T, c_grid_padded<3> > const&);
    friend void f(const_ref<T> const&);
    friend void f(const_ref<T, flex_grid<> > const&);
    friend void f(const_ref<T, c_grid<2> > const&);
    friend void f(const_ref<T, c_grid<3> > const&);
    friend void f(const_ref<T, c_grid_padded<2> > const&);
    friend void f(const_ref<T, c_grid_padded<3> > const&);
  };

  inline void
  flex_fwd_types()
  {
    flex_fwd<int>();
    flex_fwd<long>();
    flex_fwd<float>();
    flex_fwd<double>();
    flex_fwd<std::complex<double> >();
  }

}}} // namespace scitbx::af::boost_python

#endif // SCITBX_ARRAY_FAMILY_BOOST_PYTHON_FLEX_FWD_H
