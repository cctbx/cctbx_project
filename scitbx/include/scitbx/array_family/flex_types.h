/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copied from cctbx/array_family (R.W. Grosse-Kunstleve)
     2002 Aug: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_ARRAY_FAMILY_FLEX_TYPES_H
#define SCITBX_ARRAY_FAMILY_FLEX_TYPES_H

#include <complex>
#include <scitbx/array_family/flex_grid_accessor.h>
#include <scitbx/array_family/versa.h>

namespace scitbx { namespace af {

  template <typename ElementType>
  struct flex
  {
    typedef versa<ElementType, flex_grid<> > type;
  };

  template <typename ElementType>
  struct flex_const_ref
  {
    typedef const_ref<ElementType, flex_grid<> > type;
  };

  template <typename ElementType>
  struct flex_ref
  {
    typedef ref<ElementType, flex_grid<> > type;
  };

#define SCITBX_ARRAY_FAMILY_FLEX_TYPEDEFS(T, N) \
  typedef flex<T>::type flex_ ## N; \
  typedef flex_const_ref<T>::type flex_ ## N ## _const_ref; \
  typedef flex_ref<T>::type flex_ ## N ## _ref;

  SCITBX_ARRAY_FAMILY_FLEX_TYPEDEFS(bool, bool)
  SCITBX_ARRAY_FAMILY_FLEX_TYPEDEFS(int, int)
  SCITBX_ARRAY_FAMILY_FLEX_TYPEDEFS(long, long)
  SCITBX_ARRAY_FAMILY_FLEX_TYPEDEFS(std::size_t, size_t)
  SCITBX_ARRAY_FAMILY_FLEX_TYPEDEFS(float, float)
  SCITBX_ARRAY_FAMILY_FLEX_TYPEDEFS(double, double)
  SCITBX_ARRAY_FAMILY_FLEX_TYPEDEFS(std::complex<double>, complex_double)

#undef SCITBX_ARRAY_FAMILY_FLEX_TYPEDEFS

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_FLEX_TYPES_H
