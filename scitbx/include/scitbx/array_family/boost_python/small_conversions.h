/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_SMALL_CONVERSIONS_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_SMALL_CONVERSIONS_H

#include <scitbx/array_family/small.h>
#include <scitbx/boost_python/container_conversions.h>

namespace scitbx { namespace af { namespace boost_python {

  template <typename ElementType, std::size_t N>
  struct tuple_mapping_small
  {
    tuple_mapping_small()
    {
      scitbx::boost_python::container_conversions
        ::tuple_mapping_fixed_capacity<small<ElementType, N> >();
    }
  };

  void register_flex_grid_default_index_type_conversions();

}}} // namespace scitbx::af::boost_python

#endif // SCITBX_ARRAY_FAMILY_BOOST_PYTHON_SMALL_CONVERSIONS_H
