// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_BPL_UTILS_H
#define CCTBX_BPL_UTILS_H

#include <boost/python/class_builder.hpp>
#include <cctbx/array_family/flex_grid_accessor.h>
#include <cctbx/array_family/versa.h>

namespace cctbx { namespace bpl_utils {

  void raise_index_error();
  void raise_must_be_1d();
  void raise_shared_size_mismatch();
  void raise_incompatible_arrays();

  boost::python::tuple tuple_from_python_list_or_tuple(PyObject* p);

  void assert_1d(af::flex_grid<> const& grid);

  template <typename ElementType>
  af::shared_plain<ElementType>
  as_base_array(af::versa<ElementType, af::flex_grid<> > const& a)
  {
    assert_1d(a.accessor());
    af::shared_plain<ElementType> b = a.as_base_array();
    if (a.size() != b.size()) raise_shared_size_mismatch();
    return b;
  }

}}

#endif // CCTBX_BPL_UTILS_H
