/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MAPTBX_COPY_H
#define CCTBX_MAPTBX_COPY_H

#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/loops.h>
#include <cctbx/error.h>
#include <cctbx/import_scitbx_af.h>
#include <cassert>

namespace cctbx { namespace maptbx {

  template <typename FloatType,
            typename IndexType>
  af::versa<FloatType, af::flex_grid<IndexType> >
  copy(af::const_ref<FloatType, af::flex_grid<IndexType> > const& map,
       af::flex_grid<IndexType> const& result_grid)
  {
    CCTBX_ASSERT(map.accessor().origin().all_eq(result_grid.origin()));
    CCTBX_ASSERT(map.accessor().focus().all_eq(result_grid.focus()));
    typedef af::flex_grid<IndexType> f_g;
    typedef af::versa<FloatType, f_g> result_type;
    f_g m0 = map.accessor().shift_origin();
    f_g r0 = result_grid.shift_origin();
    assert(m0.focus().all_eq(r0.focus()));
    if (!m0.is_padded() && !r0.is_padded()) {
      result_type result;
      result.as_base_array().assign(map.begin(), map.end());
      result.resize(result_grid);
      return result;
    }
    af::nested_loop<IndexType> loop(m0.focus());
    result_type result(result_grid);
    FloatType* r_begin = result.begin();
    for (IndexType const& pt = loop(); !loop.over(); loop.incr()) {
      r_begin[r0(pt)] = map[m0(pt)];
    }
    return result;
  }

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_COPY_H
