#ifndef CCTBX_MAPTBX_COPY_H
#define CCTBX_MAPTBX_COPY_H

#include <cctbx/maptbx/accessors/c_grid_padded_p1.h>
#include <scitbx/array_family/accessors/c_grid_padded.h>
#include <scitbx/array_family/accessors/c_grid.h>
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

  template <typename ElementType>
  void
  copy(
    af::const_ref<ElementType, af::c_grid_padded<3> > const& source,
    af::ref<ElementType, af::c_grid<3> > const& target)
  {
    CCTBX_ASSERT(target.accessor().all_eq(source.accessor().focus()));
    typename af::c_grid_padded<3>::index_type n = target.accessor();
    typename af::c_grid_padded<3>::index_type i;
    std::size_t j = 0;
    for(i[0]=0;i[0]<n[0];i[0]++)
    for(i[1]=0;i[1]<n[1];i[1]++)
    for(i[2]=0;i[2]<n[2];i[2]++, j++) {
      target[j] = source(i);
    }
  }

  template <typename FloatType>
  af::versa<FloatType, af::flex_grid<> >
  copy(
    af::const_ref<FloatType, c_grid_padded_p1<3> > const& map_unit_cell,
    af::int3 const& first,
    af::int3 const& last)
  {
    CCTBX_ASSERT(first.all_le(last));
    af::flex_grid_default_index_type first_(af::adapt(first));
    af::flex_grid_default_index_type last_(af::adapt(last));
    af::versa<FloatType, af::flex_grid<> > result(
      af::flex_grid<>(first_, last_, false));
    FloatType* out_ptr = result.begin();
    typename c_grid_padded_p1<3>::index_type out_pt;
    for (out_pt[0] = first[0]; out_pt[0] <= last[0]; out_pt[0]++) {
    for (out_pt[1] = first[1]; out_pt[1] <= last[1]; out_pt[1]++) {
    for (out_pt[2] = first[2]; out_pt[2] <= last[2]; out_pt[2]++) {
      *out_ptr++ = static_cast<FloatType>(map_unit_cell(out_pt));
    }}}
    return result;
  }

  template <typename ElementType>
  void
  unpad_in_place(
    af::versa<ElementType, af::flex_grid<> >& map)
  {
    CCTBX_ASSERT(map.accessor().nd() == 3);
    CCTBX_ASSERT(map.accessor().is_0_based());
    CCTBX_ASSERT(map.accessor().focus()[0] == map.accessor().all()[0]);
    CCTBX_ASSERT(map.accessor().focus()[1] == map.accessor().all()[1]);
    CCTBX_ASSERT(map.accessor().focus()[2] <= map.accessor().all()[2]);
    typedef typename af::flex_grid<>::index_value_type i_v_t;
    const i_v_t n0 = map.accessor().focus()[0];
    const i_v_t n1 = map.accessor().focus()[1];
    const i_v_t n2 = map.accessor().focus()[2];
    const i_v_t d2 = map.accessor().all()[2] - n2;
    if (d2 == 0) return;
    ElementType* target = map.begin() + n2;
    const ElementType* source = target + d2;
    for(i_v_t i01=1;i01<n0*n1;i01++) {
      for(i_v_t i2=0;i2<n2;i2++) *target++ = *source++;
      source += d2;
    }
    map = af::versa<ElementType, af::flex_grid<> >(
      map, af::flex_grid<>(n0,n1,n2));
  }

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_COPY_H
