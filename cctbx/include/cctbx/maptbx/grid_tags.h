/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Moved from cctbx/maps/sym_tags.h to cctbx/maptbx (rwgk)
     2002 Jan: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MAPTBX_GRID_TAGS_H
#define CCTBX_MAPTBX_GRID_TAGS_H

#include <cctbx/maptbx/symmetry_flags.h>
#include <cctbx/maptbx/accessors/c_grid_p1.h>
#include <cctbx/maptbx/accessors/c_grid_padded_p1.h>
#include <cctbx/sgtbx/seminvariant.h>
#include <cctbx/tagged_value.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/array_family/accessors/c_grid_padded.h>
#include <scitbx/array_family/loops.h>
#include <scitbx/math/linear_regression.h>

namespace cctbx { namespace maptbx {

namespace grid_tags_detail {

  // Work-around for Visual C++ 6.0 and 7.0 overload resolution bug.
  struct space_group_symmetry_tag {};
  struct seminvariant_tag {};

  template <typename DimTupleType,
            typename IndexTupleType>
  inline tagged_value<IndexTupleType>
  multiply(DimTupleType const& n,
           sgtbx::rt_mx const& s,
           IndexTupleType const& x)
  {
    IndexTupleType result = s.r() * x;
    for(std::size_t i=0;i<3;i++) {
      result[i] = result[i] * s.t().den() + s.t()[i] * s.r().den() * n[i];
      if (result[i] % (s.r().den() * s.t().den())) {
        return tagged_value<IndexTupleType>(result, false);
      }
      result[i] /= (s.r().den() * s.t().den());
    }
    return tagged_value<IndexTupleType>(result, true);
  }

  template <typename DimensionType,
            typename GridSsType,
            typename IndexTupleType,
            typename FactorTupleType>
  inline tagged_value<IndexTupleType>
  add(DimensionType const& dim,
      GridSsType const& grid_ss,
      IndexTupleType const& pivot,
      FactorTupleType const& f)
  {
    IndexTupleType result = pivot;
    for(std::size_t i_ss=0;i_ss<grid_ss.size();i_ss++) {
      for(std::size_t i=0;i<3;i++) {
        int s = dim[i] * grid_ss[i_ss].v[i] * f[i_ss];
        if (s % grid_ss[i_ss].m) {
          return tagged_value<IndexTupleType>(result, false);
        }
        result[i] += s / grid_ss[i_ss].m;
      }
    }
    return tagged_value<IndexTupleType>(result, true);
  }

  template <typename TagVersaType,
            typename IndexTupleType>
  std::size_t
  mark_orbit(TagVersaType& p1_tags,
             sgtbx::space_group const& space_group,
             IndexTupleType const& pivot,
             space_group_symmetry_tag)
  {
    std::size_t grid_misses = 0;
    std::size_t i1d_pivot = p1_tags.accessor()(pivot);
    for(std::size_t i_smx=1;i_smx<space_group.order_z();i_smx++) {
      sgtbx::rt_mx s = space_group(i_smx);
      tagged_value<IndexTupleType>
      sym_equiv_point = multiply(p1_tags.accessor(), s, pivot);
      if (sym_equiv_point.tag) {
        std::size_t i1d_sep = p1_tags.accessor()(sym_equiv_point.value);
        while (p1_tags[i1d_sep] != -1) i1d_sep = p1_tags[i1d_sep];
        if (i1d_sep != i1d_pivot) p1_tags[i1d_sep] = i1d_pivot;
      }
      else {
        grid_misses++;
      }
    }
    return grid_misses;
  }

  template <typename TagVersaType,
            typename GridSsType,
            typename IndexTupleType>
  std::size_t
  mark_orbit(TagVersaType& p1_tags,
             GridSsType const& grid_ss,
             IndexTupleType const& pivot,
             seminvariant_tag)
  {
    std::size_t grid_misses = 0;
    std::size_t i1d_pivot = p1_tags.accessor()(pivot);
    af::small<int, 3> moduli(grid_ss.size());
    for(std::size_t i=0;i<grid_ss.size();i++) moduli[i] = grid_ss[i].m;
    af::nested_loop<af::small<int, 3> > loop(moduli);
    for (af::small<int, 3> const& f = loop(); !loop.over(); loop.incr()) {
      tagged_value<IndexTupleType>
      sym_equiv_point = add(p1_tags.accessor(), grid_ss, pivot, f);
      if (sym_equiv_point.tag) {
        std::size_t i1d_sep = p1_tags.accessor()(sym_equiv_point.value);
        while (p1_tags[i1d_sep] != -1) i1d_sep = p1_tags[i1d_sep];
        if (i1d_sep != i1d_pivot) p1_tags[i1d_sep] = i1d_pivot;
      }
      else {
        grid_misses++;
      }
    }
    return grid_misses;
  }

  template <typename TagVersaType,
            typename SymmetryType,
            typename SymmetryTypeTag>
  std::size_t
  mark_orbits(TagVersaType& p1_tags,
              SymmetryType const& symmetry,
              SymmetryTypeTag)
  {
    std::size_t grid_misses = 0;
    af::nested_loop<af::int3> loop(p1_tags.accessor());
    for (af::int3 const& pivot = loop(); !loop.over(); loop.incr()) {
      if (p1_tags(pivot) == -1) {
        grid_misses += mark_orbit(p1_tags, symmetry, pivot, SymmetryTypeTag());
      }
    }
    return grid_misses;
  }

  template <typename TagArrayType>
  std::size_t
  optimize_tags(TagArrayType tags)
  {
    std::size_t n_independent = 0;
    for(std::size_t i=0;i<tags.size();i++) {
      if (tags[i] < 0) {
        n_independent++;
      }
      else {
        std::size_t j = tags[i];
        while (tags[j] >= 0) j = tags[j];
        tags[i] = j;
      }
    }
    return n_independent;
  }

  // templated on TagType b/o Visual C++ 7 limitation
  template <typename FloatType, typename TagType>
  struct verify_grid_tags : scitbx::math::linear_regression_core<FloatType>
  {
    verify_grid_tags() {}

    verify_grid_tags(
      af::const_ref<FloatType, af::c_grid_padded<3> > const& data,
      af::const_ref<TagType, c_grid_p1<3> > const& tags,
      FloatType const& epsilon=1.e-6);

    std::size_t n_dependent;
  };

  template <typename FloatType, typename TagType>
  verify_grid_tags<FloatType, TagType>
  ::verify_grid_tags(
      af::const_ref<FloatType, af::c_grid_padded<3> > const& data,
      af::const_ref<TagType, c_grid_p1<3> > const& tags,
      FloatType const& epsilon)
  :
    n_dependent(0)
  {
    CCTBX_ASSERT(tags.accessor().all_eq(data.accessor().focus()));
    typedef af::c_grid_padded<3>::index_type index_type;
    af::nested_loop<index_type> loop(data.accessor().focus());
    af::c_grid<3> tags_accessor(tags.accessor()); // index_type conversion
    FloatType x, y, min_x, max_x, min_y, max_y;
    FloatType sum_x, sum_x2, sum_y, sum_y2, sum_xy;
    n_dependent = 0;
    std::size_t i = 0;
    for (index_type const& pt = loop(); !loop.over(); loop.incr(), i++) {
      if (tags[i] < 0) continue;
      x = data(pt);
      y = data(tags_accessor.index_nd(tags[i]));
      if (n_dependent == 0) {
        min_x = x;
        max_x = x;
        min_y = y;
        max_y = y;
        sum_x = x;
        sum_x2 = x * x;
        sum_y = y;
        sum_y2 = y * y;
        sum_xy = x * y;
      }
      if (min_x > x) min_x = x;
      if (max_x < x) max_x = x;
      if (min_y > y) min_y = y;
      if (max_y < y) max_y = y;
      sum_x += x;
      sum_x2 += x * x;
      sum_y += y;
      sum_y2 += y * y;
      sum_xy += x * y;
      n_dependent++;
    }
    if (n_dependent == 0) {
      this->reset();
    }
    else {
      this->set(n_dependent,
                min_x, max_x, min_y, max_y,
                sum_x, sum_x2, sum_y, sum_y2, sum_xy,
                epsilon);
    }
  }

} // namespace grid_tags_detail

  template <typename TagType = long>
  class grid_tags
  {
    public:
      typedef af::versa<TagType, c_grid_p1<3> > tag_array_type;
      typedef typename tag_array_type::accessor_type tag_accessor_type;
      typedef typename tag_accessor_type::index_type grid_point_type;

      grid_tags() {}

      grid_tags(af::int3 const& dim)
      :
        is_valid_(false),
        tag_array_(dim)
      {}

      bool
      is_valid() const { return is_valid_; }

      af::versa<TagType, c_grid_p1<3> >
      tag_array() const { return tag_array_; }

      void
      build(sgtbx::space_group_type const& sg_type,
            maptbx::symmetry_flags const& sym_flags);

      sgtbx::space_group_type const&
      space_group_type() const { return sg_type_; }

      maptbx::symmetry_flags const&
      symmetry_flags() const { return sym_flags_; }

      af::small<sgtbx::ss_vec_mod, 3> const&
      grid_ss() const { return grid_ss_; }

      std::size_t
      n_grid_misses() const { return n_grid_misses_; }

      std::size_t
      n_independent() const { return n_independent_; }

      std::size_t
      n_dependent() const { return tag_array_.size() - n_independent_; }

      // FUTURE: move out of class body
      template <typename FloatType>
      bool
      verify(af::const_ref<FloatType, af::c_grid_padded<3> > const& data,
             double min_correlation=0.99) const
      {
        CCTBX_ASSERT(is_valid_);
        CCTBX_ASSERT(tag_array_.accessor().all_eq(data.accessor().focus()));
        using namespace grid_tags_detail;
        verify_grid_tags<FloatType, TagType> vfy(data, tag_array_.const_ref());
        CCTBX_ASSERT(vfy.n_dependent == n_dependent());
        if (vfy.is_well_defined() && vfy.cc() < min_correlation) return false;
        return true;
      }

      // FUTURE: move out of class body
      template <typename FloatType>
      void
      sum_sym_equiv_points(
        af::ref<FloatType, c_grid_padded_p1<3> > const& data) const
      {
        CCTBX_ASSERT(is_valid_);
        CCTBX_ASSERT(tag_array_.accessor().all_eq(data.accessor().focus()));
        using namespace grid_tags_detail;
        using grid_tags_detail::multiply; // gcc 2.96 work-around
        {
          // 1. pass: accumulate contributions for sym. equiv. grid points.
          sgtbx::space_group const& space_group = sg_type_.group();
          af::nested_loop<grid_point_type> loop(tag_array_.accessor());
          for (grid_point_type const& pt = loop(); !loop.over(); loop.incr()) {
            if (tag_array_(pt) >= 0) continue;
            std::size_t i_pt = data.accessor()(pt);
            FloatType sum = data[i_pt];
            for(std::size_t i_smx=1;i_smx<space_group.order_z();i_smx++) {
              sgtbx::rt_mx m = space_group(i_smx);
              tagged_value<grid_point_type>
              sym_equiv_point = multiply(tag_array_.accessor(), m, pt);
              CCTBX_ASSERT(sym_equiv_point.tag);
              sum += data(sym_equiv_point.value);
            }
            data[i_pt] = sum;
          }
        }
        {
          // 2. pass: copy density at asym. grid point to sym. equiv. points.
          af::nested_loop<grid_point_type> loop(tag_array_.accessor());
          for (grid_point_type const& pt = loop(); !loop.over(); loop.incr()) {
            if (tag_array_(pt) < 0) continue;
            grid_point_type
              asym_pt = tag_array_.accessor().index_nd(tag_array_(pt));
            data(pt) = data(asym_pt);
          }
        }
      }

    protected:
      bool is_valid_;
      tag_array_type tag_array_;
      sgtbx::space_group_type sg_type_;
      maptbx::symmetry_flags sym_flags_;
      af::small<sgtbx::ss_vec_mod, 3> grid_ss_;
      std::size_t n_grid_misses_;
      std::size_t n_independent_;
  };

  template <typename TagType>
  void
  grid_tags<TagType>
  ::build(sgtbx::space_group_type const& sg_type,
          maptbx::symmetry_flags const& sym_flags)
  {
    using namespace grid_tags_detail;
    using grid_tags_detail::optimize_tags; // gcc 2.96 work-around
    if (   is_valid_
        && sg_type_.group() == sg_type.group()
        && sym_flags_ == sym_flags) {
      return;
    }
    sg_type_ = sg_type;
    sym_flags_ = sym_flags;
    n_grid_misses_ = 0;
    tag_array_.fill(-1);
    sgtbx::space_group sym = sym_flags.select_sub_space_group(sg_type_);
    if (mark_orbits(tag_array_, sym, space_group_symmetry_tag()) > 0) {
      throw error("Grid is not compatible with symmetry.");
    }
    if (sym_flags.use_structure_seminvariants()) {
      sgtbx::structure_seminvariant ss(sg_type.group());
      grid_ss_ = ss.grid_adapted_moduli(tag_array_.accessor());
      n_grid_misses_ = mark_orbits(tag_array_, grid_ss_, seminvariant_tag());
    }
    n_independent_ = optimize_tags(tag_array_.as_1d().ref());
    is_valid_ = true;
  }

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_GRID_TAGS_H
