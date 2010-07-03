#ifndef CCTBX_MAPTBX_GRID_TAGS_H
#define CCTBX_MAPTBX_GRID_TAGS_H

#include <cctbx/maptbx/accessors/c_grid_p1.h>
#include <cctbx/maptbx/accessors/c_grid_padded_p1.h>
#include <cctbx/sgtbx/search_symmetry.h>
#include <cctbx/tagged_value.h>
#include <scitbx/rational.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/array_family/accessors/c_grid_padded.h>
#include <scitbx/array_family/loops.h>
#include <scitbx/math/linear_correlation.h>

namespace cctbx { namespace maptbx {

namespace grid_tags_detail {

  // Work-around for Visual C++ 6.0 and 7.0 overload resolution bug.
  struct space_group_symmetry_tag {};
  struct seminvariant_tag {};

  template <typename DimTupleType>
  DimTupleType
  factors_for_common_denominator(DimTupleType const& n)
  {
    typename DimTupleType::value_type lcm = scitbx::array_lcm(n);
    DimTupleType result;
    for(std::size_t i=0;i<n.size();i++) {
      CCTBX_ASSERT(n[i] != 0);
      result[i] = lcm / n[i];
    }
    return result;
  }

  template <typename DimTupleType,
            typename IndexTupleType>
  inline tagged_value<IndexTupleType>
  multiply(DimTupleType const& n,
           DimTupleType const& f,
           sgtbx::rt_mx const& s,
           IndexTupleType const& fx)
  {
    IndexTupleType result = s.r() * fx;
    for(std::size_t i=0;i<3;i++) {
      result[i] = result[i] * s.t().den()
                + s.t()[i] * s.r().den() * f[i] * n[i];
      if (result[i] % (s.r().den() * s.t().den() * f[i])) {
        return tagged_value<IndexTupleType>(result, false);
      }
      result[i] /= (s.r().den() * s.t().den() * f[i]);
    }
    return tagged_value<IndexTupleType>(result, true);
  }

  template <typename DimensionType,
            typename GridSsType,
            typename IndexTupleType,
            typename FactorTupleType>
  inline tagged_value<IndexTupleType>
  add(DimensionType const& dim,
      GridSsType const& grid_ss_continuous,
      IndexTupleType const& pivot,
      FactorTupleType const& f)
  {
    IndexTupleType result = pivot;
    for(std::size_t i_ss=0;i_ss<grid_ss_continuous.size();i_ss++) {
      for(std::size_t i=0;i<3;i++) {
        int s = dim[i] * grid_ss_continuous[i_ss].v[i] * f[i_ss];
        if (s % grid_ss_continuous[i_ss].m) {
          return tagged_value<IndexTupleType>(result, false);
        }
        result[i] += s / grid_ss_continuous[i_ss].m;
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
    if (space_group.order_z() > 1) {
      typedef typename TagVersaType::accessor_type dim_tuple_type;
      dim_tuple_type f = factors_for_common_denominator(p1_tags.accessor());
      IndexTupleType f_pivot;
      for(std::size_t i=0;i<3;i++) f_pivot[i] = f[i] * pivot[i];
      std::size_t i1d_pivot = p1_tags.accessor()(pivot);
      for(std::size_t i_smx=1;i_smx<space_group.order_z();i_smx++) {
        sgtbx::rt_mx s = space_group(i_smx);
        tagged_value<IndexTupleType>
        sym_equiv_point = multiply(p1_tags.accessor(), f, s, f_pivot);
        if (sym_equiv_point.tag) {
          std::size_t i1d_sep = p1_tags.accessor()(sym_equiv_point.value);
          while (p1_tags[i1d_sep] != -1) i1d_sep = p1_tags[i1d_sep];
          if (i1d_sep != i1d_pivot) p1_tags[i1d_sep] = i1d_pivot;
        }
        else {
          grid_misses++;
        }
      }
    }
    return grid_misses;
  }

  template <typename TagVersaType,
            typename GridSsType,
            typename IndexTupleType>
  std::size_t
  mark_orbit(TagVersaType& p1_tags,
             GridSsType const& grid_ss_continuous,
             IndexTupleType const& pivot,
             seminvariant_tag)
  {
    std::size_t grid_misses = 0;
    std::size_t i1d_pivot = p1_tags.accessor()(pivot);
    af::small<int, 3> moduli(grid_ss_continuous.size());
    for(std::size_t i=0;i<grid_ss_continuous.size();i++) {
      moduli[i] = grid_ss_continuous[i].m;
    }
    af::nested_loop<af::small<int, 3> > loop(moduli);
    for (af::small<int, 3> const& f = loop(); !loop.over(); loop.incr()) {
      tagged_value<IndexTupleType>
      sym_equiv_point = add(p1_tags.accessor(), grid_ss_continuous, pivot, f);
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

  template <typename FloatType, typename TagType>
  scitbx::math::linear_correlation<>
  dependent_correlation(
      std::size_t n_dependent,
      af::const_ref<FloatType, af::c_grid_padded<3> > const& data,
      af::const_ref<TagType, c_grid_p1<3> > const& tags,
      double epsilon)
  {
    CCTBX_ASSERT(data.accessor().focus().all_eq(tags.accessor()));
    typedef af::c_grid_padded<3>::index_type index_type;
    af::nested_loop<index_type> loop(data.accessor().focus());
    af::c_grid<3> tags_accessor(tags.accessor()); // index_type conversion
    std::vector<FloatType> x; x.reserve(n_dependent);
    std::vector<FloatType> y; y.reserve(n_dependent);
    std::size_t i = 0;
    for (index_type const& pt = loop(); !loop.over(); loop.incr(), i++) {
      if (tags[i] < 0) continue;
      x.push_back(data(pt));
      y.push_back(data(tags_accessor.index_nd(tags[i])));
    }
    CCTBX_ASSERT(x.size() == n_dependent);
    return scitbx::math::linear_correlation<>(
      af::make_const_ref(x),
      af::make_const_ref(y),
      epsilon);
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
            sgtbx::search_symmetry_flags const& symmetry_flags);

      sgtbx::space_group_type const&
      space_group_type() const { return sg_type_; }

      sgtbx::search_symmetry_flags const&
      symmetry_flags() const { return symmetry_flags_; }

      af::small<sgtbx::ss_vec_mod, 3> const&
      grid_ss_continuous() const { return grid_ss_continuous_; }

      std::size_t
      n_grid_misses() const { return n_grid_misses_; }

      std::size_t
      n_independent() const { return n_independent_; }

      std::size_t
      n_dependent() const { return tag_array_.size() - n_independent_; }

      template <typename FloatType>
      scitbx::math::linear_correlation<>
      dependent_correlation(
        af::const_ref<FloatType, af::c_grid_padded<3> > const& data,
        double epsilon=1e-15) const
      {
        CCTBX_ASSERT(is_valid_);
        CCTBX_ASSERT(data.accessor().focus().all_eq(tag_array_.accessor()));
        return grid_tags_detail::dependent_correlation(
          n_dependent(),
          data,
          tag_array_.const_ref(),
          epsilon);
      }

      template <typename FloatType>
      bool
      verify(
        af::const_ref<FloatType, af::c_grid_padded<3> > const& data,
        double min_correlation=0.99) const
      {
        if (n_dependent() == 0) return true;
        return dependent_correlation(data).coefficient() >= min_correlation;
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
        {
          // 1. pass: accumulate contributions for sym. equiv. grid points.
          tag_accessor_type
            f = grid_tags_detail::factors_for_common_denominator(
              tag_array_.accessor());
          sgtbx::space_group const& space_group = sg_type_.group();
          af::nested_loop<grid_point_type> loop(tag_array_.accessor());
          for (grid_point_type const& pt = loop(); !loop.over(); loop.incr()) {
            if (tag_array_(pt) >= 0) continue;
            std::size_t i_pt = data.accessor()(pt);
            FloatType sum = data[i_pt];
            if (space_group.order_z() > 1) {
              grid_point_type f_pt;
              for(std::size_t i=0;i<3;i++) f_pt[i] = f[i] * pt[i];
              for(std::size_t i_smx=1;i_smx<space_group.order_z();i_smx++) {
                sgtbx::rt_mx m = space_group(i_smx);
                tagged_value<grid_point_type>
                sym_equiv_point = grid_tags_detail::multiply(
                  tag_array_.accessor(), f, m, f_pt);
                CCTBX_ASSERT(sym_equiv_point.tag);
                sum += data(sym_equiv_point.value);
              }
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

      template <typename DataType>
      std::size_t
      apply_symmetry_to_mask(
        af::ref<DataType, af::c_grid<3> > const& data) const
      {
        CCTBX_ASSERT(data.accessor().all_eq(tag_array_.accessor()));
        std::size_t n_overlap = 0;
        const TagType* tags=tag_array_.begin();
        for(std::size_t i=0;i<data.size();i++) {
          if (tags[i] < 0) continue;
          if (data[i] == 0) {
            TagType j = tags[i];
            if (data[j] == 0) {
              n_overlap++;
            }
            else {
              data[j] = 0;
            }
          }
        }
        for(std::size_t i=0;i<data.size();i++) {
          if (tags[i] < 0) continue;
          data[i] = data[tags[i]];
        }
        return n_overlap;
      }

    protected:
      bool is_valid_;
      tag_array_type tag_array_;
      sgtbx::space_group_type sg_type_;
      sgtbx::search_symmetry_flags symmetry_flags_;
      af::small<sgtbx::ss_vec_mod, 3> grid_ss_continuous_;
      std::size_t n_grid_misses_;
      std::size_t n_independent_;
  };

  template <typename TagType>
  void
  grid_tags<TagType>
  ::build(sgtbx::space_group_type const& sg_type,
          sgtbx::search_symmetry_flags const& symmetry_flags)
  {
    using namespace grid_tags_detail;
    if (   is_valid_
        && sg_type_.group() == sg_type.group()
        && symmetry_flags_ == symmetry_flags) {
      return;
    }
    sg_type_ = sg_type;
    symmetry_flags_ = symmetry_flags;
    n_grid_misses_ = 0;
    tag_array_.fill(-1);
    sgtbx::structure_seminvariants ss;
    sgtbx::space_group sym;
    if (symmetry_flags.use_seminvariants()) {
      ss = sgtbx::structure_seminvariants(sg_type.group());
      sym = sgtbx::search_symmetry(symmetry_flags_, sg_type_, ss).subgroup();
    }
    else {
      sym = sgtbx::search_symmetry(symmetry_flags_, sg_type_).subgroup();
    }
    if (mark_orbits(tag_array_, sym, space_group_symmetry_tag()) > 0) {
      throw error("Grid is not compatible with symmetry.");
    }
    if (symmetry_flags.use_seminvariants()) {
      grid_ss_continuous_ = ss.select(false).grid_adapted_moduli(
        tag_array_.accessor());
      n_grid_misses_ = mark_orbits(
        tag_array_, grid_ss_continuous_, seminvariant_tag());
    }
    n_independent_ = grid_tags_detail::optimize_tags(tag_array_.as_1d().ref());
    is_valid_ = true;
  }

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_GRID_TAGS_H
