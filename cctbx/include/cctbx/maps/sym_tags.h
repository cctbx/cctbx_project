// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MAPS_SYM_TAGS_H
#define CCTBX_MAPS_SYM_TAGS_H

#include <cctbx/array_family/tiny_types.h>
#include <cctbx/array_family/small.h>
#include <cctbx/array_family/versa.h>
#include <cctbx/array_family/loops.h>
#include <cctbx/math/linear_regression.h>
#include <cctbx/sgtbx/groups.h>
#include <cctbx/sgtbx/seminvariant.h>
#include <cctbx/maps/accessors.h>

namespace cctbx { namespace maps {

  // Work-around for Visual C++ 6.0 and 7.0 overload resolution bug.
  namespace detail {
    struct space_group_symmetry_tag {};
    struct seminvariant_tag {};
  }

  template <typename ValueType, typename TagType = bool>
  struct tagged_value
  {
    tagged_value() {}
    tagged_value(const ValueType& v)
      : value(v)
    {}
    tagged_value(const ValueType& v, const TagType& t)
      : value(v), tag(t)
    {}
    ValueType value;
    TagType tag;
  };

  template <typename DimTupleType,
            typename IndexTupleType>
  inline tagged_value<IndexTupleType>
  multiply(const DimTupleType& N,
           const sgtbx::RTMx& M,
           const IndexTupleType& X)
  {
    IndexTupleType result = M.Rpart() * X;
    for(int i=0;i<3;i++) {
      result[i] = result[i] * M.TBF() + M.Tpart()[i] * M.RBF() * N[i];
      if (result[i] % (M.RBF() * M.TBF())) {
        return tagged_value<IndexTupleType>(result, false);
      }
      result[i] /= (M.RBF() * M.TBF());
    }
    return tagged_value<IndexTupleType>(result, true);
  }

  template <typename DimensionType,
            typename GridSsType,
            typename IndexTupleType,
            typename FactorTupleType>
  inline tagged_value<IndexTupleType>
  add(const DimensionType& dim,
      const GridSsType& grid_ss,
      const IndexTupleType& pivot,
      const FactorTupleType& f)
  {
    IndexTupleType result = pivot;
    for(int i_ss=0;i_ss<grid_ss.size();i_ss++) {
      for(int i=0;i<3;i++) {
        int s = dim[i] * grid_ss[i_ss].V[i] * f[i_ss];
        if (s % grid_ss[i_ss].M) {
          return tagged_value<IndexTupleType>(result, false);
        }
        result[i] += s / grid_ss[i_ss].M;
      }
    }
    return tagged_value<IndexTupleType>(result, true);
  }

  template <typename TagVersaType,
            typename IndexTupleType>
  std::size_t
  mark_orbit(TagVersaType& p1_tags,
             const sgtbx::SpaceGroup& SgOps,
             const IndexTupleType& pivot,
             detail::space_group_symmetry_tag)
  {
    std::size_t grid_misses = 0;
    std::size_t i1d_pivot = p1_tags.accessor()(pivot);
    for(int iSMx=1;iSMx<SgOps.OrderZ();iSMx++) {
      sgtbx::RTMx M = SgOps(iSMx);
      tagged_value<IndexTupleType>
      sym_equiv_point = multiply(p1_tags.accessor(), M, pivot);
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
             const GridSsType& grid_ss,
             const IndexTupleType& pivot,
             detail::seminvariant_tag)
  {
    std::size_t grid_misses = 0;
    std::size_t i1d_pivot = p1_tags.accessor()(pivot);
    af::small<int, 3> moduli(grid_ss.size());
    for(int i=0;i<grid_ss.size();i++) moduli[i] = grid_ss[i].M;
    af::nested_loop<af::small<int, 3> > loop(moduli);
    for (const af::small<int, 3>& f = loop(); !loop.over(); loop.incr()) {
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
              const SymmetryType& symmetry,
              SymmetryTypeTag)
  {
    std::size_t grid_misses = 0;
    af::nested_loop<af::int3> loop(p1_tags.accessor());
    for (const af::int3& pivot = loop(); !loop.over(); loop.incr()) {
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

  class map_symmetry_flags
  {
    public:
      map_symmetry_flags() {}
      explicit
      map_symmetry_flags(bool use_space_group_symmetry,
                         bool use_normalizer_K2L = false,
                         bool use_structure_seminvariants = false)
        : use_space_group_symmetry_(use_space_group_symmetry),
          use_normalizer_K2L_(use_normalizer_K2L),
          use_structure_seminvariants_(use_structure_seminvariants)
      {}

      bool use_space_group_symmetry() const {
        return use_space_group_symmetry_;
      }

      bool use_normalizer_K2L() const {
        return use_normalizer_K2L_;
      }

      bool use_structure_seminvariants() const {
        return use_structure_seminvariants_;
      }

      sgtbx::SpaceGroup
      select_sub_space_group(const sgtbx::SpaceGroupInfo& SgInfo) const
      {
        sgtbx::SpaceGroup result;
        if (use_space_group_symmetry()) {
          result = SgInfo.SgOps();
        }
        else if (use_structure_seminvariants()) {
          // XXX should be method of SpaceGroup
          for(int i=1;i<SgInfo.SgOps().nLTr();i++) {
            result.expandLTr(SgInfo.SgOps()(i,0,0).Tpart());
          }
        }
        if (use_normalizer_K2L()) {
          result.expandSMxArray(
            SgInfo.getAddlGeneratorsOfEuclideanNormalizer(1, 0));
        }
        return result;
      }

      af::int3
      get_grid_factors(const sgtbx::SpaceGroupInfo& SgInfo) const
      {
        af::int3 grid_ss(1,1,1);
        if (use_structure_seminvariants()) {
          grid_ss = sgtbx::StructureSeminvariant(
            SgInfo.SgOps()).refine_gridding();
        }
        return select_sub_space_group(SgInfo).refine_gridding(grid_ss);
      }

      bool operator==(const map_symmetry_flags& rhs) const
      {
        return
             use_space_group_symmetry_ == rhs.use_space_group_symmetry_
          && use_normalizer_K2L_ == rhs.use_normalizer_K2L_
          && use_structure_seminvariants_ == rhs.use_structure_seminvariants_;
      }

    private:
      bool use_space_group_symmetry_;
      bool use_normalizer_K2L_;
      bool use_structure_seminvariants_;
  };

  template <class DataArrayType>
  struct verify_grid_tags
    : math::linear_regression_core<typename DataArrayType::size_type,
                                   typename DataArrayType::value_type>
  {
    typedef typename DataArrayType::size_type size_type;
    typedef typename DataArrayType::value_type value_type;

    verify_grid_tags() {}
    template <typename TagArrayType>
    verify_grid_tags(const DataArrayType& data,
                     const TagArrayType& tags,
                     const value_type& epsilon = 1.e-6)
    {
      cctbx_assert(data.size() <= tags.size());
      n_dependent = 0;
      size_type i;
      for(i=0;i<data.size();i++) {
        if (tags[i] >= 0) break;
      }
      if (i == data.size()) {
        reset();
        return;
      }
      n_dependent = 1;
      value_type x = data[i];
      value_type y = data[tags[i]];
      value_type min_x = x;
      value_type max_x = x;
      value_type min_y = y;
      value_type max_y = y;
      value_type sum_x = x;
      value_type sum_x2 = x * x;
      value_type sum_y = y;
      value_type sum_y2 = y * y;
      value_type sum_xy = x * y;
      for(i++;i<data.size();i++) {
        if (tags[i] < 0) continue;
        n_dependent++;
        x = data[i];
        y = data[tags[i]];
        if (min_x > x) min_x = x;
        if (max_x < x) max_x = x;
        if (min_y > y) min_y = y;
        if (max_y < y) max_y = y;
        sum_x += x;
        sum_x2 += x * x;
        sum_y += y;
        sum_y2 += y * y;
        sum_xy += x * y;
      }
      set(n_dependent,
          min_x, max_x, min_y, max_y, sum_x, sum_x2, sum_y, sum_y2, sum_xy,
          epsilon);
    }
    size_type n_dependent;
  };

  template <typename TagType>
  class grid_tags : public af::versa<TagType, grid_p1<3> >
  {
    public:
      typedef af::versa<TagType, grid_p1<3> > base_type;
      typedef typename base_type::accessor_type accessor_type;
      typedef typename accessor_type::index_type grid_point_type;

      grid_tags() {}
      grid_tags(const af::grid<3>& dim)
        : base_type(dim), m_is_valid(false)
      {}

      void build(const sgtbx::SpaceGroupInfo& SgInfo,
                 const map_symmetry_flags& sym_flags)
      {
        if (   m_is_valid
            && m_SgInfo.SgOps() == SgInfo.SgOps()
            && m_sym_flags == sym_flags) {
          return;
        }
        m_SgInfo = SgInfo;
        m_sym_flags = sym_flags;
        m_n_grid_misses = 0;
        this->fill(-1);
        sgtbx::SpaceGroup
        sym = sym_flags.select_sub_space_group(m_SgInfo);
        if (mark_orbits(*this, sym, detail::space_group_symmetry_tag()) > 0) {
          throw error("Grid is not compatible with symmetry.");
        }
        if (sym_flags.use_structure_seminvariants()) {
          sgtbx::StructureSeminvariant ss(SgInfo.SgOps());
          m_grid_ss = ss.grid_adapted_moduli(this->accessor());
          m_n_grid_misses = mark_orbits(
            *this, m_grid_ss, detail::seminvariant_tag());
        }
        m_n_independent = optimize_tags(this->as_1d().ref());
        m_is_valid = true;
      }

      bool is_valid() const { return m_is_valid; }
      const sgtbx::SpaceGroupInfo& SgInfo() const { return m_SgInfo; }
      const map_symmetry_flags& sym_flags() const { return m_sym_flags; }
      const af::small<sgtbx::ssVM, 3>& grid_ss() const {
        return m_grid_ss;
      }
      std::size_t n_grid_misses() const { return m_n_grid_misses; }
      std::size_t n_independent() const { return m_n_independent; }

      template <typename ArrayType>
      bool
      verify(const ArrayType& data, double min_correlation = 0.99) const
      {
        verify_grid_tags<ArrayType> vfy(data, this->const_ref().as_1d());
        cctbx_assert(vfy.n_dependent + m_n_independent == this->size());
        if (vfy.is_well_defined() && vfy.cc() < min_correlation) return false;
        return true;
      }

      template <typename FloatType>
      void
      sum_sym_equiv_points(af::ref<FloatType, padded_grid_p1<3> > map) const
      {
        cctbx_assert(m_is_valid);
        cctbx_assert(this->accessor() == map.accessor().n_logical());
        {
          // 1. pass: accumulate contributions for sym. equiv. grid points.
          const sgtbx::SpaceGroup& sgops = m_SgInfo.SgOps();
          af::nested_loop<grid_point_type> loop(this->accessor());
          for (const grid_point_type& pt = loop(); !loop.over(); loop.incr()) {
            if ((*this)(pt) >= 0) continue;
            std::size_t i_pt = map.accessor()(pt);
            for(int ismx=1;ismx<sgops.OrderZ();ismx++) {
              sgtbx::RTMx m = sgops(ismx);
              tagged_value<grid_point_type>
              sym_equiv_point = multiply(this->accessor(), m, pt);
              cctbx_assert(sym_equiv_point.tag);
              map[i_pt] += map(sym_equiv_point.tag);
            }
          }
        }
        {
          // 2. pass: copy density at asym. grid point to sym. equiv. points.
          af::nested_loop<grid_point_type> loop(this->accessor());
          for (const grid_point_type& pt = loop(); !loop.over(); loop.incr()) {
            if ((*this)(pt) < 0) continue;
            grid_point_type
            asym_pt = this->accessor().i_1d_as_i_nd((*this)(pt));
            map(pt) = map(asym_pt);
          }
        }
      }

    protected:
      bool m_is_valid;
      sgtbx::SpaceGroupInfo m_SgInfo;
      map_symmetry_flags m_sym_flags;
      af::small<sgtbx::ssVM, 3> m_grid_ss;
      std::size_t m_n_grid_misses;
      std::size_t m_n_independent;
  };

}} // namespace cctbx::maps

#endif // CCTBX_MAPS_SYM_TAGS_H
