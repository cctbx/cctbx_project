#ifndef SCITBX_ARRAY_FAMILY_ACCESSORS_F_GRID_H
#define SCITBX_ARRAY_FAMILY_ACCESSORS_F_GRID_H

#include <scitbx/array_family/accessors/c_grid.h>

namespace scitbx { namespace af {

  template <std::size_t Nd, typename IndexValueType=std::size_t>
  class f_grid
  {};

  template <typename IndexValueType>
  class f_grid<2, IndexValueType> : public c_grid<2, IndexValueType>
  {
    private:
      typedef c_grid<2, IndexValueType> base_t;

    public:
      typedef typename base_t::index_value_type index_value_type;
      typedef typename base_t::value_type value_type;
      typedef typename base_t::index_type index_type;

      f_grid() : base_t() {}

      f_grid(index_type const& n) : base_t(n) {}

      f_grid(index_value_type const& n0, index_value_type const& n1)
        : base_t(n0, n1)
      {}

      template <typename FlexIndexType>
      f_grid(flex_grid<FlexIndexType> const& flex_g)
        : base_t(flex_g)
      {}

      index_type
      index_nd(index_value_type const& i_1d)
      {
        return index_type(i_1d % this->elems[0], i_1d / this->elems[0]);
      }

      std::size_t
      operator()(index_type const& i) const
      {
        return i[1] * this->elems[0] + i[0];
      }

      std::size_t
      operator()(index_value_type const& i0,
                 index_value_type const& i1) const
      {
        return i1 * this->elems[0] + i0;
      }
  };

}} // scitbx::af

#endif // GUARD
