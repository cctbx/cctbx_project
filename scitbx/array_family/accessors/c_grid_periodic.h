#ifndef SCITBX_ARRAY_FAMILY_ACCESSORS_C_GRID_PERIODIC_H
#define SCITBX_ARRAY_FAMILY_ACCESSORS_C_GRID_PERIODIC_H

#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/math/modulo.h>

namespace scitbx { namespace af {

  template <std::size_t Nd>
  class c_grid_periodic;

  template <>
  class c_grid_periodic<3> : public tiny<int, 3>
  {
    public:
      typedef tiny<int, 3> index_type;
      typedef index_type::value_type index_value_type;
      typedef index_value_type value_type;

      c_grid_periodic() : index_type(0,0,0) {}

      c_grid_periodic(index_type const& n) : index_type(n) {}

      c_grid_periodic(index_value_type const& n0,
                      index_value_type const& n1,
                      index_value_type const& n2)
        : index_type(n0, n1, n2)
      {}

      template <typename FlexIndexType>
      c_grid_periodic(flex_grid<FlexIndexType> const& flex_g)
        : index_type(adapt(flex_g.all()))
      {
        SCITBX_ASSERT(flex_g.is_0_based());
      }

      flex_grid<>
      as_flex_grid() const
      {
        return flex_grid<>(adapt(*this));
      }

      std::size_t
      size_1d() const
      {
        return static_cast<std::size_t>(this->elems[0])
             * static_cast<std::size_t>(this->elems[1])
             * static_cast<std::size_t>(this->elems[2]);
      }

      index_type
      index_nd(index_value_type const& i_1d) const
      {
        index_type i_nd;
        i_nd[2] = i_1d % this->elems[2];
        i_nd[0] = i_1d / this->elems[2];
        i_nd[1] = i_nd[0] % this->elems[1];
        i_nd[0] /= this->elems[1];
        return i_nd;
      }

      std::size_t
      operator()(index_type const& i) const
      {
        const index_value_type* n = this->elems;
        return (math::mod_positive(i[0], n[0]) * n[1]
                + math::mod_positive(i[1], n[1])) * n[2]
                + math::mod_positive(i[2], n[2]);
      }

      std::size_t
      operator()(index_value_type const& i0,
                 index_value_type const& i1,
                 index_value_type const& i2) const
      {
        const index_value_type* n = this->elems;
        return (math::mod_positive(i0, n[0]) * n[1]
                + math::mod_positive(i1, n[1])) * n[2]
                + math::mod_positive(i2, n[2]);
      }
  };

}} // namespace scitbx::af

#endif // GUARD
