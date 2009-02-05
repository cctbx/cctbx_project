#ifndef SCITBX_ARRAY_FAMILY_ACCESSORS_C_GRID_PADDED_PERIODIC_H
#define SCITBX_ARRAY_FAMILY_ACCESSORS_C_GRID_PADDED_PERIODIC_H

#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/math/modulo.h>

namespace scitbx { namespace af {

  template <std::size_t Nd>
  class c_grid_padded_periodic;

  template <>
  class c_grid_padded_periodic<3>
  {
    public:
      typedef tiny<int, 3> index_type;
      typedef index_type::value_type index_value_type;

      c_grid_padded_periodic() : all_(0,0,0), focus_(0,0,0) {}

      c_grid_padded_periodic(index_type const& all)
      :
        all_(all),
        focus_(all)
      {}

      c_grid_padded_periodic(index_value_type const& all_0,
                    index_value_type const& all_1,
                    index_value_type const& all_2)
      :
        all_(all_0, all_1, all_2),
        focus_(all_0, all_1, all_2)
      {}

      c_grid_padded_periodic(index_type const& all,
                             index_type const& focus)
      :
        all_(all),
        focus_(focus)
      {}

      c_grid_padded_periodic(index_value_type const& all_0,
                             index_value_type const& all_1,
                             index_value_type const& all_2,
                             index_value_type const& focus_0,
                             index_value_type const& focus_1,
                             index_value_type const& focus_2)
      :
        all_(all_0, all_1, all_2),
        focus_(focus_0, focus_1, focus_2)
      {}

      template <typename FlexIndexType>
      c_grid_padded_periodic(af::flex_grid<FlexIndexType> const& flex_g)
      :
        all_(af::adapt(flex_g.all()))
      {
        SCITBX_ASSERT(flex_g.is_0_based());
        if (flex_g.is_padded()) {
          focus_ = index_type(af::adapt(flex_g.focus()));
        }
        else {
          focus_ = all_;
        }
      }

      flex_grid<>
      as_flex_grid() const
      {
        return flex_grid<>(adapt(all_)).set_focus(adapt(focus_));
      }

      std::size_t
      size_1d() const
      {
        return static_cast<std::size_t>(all_[0])
             * static_cast<std::size_t>(all_[1])
             * static_cast<std::size_t>(all_[2]);
      }

      index_type const&
      all() const { return all_; }

      index_type const&
      focus() const { return focus_; }

      std::size_t
      focus_size_1d() const
      {
        return static_cast<std::size_t>(focus_[0])
             * static_cast<std::size_t>(focus_[1])
             * static_cast<std::size_t>(focus_[2]);
      }

      bool is_padded() const
      {
        SCITBX_ASSERT(all_.all_ge(focus_));
        return !all_.all_eq(focus_);
      }

      std::size_t operator()(index_type const& i) const
      {
        return (math::mod_positive(i[0], focus_[0]) * all_[1]
                + math::mod_positive(i[1], focus_[1])) * all_[2]
                + math::mod_positive(i[2], focus_[2]);
      }

      std::size_t
      operator()(index_value_type const& i0,
                 index_value_type const& i1,
                 index_value_type const& i2) const
      {
        return (math::mod_positive(i0, focus_[0]) * all_[1]
                + math::mod_positive(i1, focus_[1])) * all_[2]
                + math::mod_positive(i2, focus_[2]);
      }

    protected:
      index_type all_;
      index_type focus_;
  };

}} // namespace scitbx::af

#endif // GUARD
