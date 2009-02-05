// DEPRECATED FILE!
// replaced by scitbx/array_family/accessors/c_grid_padded_periodic.h

#ifndef CCTBX_MAPTBX_ACCESSORS_C_GRID_PADDED_P1_H
#define CCTBX_MAPTBX_ACCESSORS_C_GRID_PADDED_P1_H

#include <scitbx/array_family/accessors/flex_grid.h>
#include <cctbx/import_scitbx_af.h>
#include <scitbx/math/modulo.h>

namespace cctbx { namespace maptbx {

  template <std::size_t Nd>
  class c_grid_padded_p1;

  template <>
  class c_grid_padded_p1<3>
  {
    public:
      typedef af::tiny<int, 3> index_type;
      typedef index_type::value_type index_value_type;

      c_grid_padded_p1()
      : all_(0,0,0),
        focus_(0,0,0)
      {}

      c_grid_padded_p1(index_type const& all,
                       index_type const& focus)
      : all_(all),
        focus_(focus)
      {}

      template <typename FlexIndexType>
      c_grid_padded_p1(af::flex_grid<FlexIndexType> const& flex_g)
      :
        all_(af::adapt(flex_g.all()))
      {
        SCITBX_ASSERT(flex_g.is_0_based());
        focus_ = index_type(af::adapt(flex_g.focus()));
      }

      af::flex_grid<>
      as_flex_grid() const
      {
        return af::flex_grid<>(af::adapt(all_))
                   .set_focus(af::adapt(focus_));
      }

      std::size_t
      size_1d() const
      {
        return all_[0] * all_[1] * all_[2];
      }

      index_type const&
      all() const { return all_; }

      index_type const&
      focus() const { return focus_; }

      std::size_t
      focus_size_1d() const
      {
        return focus_[0] * focus_[1] * focus_[2];
      }

      bool is_padded() const
      {
        SCITBX_ASSERT(all_.all_ge(focus_));
        return !all_.all_eq(focus_);
      }

      std::size_t operator()(index_type const& i) const
      {
        return (scitbx::math::mod_positive(i[0], focus_[0]) * all_[1]
              + scitbx::math::mod_positive(i[1], focus_[1])) * all_[2]
              + scitbx::math::mod_positive(i[2], focus_[2]);
      }

      std::size_t
      operator()(index_value_type const& i0,
                 index_value_type const& i1,
                 index_value_type const& i2) const
      {
        return (scitbx::math::mod_positive(i0, focus_[0]) * all_[1]
              + scitbx::math::mod_positive(i1, focus_[1])) * all_[2]
              + scitbx::math::mod_positive(i2, focus_[2]);
      }

    protected:
      index_type all_;
      index_type focus_;
  };

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_ACCESSORS_C_GRID_PADDED_P1_H
