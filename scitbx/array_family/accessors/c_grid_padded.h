// The implementations in this file are highly redundant in order
// to help the optimizer.

#ifndef SCITBX_ARRAY_FAMILY_ACCESSORS_C_GRID_PADDED_H
#define SCITBX_ARRAY_FAMILY_ACCESSORS_C_GRID_PADDED_H

#include <scitbx/array_family/accessors/c_index_1d_calculator.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/math/compile_time.h>

namespace scitbx { namespace af {

  template <std::size_t Nd, typename IndexValueType=std::size_t>
  class c_grid_padded
  {
    public:
      typedef tiny<IndexValueType, Nd> index_type;
      typedef typename index_type::value_type index_value_type;

      c_grid_padded()
      {
        for(std::size_t i=0;i<Nd;i++) all_[i] = 0;
        for(std::size_t i=0;i<Nd;i++) focus_[i] = 0;
      }

      c_grid_padded(index_type const& all)
      :
        all_(all),
        focus_(all)
      {}

      c_grid_padded(index_type const& all,
                    index_type const& focus)
      :
        all_(all),
        focus_(focus)
      {}

      template <typename FlexIndexType>
      c_grid_padded(af::flex_grid<FlexIndexType> const& flex_g)
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

      af::flex_grid<>
      as_flex_grid() const
      {
        return af::flex_grid<>(af::adapt(all_))
                    .set_focus(af::adapt(focus_));
      }

      std::size_t
      size_1d() const
      {
        return math::compile_time::product<Nd>::get(all_);
      }

      index_type const&
      all() const { return all_; }

      index_type const&
      focus() const { return focus_; }

      std::size_t
      focus_size_1d() const
      {
        return math::compile_time::product<Nd>::get(focus_);
      }

      bool is_padded() const
      {
        SCITBX_ASSERT(all_.all_ge(focus_));
        return !all_.all_eq(focus_);
      }

      std::size_t
      operator()(index_type const& i) const
      {
        return c_index_1d_calculator<Nd>::get(all_, i);
      }

    protected:
      index_type all_;
      index_type focus_;
  };

  template <typename IndexValueType>
  class c_grid_padded<2, IndexValueType>
  {
    public:
      typedef IndexValueType index_value_type;
      typedef tiny<IndexValueType, 2> index_type;

      c_grid_padded() : all_(0,0), focus_(0,0) {}

      c_grid_padded(index_type const& all)
      :
        all_(all),
        focus_(all)
      {}

      c_grid_padded(index_value_type const& all_0,
                    index_value_type const& all_1)
      :
        all_(all_0, all_1),
        focus_(all_0, all_1)
      {}

      c_grid_padded(index_type const& all,
                    index_type const& focus)
      :
        all_(all),
        focus_(focus)
      {}

      c_grid_padded(index_value_type const& all_0,
                    index_value_type const& all_1,
                    index_value_type const& focus_0,
                    index_value_type const& focus_1)
      :
        all_(all_0, all_1),
        focus_(focus_0, focus_1)
      {}

      template <typename FlexIndexType>
      c_grid_padded(af::flex_grid<FlexIndexType> const& flex_g)
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

      af::flex_grid<>
      as_flex_grid() const
      {
        return af::flex_grid<>(af::adapt(all_))
                    .set_focus(af::adapt(focus_));
      }

      std::size_t
      size_1d() const
      {
        return static_cast<std::size_t>(all_[0])
             * static_cast<std::size_t>(all_[1]);
      }

      index_type const&
      all() const { return all_; }

      index_type const&
      focus() const { return focus_; }

      std::size_t
      focus_size_1d() const
      {
        return static_cast<std::size_t>(focus_[0])
             * static_cast<std::size_t>(focus_[1]);
      }

      bool is_padded() const
      {
        SCITBX_ASSERT(all_[0] >= focus_[0] && all_[1] >= focus_[1]);
        return all_[0] != focus_[0] || all_[1] != focus_[1];
      }

      std::size_t
      operator()(index_type const& i) const
      {
        return i[0] * all_[1] + i[1];
      }

      std::size_t
      operator()(index_value_type const& i0,
                 index_value_type const& i1) const
      {
        return i0 * all_[1] + i1;
      }

    protected:
      index_type all_;
      index_type focus_;
  };

  template <typename IndexValueType>
  class c_grid_padded<3, IndexValueType>
  {
    public:
      typedef IndexValueType index_value_type;
      typedef tiny<IndexValueType, 3> index_type;

      c_grid_padded() : all_(0,0,0), focus_(0,0,0) {}

      c_grid_padded(index_type const& all)
      :
        all_(all),
        focus_(all)
      {}

      c_grid_padded(index_value_type const& all_0,
                    index_value_type const& all_1,
                    index_value_type const& all_2)
      :
        all_(all_0, all_1, all_2),
        focus_(all_0, all_1, all_2)
      {}

      c_grid_padded(index_type const& all,
                    index_type const& focus)
      :
        all_(all),
        focus_(focus)
      {}

      c_grid_padded(index_value_type const& all_0,
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
      c_grid_padded(af::flex_grid<FlexIndexType> const& flex_g)
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

      af::flex_grid<>
      as_flex_grid() const
      {
        return af::flex_grid<>(af::adapt(all_))
                    .set_focus(af::adapt(focus_));
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
        SCITBX_ASSERT(   all_[0] >= focus_[0]
                      && all_[1] >= focus_[1]
                      && all_[2] >= focus_[2]);
        return    all_[0] != focus_[0]
               || all_[1] != focus_[1]
               || all_[2] != focus_[2];
      }

      std::size_t
      operator()(index_type const& i) const
      {
        return (i[0] * all_[1] + i[1]) * all_[2] + i[2];
      }

      std::size_t
      operator()(index_value_type const& i0,
                 index_value_type const& i1,
                 index_value_type const& i2) const
      {
        return (i0 * all_[1] + i1) * all_[2] + i2;
      }

    protected:
      index_type all_;
      index_type focus_;
  };

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_ACCESSORS_C_GRID_PADDED_H
