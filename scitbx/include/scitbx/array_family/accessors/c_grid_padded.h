/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Created (R.W. Grosse-Kunstleve)
 */

// The implementations in this file are highly redundant in order
// to help the optimizer.

#ifndef SCITBX_ARRAY_FAMILY_ACCESSORS_C_GRID_PADDED_H
#define SCITBX_ARRAY_FAMILY_ACCESSORS_C_GRID_PADDED_H

#include <scitbx/array_family/accessors/c_index_1d_calculator.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/math/compile_time.h>

namespace scitbx { namespace af {

  template <std::size_t Nd>
  class c_grid_padded
  {
    public:
      typedef tiny<std::size_t, Nd> index_type;
      typedef typename index_type::value_type index_value_type;

      c_grid_padded()
      {
        for(std::size_t i=0;i<Nd;i++) grid_[i] = 0;
        for(std::size_t i=0;i<Nd;i++) layout_[i] = 0;
      }

      c_grid_padded(index_type const& grid)
      :
        grid_(grid),
        layout_(grid)
      {}

      c_grid_padded(index_type const& grid,
                    index_type const& layout)
      :
        grid_(grid),
        layout_(layout)
      {}

      template <typename FlexIndexType>
      c_grid_padded(af::flex_grid<FlexIndexType> const& flex_g)
      :
        grid_(af::adapt(flex_g.grid()))
      {
        SCITBX_ASSERT(flex_g.is_0_based());
        if (flex_g.is_padded()) {
          layout_ = index_type(af::adapt(flex_g.layout()));
        }
        else {
          layout_ = grid_;
        }
      }

      af::flex_grid<>
      as_flex_grid() const
      {
        return af::flex_grid<>(af::adapt(grid_))
                   .set_layout(af::adapt(layout_));
      }

      std::size_t
      size_1d() const
      {
        return math::compile_time::product<Nd>::get(grid_);
      }

      index_type const&
      grid() const { return grid_; }

      index_type const&
      layout() const { return layout_; }

      std::size_t
      layout_size_1d() const
      {
        return math::compile_time::product<Nd>::get(layout_);
      }

      bool is_padded() const
      {
        SCITBX_ASSERT(grid_.all_ge(layout_));
        return !grid_.all_eq(layout_);
      }

      std::size_t
      operator()(index_type const& i) const
      {
        return c_index_1d_calculator<Nd>::get(grid_, i);
      }

    protected:
      index_type grid_;
      index_type layout_;
  };

  template <>
  class c_grid_padded<2>
  {
    public:
      typedef tiny<std::size_t, 2> index_type;
      typedef index_type::value_type index_value_type;

      c_grid_padded() : grid_(0,0), layout_(0,0) {}

      c_grid_padded(index_type const& grid)
      :
        grid_(grid),
        layout_(grid)
      {}

      c_grid_padded(index_value_type const& g0,
                    index_value_type const& g1)
      :
        grid_(g0, g1),
        layout_(g0, g1)
      {}

      c_grid_padded(index_type const& grid,
                    index_type const& layout)
      :
        grid_(grid),
        layout_(layout)
      {}

      c_grid_padded(index_value_type const& g0,
                    index_value_type const& g1,
                    index_value_type const& l0,
                    index_value_type const& l1)
      :
        grid_(g0, g1),
        layout_(l0, l1)
      {}

      template <typename FlexIndexType>
      c_grid_padded(af::flex_grid<FlexIndexType> const& flex_g)
      :
        grid_(af::adapt(flex_g.grid()))
      {
        SCITBX_ASSERT(flex_g.is_0_based());
        if (flex_g.is_padded()) {
          layout_ = index_type(af::adapt(flex_g.layout()));
        }
        else {
          layout_ = grid_;
        }
      }

      af::flex_grid<>
      as_flex_grid() const
      {
        return af::flex_grid<>(af::adapt(grid_))
                   .set_layout(af::adapt(layout_));
      }

      std::size_t
      size_1d() const
      {
        return grid_[0] * grid_[1];
      }

      index_type const&
      grid() const { return grid_; }

      index_type const&
      layout() const { return layout_; }

      std::size_t
      layout_size_1d() const
      {
        return layout_[0] * layout_[1];
      }

      bool is_padded() const
      {
        SCITBX_ASSERT(grid_[0] >= layout_[0] && grid_[1] >= layout_[1]);
        return grid_[0] != layout_[0] || grid_[1] != layout_[1];
      }

      std::size_t
      operator()(index_type const& i) const
      {
        return i[0] * grid_[1] + i[1];
      }

      std::size_t
      operator()(index_value_type const& i0,
                 index_value_type const& i1) const
      {
        return i0 * grid_[1] + i1;
      }

    protected:
      index_type grid_;
      index_type layout_;
  };

  template <>
  class c_grid_padded<3>
  {
    public:
      typedef tiny<std::size_t, 3> index_type;
      typedef index_type::value_type index_value_type;

      c_grid_padded() : grid_(0,0,0), layout_(0,0,0) {}

      c_grid_padded(index_type const& grid)
      :
        grid_(grid),
        layout_(grid)
      {}

      c_grid_padded(index_value_type const& g0,
                    index_value_type const& g1,
                    index_value_type const& g2)
      :
        grid_(g0, g1, g2),
        layout_(g0, g1, g2)
      {}

      c_grid_padded(index_type const& grid,
                    index_type const& layout)
      :
        grid_(grid),
        layout_(layout)
      {}

      c_grid_padded(index_value_type const& g0,
                    index_value_type const& g1,
                    index_value_type const& g2,
                    index_value_type const& l0,
                    index_value_type const& l1,
                    index_value_type const& l2)
      :
        grid_(g0, g1, g2),
        layout_(l0, l1, l2)
      {}

      template <typename FlexIndexType>
      c_grid_padded(af::flex_grid<FlexIndexType> const& flex_g)
      :
        grid_(af::adapt(flex_g.grid()))
      {
        SCITBX_ASSERT(flex_g.is_0_based());
        if (flex_g.is_padded()) {
          layout_ = index_type(af::adapt(flex_g.layout()));
        }
        else {
          layout_ = grid_;
        }
      }

      af::flex_grid<>
      as_flex_grid() const
      {
        return af::flex_grid<>(af::adapt(grid_))
                   .set_layout(af::adapt(layout_));
      }

      std::size_t
      size_1d() const
      {
        return grid_[0] * grid_[1] * grid_[2];
      }

      index_type const&
      grid() const { return grid_; }

      index_type const&
      layout() const { return layout_; }

      std::size_t
      layout_size_1d() const
      {
        return layout_[0] * layout_[1] * layout_[2];
      }

      bool is_padded() const
      {
        SCITBX_ASSERT(   grid_[0] >= layout_[0]
                      && grid_[1] >= layout_[1]
                      && grid_[2] >= layout_[2]);
        return    grid_[0] != layout_[0]
               || grid_[1] != layout_[1]
               || grid_[2] != layout_[2];
      }

      std::size_t
      operator()(index_type const& i) const
      {
        return (i[0] * grid_[1] + i[1]) * grid_[2] + i[2];
      }

      std::size_t
      operator()(index_value_type const& i0,
                 index_value_type const& i1,
                 index_value_type const& i2) const
      {
        return (i0 * grid_[1] + i1) * grid_[2] + i2;
      }

    protected:
      index_type grid_;
      index_type layout_;
  };

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_ACCESSORS_C_GRID_PADDED_H
