/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Created (R.W. Grosse-Kunstleve)
 */

// The implementations in this file are highly redundant in order
// to help the optimizer.

#ifndef SCITBX_ARRAY_FAMILY_ACCESSORS_C_GRID_H
#define SCITBX_ARRAY_FAMILY_ACCESSORS_C_GRID_H

#include <scitbx/array_family/accessors/c_index_1d_calculator.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/math/compile_time.h>

namespace scitbx { namespace af {

  template <std::size_t Nd>
  class c_grid : public tiny<std::size_t, Nd>
  {
    public:
      typedef tiny<std::size_t, Nd> index_type;
      typedef typename index_type::value_type index_value_type;
      typedef index_value_type value_type;

      c_grid() { for(std::size_t i=0;i<Nd;i++) this->elems[i] = 0; }

      c_grid(index_type const& n) : index_type(n) {}

      template <typename FlexIndexType>
      c_grid(flex_grid<FlexIndexType> const& flex_g)
      :
        index_type(af::adapt(flex_g.all()))
      {
        SCITBX_ASSERT(flex_g.is_0_based());
      }

      af::flex_grid<>
      as_flex_grid() const
      {
        return af::flex_grid<>(af::adapt(*this));
      }

      std::size_t
      size_1d() const
      {
        return math::compile_time::product<Nd>::get(*this);
      }

      index_type
      index_nd(index_value_type const& i_1d) const
      {
        index_type i_nd;
        i_nd[0] = i_1d;
        for(std::size_t j=this->size()-1;j;j--) {
          i_nd[j] = i_nd[0] % this->elems[j];
          i_nd[0] /= this->elems[j];
        }
        return i_nd;
      }

      bool
      is_valid_index(index_type const& i) const
      {
        for(std::size_t j=0;j<this->size();j++) {
          if (i[j] >= this->elems[j]) return false;
        }
        return true;
      }

      std::size_t
      operator()(index_type const& i) const
      {
        return c_index_1d_calculator<Nd>::get(*this, i);
      }
  };

  template <>
  class c_grid<2> : public tiny<std::size_t, 2>
  {
    public:
      typedef tiny<std::size_t, 2> index_type;
      typedef index_type::value_type index_value_type;
      typedef index_value_type value_type;

      c_grid() : index_type(0,0) {}

      c_grid(index_type const& n) : index_type(n) {}

      c_grid(index_value_type const& n0,
             index_value_type const& n1)
      : index_type(n0, n1)
      {}

      template <typename FlexIndexType>
      c_grid(flex_grid<FlexIndexType> const& flex_g)
      :
        index_type(af::adapt(flex_g.all()))
      {
        SCITBX_ASSERT(flex_g.is_0_based());
      }

      af::flex_grid<>
      as_flex_grid() const
      {
        return af::flex_grid<>(af::adapt(*this));
      }

      std::size_t
      size_1d() const
      {
        return this->elems[0] * this->elems[1];
      }

      index_type
      index_nd(index_value_type const& i_1d)
      {
        return index_type(i_1d / this->elems[1], i_1d % this->elems[1]);
      }

      bool
      is_valid_index(index_type const& i) const
      {
        if (i[0] >= this->elems[0]) return false;
        if (i[1] >= this->elems[1]) return false;
        return true;
      }

      std::size_t
      operator()(index_type const& i) const
      {
        return i[0] * this->elems[1] + i[1];
      }

      std::size_t
      operator()(index_value_type const& i0,
                 index_value_type const& i1) const
      {
        return i0 * this->elems[1] + i1;
      }
  };

  template <>
  class c_grid<3> : public tiny<std::size_t, 3>
  {
    public:
      typedef tiny<std::size_t, 3> index_type;
      typedef index_type::value_type index_value_type;
      typedef index_value_type value_type;

      c_grid() : index_type(0,0,0) {}

      c_grid(index_type const& n) : index_type(n) {}

      c_grid(index_value_type const& n0,
             index_value_type const& n1,
             index_value_type const& n2)
      : index_type(n0, n1, n2)
      {}

      template <typename FlexIndexType>
      c_grid(flex_grid<FlexIndexType> const& flex_g)
      :
        index_type(af::adapt(flex_g.all()))
      {
        SCITBX_ASSERT(flex_g.is_0_based());
      }

      af::flex_grid<>
      as_flex_grid() const
      {
        return af::flex_grid<>(af::adapt(*this));
      }

      std::size_t
      size_1d() const
      {
        return this->elems[0] * this->elems[1] * this->elems[2];
      }

      index_type
      index_nd(index_value_type const& i_1d)
      {
        index_type i_nd;
        i_nd[2] = i_1d % this->elems[2];
        i_nd[0] = i_1d / this->elems[2];
        i_nd[1] = i_nd[0] % this->elems[1];
        i_nd[0] /= this->elems[1];
        return i_nd;
      }

      bool
      is_valid_index(index_type const& i) const
      {
        if (i[0] >= this->elems[0]) return false;
        if (i[1] >= this->elems[1]) return false;
        if (i[2] >= this->elems[2]) return false;
        return true;
      }

      std::size_t
      operator()(index_type const& i) const
      {
        return (i[0] * this->elems[1] + i[1]) * this->elems[2] + i[2];
      }

      std::size_t
      operator()(index_value_type const& i0,
                 index_value_type const& i1,
                 index_value_type const& i2) const
      {
        return (i0 * this->elems[1] + i1) * this->elems[2] + i2;
      }
  };

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_ACCESSORS_C_GRID_H
