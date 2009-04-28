// The implementations in this file are highly redundant in order
// to help the optimizer.

#ifndef SCITBX_ARRAY_FAMILY_ACCESSORS_C_GRID_H
#define SCITBX_ARRAY_FAMILY_ACCESSORS_C_GRID_H

#include <scitbx/array_family/accessors/c_index_1d_calculator.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/math/compile_time.h>

namespace scitbx { namespace af {

  template <std::size_t Nd, typename IndexValueType=std::size_t>
  class c_grid : public tiny<IndexValueType, Nd>
  {
    public:
      typedef IndexValueType index_value_type;
      typedef index_value_type value_type;
      typedef tiny<IndexValueType, Nd> index_type;

      c_grid() { for(std::size_t i=0;i<Nd;i++) this->elems[i] = 0; }

      c_grid(index_type const& n) : index_type(n) {}

      template <typename FlexIndexType>
      c_grid(flex_grid<FlexIndexType> const& flex_g)
      :
        index_type(af::adapt(flex_g.all()))
      {
        SCITBX_ASSERT(flex_g.is_0_based());
        SCITBX_ASSERT(!flex_g.is_padded());
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

  template <typename IndexValueType>
  class c_grid<2, IndexValueType> : public tiny<IndexValueType, 2>
  {
    public:
      typedef IndexValueType index_value_type;
      typedef index_value_type value_type;
      typedef tiny<IndexValueType, 2> index_type;

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
        SCITBX_ASSERT(!flex_g.is_padded());
      }

      af::flex_grid<>
      as_flex_grid() const
      {
        return af::flex_grid<>(af::adapt(*this));
      }

      std::size_t
      size_1d() const
      {
        return static_cast<std::size_t>(this->elems[0])
             * static_cast<std::size_t>(this->elems[1]);
      }

      bool
      is_square() const { return this->elems[0] == this->elems[1]; }

      std::size_t n_rows() const { return this->elems[0]; }
      std::size_t n_columns() const { return this->elems[1]; }

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

      bool
      is_valid_index(
        index_value_type const& i0,
        index_value_type const& i1) const
      {
        if (i0 >= this->elems[0]) return false;
        if (i1 >= this->elems[1]) return false;
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

  template <typename IndexValueType>
  class c_grid<3, IndexValueType> : public tiny<IndexValueType, 3>
  {
    public:
      typedef IndexValueType index_value_type;
      typedef index_value_type value_type;
      typedef tiny<IndexValueType, 3> index_type;

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
        SCITBX_ASSERT(!flex_g.is_padded());
      }

      af::flex_grid<>
      as_flex_grid() const
      {
        return af::flex_grid<>(af::adapt(*this));
      }

      std::size_t
      size_1d() const
      {
        return static_cast<std::size_t>(this->elems[0])
             * static_cast<std::size_t>(this->elems[1])
             * static_cast<std::size_t>(this->elems[2]);
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

      bool
      is_valid_index(
        index_value_type const& i0,
        index_value_type const& i1,
        index_value_type const& i2) const
      {
        if (i0 >= this->elems[0]) return false;
        if (i1 >= this->elems[1]) return false;
        if (i2 >= this->elems[2]) return false;
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
