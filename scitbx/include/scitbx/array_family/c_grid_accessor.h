/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_ARRAY_FAMILY_C_GRID_ACCESSOR_H
#define SCITBX_ARRAY_FAMILY_C_GRID_ACCESSOR_H

#include <scitbx/array_family/tiny.h>

namespace scitbx { namespace af {

  template <std::size_t Nd>
  class c_grid_base : public tiny<std::size_t, Nd>
  {
    public:
      typedef tiny<std::size_t, Nd> index_type;
      typedef typename index_type::value_type index_value_type;
      typedef index_value_type value_type;

      c_grid_base(index_type const& n) : index_type(n) {}

      bool
      is_valid_index(index_type const& i) const
      {
        for(std::size_t j=0;j<this->size();j++) {
          if (i[j] >= this->elems[j]) return false;
        }
        return true;
      }

      index_type
      index_nd(index_value_type i_1d)
      {
        index_type i_nd;
        for(std::size_t j=this->size()-1;j;j--) {
          i_nd[j] = i_1d % this->elems[j];
          i_1d /= this->elems[j];
        }
        i_nd[0] = i_1d;
        return i_nd;
      }
  };

  template <std::size_t Nd>
  class c_grid;

  template <>
  class c_grid<2> : public c_grid_base<2>
  {
    public:
      typedef c_grid_base<2> base_type;

      c_grid() : base_type(index_type(0,0)) {}

      c_grid(index_type const& n) : base_type(n) {}

      c_grid(index_value_type const& n0,
             index_value_type const& n1)
      : base_type(index_type(n0, n1))
      {}

      template <typename OtherArrayType>
      c_grid(af::array_adaptor<OtherArrayType> const& a_a)
      : base_type(index_type(a_a))
      {}

      std::size_t
      size_1d() const
      {
        return this->elems[0] * this->elems[1];
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
  class c_grid<3> : public c_grid_base<3>
  {
    public:
      typedef c_grid_base<3> base_type;

      c_grid() : base_type(index_type(0,0,0)) {}

      c_grid(index_type const& n) : base_type(n) {}

      c_grid(index_value_type const& n0,
             index_value_type const& n1,
             index_value_type const& n2)
      : base_type(index_type(n0, n1, n2))
      {}

      template <typename OtherArrayType>
      c_grid(af::array_adaptor<OtherArrayType> const& a_a)
      : base_type(index_type(a_a))
      {}

      std::size_t
      size_1d() const
      {
        return this->elems[0] * this->elems[1] * this->elems[2];
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

#endif // SCITBX_ARRAY_FAMILY_C_GRID_ACCESSOR_H
