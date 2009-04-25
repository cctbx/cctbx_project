// The implementations in this file are highly redundant in order
// to help the optimizer.

#ifndef SCITBX_ARRAY_FAMILY_ACCESSORS_INTERVAL_GRID_H
#define SCITBX_ARRAY_FAMILY_ACCESSORS_INTERVAL_GRID_H

#include <scitbx/error.h>
#include <scitbx/array_family/accessors/c_index_1d_calculator.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/math/compile_time.h>

namespace scitbx { namespace af {

  template <std::size_t N>
  struct interval_calculator
  {
    template <typename ExtendArrayType, typename IndexType>
      static std::size_t
    get_1d_index(ExtendArrayType const& e, IndexType const& i)
    {
      return interval_calculator<N-1>::get_1d_index(e, i) * e.size()[N-1] + i[N-1]-e.first()[N-1];
    }

    template <typename IndexType>
      static void get_size(IndexType &size, const IndexType &first, const IndexType &last)
    {
      interval_calculator<N-1>::get_size(size, first, last);
      size[N-1] = last[N-1]-first[N-1];
    }


    template <typename ExtendArrayType, typename IndexType>
      static bool is_index_valid(const ExtendArrayType &e, const IndexType &ind_nd)
    {
      return interval_calculator<N-1>::is_index_valid(e,ind_nd) &&
        (ind_nd[N-1] >=e.first()[N-1]) && (ind_nd[N-1]<(e.first()[N-1] + e.size()[N-1]));
    }
  };

  template<>
  struct interval_calculator<1>
  {
    template <typename ExtendArrayType, typename IndexType>
    static std::size_t
    get_1d_index(ExtendArrayType const& e, IndexType const& i)
    {
      return i[0] - e.first()[0];
    }

    template <typename ExtendArrayType, typename IndexType>
      static bool is_index_valid(const ExtendArrayType &e, const IndexType &ind_nd)
    {
      return (ind_nd[0] >= e.first()[0]) && (ind_nd[0] < (e.first()[0] + e.size()[0]));
    }

    template <typename IndexType>
      static void get_size(IndexType &size, const IndexType &first, IndexType const& last)
    {
      size[0] = last[0]-first[0];
    }
  };


  template <std::size_t Nd, typename IndexValueType=long>
  class c_interval_grid
  {
    public:
      typedef tiny<IndexValueType, Nd> index_type;
      typedef IndexValueType index_value_type;
      typedef index_value_type value_type;

      c_interval_grid() { for(std::size_t i=0;i<Nd;i++) i3first[i] = i3size[i] = 0; }

      c_interval_grid(index_type const& first_, const index_type &last_) : i3first(first_)
      {
        interval_calculator<Nd>::get_size(this->i3size, this->first(), last_);
        for(unsigned i=0; i<Nd; ++i)
          if( this->size()[i]<=0 )
            throw scitbx::error("Invalid interval");
      }

      const index_type& size() const { return i3size; }
      const index_type& first() const { return i3first; }
      index_type last() const // first + size
      {
        index_type last_;
        for(register unsigned short i=0; i<Nd; ++i)
          last_[i] = this->first()[i] + this->size()[i];
        return last_;
      }

      std::size_t size_1d() const
      {
        return math::compile_time::product<Nd>::get(this->size());
      }

      index_type index_nd(index_value_type const& i_1d) const
      {
        index_type i_nd;
        i_nd[0] = i_1d;
        for(std::size_t j=this->size().size()-1;j;j--) {
          i_nd[j] = i_nd[0] % this->size()[j];
          i_nd[0] /= this->size()[j];
        }
        return i_nd + this->first();
      }

      bool is_valid_index(index_type const& i) const
      {
        return interval_calculator<Nd>::is_index_valid(*this, i);
      }

      std::size_t operator()(index_type const& index_nd) const
      {
        return interval_calculator<Nd>::get_1d_index(*this, index_nd);
      }

    private:
      index_type i3first, i3size;
  };

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_ACCESSORS_INTERVAL_H
