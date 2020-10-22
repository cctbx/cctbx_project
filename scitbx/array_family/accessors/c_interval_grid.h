// The implementations in this file are highly redundant in order
// to help the optimizer.

#ifndef SCITBX_ARRAY_FAMILY_ACCESSORS_INTERVAL_GRID_H
#define SCITBX_ARRAY_FAMILY_ACCESSORS_INTERVAL_GRID_H

#include <scitbx/error.h>
#include <scitbx/array_family/accessors/c_index_1d_calculator.h>
#include <scitbx/math/compile_time.h>

namespace scitbx { namespace af {

  template <std::size_t N>
  struct interval_calculator
  {
    template <typename ExtendArrayType, typename IndexType>
      static std::size_t
    get_1d_index(ExtendArrayType const& e, IndexType const& i)
    {
      return interval_calculator<N-1>::get_1d_index(e, i) * e.all()[N-1]
        + i[N-1]-e.origin()[N-1];
    }

    template <typename IndexType>
      static void get_size(IndexType &size, const IndexType &origin,
        const IndexType &last)
    {
      interval_calculator<N-1>::get_size(size, origin, last);
      size[N-1] = last[N-1]-origin[N-1];
    }

    template <typename ExtendArrayType, typename IndexType>
      static bool is_index_valid(const ExtendArrayType &e,
        const IndexType &ind_nd)
    {
      return interval_calculator<N-1>::is_index_valid(e,ind_nd) &&
        (ind_nd[N-1] >=e.origin()[N-1]) && (ind_nd[N-1]<(e.origin()[N-1]
        + e.size()[N-1]));
    }
  };

  template<>
  struct interval_calculator<1>
  {
    template <typename ExtendArrayType, typename IndexType>
    static std::size_t
    get_1d_index(ExtendArrayType const& e, IndexType const& i)
    {
      return i[0] - e.origin()[0];
    }

    template <typename ExtendArrayType, typename IndexType>
      static bool is_index_valid(const ExtendArrayType &e,
        const IndexType &ind_nd)
    {
      return (ind_nd[0] >= e.origin()[0]) && (ind_nd[0] < (e.origin()[0]
        + e.size()[0]));
    }

    template <typename IndexType>
      static void get_size(IndexType &size, const IndexType &origin,
        IndexType const& last)
    {
      size[0] = last[0]-origin[0];
    }
  };


  template <std::size_t Nd, typename IndexValueType=long>
  class c_interval_grid
  {
    public:
      typedef tiny<IndexValueType, Nd> index_type;
      typedef IndexValueType index_value_type;
      typedef index_value_type value_type;

      c_interval_grid()
      {
        for(std::size_t i=0;i<Nd;i++) origin_[i] = all_[i] = 0;
      }

      c_interval_grid(index_type const& origin, const index_type &last)
        : origin_(origin)
      {
        interval_calculator<Nd>::get_size(this->all_, this->origin(), last);
        for(unsigned i=0; i<Nd; ++i)
          if( this->all()[i]<=0 )
            throw scitbx::error("Invalid interval");
      }

      template <typename FlexIndexType>
      c_interval_grid(af::flex_grid<FlexIndexType> const& flex_g)
        : origin_(af::adapt(flex_g.origin())), all_(af::adapt(flex_g.all()))
      {
        SCITBX_ASSERT(flex_g.is_0_based());
        SCITBX_ASSERT(!flex_g.is_padded());
      }

      af::flex_grid<> as_flex_grid() const
      {
        return af::flex_grid<>(af::adapt(origin_), af::adapt(this->last()));
      }

      const index_type& all() const { return all_; }
      const index_type& origin() const { return origin_; }
      index_type last(bool open_range=true) const // origin + size
      {
        index_type result;
        for(unsigned short i=0; i<Nd; ++i)
          result[i] = this->origin()[i] + this->all()[i];
        if( !open_range )
          result -= static_cast<index_value_type> (1);
        return result;
      }

      std::size_t size_1d() const
      {
        return math::compile_time::product<Nd>::get(this->all());
      }

      index_type index_nd(index_value_type const& i_1d) const
      {
        index_type i_nd;
        i_nd[0] = i_1d;
        for(std::size_t j=this->all().size()-1;j;j--) {
          i_nd[j] = i_nd[0] % this->all()[j];
          i_nd[0] /= this->all()[j];
        }
        return i_nd + this->origin();
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
      index_type origin_, all_;
  };

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_ACCESSORS_INTERVAL_H
