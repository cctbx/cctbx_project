#ifndef CCTBX_SGTBX_ASU_REFERENCE_TABLE_H
#define CCTBX_SGTBX_ASU_REFERENCE_TABLE_H

// \cond

///
/// This file is for internal use only
///

#include <boost/static_assert.hpp>

#include "facet_collection.h"
#include "expressions.h"
#include "shortcuts.h"

namespace cctbx { namespace sgtbx { namespace asu {

  template< typename T> class expression_adaptor : public facet_collection
  {
      BOOST_STATIC_ASSERT( is_facet_expression<T>::value );

    public:
      T obj;

    void change_basis(const change_of_basis_op &o)
    {
      obj.change_basis(o);
    }

    void print(std::ostream &os) const
    {
      print_lines<T,true>::execute(obj,os);
    }

    bool is_inside(const rvector3_t &p) const
    {
      return obj.is_inside(p);
    }

    bool is_inside(const scitbx::af::int3 &num, const scitbx::af::int3 &den) const
    {
      return obj.is_inside(num,den);
    }

    short where_is(const scitbx::af::int3 &num, const scitbx::af::int3 &den) const
    {
      return obj.where_is(num,den);
    }


    size_type size() const
    {
      return n_faces<T>::value;
    }

    void get_nth_plane(size_type i, cut &plane) const
    {
      cctbx::sgtbx::asu::get_nth_plane(obj, i, plane);  // obj must be and_expression
    }

    expression_adaptor(const T &o) : obj(o) { }

    facet_collection::pointer new_volume_only() const
    {
      typedef typename strip<T>::return_type return_type;
      return facet_collection::pointer( new expression_adaptor< return_type >( strip<T>::execute(obj) ) );
    }

    facet_collection::pointer new_copy() const
    {
      return facet_collection::pointer( new expression_adaptor<T>(*this) );
    }

  };

  template<typename TL, typename TR>
    facet_collection::pointer facet_collection_asu(const and_expression<TL,TR> &expr)
  {
      // BOOST_STATIC_ASSERT( is_facet_expression< and_expression<TL,TR> >::value );
      return expression_adaptor< and_expression<TL,TR> > ( expr ).new_copy();
  }

}}}
// \endcond
#endif

