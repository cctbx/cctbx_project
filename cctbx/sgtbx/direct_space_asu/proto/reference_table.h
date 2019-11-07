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

    void print_as_xyz(std::ostream &os) const
    {
      print_lines_as_xyz<T,true>::execute(obj,os);
    }

    bool is_inside(const rvector3_t &p) const
    {
      return obj.is_inside(p);
    }

    bool is_inside(const scitbx::af::int3 &num, const scitbx::af::int3 &den)
      const
    {
      return obj.is_inside(num,den);
    }

    bool is_inside(const scitbx::af::int3 &num) const
    {
      return obj.is_inside(num);
    }

    bool is_inside_shape_only(const scitbx::af::double3 &point, double tol)
      const
    {
      return obj.is_inside_shape_only(point,tol);
    }


    short where_is(const scitbx::af::int3 &num, const scitbx::af::int3 &den)
      const
    {
      return obj.where_is(num,den);
    }

    short where_is(const scitbx::af::int3 &num) const
    {
      return obj.where_is(num);
    }

    void optimize_for_grid(const scitbx::af::int3 &grid_size)
    {
      cctbx::sgtbx::asu::optimize_for_grid(obj, grid_size);
    }

    void get_optimized_grid_limits(scitbx::af::long3 &max_p) const
    {
      cctbx::sgtbx::asu::get_optimized_grid_limits(obj, max_p);
    }

    size_type size() const
    {
      return n_faces<T>::value;
    }

    void get_nth_plane(size_type i, cut &plane) const
    {
      // obj must be and_expression
      cctbx::sgtbx::asu::get_nth_plane(obj, i, plane);
    }

    double get_tolerance(const scitbx::af::double3 &tol3d) const
    {
      return cctbx::sgtbx::asu::get_tolerance(obj, tol3d);
    }

    expression_adaptor(const T &o) : obj(o) { }

    facet_collection::pointer new_shape_only() const
    {
      typedef typename strip<T>::return_type return_type;
      return facet_collection::pointer( new expression_adaptor< return_type >( strip<T>::execute(obj) ) );
    }

    // DO NOT USE!! Experimental
    facet_collection::pointer new_shape_only_keep_inclusive_flag() const
    {
      typedef typename strip_keep_inclusive_flag<T>::return_type return_type;
      return facet_collection::pointer( new expression_adaptor< return_type >(
        strip_keep_inclusive_flag<T>::execute(obj) ) );
    }

    // DO NOT USE!!!
    facet_collection::pointer add_face(const cut &) const
    {
      /* and_expression<T,cut> r = (this->obj) & face;
      expression_adaptor< and_expression<T,cut> > rr(r);
      return rr.new_copy(); */
      return this->new_copy();
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
