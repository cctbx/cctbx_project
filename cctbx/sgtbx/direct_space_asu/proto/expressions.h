#ifndef CCTBX_SGTBX_ASU_EXPRESSIONS_H
#define CCTBX_SGTBX_ASU_EXPRESSIONS_H

// \cond

#include <boost/static_assert.hpp>
#include "cut.h"

namespace cctbx { namespace sgtbx { namespace asu {

  template<typename TR> class cut_expression
  {
  public:
    cut lhs;
    TR rhs;

    void change_basis(const change_of_basis_op &o)
    {
      lhs.change_basis(o);
      rhs.change_basis(o);
    }

    void print(std::ostream &os, bool x=false) const
    {
      lhs.print(os);
      os << "[";
      rhs.print(os); 
      os << "]";
    }

    bool is_inside(const rvector3_t &p) const
    {
      return lhs.is_inside(p, rhs);
    }

    cut_expression<TR> operator- () const
    {
      cut_expression<TR> r(*this);
      r.lhs = -lhs;
      return r;
    }

    cut_expression<TR> operator~ () const
    {
      cut_expression<TR> r(*this);
      r.lhs = ~lhs;
      return r;
    }
 
    cut_expression(const cut &l, const TR &r) : lhs(l), rhs(r) { }

  }; // class cut_expression


  template<typename T> struct is_facet
  {
    static const bool value = false;
  };

  template<> struct is_facet< cut >
  {
    static const bool value = true;
  };

  template< typename TR> struct is_facet< cut_expression< TR>  >
  {
    static const bool value = true;
  };


  template<typename TL, typename TR> 
  class and_expression
  {
  public:
    TL lhs;
    TR rhs;

    void change_basis(const change_of_basis_op &o)
    {
      lhs.change_basis(o);
      rhs.change_basis(o);
    }

    void print(std::ostream &os, bool brk=false) const
    {
      if( brk )
      {
        lhs.print(os, true);
        os << "\n & ";
        rhs.print(os);
      }
      else
      {
        os << "(";
        lhs.print(os);
        os << " & ";
        rhs.print(os); 
        os << ")";
      }
    }

    bool is_inside(const rvector3_t &p) const
    {
      return lhs.is_inside(p) && rhs.is_inside(p);
    }

    and_expression(const TL &l, const TR &r) : lhs(l), rhs(r) { }

  }; // class and_expression


  // for cut and cut_epxression
  template<typename T>
    struct n_faces
  {
    static const size_type value =  1;
  };

  // specialization for and_expression
  template<typename TL, typename TR>
    struct n_faces< and_expression<TL,TR> >
  {
      static const size_type value = n_faces<TL>::value + n_faces<TR>::value;
  };


  inline size_type get_nth_plane(const cut &expr, size_type , cut &plane)
  {
    plane = expr;
    return 1;
  }

  template< typename TR>
    inline size_type get_nth_plane(const cut_expression<TR> &expr, size_type , cut &plane)
  {
    plane = expr.lhs;
    return 1;
  }

  template<typename TL, typename TR>  
  inline size_type get_nth_plane(const and_expression<TL,TR> &expr, size_type i, cut &plane)
  {
    BOOST_STATIC_ASSERT( is_facet< TR >::value );  // right-hand-side must be a facet
    if( i==0 )
    {
      get_nth_plane(expr.rhs, 0, plane);
      return 0;
    }
    return  1 + get_nth_plane(expr.lhs, i-1, plane);
  }


  template< typename T >
    struct strip
  {
    typedef int return_type;
    static return_type execute(const T &epxr)
    {
      return 0;
    }
  };

  template< >
    struct strip<cut>
  {
    typedef cut return_type;
    static return_type execute(const cut &expr)
    {
      return expr;
    }
  };


  template< typename TR >
    struct strip< cut_expression<TR> >
  {
    typedef cut return_type;
    static return_type execute(const cut_expression<TR> &expr)
    {
      return expr.lhs;
    }
  };


  template< typename TL, typename TR >
    struct strip< and_expression<TL,TR> >
  {
    typedef typename strip<TL>::return_type left_type;
    typedef typename strip<TR>::return_type right_type;

    typedef and_expression< left_type, right_type > return_type;

    static return_type execute(const and_expression<TL,TR> &expr)
    {
      return return_type( strip<TL>::execute(expr.lhs), strip<TR>::execute(expr.rhs) );
    }
  };

  template<typename TL, typename TR> class or_expression
  {
  public:
    TL lhs;
    TR rhs;

    void change_basis(const change_of_basis_op &o)
    {
      lhs.change_basis(o);
      rhs.change_basis(o);
    }

    void print(std::ostream &os, bool x=false) const
    {
      os << "(";
      lhs.print(os);
      os << " | ";
      rhs.print(os); 
      os << ")";
    }

    bool is_inside(const rvector3_t &p) const
    {
      return lhs.is_inside(p) || rhs.is_inside(p);
    }

    or_expression(const TL &l, const TR &r) : lhs(l), rhs(r) { }

  }; // or_expression


  // operations
  template< typename TL, typename TR> 
  inline or_expression<TL,TR> operator| (const TL &a, const TR &b) 
  {
    return or_expression<TL,TR>(a, b);
  }

  template< typename TL, typename TR>
  inline and_expression<TL,TR> operator& (const TL &a, const TR &b)
  {
    return and_expression<TL,TR>(a, b);
  }


  template<typename T> struct is_facet_expression
  {
    static const bool value = false;
  };

  template<typename TL, typename TR> struct is_facet_expression< and_expression< TL, cut_expression<TR> > >
  {
    static const bool value = true;
  };

  template< typename TL > struct is_facet_expression< and_expression< TL, cut > >
  {
    static const bool value = true;
  };



}}}
// \endcond
#endif

