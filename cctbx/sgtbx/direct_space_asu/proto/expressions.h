#ifndef CCTBX_SGTBX_ASU_EXPRESSIONS_H
#define CCTBX_SGTBX_ASU_EXPRESSIONS_H

// \cond

#include <boost/static_assert.hpp>
#include "cut.h"

namespace cctbx { namespace sgtbx { namespace asu {

  template<typename TR>
    class cut_expression
  {
  public:
    cut lhs;
    TR rhs;

    void change_basis(const change_of_basis_op &o)
    {
      lhs.change_basis(o);
      rhs.change_basis(o);
    }

    void print(std::ostream &os ) const
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

    bool is_inside(const scitbx::af::int3 &num, const scitbx::af::int3 &den)
      const
    {
      return lhs.is_inside(num,den,rhs);
    }

    bool is_inside(const scitbx::af::int3 &num) const
    {
      return lhs.is_inside(num,rhs);
    }

    short where_is(const scitbx::af::int3 &num, const scitbx::af::int3 &den)
      const
    {
      return lhs.where_is(num,den,rhs);
    }

    short where_is(const scitbx::af::int3 &num) const
    {
      return lhs.where_is(num,rhs);
    }


    bool is_inside_shape_only(const scitbx::af::double3 &p, double tol) const
    {
      return lhs.is_inside_shape_only(p,tol);
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


  template<typename T>
    struct is_facet
  {
    static const bool value = false;
  };

  template<>
    struct is_facet< cut >
  {
    static const bool value = true;
  };

  template< typename TR>
    struct is_facet< cut_expression< TR>  >
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

    void print(std::ostream &os) const
    {
        os << "(";
        lhs.print(os);
        os << " & ";
        rhs.print(os);
        os << ")";
    }

    bool is_inside(const rvector3_t &p) const
    {
      return lhs.is_inside(p) && rhs.is_inside(p);
    }

    bool is_inside(const scitbx::af::int3 &num, const scitbx::af::int3 &den) const
    {
      return lhs.is_inside(num,den) && rhs.is_inside(num,den);
    }

    bool is_inside(const scitbx::af::int3 &num) const
    {
      return lhs.is_inside(num) && rhs.is_inside(num);
    }


    short where_is(const scitbx::af::int3 &num, const scitbx::af::int3 &den) const
    {
      const short l = lhs.where_is(num,den),
        r = rhs.where_is(num,den);
      if( r==1 && l==1 ) // fully inside
        return 1;
      if( r==0 || l==0 ) // fully outside
        return 0;
      return -1; // on the face
    }

    short where_is(const scitbx::af::int3 &num) const
    {
      const short l = lhs.where_is(num),
        r = rhs.where_is(num);
      if( r==1 && l==1 ) // fully inside
        return 1;
      if( r==0 || l==0 ) // fully outside
        return 0;
      return -1; // on the face
    }


    bool is_inside_shape_only(const scitbx::af::double3 &p, double tol) const
    {
      return lhs.is_inside_shape_only(p,tol) && rhs.is_inside_shape_only(p,tol);
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

  inline double get_tolerance(const cut &c, const scitbx::af::double3 &tol)
  {
    return c.get_tolerance(tol);
  }

  template<typename Expr>
    inline double get_tolerance(const Expr &expr, const scitbx::af::double3 &tol)
  {
    return std::max(get_tolerance(expr.lhs, tol), get_tolerance(expr.rhs, tol));
  }

  inline void optimize_for_grid(cut &c, const scitbx::af::int3 &grid_size)
  {
    c.optimize_for_grid(grid_size);
  }

  template<typename Expr>
    inline void optimize_for_grid(Expr &expr, const scitbx::af::int3 &grid_size)
  {
    optimize_for_grid(expr.lhs, grid_size);
    optimize_for_grid(expr.rhs, grid_size);
  }

  inline void get_optimized_grid_limits(const cut &c, scitbx::af::long3 &max_p)
  {
    c.get_optimized_grid_limits(max_p);
  }

  template<typename Expr> inline
    void get_optimized_grid_limits(const Expr &expr, scitbx::af::long3 &max_p)
  {
    scitbx::af::long3 m1, m2;
    get_optimized_grid_limits(expr.lhs, m1);
    get_optimized_grid_limits(expr.rhs, m2);
    for(unsigned char i=0; i<3U; ++i)
      max_p[i] = std::min(m1[i],m2[i]);
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
      cut cut_inc(expr);
      cut_inc.inclusive = true;
      return cut_inc;
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

  template< typename T, bool Break = false >
    struct print_lines
  {
    static void execute(const T &expr, std::ostream &os)
    {
      expr.print(os);
    }
  };


  template< typename TL, typename TR >
    struct print_lines< and_expression<TL,TR>, true >
  {
    static void execute(const and_expression<TL,TR> &expr, std::ostream &os)
    {
      print_lines<TL,true>::execute(expr.lhs,os);
      os << "\n & ";
      print_lines<TR>::execute(expr.rhs,os);
    }
  };


  template<typename TL, typename TR>
    class or_expression
  {
  public:
    TL lhs;
    TR rhs;

    void change_basis(const change_of_basis_op &o)
    {
      lhs.change_basis(o);
      rhs.change_basis(o);
    }

    void print(std::ostream &os) const
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

    bool is_inside(const scitbx::af::int3 &num, const scitbx::af::int3 &den) const
    {
      return lhs.is_inside(num,den) || rhs.is_inside(num,den);
    }

    bool is_inside(const scitbx::af::int3 &num) const
    {
      return lhs.is_inside(num) || rhs.is_inside(num);
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


  template<typename T>
    struct is_facet_expression
  {
    static const bool value = false;
  };

  template<typename TL, typename TR>
    struct is_facet_expression< and_expression< TL, cut_expression<TR> > >
  {
    static const bool value = true;
  };

  template< typename TL >
    struct is_facet_expression< and_expression< TL, cut > >
  {
    static const bool value = true;
  };


  // These are experimental
  template< typename T >
    struct strip_keep_inclusive_flag
  {
    typedef int return_type;
    static return_type execute(const T &epxr)
    {
      return 0;
    }
  };

  template< >
    struct strip_keep_inclusive_flag<cut>
  {
    typedef cut return_type;
    static return_type execute(const cut &expr)
    {
      cut cut_inc(expr);
      return cut_inc;
    }
  };


  template< typename TR >
    struct strip_keep_inclusive_flag< cut_expression<TR> >
  {
    typedef cut return_type;
    static return_type execute(const cut_expression<TR> &expr)
    {
      return expr.lhs;
    }
  };


  template< typename TL, typename TR >
    struct strip_keep_inclusive_flag< and_expression<TL,TR> >
  {
    typedef typename strip_keep_inclusive_flag<TL>::return_type left_type;
    typedef typename strip_keep_inclusive_flag<TR>::return_type right_type;

    typedef and_expression< left_type, right_type > return_type;

    static return_type execute(const and_expression<TL,TR> &expr)
    {
      return return_type( strip_keep_inclusive_flag<TL>::execute(expr.lhs),
        strip_keep_inclusive_flag<TR>::execute(expr.rhs) );
    }
  };



}}}
// \endcond
#endif

