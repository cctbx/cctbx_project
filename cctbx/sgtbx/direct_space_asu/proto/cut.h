#ifndef CCTBX_SGTBX_ASU_CUT_H
#define CCTBX_SGTBX_ASU_CUT_H

#include <iosfwd>
#include <vector>
#include <set>

#include <cmath>

#include <boost/rational.hpp>
#include <scitbx/mat3.h>
#include <cctbx/sgtbx/space_group.h>
#include <cctbx/sgtbx/symbols.h>
#include <cctbx/sgtbx/change_of_basis_op.h>
#include <cctbx/sgtbx/space_group_type.h>


namespace cctbx { namespace sgtbx { 
  
//! Direct space asymmetric unit classes
namespace asu {

  typedef unsigned short size_type;
  // typedef short int_type;
  typedef int int_type;
  typedef boost::rational<int_type> rational_t;

  typedef scitbx::mat3<int_type> imatrix3_t;
  typedef scitbx::vec3<int_type> ivector3_t;
  typedef scitbx::vec3<double> dvector3_t;
  typedef scitbx::vec3< rational_t > rvector3_t;
  typedef scitbx::mat3< rational_t > rmatrix3_t;
  typedef std::set< rvector3_t > set_rvector3_t;

  using cctbx::sgtbx::change_of_basis_op;

  inline std::ostream &operator<< (std::ostream &s, const ivector3_t &v)
  {
    s << "("<< v[0] <<", "<< v[1] << ", " << v[2] << ")";
    return s;
  }

  inline std::ostream &operator<< (std::ostream &s, const rvector3_t &v)
  {
    s << "("<< v[0] <<", "<< v[1] << ", " << v[2] << ")";
    return s;
  }

}}} // namespace cctbx::sgtbx::asu


namespace boost {

  // mrt: WHY DO I HAVE TO DEFINE this operator in boost:: { } ?????????????????
  //
  inline bool operator< (const cctbx::sgtbx::asu::rvector3_t &a, const cctbx::sgtbx::asu::rvector3_t &b)
  {
    return (a[0]<b[0]) || (a[0]==b[0] && a[1]<b[1]) || (a[0]==b[0] && a[1]==b[1] && a[2]<b[2]);
  }
} // namespace boost

namespace cctbx { namespace sgtbx { namespace asu {

  template<typename T> class cut_expression;

  //! Plane in 3D for representation of faces of the asymmetric unit polyhedron
  /*! \class cut cut.h cctbx/sgtbx/direct_space_asu/proto/cut.h
   */
  class cut
  {
  public:
    //! Plane normal
    ivector3_t n;
    //! Plane constant
    int_type c;
    //! Flag indicating if plane points belong to the asymmetric unit
    bool inclusive;

    void print(std::ostream &os, bool x=false ) const;

    //! Returns an arbitrary point, which belongs to the plane
    void get_point_in_plane(rvector3_t &r) const;

    //! Applies basis change to the plane
    void change_basis(const change_of_basis_op &cb_op);

    //! Evaluates plane equation
    rational_t evaluate(const rvector3_t &p) const
    {
      return  n[0]*p[0] + n[1]*p[1] + n[2]*p[2] + c;
    }

    template<typename TR> bool is_inside(const rvector3_t &p, const TR &expr) const
    {
      rational_t v = evaluate(p);
      if( v>0 )
        return true;
      if( v<0 )
        return false;
      return expr.is_inside(p);
    }

    //! Tests if point is on the inside part of the space
    bool is_inside(const rvector3_t &p) const
    {
      rational_t v = evaluate(p);
      if( v>0 )
        return true;
      if( v<0 )
        return false;
      return inclusive;
    }

    template<typename TR> 
    cut_expression<TR> operator() (const TR &o) const // better:  operator[]
    {
      return cut_expression<TR>(*this, o );
    }

    //! Unsets inclusive flag
    cut operator+ () const
    {
      CCTBX_ASSERT( inclusive );
      cut r(*this);
      r.inclusive = false;
      return r;
    }

    //! Returns -n, -c plane
    cut operator- () const
    {
      return cut(-n, -c, inclusive);
    }

    //!  Returns -n plane
    cut operator~ () const
    {
      cut r(*this);
      r.n = -n;
      return r;
    }

    //! Return c*x plane
    cut operator* (int_type x) const
    {
      return cut(n, c*x, inclusive);
    }

    //! Returns c/x plane
    cut operator/ (int_type x) const
    {
      CCTBX_ASSERT( x!=0 && c!=0 );
      return cut(n*std::abs(x), c*x/std::abs(x), inclusive);
    }

    //! Returns c=1  plane
    cut one() const
    {
      return cut(n, 1, inclusive);
    }

    //! Constructs a plane from a normal vector and a rational constant
    cut(const ivector3_t &n_, rational_t c_, bool inc_=true) : inclusive(inc_)
    { 
      CCTBX_ASSERT( c_.denominator() > 0 );
      n = n_ * c_.denominator();
      c = c_.numerator();
      int_type g = boost::gcd(boost::gcd(n[0],n[1]),boost::gcd(n[2],c)); // is that right?
      CCTBX_ASSERT( g>0 );
      CCTBX_ASSERT( c%g == 0 && n[0]%g==0 && n[1]%g==0 && n[2]%g==0 );
      n /= g;
      c /= g;
    }

    //! Constructs a from normal coordinates and a rational constant
    cut(int_type x, int_type y, int_type z, rational_t c_)
    { 
      new(this) cut(ivector3_t(x,y,z), c_);
    }

    cut() {}

  }; // class cut


}}} // namespace cctbx::sgtbx::asu
#endif

