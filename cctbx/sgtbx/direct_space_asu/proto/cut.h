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


namespace cctbx { namespace sgtbx { namespace asu {

  typedef unsigned short size_type;
  // typedef short int_type;
  typedef int int_type;
  typedef boost::rational<int_type> rational_t;

  // increase this to raise number of points on the facet
  // const rational_t zero = 0; // 1.0/( 1000.0*100.0 );

  typedef scitbx::mat3<int_type> imatrix3_t;
  typedef scitbx::vec3<int_type> ivector3_t;
  typedef scitbx::vec3<double> dvector3_t;
  typedef scitbx::vec3< rational_t > rvector3_t;
  typedef scitbx::mat3< rational_t > rmatrix3_t;
  typedef std::set< rvector3_t > set_rvector3_t;

  using cctbx::sgtbx::change_of_basis_op;

  inline int_type get_int(rational_t r)
  {
    return r.numerator()/r.denominator();
  }

  inline ivector3_t get_int(rvector3_t r)
  {
    return ivector3_t(get_int(r[0]), get_int(r[1]), get_int(r[2]));
  }

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


}}} // namespace 

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
  class cut
  {
  public:
    ivector3_t n;
    int_type c;
    bool inclusive;

    void print(std::ostream &os, bool x=false ) const;

    //! Returns an arbitrary point, which belongs to the plane
    void get_point_in_plane(rvector3_t &r) const;

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

    cut operator+ () const
    {
      CCTBX_ASSERT( inclusive );
      cut r(*this);
      r.inclusive = false;
      return r;
    }

    cut operator- () const
    {
      return cut(-n, -c, inclusive);
    }

    cut operator~ () const
    {
      cut r(*this);
      r.n = -n;
      return r;
    }

    cut operator* (int_type x) const
    {
      return cut(n, c*x, inclusive);
    }

    cut operator/ (int_type x) const
    {
      CCTBX_ASSERT( x!=0 && c!=0 );
      return cut(n*std::abs(x), c*x/std::abs(x), inclusive);
    }

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

    cut(int_type x, int_type y, int_type z, rational_t c_)
    { 
      new(this) cut(ivector3_t(x,y,z), c_);
    }

    cut() {}

  }; // class cut


}}} // namespace cctbx::sgtbx::asu
#endif

