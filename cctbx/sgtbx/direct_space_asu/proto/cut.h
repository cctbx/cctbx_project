#ifndef CCTBX_SGTBX_ASU_CUT_H
#define CCTBX_SGTBX_ASU_CUT_H

#include <iosfwd>
#include <vector>
#include <set>
#include <limits>

#include <cmath>

#include <boost/rational.hpp>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/simple_tiny_io.h>
#include <scitbx/array_family/simple_io.h>
#include <scitbx/mat3.h>
#include <cctbx/sgtbx/space_group.h>
#include <cctbx/sgtbx/symbols.h>
#include <cctbx/sgtbx/change_of_basis_op.h>
#include <cctbx/sgtbx/space_group_type.h>


namespace cctbx { namespace sgtbx {

//! Direct space asymmetric unit classes
namespace asu {

  typedef long int_type;
  typedef scitbx::af::tiny<int_type,4> int4_t;
  typedef scitbx::vec3<int_type> int3_t;
  typedef scitbx::mat3<int_type> mat3_t;
  // typedef sg_vec3::value_type int_type;

  namespace {
    // sanity checks
    BOOST_STATIC_ASSERT( std::numeric_limits<int_type>::is_integer );
    BOOST_STATIC_ASSERT( std::numeric_limits<int_type>::digits >= 31 );
    BOOST_STATIC_ASSERT( std::numeric_limits<mat3_t::value_type>::is_integer );
    BOOST_STATIC_ASSERT( sizeof(mat3_t::value_type) == sizeof(int_type) );
  }

  typedef unsigned short size_type;
  typedef boost::rational<int> rational_t;

  typedef scitbx::vec3< rational_t > rvector3_t;
  typedef std::set< rvector3_t > set_rvector3_t;

  using cctbx::sgtbx::change_of_basis_op;

}}} // namespace cctbx::sgtbx::asu


namespace boost {

  // mrt: WHY DO I HAVE TO DEFINE this operator in boost:: { } ?????????????????
  //
  inline bool operator< (const cctbx::sgtbx::asu::rvector3_t &a,
    const cctbx::sgtbx::asu::rvector3_t &b)
  {
    return (a[0]<b[0]) || (a[0]==b[0] && a[1]<b[1]) || (a[0]==b[0] && a[1]==b[1]
      && a[2]<b[2]);
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
    int3_t n;
    //! Plane constant
    int_type c;
    //! Flag indicating if plane points belong to the asymmetric unit
    bool inclusive;

    void print(std::ostream &os ) const;
    std::string as_string() const
    {
      std::stringstream buf;
      this->print(buf);
      return buf.str();
    }

    //! Returns an arbitrary point, which belongs to the plane
    void get_point_in_plane(rvector3_t &r) const;

    //! Applies basis change to the plane
    void change_basis(const change_of_basis_op &cb_op);

    // TODO: do not use. Is it the same as change_basis?
    void apply_symop(const rt_mx &symop);

    double get_tolerance(const scitbx::af::double3 &tol3) const
    {
      return std::fabs(n[0]*tol3[0]) + std::fabs(n[1]*tol3[1])
        + std::fabs(n[2]*tol3[2]);
    }

    //! Evaluates plane equation
    rational_t evaluate(const rvector3_t &p) const
    {
      return  n[0]*p[0] + n[1]*p[1] + n[2]*p[2] + c;
    }

    long evaluate_int(const scitbx::af::int3 &num,
      const scitbx::af::int3 &den) const
    {
      return // this limits grid size to about (max(long)/4)^(1/3)
        // 812 for 32 bit systems, 1,321,122 for 64bit systems
        static_cast<long>(num[0])*n[0]*static_cast<long>(den[1])*den[2]
        + static_cast<long>(num[1])*n[1]*static_cast<long>(den[0])*den[2]
        + static_cast<long>(num[2])*n[2]*static_cast<long>(den[0])*den[1]
        + static_cast<long>(c)*den[0]*static_cast<long>(den[1])*den[2];
    }

    //! Evaluate plane. Use with optimized plane only
    long evaluate_int(const scitbx::af::int3 &num) const
    {
      return
          static_cast<long>(num[0])*n[0]
        + static_cast<long>(num[1])*n[1]
        + static_cast<long>(num[2])*n[2]
        + static_cast<long>(c);
    }

    //! Optimize plane for use with grid of size grid_size
    void optimize_for_grid(const scitbx::af::int3 &grid_size)
    {
      // Actual transformation:
      //   n[0]*= (static_cast<long>(grid_size[1])*grid_size[2]);
      //   n[1]*= (static_cast<long>(grid_size[2])*grid_size[0]);
      //   n[2]*= (static_cast<long>(grid_size[0])*grid_size[1]);
      //   c*= (grid_size[0]*static_cast<long>(grid_size[1])*grid_size[2]);
      // It is applied in a way reducing the risk of overflows
      BOOST_STATIC_ASSERT( sizeof(int_type) >= sizeof(long) );
      std::ostringstream errmsg;
      errmsg << "Integer overflow. Grid: " << grid_size << ",  asu cut: ";
      this->print(errmsg);
      int_type g = boost::gcd(boost::gcd(grid_size[0],grid_size[1]),
          grid_size[2]);
      CCTBX_ASSERT( g>0 );
      long gsz[3];
      for(unsigned char i=0; i<3U; ++i)
      {
        CCTBX_ASSERT( grid_size[i]%g == 0 );
        gsz[i] = grid_size[i]/g;
      }
      long szf[4];
      const double max_int = std::numeric_limits<long>::max()-3;
      if( (static_cast<double>(gsz[1])*gsz[2] > max_int)
          ||(static_cast<double>(gsz[0])*gsz[2] > max_int)
          ||(static_cast<double>(gsz[1])*gsz[0] > max_int)
        )
        throw error(errmsg.str());
      szf[0] = gsz[1]*gsz[2];
      szf[1] = gsz[0]*gsz[2];
      szf[2] = gsz[0]*gsz[1];
      g = boost::gcd(boost::gcd(szf[0],szf[1]),szf[2]);
      CCTBX_ASSERT( g>0 );
      for(unsigned char i=0; i<3U; ++i)
      {
        CCTBX_ASSERT( szf[i]%g == 0 );
        szf[i] /= g;
      }
      if( static_cast<double>(szf[2])*grid_size[2] > max_int )
        throw error(errmsg.str());
      szf[3] = szf[2]*grid_size[2];

      for(unsigned char i=0; i<3U; ++i)
      {
        if( static_cast<double>(n[i])*szf[i] > max_int )
          throw error(errmsg.str());
        n[i] *= szf[i];
      }
      if( static_cast<double>(c)*szf[3] > max_int )
        throw error(errmsg.str());
      c *= szf[3];
      this->normalize();
    }

    void get_optimized_grid_limits(scitbx::af::long3 &max_p) const
    {
      const long C = std::numeric_limits<long>::max()-3;
      const long Cc = C - std::abs(c);
      CCTBX_ASSERT( C>0 && Cc>0 );
      unsigned char nnz=0;
      for(unsigned char i=0; i<3U; ++i)
      {
        if( n[i]!=0 )
          ++nnz;
      }
      CCTBX_ASSERT( nnz>0U && nnz<=3U );
      for(unsigned char i=0; i<3U; ++i)
      {
        max_p[i] = (n[i]!=0) ? (Cc/nnz/std::abs(n[i])) : C;
        CCTBX_ASSERT( max_p[i]>=0 );
      }
    }

    double evaluate_double(const scitbx::af::double3 &p) const
    {
      return n[0]*p[0] + n[1]*p[1] + n[2]*p[2] + c;
    }

    template<typename TR>
      bool is_inside(const rvector3_t &p, const TR &expr) const
    {
      rational_t v = evaluate(p);
      if( v>0 )
        return true;
      if( v<0 )
        return false;
      return expr.is_inside(p);
    }

    template<typename TR>
      bool is_inside(const scitbx::af::int3 &num, const scitbx::af::int3 &den,
        const TR &expr) const
    {
      long v = evaluate_int(num,den);
      if( v>0 )
        return true;
      if( v<0 )
        return false;
      return expr.is_inside(num,den);
    }

    template<typename TR>
      bool is_inside(const scitbx::af::int3 &num, const TR &expr) const
    {
      long v = evaluate_int(num);
      if( v>0 )
        return true;
      if( v<0 )
        return false;
      return expr.is_inside(num);
    }


    /*
    template<typename TR>
      bool is_inside(const scitbx::af::double3 &p, double dtol, const TR &expr) const
    {
      double v = evaluate_double(p);
      // double dtol = std::fabs( tol[0]*n[0]+tol[1]*n[1]+tol[2]*n[2] );
      if( v>dtol )
        return true;
      if( v<-dtol )
        return false;
      return expr.is_inside(p,dtol);
    } */


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

    bool is_inside(const scitbx::af::int3 &num,
      const scitbx::af::int3 &den) const
    {
      long v = evaluate_int(num,den);
      if( v>0 )
        return true;
      if( v<0 )
        return false;
      return inclusive;
    }

    bool is_inside(const scitbx::af::int3 &num) const
    {
      long v = evaluate_int(num);
      if( v>0 )
        return true;
      if( v<0 )
        return false;
      return inclusive;
    }

    //! Returns 1 if point is inside, -1 if inside on the face, 0 otherwise
    short where_is(const scitbx::af::int3 &num,
      const scitbx::af::int3 &den) const
    {
      long v = evaluate_int(num,den);
      if( v>0 )
        return 1;
      if( v<0 )
        return 0;
      return inclusive?-1:0;
    }

    template<typename TR>
      short where_is(const scitbx::af::int3 &num, const scitbx::af::int3 &den,
        const TR &expr) const
    {
      long v = evaluate_int(num,den);
      if( v>0 )
        return 1;
      if( v<0 )
        return 0;
      return expr.is_inside(num,den) ? -1:0;
    }

    short where_is(const scitbx::af::int3 &num) const
    {
      long v = evaluate_int(num);
      if( v>0 )
        return 1;
      if( v<0 )
        return 0;
      return inclusive?-1:0;
    }

    template<typename TR>
      short where_is(const scitbx::af::int3 &num, const TR &expr) const
    {
      long v = evaluate_int(num);
      if( v>0 )
        return 1;
      if( v<0 )
        return 0;
      return expr.is_inside(num) ? -1:0;
    }


    bool is_inside_volume_only(const scitbx::af::double3 &p, double dtol) const
    {
      return evaluate_double(p) >= -dtol;
    }


    template<typename TR>
    cut_expression<TR> operator() (const TR &o) const
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
    cut(const int3_t &n_, rational_t c_, bool inc_=true) : inclusive(inc_)
    {
      CCTBX_ASSERT( c_.denominator() > 0 );
      n = n_ * static_cast<int_type>( c_.denominator() );
      c = c_.numerator();
      this->normalize();
    }

    //! Constructs a plane from normal coordinates and a rational constant
    cut(int_type x, int_type y, int_type z, rational_t c_, bool inc_ = true)
    {
      new(this) cut(int3_t(x,y,z), c_, inc_);
    }

    //! Constructs a plane from a normal vector and a rational constant
    cut(const scitbx::vec3<int> &n_, rational_t c_, bool inc_=true)
      : inclusive(inc_)
    {
      new(this) cut(int3_t(n_), c_, inc_);
    }

    cut() {}

  private:

    void normalize()
    {
      // is that right?
      const int_type g = boost::gcd(boost::gcd(n[0],n[1]),boost::gcd(n[2],c));
      CCTBX_ASSERT( g>0 );
      CCTBX_ASSERT( c%g == 0 && n[0]%g==0 && n[1]%g==0 && n[2]%g==0 );
      n /= g;
      c /= g;
    }

  }; // class cut


}}} // namespace cctbx::sgtbx::asu
#endif

