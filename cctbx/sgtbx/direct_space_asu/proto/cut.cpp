#include <iostream>

#include "cut.h"

namespace cctbx { namespace sgtbx { namespace asu {

  // translated from cut:as_xyz in cctbx/sgtbx/direc_space_asu/cut_play.py
  void cut::print_as_xyz(std::ostream &ostr ) const
  {
    int n_negative = 0;
    int n_non_zero = 0;
    for(short j=0; j<3; ++j)
    {
      if (n[j] < 0)
        n_negative += 1;
      if (n[j] != 0)
        n_non_zero += 1;
    }
    int f = 0;
    if (n_non_zero == 1)
    {
      if (n_negative == 0)
        f = 1;
      else
        f = -1;
    }
    else
    {
      if (this->c > 0)
        n_negative += 1;
      if (n_negative <= n_non_zero/2)
        f = 1;
      else
        f = -1;
    }
    int_type g = boost::gcd(n[0], boost::gcd(n[1],n[2]));
    std::ostringstream s;
    char x[]="xyz";
    for(short j=0; j<3; ++j)
    {
      rational_t nf(this->n[j]*f,g);
      if (nf == 0)
        continue;
      if (nf > 0)
        s << '+';
      if (abs(nf) != 1)
      {
        if( nf.denominator()==1 )
          s << nf.numerator();
        else
          s << nf;
        s << '*' << x[j];
      }
      else
      {
        if (nf < 0)
          s << '-';
        s << x[j];
      }
    }
    if (f > 0)
      s << '>';
    else
      s << '<';
    if (this->inclusive)
      s << '=';
    rational_t cc(-this->c*f,g);
    if( cc.denominator()==1 )
      s << cc.numerator();
    else
      s << cc;
    // s << (-rational_t(this->c*f,g));
    std::string str = s.str();
    if (str[0] == '+')
      str = str.substr(1);
    // if (this->has_cuts()):
    //  s += " [" + self.cut_expr.as_xyz() + "]"
    //return s
    ostr << str;
  }

  void cut::print(std::ostream &os ) const
  {
    int_type g = boost::gcd(n[0], boost::gcd(n[1],n[2]));
    CCTBX_ASSERT( g>0 );
    if( !inclusive )
      os << "+";
    os << "cut((" << n[0]/g << "," << n[1]/g << "," << n[2]/g << "), "
      <<  rational_t(c,g) << ")";
  }

  void cut::get_point_in_plane(rvector3_t &r) const
  {
    r = rvector3_t(0,0,0);
    for(short i=0; i<3; ++i)
    {
      if( n[i]!=0 )
      {
        r[i] = rational_t(-c,n[i]);
        return;
      }
    }
    throw cctbx::error("cut_plane normal vector is null vector");
  }

  namespace {
    inline rational_t dot(const tr_vec &lhs, const tr_vec &rhs)
    {
      return rational_t( lhs.num()*rhs.num(), lhs.den()*rhs.den() );
    }
  }

  void cut::change_basis(const change_of_basis_op &cb_op)
  {
    CCTBX_ASSERT( this->n.length_sq()!= 0 );
    rot_mx r_inv_transpose( cb_op.c_inv().r().transpose() );
    tr_vec np = r_inv_transpose * tr_vec(cctbx::sg_vec3(this->n), 1);
    tr_vec t = cb_op.c().t();
    rational_t cp = rational_t(this->c) - dot(np,t);
    CCTBX_ASSERT( np.den()>0 );
    new(this) cut(int3_t(np.num()), cp*np.den(), inclusive);
  }

  // TODO: is it the same as change_basis ?
  void cut::apply_symop(const rt_mx &symop)
  {
    CCTBX_ASSERT( this->n.length_sq()!= 0 );
    rot_mx r_trans_inv( symop.r().transpose().inverse() );
    tr_vec np = r_trans_inv * tr_vec(cctbx::sg_vec3(this->n), 1);
    tr_vec t = symop.t();
    rational_t cp = rational_t(this->c) - dot(np,t);
    CCTBX_ASSERT( np.den()>0 );
    new(this) cut(int3_t(np.num()), cp*np.den(), inclusive);
  }

}}}
