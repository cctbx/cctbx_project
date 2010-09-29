#include <iostream>

#include "cut.h"

namespace cctbx { namespace sgtbx { namespace asu {

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

