#include <iostream>

#include "cut.h"

namespace cctbx { namespace sgtbx { namespace asu {

  void cut::print(std::ostream &os) const
  {
    int_type g = boost::gcd(n[0], boost::gcd(n[1],n[2]));
    os << "n = (" << n[0]/g << ", " << n[1]/g << ", " << n[2]/g << ")   c= "
      <<  rational_t(c,g);
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
    throw "cut_plane normal vector is null vector";
  }

  void cut::change_basis(const change_of_basis_op &cb_op)
  {
    int_type n_norm_sq = n.length_sq();
    CCTBX_ASSERT( n_norm_sq != 0 );
    cctbx::sgtbx::tr_vec x( n*( -this->c) , n_norm_sq );
    rvector3_t v = rvector3_t(n);
    CCTBX_ASSERT( n*x.num()+this->c*x.den() == 0 );
    cctbx::sgtbx::rot_mx r_inv_transpose( cb_op.c_inv().r().transpose() );
    cctbx::sgtbx::tr_vec np = r_inv_transpose * cctbx::sgtbx::tr_vec(n, 1);
    cctbx::sgtbx::tr_vec xp = cb_op.c().r() * x + cb_op.c().t(); // this should fail?
    rational_t cp( -(np.num() * xp.num()), np.den()*xp.den() );
    new(this) cut(np.num(), cp*np.den());
  }

}}}

