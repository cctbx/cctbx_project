#ifndef CCTBX_SGTBX_GRID_SYMOP_H
#define CCTBX_SGTBX_GRID_SYMOP_H

#include <scitbx/array_family/tiny_plain.h>
#include <scitbx/array_family/tiny_types.h>
#include <cctbx/sgtbx/rt_mx.h>

namespace cctbx { namespace sgtbx {

  class grid_symop : public af::tiny_plain<int, 3*4>
  {
    public:

    int& operator()(unsigned char r, unsigned char c)
    {
      return this->elems[r*4+c];
    }

    const int& operator()(unsigned char r, unsigned char c) const
    {
      return this->elems[r*4+c];
    }

    grid_symop(const cctbx::sgtbx::rt_mx &symop, const scitbx::af::int3 &grid)
    {
      const cctbx::sgtbx::rot_mx &rot = symop.r();
      const cctbx::sgtbx::tr_vec &tr = symop.t();
      const int rd = rot.den();
      const int td = tr.den();
      for(unsigned char r=0; r<3; ++r)
      {
        for(unsigned char c=0; c<3; ++c)
        {
          int tmp = rot(r,c);
          SCITBX_ASSERT( tmp % rd == 0 );
          tmp = (tmp/rd)*grid[r];
          if( tmp % grid[c] != 0 )
            throw cctbx::error("The grid is not compatible with the spacegroup");
          (*this)(r,c) = tmp/grid[c];
        }
        int tmp = tr[r]*grid[r];
        if( tmp % td != 0 )
          throw cctbx::error("The grid is not compatible with the spacegroup");
        (*this)(r,3) = tmp/td;
      }
    }

    af::int3 apply_to(const af::int3 &rhs) const
    {
      af::int3 result;
      for(unsigned char r=0; r<3; ++r )
      {
        result[r] = (*this)(r,0) * rhs[0] + (*this)(r,1) * rhs[1] + (*this)(r,2) * rhs[2] + (*this)(r,3);
      }
      return result;
    }

  };

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_GRID_SYMOP_H

