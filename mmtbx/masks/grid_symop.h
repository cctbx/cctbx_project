#ifndef CCTBX_SGTBX_GRID_SYMOP_H
#define CCTBX_SGTBX_GRID_SYMOP_H

#include <scitbx/array_family/tiny_plain.h>
#include <scitbx/array_family/tiny_types.h>
#include <cctbx/sgtbx/rt_mx.h>

namespace cctbx { namespace sgtbx {

  //! Symmetry operator optimized for grid
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

    //! Construct symmetry operator from symop and grid of size grid
    grid_symop(const cctbx::sgtbx::rt_mx &symop, const scitbx::af::int3 &grid)
    {
      const cctbx::sgtbx::rot_mx &rot = symop.r();
      const cctbx::sgtbx::tr_vec &tr = symop.t();
      const int rd = rot.den();
      const int td = tr.den();
      const double max_int = std::numeric_limits<int>::max()-3;
      std::string err1("Integer overflow. "),
        err2("The grid is not compatible with the spacegroup. ");
      {
        std::ostringstream str;
        str << "Symop: " << symop.as_xyz()
          << " on the grid: " << grid;
        err2 += str.str();
        str << ". Max int: " << max_int << ". May be grid is too large.";
        err1 += str.str();
      }
      for(unsigned char r=0; r<3; ++r)
      {
        for(unsigned char c=0; c<3; ++c)
        {
          int tmp = rot(r,c);
          SCITBX_ASSERT( tmp % rd == 0 );
          tmp /= rd;
          if( static_cast<double>(tmp)*grid[r] > max_int )
            throw error(err1);
          tmp = tmp*grid[r];
          if( tmp % grid[c] != 0 )
            throw error(err2);
          (*this)(r,c) = tmp/grid[c];
        }
        if( static_cast<double>(tr[r])*grid[r] > max_int )
          throw error(err1);
        int tmp = tr[r]*grid[r];
        if( tmp % td != 0 )
          throw error(err2);
        (*this)(r,3) = tmp/td;
      }
    }

    //! Apply symmetry operator to poit rhs/grid
    af::int3 apply_to(const af::int3 &rhs) const
    {
      af::int3 result;
      for(unsigned char r=0; r<3; ++r )
      {
        result[r] = (*this)(r,0) * rhs[0] + (*this)(r,1) * rhs[1] + (*this)(r,2)
          * rhs[2] + (*this)(r,3);
      }
      return result;
    }

    void get_grid_limits(scitbx::af::int3 &max_p) const
    {
      const int C = std::numeric_limits<int>::max()-3;
      CCTBX_ASSERT( C>0 );
      max_p[0] = max_p[1] = max_p[2] = C;
      const grid_symop &m = *this;
      for(unsigned char r=0; r<3U; ++r)
      {
        unsigned char nnz=0;
        for(unsigned char c=0; c<3U; ++c)
        {
          if( m(r,c)!=0 )
            ++nnz;
        }
        CCTBX_ASSERT( nnz>0U && nnz<=3U );
        const int Cc = C - std::abs(m(r,3));
        CCTBX_ASSERT( Cc>0 );
        for(unsigned char c=0; c<3U; ++c)
        {
          const int t = std::abs(m(r,c));
          const int mx = (t!=0) ? (Cc / nnz / t) : C;
          max_p[c] = std::min(max_p[c],mx);
          CCTBX_ASSERT( max_p[c]>=0 );
        }
      }
    }


  };

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_GRID_SYMOP_H
