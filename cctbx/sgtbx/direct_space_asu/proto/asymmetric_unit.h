#ifndef CCTBX_SGTBX_ASU_ASYMMETRIC_UNIT_H
#define CCTBX_SGTBX_ASU_ASYMMETRIC_UNIT_H

#include <boost/config.hpp>

#include <cctbx/sgtbx/direct_space_asu/proto/direct_space_asu.h>

namespace cctbx { namespace sgtbx { namespace asu {

#if (2*2==4) || defined(BOOST_NO_CXX11_SCOPED_ENUMS)
namespace space { enum id {direct, reciprocal}; }
namespace optimization { enum id {unoptimized, optimized }; }
using namespace space;
using namespace optimization;
template<space::id S, optimization::id O=unoptimized> class asymmetric_unit;
#else
enum class optimization : short {unoptimized, optimized };
enum class space : short {direct, reciprocal };

using space::direct;  this fails with gcc on linux
using space::reciprocal;
using optimization::unoptimized;
using optimization::optimized;
template<space S, optimization O=unoptimized> class asymmetric_unit;
#endif

template<> class asymmetric_unit<direct,optimized> : private detail::faces
{
  typedef detail::faces base_t;

public:
  typedef scitbx::int3 grid_size_t;

  asymmetric_unit(const direct_space_asu &a, const scitbx::af::int3 &g)
    : base_t(a.get_faces()), grid_size_(g)
  {
    faces_->optimize_for_grid(grid_size_);
  }

  const grid_size_t &grid_size() const {return grid_size_;}

  scitbx::af::long3 get_optimized_grid_limits() const
  {
    scitbx::af::long3 max_p;
    faces_->get_optimized_grid_limits(max_p);
    return max_p;
  }

  //! Tests where num/grid_size is: fully inside, fully outside or on the face
  /*! Returns 1 : inside, 0 : outside, -1 : on the face
   * Must be used only with asymmetric unit optimized_for_grid grid_size.
   */
  short where_is(const scitbx::int3 &num) const
  {
    return faces_->where_is(num);
  }

private:

  const grid_size_t grid_size_;
};


}}}
#endif
