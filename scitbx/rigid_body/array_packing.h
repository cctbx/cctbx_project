#ifndef SCITBX_RIGID_BODY_ARRAY_PACKING_H
#define SCITBX_RIGID_BODY_ARRAY_PACKING_H

#include <scitbx/rigid_body/body_t.h>
#include <scitbx/array_family/shared.h>
#include <boost/numeric/conversion/cast.hpp>

namespace scitbx { namespace rigid_body {

//! Array packing and unpacking, mainly to aid Python interfaces.
namespace array_packing {

  template <typename ElementType, std::size_t N>
  af::shared<af::tiny<ElementType, N> >
  unpack_ref_tiny(
    af::const_ref<ElementType> const& packed,
    std::size_t result_size)
  {
    SCITBX_ASSERT(
      packed.size() == (packed.begin() == 0 ? 0 : result_size * N));
    af::shared<af::tiny<ElementType, N> > result;
    if (packed.begin() != 0) {
      result.resize(result_size);
      unsigned j = 0;
      for(std::size_t i=0;i<result_size;i++,j+=N) {
        std::copy(
          &packed[j],
          &packed[j+N], result[i].begin());
      }
    }
    return result;
  }

  template <typename FloatType>
  af::shared<af::small<FloatType, 6> >
  unpack_ref_small_6(
    af::const_ref<shared_ptr<body_t<FloatType> > > const& bodies,
    unsigned degrees_of_freedom,
    af::const_ref<FloatType> const& packed)
  {
    SCITBX_ASSERT(
      packed.size() == (packed.begin() == 0 ? 0 : degrees_of_freedom));
    typedef FloatType ft;
    af::shared<af::small<ft, 6> > result;
    if (packed.begin() != 0) {
      unsigned nb = boost::numeric_cast<unsigned>(bodies.size());
      result.reserve(nb);
      unsigned j = 0;
      for(unsigned ib=0;ib<nb;ib++) {
        body_t<ft> const* body = bodies[ib].get();
        unsigned n = body->joint->degrees_of_freedom;
        result.push_back(
          af::small<ft, 6>(af::adapt(
            af::const_ref<ft>(&packed[j], n))));
        j += n;
      }
      SCITBX_ASSERT(j == degrees_of_freedom);
    }
    return result;
  }

}}} // namespace scitbx::rigid_body::array_packing

#endif // GUARD
