#ifndef CCTBX_SGTBX_DIRECT_SPACE_ASU_H
#define CCTBX_SGTBX_DIRECT_SPACE_ASU_H

#include <scitbx/vec3.h>
#include <scitbx/array_family/small.h>
#include <cctbx/import_scitbx_af.h>

namespace cctbx { namespace sgtbx { namespace direct_space_asu {

  template <typename FloatType=double>
  class float_cut_plane
  {
    public:
      float_cut_plane() {}

      float_cut_plane(
        scitbx::vec3<FloatType> const& n_,
        FloatType const& c_)
      :
        n(n_),
        c(c_)
      {}

      scitbx::vec3<FloatType> n;
      FloatType c;

      FloatType
      evaluate(scitbx::vec3<FloatType> const& point)
      {
        return n * point + c;
      }

      bool
      is_inside(scitbx::vec3<FloatType> const& point)
      {
        if (evaluate(point) < 0) return false;
        return true;
      }
  };

  template <typename FloatType=double>
  class float_asu
  {
    public:
      typedef af::small<float_cut_plane<FloatType>, 9> facets_t;

      float_asu() {}

      float_asu(facets_t const& facets_)
      :
        facets(facets_)
      {}

      facets_t facets;

      bool
      is_inside(scitbx::vec3<FloatType> const& point)
      {
        for(std::size_t i=0;i<facets.size();i++) {
          if (!facets[i].is_inside(point)) return false;
        }
        return true;
      }
  };

}}} // namespace cctbx::sgtbx::direct_space_asu

#endif // CCTBX_SGTBX_DIRECT_SPACE_ASU_H
