#ifndef SCITBX_MATH_DIHEDRAL_H
#define SCITBX_MATH_DIHEDRAL_H

#include <scitbx/vec3.h>
#include <scitbx/constants.h>
#include <tbxx/error_utils.hpp>

namespace scitbx { namespace math {

  struct dihedral
  {
    vec3<double> d_01;
    vec3<double> d_21;
    vec3<double> d_23;
    vec3<double> n_0121;
    vec3<double> n_2123;
    double n_0121_norm;
    double n_2123_norm;

    //! Default constructor. Some data members are not initialized!
    dihedral() {}

    dihedral(
      af::tiny<vec3<double>, 4> const& sites)
    {
      init(sites.begin());
    }

    dihedral(
      af::const_ref<vec3<double> > const& sites)
    {
      TBXX_ASSERT(sites.size() == 4);
      init(sites.begin());
    }

    void
    init(
      vec3<double> const* sites)
    {
      d_01 = sites[0] - sites[1];
      d_21 = sites[2] - sites[1];
      d_23 = sites[2] - sites[3];
      n_0121 = d_01.cross(d_21);
      n_0121_norm = n_0121.length_sq();
      n_2123 = d_21.cross(d_23);
      n_2123_norm = n_2123.length_sq();
    }

    boost::optional<double>
    angle(
      bool deg=false)
    {
      if (n_0121_norm == 0 || n_2123_norm == 0) {
        return boost::optional<double>();
      }
      double cos_angle = std::max(-1.,std::min(1.,
        n_0121 * n_2123 / std::sqrt(n_0121_norm * n_2123_norm)));
      double result = std::acos(cos_angle);
      if (d_21 * (n_0121.cross(n_2123)) < 0) {
        result *= -1;
      }
      if (deg) result /= constants::pi_180;
      return boost::optional<double>(result);
    }
  };

}} // namespace scitbx::math

#endif // GUARD
