#ifndef MMTBX_DYNAMICS_DYNAMICS_H
#define MMTBX_DYNAMICS_DYNAMICS_H

#include <mmtbx/import_scitbx_af.h>
#include <mmtbx/error.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/vec3.h>

namespace mmtbx { namespace dynamics {

using scitbx::vec3;

template <typename FloatType>
FloatType
kinetic_energy(
  af::const_ref<vec3<FloatType> > const& velocities,
  af::const_ref<FloatType> const& masses)
{
  MMTBX_ASSERT(masses.size() == velocities.size());
  FloatType result = 0;
  for (std::size_t i=0;i<velocities.size();i++) {
    result += masses[i] * velocities[i].length_sq();
  }
  result *= 0.5;
  return result;
}

class center_of_mass_info {
public:
   center_of_mass_info(vec3<double> center_of_mass,
                       af::shared<vec3<double> > sites_cart,
                       af::shared<vec3<double> > velocities,
                       af::shared<double> const& weights);
   double ekcm() const { return ekcm_; }
   vec3<double> vcm() const { return vcm_; }
   vec3<double> acm() const { return acm_; }
private:
   double ekcm_;
   vec3<double> vcm_;
   vec3<double> acm_;
};

void vxyz_at_t_plus_dt_over_2(af::shared<vec3<double> > vxyz,
                            af::shared<double> const& weights,
                            af::shared<vec3<double> > const& grad,
                            double tstep);

af::shared<vec3<double> > stop_center_of_mass_motion(
                            vec3<double> center_of_mass,
                            vec3<double> acm,
                            vec3<double> vcm,
                            af::shared<vec3<double> > sites_cart,
                            af::shared<vec3<double> > velocities,
                            af::shared<double> const& weights);

}} // namespace mmtbx::dynamics

#endif // MMTBX_DYNAMICS_DYNAMICS_H
