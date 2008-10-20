#ifndef MMTBX_DYNAMICS_DYNAMICS_H
#define MMTBX_DYNAMICS_DYNAMICS_H

#include <cctbx/sgtbx/space_group.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <vector>

using scitbx::vec3;
namespace af=scitbx::af;

namespace mmtbx { namespace dynamics {

class kinetic_energy_and_temperature {
public:
   kinetic_energy_and_temperature(af::shared<vec3<double> > const& vxyz,
                                  af::shared<double> const& weights);
   double kinetic_energy() const { return ekin; }
   double temperature() const { return temp; }
private:
   double ekin, temp;
};

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
