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

void vxyz_at_t_plus_dt_over_2(af::shared<vec3<double> > vxyz,
                            af::shared<double> const& weights,
                            af::shared<vec3<double> > const& grad,
                            double tstep);

}} // namespace mmtbx::dynamics

#endif // MMTBX_DYNAMICS_DYNAMICS_H
