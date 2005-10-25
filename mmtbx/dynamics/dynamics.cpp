#include <assert.h>
#include <math.h>
#include <iostream>
#include <mmtbx/dynamics/dynamics.h>
#include <mmtbx/error.h>

using namespace std;
namespace mmtbx { namespace dynamics {

kinetic_energy_and_temperature::kinetic_energy_and_temperature(
                                  af::shared<vec3<double> > const& vxyz,
                                  af::shared<double> const& weights)
{
    double k_boltz = 1.380662e-03;
    ekin = 0.0;
    int ndegf = 0;
    for (std::size_t i=0; i < vxyz.size(); i++) {
      ndegf += 3;
      vec3<double> const& v = vxyz[i];
      ekin += weights[i] * (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
     }
    ekin *= 0.5;
    temp = 2.0 * ekin / (ndegf * k_boltz);
}

void vxyz_at_t_plus_dt_over_2(af::shared<vec3<double> > vxyz,
                              af::shared<double> const& weights,
                              af::shared<vec3<double> > const& grad,
                              double tstep)
{
    for (std::size_t i=0; i < weights.size(); i++) {
      double factor = tstep / weights[i];
      vec3<double> v = vxyz[i];
      vec3<double> const& g = grad[i];
      v[0] -= g[0] * factor;
      v[1] -= g[1] * factor;
      v[2] -= g[2] * factor;
      vxyz[i] = v;
    }
}


}} // namespace mmtbx::dynamics
