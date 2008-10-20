#include <assert.h>
#include <math.h>
#include <iostream>
#include <mmtbx/dynamics/dynamics.h>
#include <mmtbx/error.h>

using namespace std;
namespace mmtbx { namespace dynamics {

namespace af=scitbx::af;
using scitbx::vec3;
using scitbx::mat3;
using scitbx::sym_mat3;

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

center_of_mass_info::center_of_mass_info(
                       vec3<double> center_of_mass,
                       af::shared<vec3<double> > sites_cart,
                       af::shared<vec3<double> > velocities,
                       af::shared<double> const& weights)
{
    double timfac = 0.04888821;
    double vxcm = 0.0;
    double vycm = 0.0;
    double vzcm = 0.0;
    double axcm = 0.0;
    double aycm = 0.0;
    double azcm = 0.0;
    double xcm = 0.0;
    double ycm = 0.0;
    double zcm = 0.0;
    double tmass = 0;
    for (std::size_t i=0; i < weights.size(); i++) {
      double weight = weights[i];
      vec3<double> velocity = velocities[i];
      vec3<double> site = sites_cart[i];
      tmass += weight;
      vxcm += velocity[0] * weight;
      vycm += velocity[1] * weight;
      vzcm += velocity[2] * weight;
      xcm += site[0] * weight;
      ycm += site[1] * weight;
      zcm += site[2] * weight;
      axcm += (site[1] * velocity[2] - site[2] * velocity[1]) * weight;
      aycm += (site[2] * velocity[0] - site[0] * velocity[2]) * weight;
      azcm += (site[0] * velocity[1] - site[1] * velocity[0]) * weight;
    }
    axcm -= (ycm * vzcm - zcm * vycm) / tmass;
    aycm -= (zcm * vxcm - xcm * vzcm) / tmass;
    azcm -= (xcm * vycm - ycm * vxcm) / tmass;
    vxcm /= tmass;
    vycm /= tmass;
    vzcm /= tmass;
    vcm_ = vec3<double> (vxcm,vycm,vzcm);
    acm_ = vec3<double> (axcm,aycm,azcm);
    ekcm_ = (vxcm*vxcm + vycm*vycm + vzcm*vzcm) * tmass * 0.5;
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

af::shared<vec3<double> > stop_center_of_mass_motion(
                            vec3<double> center_of_mass,
                            vec3<double> acm,
                            vec3<double> vcm,
                            af::shared<vec3<double> > sites_cart,
                            af::shared<vec3<double> > velocities,
                            af::shared<double> const& weights)
{
  double xx = 0.0;
  double xy = 0.0;
  double xz = 0.0;
  double yy = 0.0;
  double yz = 0.0;
  double zz = 0.0;
  mat3<double> tcm_inv;
  af::shared<vec3<double> > result(weights.size()) ;
  for (std::size_t i=0; i < weights.size(); i++) {
    vec3<double> ri = sites_cart[i] - center_of_mass;
    xx += ri[0]*ri[0] * weights[i];
    xy += ri[0]*ri[1] * weights[i];
    xz += ri[0]*ri[2] * weights[i];
    yy += ri[1]*ri[1] * weights[i];
    yz += ri[1]*ri[2] * weights[i];
    zz += ri[2]*ri[2] * weights[i];
  }
  tcm_inv = mat3<double>(yy+zz,-xy,-xz,  -xy,xx+zz,-yz,  -xz,-yz,xx+yy);
  double det = tcm_inv.determinant();
  if(det > 1.e-4) {
    tcm_inv = tcm_inv.inverse();
    // get angular velocity OXCM, OYCM, OZCM
    double oxcm = acm[0]*tcm_inv[0] + acm[1]*tcm_inv[3] + acm[2]*tcm_inv[6];
    double oycm = acm[0]*tcm_inv[1] + acm[1]*tcm_inv[4] + acm[2]*tcm_inv[7];
    double ozcm = acm[0]*tcm_inv[2] + acm[1]*tcm_inv[5] + acm[2]*tcm_inv[8];
    // remove CM translational and rotational motion from velocities
    for (std::size_t i=0; i < weights.size(); i++) {
      vec3<double> ri = sites_cart[i] - center_of_mass;
      double vx = velocities[i][0];
      vx += -vcm[0] - oycm*ri[2] + ozcm*ri[1];
      double vy = velocities[i][1];
      vy += -vcm[1] - ozcm*ri[0] + oxcm*ri[2];
      double vz = velocities[i][2];
      vz += -vcm[2] - oxcm*ri[1] + oycm*ri[0];
      result[i] = vec3<double> (vx,vy,vz);
    }
  }
  return result;
}

}} // namespace mmtbx::dynamics
