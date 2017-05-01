#ifndef MMTBX_BUILDING_LOOP_CLOSURE_CCD_H
#define MMTBX_BUILDING_LOOP_CLOSURE_CCD_H

#include <mmtbx/validation/ramachandran/rama_eval.h>
#include <iotbx/pdb/hierarchy.h>
#include <iostream>
#include <math.h>
#include <scitbx/constants.h>
#include <stdlib.h>

namespace mmtbx { namespace building { namespace loop_closure {

// using namespace mmtbx::validation::ramachandran;
// using namespace scitbx::array_family;
// using namespace scitbx::af;
// using namespace iotbx::pdb::hierarchy;


class ccd_cpp {
  private:
    // constructor arguments
    scitbx::af::shared<iotbx::pdb::hierarchy::atom> moving_h_atoms;


  template <typename FloatType>
  scitbx::vec3<FloatType>
  project_point_on_axis(
    scitbx::vec3<FloatType> axis_point_1,
    scitbx::vec3<FloatType> axis_point_2,
    scitbx::vec3<FloatType> point)
  {
    scitbx::vec3<FloatType> ab = axis_point_2 - axis_point_1;
    scitbx::vec3<FloatType> ap = point - axis_point_1;
    scitbx::vec3<FloatType> result = axis_point_1 + ((ap*ab) / (ab*ab)) * ab;
    return result;
  }

  public:
  double
  _find_angle(
      scitbx::vec3<double> axis_point_1,
      scitbx::vec3<double> axis_point_2)
  {
    scitbx::af::tiny<scitbx::vec3<double>, 3 > f_all, s_home_all, r_all, r_home_all;
    double b=0, c=0;
    for (int i=0; i<3; i++) {
      scitbx::vec3<double> fixed_coor = fixed_ref_atoms[i];
      scitbx::vec3<double> moving_coor = moving_h_atoms[moving_ref_atoms_iseqs[i]].data->xyz;
      // _get_f_r_s
      scitbx::vec3<double> f, s_home(0,0,0), r_home(0,0,0);
      double r_norm;
      scitbx::vec3<double> fc_proj, mc_proj, r, ap_21;
      fc_proj = project_point_on_axis(axis_point_1, axis_point_2, fixed_coor);
      mc_proj = project_point_on_axis(axis_point_1, axis_point_2, moving_coor);
      f = fixed_coor - fc_proj;
      r = moving_coor - mc_proj;
      ap_21 = axis_point_2 - axis_point_1;
      r_norm = r.length();
      // std::cout << "r.length() is godd";
      if (fabs(r_norm) > 1e-10) r_home = r.normalize();
      // double ap_21_norm = ap_21.length(); // not needed?
      scitbx::vec3<double> theta_home(0,0,0);
      if (fabs(ap_21.length()) > 1e-10) theta_home = ap_21.normalize();
      scitbx::vec3<double> tt = theta_home.cross(r_home);
      if (fabs(tt.length()) > 1e-10) s_home = tt.normalize();


      b += 2*r_norm*(f*r_home);
      c += 2*r_norm*(f*s_home);

      // f_all[i] = f;
      // s_home_all[i] = s_home;
      // r_all[i] = r_norm;
      // r_home_all[i] = r_home;
    }
    double znam = std::sqrt(b*b+c*c);
    double sin_alpha = c/znam;
    double cos_alpha = b/znam;
    double alpha = atan2(sin_alpha, cos_alpha);

    return scitbx::rad_as_deg(alpha);
  }

  double _modify_angle(double angle)
  {
    double threshold = 1;
    if (std::abs(angle) > threshold) {
      if (angle > 0) return threshold;
      else return -threshold;
    }
    else return angle;
  }

    // others
    scitbx::af::tiny<size_t, 3> moving_ref_atoms_iseqs;
    scitbx::af::tiny<scitbx::vec3<double>, 3 > fixed_ref_atoms;
    iotbx::pdb::hierarchy::root moving_h;
    mmtbx::validation::ramachandran::rama_eval r;
    double convergence_diff;
    double resulting_rmsd;
    double needed_rmsd;
    bool early_exit;
    int max_number_of_iterations;


  ccd_cpp(
      scitbx::af::tiny<scitbx::vec3<double>, 3 > fixed_ref_atoms_,
      iotbx::pdb::hierarchy::root & moving_h_,
      scitbx::af::tiny<size_t, 3> moving_ref_atoms_iseqs_,
      const int& max_number_of_iterations_=500,
      const double& needed_rmsd_=0.1)
  :
    fixed_ref_atoms(fixed_ref_atoms_),
    moving_h(moving_h_),
    moving_ref_atoms_iseqs(moving_ref_atoms_iseqs_),
    max_number_of_iterations(max_number_of_iterations_),
    needed_rmsd(needed_rmsd_)
  {
    r = mmtbx::validation::ramachandran::rama_eval();
    moving_h_atoms = moving_h.atoms();
    convergence_diff = 1e-5;

  };


};

}}} // namespace mmtbx::building::loop_closure

#endif // MMTBX_BUILDING_LOOP_CLOSURE_CCD_H
