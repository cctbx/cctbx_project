#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
//#include <boost/python/return_value_policy.hpp>
//#include <boost/python/return_by_value.hpp>
#include <boost/optional.hpp>

#include <mmtbx/error.h>
#include <cctbx/geometry_restraints/dihedral.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/accessors/c_grid.h>

#include <cmath>
#include <iostream>

namespace mmtbx { namespace ramachandran {
  namespace af = scitbx::af;

  double convert_angle (double theta) {
    if (theta > 180) {
      theta = -360 + theta;
    } else if (theta < -180) {
      theta = 360 + theta;
    }
    return theta;
  }

  class lookup_table
  {
    public:
      af::versa<double, af::flex_grid<> > plot;
      double values_max;

      lookup_table (
        af::const_ref< double > values,
        int n_angles,
        double scale_allowed=1.0)
      {
        MMTBX_ASSERT(values.size() == (n_angles * n_angles));
        MMTBX_ASSERT(scale_allowed > 0.0);
        af::flex_grid<>::index_type fg_origin;
        af::flex_grid<>::index_type fg_last;
        for (unsigned i = 0; i < 2; i++) {
          fg_origin.push_back(0);
          fg_last.push_back(n_angles);
        }
        plot = af::versa<double, af::flex_grid<> >(
          af::flex_grid<>(fg_origin, fg_last, true));
        values_max = 0.0;
        for (unsigned i = 0; i < values.size(); i++) {
          if (values[i] > 0.0) {
            plot[i] = values[i] * scale_allowed;
          } else {
            plot[i] = values[i];
          }
          if (values[i] > values_max) {
            values_max = values[i];
          }
        }
      }

      double get_score (
        double phi,
        double psi)
      {
        phi = convert_angle(phi);
        psi = convert_angle(psi);
        //std::cout << phi << " " << psi << "\n";
        MMTBX_ASSERT((phi <= 180.0) && (phi >= -180.0));
        MMTBX_ASSERT((psi <= 180.0) && (psi >= -180.0));
        int phi_1 = (int) floor(phi);
        int phi_2 = (int) ceil(phi);
        int psi_1 = (int) floor(psi);
        int psi_2 = (int) ceil(psi);
        if ((phi_1 % 2) == 0) {
          if (phi_2 == phi_1)
            phi_2 += 1;
          phi_1 -= 1;
        } else if ((phi_2 % 2) == 0) {
          phi_2 += 1;
        }
        if ((psi_1 % 2) == 0) {
          if (psi_2 == psi_1)
            psi_2 += 1;
          psi_1 -= 1;
        } else if ((psi_2 % 2) == 0) {
          psi_2 += 1;
        }
        double r_phi_psi = 0;
        if (phi_2 == phi_1) {
          if (psi_2 == psi_1) {
            r_phi_psi = get_point(phi_1, psi_1);
          } else {
            double r11 = get_point(phi_1, psi_1);
            double r12 = get_point(phi_1, psi_2);
            r_phi_psi = (((psi-psi_1)*r12)+((psi_2-psi)*r11)) / (psi_2-psi_1);
          }
        } else if (psi_2 == psi_1) {
          double r11 = get_point(phi_1, psi_1);
          double r21 = get_point(phi_2, psi_1);
          r_phi_psi = (((phi-phi_1)*r21)+((phi_2-phi)*r11)) / (phi_2-phi_1);
        } else {
          double r11 = get_point(phi_1, psi_1);
          double r12 = get_point(phi_1, psi_2);
          double r21 = get_point(phi_2, psi_1);
          double r22 = get_point(phi_2, psi_2);
          double d_phi_d_psi = (double) (phi_2 - phi_1) * (psi_2 - psi_1);
          MMTBX_ASSERT(d_phi_d_psi != 0);
          r_phi_psi = ((r11/d_phi_d_psi) * (phi_2-phi) * (psi_2-psi)) +\
                      ((r21/d_phi_d_psi) * (phi-phi_1) * (psi_2-psi)) +\
                      ((r12/d_phi_d_psi) * (phi_2-phi) * (psi-psi_1)) +\
                      ((r22/d_phi_d_psi) * (phi-phi_1) * (psi-psi_1));
        }
        return r_phi_psi;
      }

      // score inverted and adjusted to have a minimum value of 0
      double get_energy (
        double phi,
        double psi)
      {
        double score = get_score(phi, psi);
        return - (score - values_max);
      }

      double
      compute_gradients (
        af::ref<scitbx::vec3<double> > const& gradient_array,
        af::const_ref<scitbx::vec3<double> > const& sites_cart,
        af::tiny<unsigned, 5> const& i_seqs,
        double weight=1.0,
        double epsilon=0.1)
      {
        MMTBX_ASSERT(gradient_array.size() == sites_cart.size());
        MMTBX_ASSERT(epsilon > 0.0);
        using cctbx::geometry_restraints::dihedral;
        af::tiny<scitbx::vec3<double>, 4> phi_sites;
        af::tiny<scitbx::vec3<double>, 4> psi_sites;
        for (unsigned i = 0; i < 4; i++) {
          phi_sites[i] = sites_cart[i_seqs[i]];
          psi_sites[i] = sites_cart[i_seqs[i+1]];
        }
        dihedral phi(phi_sites, 0, 1.0);
        dihedral psi(psi_sites, 0, 1.0);
        double phi_deg = phi.angle_model;
        double psi_deg = psi.angle_model;
        double residual = get_energy(phi_deg, psi_deg);
        double r_phi_1 = get_energy(phi_deg - epsilon, psi_deg);
        double r_phi_2 = get_energy(phi_deg + epsilon, psi_deg);
        double d_r_d_phi = (r_phi_2 - r_phi_1) / (epsilon * 2);
        double r_psi_1 = get_energy(phi_deg, psi_deg - epsilon);
        double r_psi_2 = get_energy(phi_deg, psi_deg + epsilon);
        double d_r_d_psi = (r_psi_2 - r_psi_1) / (epsilon * 2);
        af::tiny<scitbx::vec3<double>, 4> d_phi_d_xyz = - phi.grad_delta();
        af::tiny<scitbx::vec3<double>, 4> d_psi_d_xyz = - psi.grad_delta();
        for (unsigned k = 0; k < 5; k++) {
          std::size_t i_seq = i_seqs[k];
          if (k < 4) {
            gradient_array[i_seq] += d_r_d_phi * d_phi_d_xyz[k] * weight;
          }
          if (k > 0) {
            gradient_array[i_seq] += d_r_d_psi * d_psi_d_xyz[k-1] * weight;
          }
        }
        return residual;
      }

      double
      compute_gradients_finite_differences (
        af::ref<scitbx::vec3<double> > const& gradient_array,
        af::const_ref<scitbx::vec3<double> > const& sites_cart,
        af::tiny<unsigned, 5> const& i_seqs,
        double weight=1.0,
        double epsilon=0.01)
      {
        MMTBX_ASSERT(gradient_array.size() == sites_cart.size());
        MMTBX_ASSERT(epsilon > 0.0);
        using cctbx::geometry_restraints::dihedral;
        af::tiny<scitbx::vec3<double>, 4> phi_sites;
        af::tiny<scitbx::vec3<double>, 4> psi_sites;
        for (unsigned k = 0; k < 4; k++) {
          phi_sites[k] = sites_cart[i_seqs[k]];
          psi_sites[k] = sites_cart[i_seqs[k+1]];
        }
        dihedral phi(phi_sites, 0, 1.0);
        dihedral psi(psi_sites, 0, 1.0);
        double residual = get_energy(phi.angle_model, psi.angle_model);
        for (unsigned k = 0; k < 5; k++) {
          std::size_t i_seq = i_seqs[k];
          for (unsigned u = 0; u < 3; u++) {
            if (k < 4) {
              phi_sites[k][u] -= epsilon;
            }
            if (k > 0) {
              psi_sites[k-1][u] -= epsilon;
            }
            dihedral phi1(phi_sites, 0, 1.0);
            dihedral psi1(psi_sites, 0, 1.0);
            if (k < 4) {
              phi_sites[k][u] += epsilon * 2;
            }
            if (k > 0) {
              psi_sites[k-1][u] += epsilon * 2;
            }
            dihedral phi2(phi_sites, 0, 1.0);
            dihedral psi2(psi_sites, 0, 1.0);
            if (k < 4) {
              phi_sites[k][u] -= epsilon;
            }
            if (k > 0) {
              psi_sites[k-1][u] -= epsilon;
            }
            double r1 = get_energy(phi1.angle_model, psi1.angle_model);
            double r2 = get_energy(phi2.angle_model, psi2.angle_model);
            gradient_array[i_seq][u] += weight * ((r2-r1) / (epsilon*2));
          }
        }
        return residual;
      }

    private :
      double get_point (
        int phi,
        int psi)
      {
        phi = (int) convert_angle(phi);
        psi = (int) convert_angle(psi);
        MMTBX_ASSERT((phi < 180) && (phi > -180));
        MMTBX_ASSERT((psi < 180) && (psi > -180));
        MMTBX_ASSERT((abs(phi % 2) == 1) && (abs(psi % 2) == 1));
        int i = (int) (phi + 179) / 2;
        int j = (int) (psi + 179) / 2;
        return plot(i,j);
      }
  };

  void init_module()
  {
    using namespace boost::python;
    class_<lookup_table>("lookup_table", no_init)
      .def(init<af::const_ref< double >,
                int,
                double>((
        arg("values"),
        arg("n_angles"),
        arg("scale_allowed")=1.0)))
      .def("get_score", &lookup_table::get_score, (
        arg("phi"),
        arg("psi")))
      .def("get_energy", &lookup_table::get_energy, (
        arg("phi"),
        arg("psi")))
      .def("compute_gradients", &lookup_table::compute_gradients, (
        arg("gradient_array"),
        arg("sites_cart"),
        arg("i_seqs"),
        arg("weight")=1.0,
        arg("epsilon")=0.1))
      .def("compute_gradients_finite_differences",
        &lookup_table::compute_gradients_finite_differences, (
        arg("gradient_array"),
        arg("sites_cart"),
        arg("i_seqs"),
        arg("weight")=1.0,
        arg("epsilon")=0.01));
  }

}} // namespace mmtbx::ramachandran

BOOST_PYTHON_MODULE(mmtbx_ramachandran_ext)
{
  mmtbx::ramachandran::init_module();
}
