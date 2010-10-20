
#include <mmtbx/error.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/accessors/c_grid.h>

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
//#include <boost/python/return_value_policy.hpp>
//#include <boost/python/return_by_value.hpp>
#include <boost/optional.hpp>

#include <cmath>
#include <iostream>

namespace mmtbx { namespace ramachandran {
  namespace af = scitbx::af;

  class lookup_table
  {
    public:
      af::versa<double, af::flex_grid<> > plot;
      double max;

      lookup_table (
        af::const_ref< double > values,
        int n_angles)
      {
        MMTBX_ASSERT(values.size() == (n_angles * n_angles));
        af::flex_grid<>::index_type fg_origin;
        af::flex_grid<>::index_type fg_last;
        for (unsigned i = 0; i < 2; i++) {
          fg_origin.push_back(0.0);
          fg_last.push_back((double) n_angles);
        }
        plot = af::versa<double, af::flex_grid<> >(
          af::flex_grid<>(fg_origin, fg_last, true));
        max = 0.0;
        for (unsigned i = 0; i < values.size(); i++) {
          plot[i] = values[i];
          if (values[i] > max) {
            max = values[i];
          }
        }
      }

      double get_score (
        double phi,
        double psi)
      {
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
        double r11 = get_point(phi_1, psi_1);
        double r12 = get_point(phi_1, psi_2);
        double r21 = get_point(phi_2, psi_1);
        double r22 = get_point(phi_2, psi_2);
        double d_phi_d_psi = (double) (phi_2 - phi_1) * (psi_2 - psi_1);
        double r_phi_psi = ((r11/d_phi_d_psi) * (phi_2-phi) * (psi_2-psi)) +\
                           ((r21/d_phi_d_psi) * (phi-phi_1) * (psi_2-psi)) +\
                           ((r12/d_phi_d_psi) * (phi_2-phi) * (psi-psi_1)) +\
                           ((r22/d_phi_d_psi) * (phi-phi_1) * (psi-psi_1));
        return r_phi_psi;
      }

      double get_energy (
        double phi,
        double psi)
      {
        double score = get_score(phi, psi);
        return - (score - max);
      }

    private :
      double get_point (
        int phi,
        int psi)
      {
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
                int>((
        arg("values"),
        arg("n_angles"))))
      .def("get_score", &lookup_table::get_score, (
        arg("phi"),
        arg("psi")))
      .def("get_energy", &lookup_table::get_energy, (
        arg("phi"),
        arg("psi")));
  }
}} // namespace mmtbx::ramachandran

BOOST_PYTHON_MODULE(mmtbx_ramachandran_ext)
{
  mmtbx::ramachandran::init_module();
}
