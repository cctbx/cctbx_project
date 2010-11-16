#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/optional.hpp>

#include <mmtbx/error.h>
#include <cctbx/geometry_restraints/dihedral.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>

#include <cmath>
#include <iostream>

namespace mmtbx { namespace ramachandran { namespace restraints {
  namespace af = scitbx::af;
  using cctbx::geometry_restraints::dihedral;
  namespace gr = cctbx::geometry_restraints;

  class target_and_gradients
  {
    public:
      target_and_gradients (
        af::ref<scitbx::vec3<double> > const& gradient_array,
        double const& phi_target,
        double const& psi_target,
        double const& weight,
        af::const_ref<scitbx::vec3<double> > const& rama_table,
        af::const_ref<scitbx::vec3<double> > const& sites_cart,
        af::tiny<unsigned, 5> const& i_seqs)
      {
        MMTBX_ASSERT(gradient_array.size() == sites_cart.size());
        gradients_.resize(sites_cart.size(), scitbx::vec3<double>(0,0,0));
        using cctbx::geometry_restraints::dihedral;
        af::tiny<scitbx::vec3<double>, 4> phi_sites;
        af::tiny<scitbx::vec3<double>, 4> psi_sites;
        for (unsigned i = 0; i < 4; i++) {
          phi_sites[i] = sites_cart[i_seqs[i]];
          psi_sites[i] = sites_cart[i_seqs[i+1]];
        }
        double esd = 10.0;
        double w = 1./(esd * esd) * weight;
        dihedral phi(phi_sites, phi_target, w);
        dihedral psi(psi_sites, psi_target, w);
        target_ = phi.residual()+psi.residual();
        af::tiny<scitbx::vec3<double>, 4> d_phi_d_xyz = phi.gradients();
        af::tiny<scitbx::vec3<double>, 4> d_psi_d_xyz = psi.gradients();
        for (unsigned k = 0; k < 5; k++) {
          std::size_t i_seq = i_seqs[k];
          if(k < 4) gradient_array[i_seq] += d_phi_d_xyz[k] ;
          if(k > 0) gradient_array[i_seq] += d_psi_d_xyz[k-1] ;
        }
      }

      double target() { return target_; }

      af::shared<scitbx::vec3<double> > gradients() { return gradients_; }

    private :
      double target_;
      af::shared<scitbx::vec3<double> > gradients_;
  };

template <typename FloatType>
af::tiny<FloatType, 3>
  target_phi_psi(af::const_ref<scitbx::vec3<double> > const& rama_table,
                 af::const_ref<scitbx::vec3<double> > const& sites_cart,
                 af::tiny<unsigned, 5> const& i_seqs)
{
  af::tiny<scitbx::vec3<double>, 4> phi_sites;
  af::tiny<scitbx::vec3<double>, 4> psi_sites;
  for (unsigned i = 0; i < 4; i++) {
    phi_sites[i] = sites_cart[i_seqs[i]];
    psi_sites[i] = sites_cart[i_seqs[i+1]];
  }
  dihedral phi1(phi_sites, 0, 1.0);
  dihedral psi1(psi_sites, 0, 1.0);
  double phi_deg = phi1.angle_model;
  double psi_deg = psi1.angle_model;
  double phi_t=0, psi_t=0;
  double dist_to_current = 1.e+9;
  double score_current = 0;
  for(int i=0; i<rama_table.size(); i++) {
    scitbx::vec3<double> point = rama_table[i];
    double d1 = gr::angle_delta_deg(point[0],phi_deg);
    double d2 = gr::angle_delta_deg(point[1],psi_deg);
    double d = std::sqrt(d1 * d1 + d2 * d2);
    if(d<dist_to_current) {
      dist_to_current = d;
      score_current = point[2];
    }
  }
  double dist_to_allowed = 1.e+9;
  for(int i=0; i<rama_table.size(); i++) {
    scitbx::vec3<double> point = rama_table[i];
    double d1 = gr::angle_delta_deg(point[0],phi_deg);
    double d2 = gr::angle_delta_deg(point[1],psi_deg);
    double d = std::sqrt(d1 * d1 + d2 * d2);
    if(point[2] >= score_current && d<dist_to_allowed) {
      dist_to_allowed = d;
      phi_t = point[0];
      psi_t = point[1];
    }
  }
  return af::tiny<double, 3> (phi_t,psi_t,dist_to_allowed) ;
};


  void init_module()
  {
    using namespace boost::python;
    class_<target_and_gradients>("target_and_gradients", no_init)
      .def(init<af::ref<scitbx::vec3<double> > const&,
                double const&,
                double const&,
                double const&,
                af::const_ref<scitbx::vec3<double> > const&,
                af::const_ref<scitbx::vec3<double> > const&,
                af::tiny<unsigned, 5> const&>((
        arg("gradient_array"),
        arg("phi_target"),
        arg("psi_target"),
        arg("weight"),
        arg("rama_table"),
        arg("sites_cart"),
        arg("i_seqs"))))
      .def("target", &target_and_gradients::target)
      .def("gradients", &target_and_gradients::gradients)
    ;
    def("target_phi_psi",
         (af::tiny<double, 3>(*)
           (af::const_ref<scitbx::vec3<double> > const&,
            af::const_ref<scitbx::vec3<double> > const&,
            af::tiny<unsigned, 5> const&)) target_phi_psi,
              (arg("rama_table"),
               arg("sites_cart"),
               arg("i_seqs")))
   ;
  }

  }}} // namespace mmtbx::ramachandran::restraints

BOOST_PYTHON_MODULE(mmtbx_ramachandran_restraints_ext)
{
  mmtbx::ramachandran::restraints::init_module();
}
