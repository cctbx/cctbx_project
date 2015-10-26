#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/optional.hpp>

#include <mmtbx/geometry_restraints/ramachandran.h>
#include <cctbx/geometry_restraints/proxy_select.h>
#include <scitbx/array_family/boost_python/shared_wrapper.h>

namespace mmtbx { namespace geometry_restraints {
namespace boost_python {

  void wrap_ramachandran_proxies ()
  {
    using namespace boost::python;
    typedef phi_psi_proxy w_t;
    class_<w_t>("phi_psi_proxy", no_init)
      .def(init<
        w_t::i_seqs_type const&,
        std::string const&,
        std::string const&,
        optional<size_t> >((
          arg("i_seqs"),
          arg("residue_name"),
          arg("residue_type"),
          arg("residue_index")=1)))
      .def_readonly("residue_name", &w_t::residue_name)
      .def_readonly("residue_type", &w_t::residue_type)
      .def("get_i_seqs", &w_t::get_i_seqs)
    //  .def_readonly("phi_psi", &w_t::phi_psi_i_seqs)
    ;
    {
      typedef return_internal_reference<> rir;
      scitbx::af::boost_python::shared_wrapper<phi_psi_proxy, rir>::wrap(
        "shared_phi_psi_proxy")
        .def("proxy_select",
          (af::shared<w_t>(*)(
           af::const_ref<w_t> const&,
           std::size_t,
           af::const_ref<std::size_t> const&))
           cctbx::geometry_restraints::shared_proxy_select, (
         arg("n_seq"), arg("iselection")));
    }
  }

  void wrap_ramachandran_targets ()
  {
    using namespace boost::python;
    // COOT-like restraints
    class_<lookup_table>("lookup_table", no_init)
      .def(init<af::const_ref< double >,
                int,
                double>((
        arg("values"),
        arg("n_angles"),
        arg("scale_allowed")=1.0)))
      .def("get_score", &lookup_table::get_score, (
        arg("phi"),
        arg("psi"),
        arg("use_splines")=false))
      .def("get_energy", &lookup_table::get_energy, (
        arg("phi"),
        arg("psi")))
      .def("compute_gradients", &lookup_table::compute_gradients, (
        arg("gradient_array"),
        arg("sites_cart"),
        arg("proxy"),
        arg("weight")=1.0,
        arg("epsilon")=0.1));

    def("target_phi_psi",
         (af::tiny<double, 3>(*)
           (af::const_ref<scitbx::vec3<double> > const&,
            af::const_ref<scitbx::vec3<double> > const&,
            phi_psi_proxy const&)) target_phi_psi,
              (arg("rama_table"),
               arg("sites_cart"),
               arg("proxy")));
    // QUANTA-style harmonic restraints
    def("ramachandran_residual_sum",
        (double(*) (
          af::const_ref<scitbx::vec3<double> > const&,
          af::const_ref<phi_psi_proxy> const&,
          af::ref<scitbx::vec3<double> > const&,
          af::const_ref<scitbx::vec3<double> > const&,
          af::const_ref<scitbx::vec3<double> > const&,
          af::const_ref<scitbx::vec3<double> > const&,
          af::const_ref<scitbx::vec3<double> > const&,
          af::const_ref<scitbx::vec3<double> > const&,
          af::const_ref<scitbx::vec3<double> > const&,
          af::tiny<double, 4>,
          af::ref<double> const&))
        ramachandran_residual_sum,
        (arg("sites_cart"), arg("proxies"), arg("gradient_array"),
        arg("general_table"), arg("gly_table"), arg("cispro_table"),
        arg("transpro_table"), arg("prepro_table"), arg("ileval_table"),
        arg("weights"), arg("residuals_array")));
        //sites_cart, proxies,gr_arr, 6 tables
  }

  void init_module ()
  {
    wrap_ramachandran_proxies();
    wrap_ramachandran_targets();
  }

}}} // namespace mmtbx::geometry_restraints::boost_python

BOOST_PYTHON_MODULE(mmtbx_ramachandran_restraints_ext)
{
  mmtbx::geometry_restraints::boost_python::init_module();
}
