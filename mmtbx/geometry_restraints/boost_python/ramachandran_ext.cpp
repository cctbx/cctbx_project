
#include <mmtbx/geometry_restraints/ramachandran.h>
#include <cctbx/geometry_restraints/proxy_select.h>
#include <scitbx/array_family/boost_python/shared_wrapper.h>

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/optional.hpp>

namespace mmtbx { namespace geometry_restraints {
namespace {
  using namespace boost::python;

  void wrap_ramachandran_proxies ()
  {
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

    // QUANTA-style harmonic restraints
    typedef rama_target_and_gradients w_t;
    class_<w_t>("rama_target_and_gradients", no_init)
      .def(init<af::ref<scitbx::vec3<double> > const&,
                double const&,
                double const&,
                double const&,
                af::const_ref<scitbx::vec3<double> > const&,
                af::const_ref<scitbx::vec3<double> > const&,
                phi_psi_proxy const&>((
        arg("gradient_array"),
        arg("phi_target"),
        arg("psi_target"),
        arg("weight"),
        arg("rama_table"),
        arg("sites_cart"),
        arg("proxy"))))
      .def("target", &w_t::target)
      .def("gradients", &w_t::gradients);
    def("target_phi_psi",
         (af::tiny<double, 3>(*)
           (af::const_ref<scitbx::vec3<double> > const&,
            af::const_ref<scitbx::vec3<double> > const&,
            phi_psi_proxy const&)) target_phi_psi,
              (arg("rama_table"),
               arg("sites_cart"),
               arg("proxy")));
  }

} // namespace anonymous

namespace boost_python {
  void init_module ()
  {
    wrap_ramachandran_proxies();
    wrap_ramachandran_targets();
  }

} // namespace boost_python
}} // namespace mmtbx::geometry_restraints

BOOST_PYTHON_MODULE(mmtbx_ramachandran_restraints_ext)
{
  mmtbx::geometry_restraints::boost_python::init_module();
}
