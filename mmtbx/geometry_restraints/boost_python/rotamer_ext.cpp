
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/optional.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>

#include <mmtbx/geometry_restraints/ramachandran.h>
#include <mmtbx/geometry_restraints/rotamer.h>

namespace mmtbx { namespace geometry_restraints {
namespace {

  void wrap_rotamer_proxy ()
  {
    using namespace boost::python;
    typedef rotamer_proxy w_t;
    class_<w_t>("rotamer_proxy", no_init)
      .def(init<
        std::string const&,
        std::string const&,
        size_t,
        af::tiny<unsigned, 5> const&,
        af::tiny<unsigned, 4> const&,
        af::tiny<unsigned, 4> const&,
        af::tiny<unsigned, 4> const&,
        af::tiny<unsigned, 4> const&,
        optional< size_t const& > >((
          arg("residue_name"),
          arg("residue_type"),
          arg("n_angles"),
          arg("phi_psi"),
          arg("chi1"),
          arg("chi2"),
          arg("chi3"),
          arg("chi4"),
          arg("residue_index")=1)))
      .def(init<
        std::string const&,
        size_t,
        af::tiny<unsigned, 4> const&,
        af::tiny<unsigned, 4> const&,
        af::tiny<unsigned, 4> const&,
        af::tiny<unsigned, 4> const&,
        optional< size_t const& > >((
          arg("residue_name"),
          arg("n_angles"),
          arg("chi1"),
          arg("chi2"),
          arg("chi3"),
          arg("chi4"),
          arg("residue_index")=1)))
      .def(init<
        std::string const&,
        std::string const&,
        af::tiny<unsigned, 5> const&,
        optional<size_t> >((
          arg("residue_name"),
          arg("residue_type"),
          arg("phi_psi"),
          arg("residue_index")=1)))
      .def_readonly("have_phi_psi", &w_t::have_phi_psi)
      .def_readonly("n_angles", &w_t::n_angles)
      .def_readonly("residue_name", &w_t::residue_name)
      .def_readonly("residue_type", &w_t::residue_type)
      .def("get_rotamer_rmsd", &w_t::get_rotamer_rmsd, (
        arg("angles"),
        arg("sites_cart")))
      .def("find_dihedral_proxy", &w_t::find_dihedral_proxy, (
        arg("dihedral_proxy")));
    //  .def_readonly("phi_psi", &w_t::phi_psi_i_seqs)
    ;
    {
      typedef return_internal_reference<> rir;
      scitbx::af::boost_python::shared_wrapper<rotamer_proxy, rir>::wrap(
        "shared_rotamer_proxy");
    }
  }

  void wrap_ramachandran ()
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
        arg("psi")))
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
                rotamer_proxy const&>((
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
            rotamer_proxy const&)) target_phi_psi,
              (arg("rama_table"),
               arg("sites_cart"),
               arg("proxy")));
  }

} // namespace anonymous

namespace boost_python {
  void init_module ()
  {
    wrap_ramachandran();
    wrap_rotamer_proxy();
  }

} // namespace boost_python
}} // namespace mmtbx::geometry_restraints

BOOST_PYTHON_MODULE(mmtbx_rotamer_restraints_ext)
{
  mmtbx::geometry_restraints::boost_python::init_module();
}
