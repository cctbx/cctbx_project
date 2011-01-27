#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
//#include <boost/python/return_value_policy.hpp>
//#include <boost/python/return_by_value.hpp>
#include <boost/optional.hpp>

#include <mmtbx/geometry_restraints/ramachandran.h>

namespace mmtbx { namespace ramachandran {
  namespace boost_python {

  void init_module()
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

    // QUANTA-style harmonic restraints
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
      .def("gradients", &target_and_gradients::gradients);
    def("target_phi_psi",
         (af::tiny<double, 3>(*)
           (af::const_ref<scitbx::vec3<double> > const&,
            af::const_ref<scitbx::vec3<double> > const&,
            af::tiny<unsigned, 5> const&)) target_phi_psi,
              (arg("rama_table"),
               arg("sites_cart"),
               arg("i_seqs")));
  }

}}} // namespace mmtbx::ramachandran::boost_python

BOOST_PYTHON_MODULE(mmtbx_ramachandran_restraints_ext)
{
  mmtbx::ramachandran::boost_python::init_module();
}
