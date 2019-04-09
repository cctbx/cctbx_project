#include <cctbx/boost_python/flex_fwd.h>

#include <scitbx/boost_python/is_polymorphic_workaround.h>

#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/shared_wrapper.h>

#include <scitbx/array_family/shared.h>

#include <mmtbx/tls/optimise_amplitudes.h>

SCITBX_BOOST_IS_POLYMORPHIC_WORKAROUND(mmtbx::tls::common)

namespace mmtbx { namespace tls { namespace optimise {
  namespace bp = boost::python;
  namespace af = scitbx::af;

namespace {
  void init_module()
  {
    using namespace boost::python;
    using boost::python::arg;

    class_<MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator>(
        "MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator",
        init< symArrNd const&,
              dblArrNd const&,
              dblArrNd const&,
              bp::list const&,
              bp::list const&,
              selArr1d const&,
              symArr1d const& >(
                (arg("target_uijs"),
                 arg("target_weights"),
                 arg("base_amplitudes"),
                 arg("base_uijs"),
                 arg("base_atom_indices"),
                 arg("dataset_hash"),
                 arg("residual_uijs"))
                )
              )
      .def("set_current_amplitudes", &MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator::setCurrentAmplitudes)
      .def("get_current_amplitudes", &MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator::getCurrentAmplitudes)
      .def("print_current_amplitudes", &MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator::printCurrentAmplitudes)
      .add_property("x",
          &MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator::getCurrentAmplitudes,
          &MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator::setCurrentAmplitudes
          )
      .def("set_residual_mask",
          &MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator::setResidualMask,
          ( args("mask") ),
          "Select which datasets are used to optimise the residual levels (if any)"
          )
      .def("compute_functional_and_gradients",
          &MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator::computeFunctionalAndGradients
          //return_value_policy<manage_new_object>()
          )
    ;
  }

} // Close unnamed

}}} // close mmtbx::tls::optimise

BOOST_PYTHON_MODULE(mmtbx_tls_optimise_amplitudes_ext)
{
  mmtbx::tls::optimise::init_module();
}
