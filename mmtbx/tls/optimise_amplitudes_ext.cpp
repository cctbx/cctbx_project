#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python.hpp>
#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/shared_wrapper.h>

#include <mmtbx/tls/optimise_amplitudes.h>

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
              symArr1d const&,
              double,
              double,
              double >(
                (arg("target_uijs"),
                 arg("target_weights"),
                 arg("base_amplitudes"),
                 arg("base_uijs"),
                 arg("base_atom_indices"),
                 arg("base_dataset_hash"),
                 arg("atomic_uijs"),
                 arg("weight_sum_of_amplitudes"),
                 arg("weight_sum_of_squared_amplitudes"),
                 arg("weight_sum_of_amplitudes_squared"))
                )
              )
      .def("set_current_amplitudes", &MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator::setCurrentAmplitudes)
      .def("get_current_amplitudes", &MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator::getCurrentAmplitudes)
      .def("print_current_amplitudes", &MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator::printCurrentAmplitudes)
      .add_property("x",
          &MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator::getCurrentAmplitudes,
          &MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator::setCurrentAmplitudes
          )
      .def("set_atomic_optimisation_mask",
          &MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator::setAtomicOptimisationMask,
          ( args("mask") ),
          "Select which datasets are used to optimise the atomic level"
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
