#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <cctbx/miller/randomize_amplitude_and_phase.h>
#include <boost/python/def.hpp>

namespace cctbx { namespace miller { namespace boost_python {

  template <typename ComplexType, typename FloatType>
  struct randomize_amplitude_and_phase_wrapper
  {
    static void wrap() {
      using namespace boost::python;
      return_value_policy<return_by_value> rbv;

      def("randomize_amplitude_and_phase",
        (af::shared<ComplexType>(*)
          (af::const_ref<ComplexType> const&,
           af::const_ref<bool> const&,
           FloatType const&,
           FloatType const&,
           int))
             randomize_amplitude_and_phase, (
               arg("data"),
               arg("selection"),
               arg("amplitude_error"),
               arg("phase_error_deg"),
               arg("random_seed")))
      ;
    }
  };

  void wrap_randomize_amplitude_and_phase() {
    randomize_amplitude_and_phase_wrapper<std::complex<double>, double>::wrap();
  }
}}} // cctbx::miller::boostpython
