#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/module.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <cctbx/anharmonic.h>

namespace cctbx { namespace boost_python { namespace anharmonic {
  using namespace cctbx::adptbx::anharmonic;

  template <typename FloatType>
  struct anharmonic_adp_wrapper {
    typedef GramCharlier<FloatType> wt;

    static void wrap() {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;

      class_<wt>("gram_charlier", no_init)
        .def(init<af::shared<FloatType> const &>(
              (arg("Cijk"))))
        .def(init<af::shared<FloatType> const&,
          af::shared<FloatType> const&>(
            (arg("Cijk"),
              arg("Dijkl"))))
        .def("calculate", &wt::calculate,
          (arg("miller_index")))
        .def("get_order", &wt::get_order)
        .def("gradient_coefficients", &wt::gradient_coefficients,
          (arg("miller_index")))
        .def("gradient_coefficients_in_place", &wt::gradient_coefficients_in_place,
          (arg("miller_index"), arg("result_vector")))
        .def("data", &wt::data)
        ;
    }
  };

  void wrap_anharmonic_adp() {
    anharmonic_adp_wrapper<double>::wrap();
  }

  void init_module()
  {
    using namespace boost::python;

    wrap_anharmonic_adp();
  }
}}} // cctbx::boost_python::anharmonic

BOOST_PYTHON_MODULE(cctbx_anharmonic_ext)
{
  cctbx::boost_python::anharmonic::init_module();
}
