#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/lbfgsb.h>
#include <boost/python/module.hpp>
#include <boost/python/class.hpp>

namespace scitbx { namespace lbfgsb { namespace {

  template <typename FloatType>
  struct minimizer_wrappers
  {
    typedef minimizer<FloatType> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("minimizer", no_init)
        .def(init<int const&,
                  int const&,
                  af::shared<FloatType>,
                  af::shared<FloatType>,
                  af::shared<int>,
                  FloatType const&,
                  FloatType const&,
                  int const&>())
        .def("task", &w_t::task)
        .def("request_stop", &w_t::request_stop)
        .def("request_stop_with_restore", &w_t::request_stop_with_restore)
        .def("process", &w_t::process)
      ;
    }
  };

  void init_module()
  {
    using namespace boost::python;

    minimizer_wrappers<double>::wrap();
  }

}}} // namespace scitbx::lbfgs::<anonymous>

BOOST_PYTHON_MODULE(scitbx_lbfgsb_ext)
{
  scitbx::lbfgsb::init_module();
}
