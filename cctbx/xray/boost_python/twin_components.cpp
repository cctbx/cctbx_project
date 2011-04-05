#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/twin_component.h>

#include <scitbx/array_family/shared.h>
#include <scitbx/stl/vector_wrapper.h>

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>


namespace cctbx { namespace xray { namespace boost_python {

namespace {

  template <typename FloatType>
  struct twin_fraction_wrapper {
    static void wrap() {
      typedef twin_fraction<FloatType> wt;
      using namespace boost::python;
      class_<wt>("twin_fraction", no_init)
        .def(init<FloatType,
                  optional<bool> >
             ((arg("value"),
               arg("grad"))))
        .def(init<twin_fraction<FloatType> const&>
             ((arg("source"))))
        .def_readwrite("grad_index", &wt::grad_index)
        .def_readwrite("grad", &wt::grad)
        .def_readwrite("value", &wt::value)
        .def("deep_copy", &wt::deep_copy)
        ;
    }
  };

  template <typename FloatType>
  struct twin_component_wrappers : boost::python::pickle_suite
  {
    typedef twin_component<FloatType> wt;

    static boost::python::tuple
    getinitargs(wt const &self) {
      return boost::python::make_tuple(self.twin_law,
        self.value, self.grad);
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef default_call_policies dcp;
      class_<wt, bases<twin_fraction<FloatType> > >("twin_component", no_init)
        .def(init<sgtbx::rot_mx const &,
                  FloatType,
                  bool>
                  ((arg("twin_law"),
                    arg("value"),
                    arg("grad"))))
        .def("deep_copy", &wt::deep_copy)
        .def_pickle(twin_component_wrappers())
        .def_readonly("twin_law", &wt::twin_law)
      ;
    }
  };
} // namespace anonymous

  void wrap_twin_component()
  {
    using namespace boost::python;

    twin_fraction_wrapper<double>::wrap();
    twin_component_wrappers<double>::wrap();

    def("set_grad_twin_fraction",
      (void(*)(
        af::shared<twin_component<double> *>, bool))
        set_grad_twin_fraction, (
          arg("twin_components"),
          arg("grad_twin_fraction")=true));

    def("sum_twin_fractions",
      (double(*)(
        af::shared<twin_component<double> *>))
        sum_twin_fractions, (arg("twin_components")));

    {
      using namespace scitbx::boost_python::container_conversions;
      tuple_mapping_variable_capacity<af::shared<twin_fraction<double> *> >();
      tuple_mapping_variable_capacity<af::shared<twin_component<double> *> >();
    }
  }

}}} // namespace cctbx::xray::boost_python
