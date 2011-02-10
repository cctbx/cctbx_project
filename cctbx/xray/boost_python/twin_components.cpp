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
  struct twin_component_wrappers : boost::python::pickle_suite
  {
    typedef twin_component<FloatType> wt;

    static boost::python::tuple
    getinitargs(wt const &self) {
      return boost::python::make_tuple(self.twin_law, self.twin_fraction,
                                       self.grad_twin_fraction);
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      typedef default_call_policies dcp;
      class_<wt>("twin_component", no_init)
        .def(init<sgtbx::rot_mx const &,
                  FloatType const &,
                  optional<bool> >
                  ((arg("twin_law"),
                    arg("twin_fraction"),
                    arg("grad_twin_fraction"))))
        .def_pickle(twin_component_wrappers())
        .def("set_grad_twin_fraction", &wt::set_grad_twin_fraction)
        .def_readonly("twin_law", &wt::twin_law)
        .add_property("grad_twin_fraction",
                      make_getter(&wt::grad_twin_fraction, rbv()),
                      make_setter(&wt::grad_twin_fraction, dcp()))
        .add_property("twin_fraction", make_getter(&wt::twin_fraction, rbv()),
                                       make_setter(&wt::twin_fraction, dcp()))
      ;
    }
  };
} // namespace anonymous

  void wrap_twin_component()
  {
    using namespace boost::python;

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
      tuple_mapping_variable_capacity<af::shared<twin_component<double> *> >();
    }
  }

}}} // namespace cctbx::xray::boost_python
