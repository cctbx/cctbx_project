#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/iterator.hpp>

#include <cctbx/xray/parameter_map.h>
#include <cctbx/xray/scatterer.h>

namespace cctbx { namespace xray { namespace boost_python {

struct parameter_indices_wrapper
{
  typedef parameter_indices wt;

  static int invariable() {
    return wt::invariable;
  }

  static void wrap() {
    using namespace boost::python;
    class_<wt>("parameter_indices", no_init)
      .add_static_property("invariable", invariable)
      .def_readonly("site", &wt::site)
      .def_readonly("u_iso", &wt::u_iso)
      .def_readonly("u_aniso", &wt::u_aniso)
      .def_readonly("occupancy", &wt::occupancy)
      .def_readonly("fp", &wt::fp)
      .def_readonly("fdp", &wt::fdp)
      ;
  }
};

template <class XRayScattererType>
struct parameter_map_wrapper
{
  typedef parameter_map<XRayScattererType> wt;

  static void wrap(char const *name) {
    using namespace boost::python;
    typedef return_value_policy<copy_const_reference> ccr;
    typedef return_value_policy<return_by_value> rbv;
    class_<wt>("parameter_map", no_init)
      .def(init<af::const_ref<typename wt::xray_scatterer_type> const &>((
            arg("scatterers"))))
      .def(init<af::const_ref<typename wt::xray_scatterer_type> const &,
                af::shared<twin_component<double> *> const &>((
            arg("scatterers"), arg("twin_components"))))
      .def("__len__", &wt::size)
      .def("__getitem__", &wt::operator[], ccr())
      .def("__iter__", iterator<wt, ccr>())
      .add_property("n_parameters", &wt::n_parameters)
      .add_property("n_scatterers", &wt::n_scatterers)
      .add_property("twin_fractions", make_getter(&wt::twin_fractions, rbv()))
      ;
  }
};

void wrap_parameter_map() {
  parameter_indices_wrapper::wrap();
  parameter_map_wrapper<xray::scatterer<> >::wrap("parameter_map");
}

}}}
