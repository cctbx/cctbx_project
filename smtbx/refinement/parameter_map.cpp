#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>

#include <cctbx/xray/scatterer.h>
#include <smtbx/import_cctbx.h>

#include <smtbx/refinement/parameter_map.h>

namespace smtbx { namespace refinement { namespace boost_python {

struct parameter_indices_wrapper
{
  typedef parameter_indices wt;

  static void wrap() {
    using namespace boost::python;
    class_<wt>("parameter_indices", no_init)
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
    return_value_policy<copy_const_reference> ccr;
    class_<wt>("parameter_map", no_init)
      .def(init<af::const_ref<typename wt::xray_scatterer_type> const &>(
            arg_("scatterers")))
      .def("__len__", &wt::size)
      .def("__getitem__", &wt::operator[], ccr)
      .def("n_parameters", &wt::n_parameters)
      ;
  }
};

void wrap_parameter_map() {
  parameter_indices_wrapper::wrap();
  parameter_map_wrapper<xray::scatterer<> >::wrap("parameter_map");
}

}}}
