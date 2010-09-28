#include <smtbx/refinement/constraints/scatterer_parameters.h>

#include <boost/python/class.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <scitbx/array_family/boost_python/shared_wrapper.h>

namespace smtbx { namespace refinement { namespace constraints {
namespace boost_python {

struct scatterer_parameters_wrapper
{
  typedef scatterer_parameters wt;

  static
  af::shared<wt>
  *init_shared_scatterer_parameters(af::const_ref<wt::scatterer_type> const
                                   &scatterers)
  {
    af::shared<wt>
    *result = new af::shared<wt>((af::reserve(scatterers.size())));
    for (std::size_t i_sc=0; i_sc<scatterers.size(); ++i_sc) {
      result->push_back(wt(&scatterers[i_sc]));
    }
    return result;
  }

  static void wrap() {
    using namespace boost::python;
    typedef return_internal_reference<> rir_t;
    rir_t rir;
    class_<wt>("scatterer_parameters", no_init)
      .def(init<wt::scatterer_type *>(arg("scatterer")))
      .def(init<wt::scatterer_type *,
                crystallographic_parameter *,
                crystallographic_parameter *,
                crystallographic_parameter *>(
           (arg("scatterer"), arg("site"), arg("occupancy"), arg("u"))))
      .add_property("scatterer", make_getter(&wt::scatterer, rir))
      .add_property("site"     , make_getter(&wt::site, rir)
                               , make_setter(&wt::site))
      .add_property("occupancy", make_getter(&wt::occupancy, rir)
                               , make_setter(&wt::occupancy))
      .add_property("u"        , make_getter(&wt::u, rir)
                               , make_setter(&wt::u))
      ;

    scitbx::af::boost_python::shared_wrapper<wt, rir_t>
    ::wrap("shared_scatterer_parameters")
      .def("__init__", make_constructor(init_shared_scatterer_parameters))
      .def("mapping_to_grad_fc", mapping_to_grad_fc)
      ;
  }
};

void wrap_scatterer_parameters() {
  scatterer_parameters_wrapper::wrap();
}


}}}}
