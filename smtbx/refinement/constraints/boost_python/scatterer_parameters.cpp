#include <smtbx/refinement/constraints/scatterer_parameters.h>

#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <scitbx/array_family/boost_python/shared_wrapper.h>

namespace smtbx { namespace refinement { namespace constraints {
namespace boost_python {

struct scatterer_parameters_wrapper
{
  typedef scatterer_parameters wt;

  static void wrap() {
    using namespace boost::python;
    typedef return_internal_reference<> rir_t;
    rir_t rir;
    class_<wt>("scatterer_parameters", no_init)
      .def(init<>())
      .def(init<parameter *, parameter *, parameter *>(
           (arg("site"), arg("occupancy"), arg("u"))))
      .add_property("site"     , make_getter(&wt::site, rir)
                               , make_setter(&wt::site))
      .add_property("occupancy", make_getter(&wt::occupancy, rir)
                               , make_setter(&wt::occupancy))
      .add_property("u"        , make_getter(&wt::u, rir)
                               , make_setter(&wt::u))
      ;
    scitbx::af::boost_python::shared_wrapper<scatterer_parameters, rir_t>
    ::wrap("shared_scatterer_parameters");
  }
};

void wrap_scatterer_parameters() {
  scatterer_parameters_wrapper::wrap();
}


}}}}
