#include <smtbx/refinement/constraints/reparametrisation.h>

#include <boost/python/class.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/str.hpp>

#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <sstream>
#include <iostream>

namespace smtbx { namespace refinement { namespace constraints {
namespace boost_python {

af::shared<std::size_t>
mapping_to_grad_fc(af::const_ref<independent_scalar_parameter *> const &params) {
  af::shared<std::size_t> result((af::reserve(params.size())));
  for (std::size_t i=0; i<params.size(); ++i) {
    if (params[i]->is_variable()) { result.push_back(params[i]->index()); }
  }
  return result;
}

struct independent_scalar_parameters_wrapper
{
  typedef independent_scalar_parameter wt;

  static void wrap() {
    using namespace boost::python;
    typedef return_internal_reference<> rir_t;
    rir_t rir;

    scitbx::af::boost_python::shared_wrapper<wt *, rir_t>
    ::wrap("shared_independent_shared_parameters")
      .def("mapping_to_grad_fc", mapping_to_grad_fc)
      ;
  }
};

void wrap_independent_scalar_parameters() {
  independent_scalar_parameters_wrapper::wrap();
}

}}}}
