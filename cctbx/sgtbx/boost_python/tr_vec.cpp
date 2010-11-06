#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/operators.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <cctbx/sgtbx/tr_vec.h>
#include <cctbx/sgtbx/tr_vec_hash.h>
#include <boost_adaptbx/hash.h>

namespace cctbx { namespace sgtbx { namespace boost_python {

namespace {

  struct tr_vec_wrappers : boost_adaptbx::py_hashable<tr_vec>
  {
    typedef tr_vec w_t;

    static std::string
    str(w_t const& o) { return o.as_string(); }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("tr_vec")
        .def(init<optional<int> >((arg("tr_den")=sg_t_den)))
        .def(init<sg_vec3 const&, optional<int> >((
          arg("v"),
          arg("tr_den")=sg_t_den)))
        .def("num", (sg_vec3 const&(w_t::*)() const) &w_t::num, ccr())
        .def("den", (int const&(w_t::*)() const) &w_t::den, ccr())
        .def(self == self)
        .def(self != self)
        .def(-self)
        .def("is_valid", &w_t::is_valid)
        .def("is_zero", &w_t::is_zero)
        .def("new_denominator", &w_t::new_denominator, (arg("new_den")))
        .def("scale", &w_t::scale, (arg("factor")))
        .def("mod_positive", &w_t::mod_positive)
        .def("mod_short", &w_t::mod_short)
        .def("cancel", &w_t::cancel)
        .def("as_double", &w_t::as_double)
        .def("plus", &w_t::plus, (arg("rhs")))
        .def("minus", &w_t::minus, (arg("rhs")))
        .def("as_string", &w_t::as_string, (
          arg("decimal")=false,
          arg("separator")=","))
        .def("__str__", str)
        .def("__hash__", py_hash)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_tr_vec()
  {
    tr_vec_wrappers::wrap();
  }

}}} // namespace cctbx::sgtbx::boost_python
