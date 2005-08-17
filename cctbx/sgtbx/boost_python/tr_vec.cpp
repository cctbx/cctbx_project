#include <boost/python/class.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <cctbx/sgtbx/tr_vec.h>

namespace cctbx { namespace sgtbx { namespace boost_python {

namespace {

  struct tr_vec_wrappers
  {
    typedef tr_vec w_t;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      as_string_overloads, as_string, 0, 2)

    static std::string
    str(w_t const& o) { return o.as_string(); }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("tr_vec")
        .def(init<optional<int> >((arg_("tr_den")=sg_t_den)))
        .def(init<sg_vec3 const&, optional<int> >((
          arg_("v"),
          arg_("tr_den")=sg_t_den)))
        .def("num", (sg_vec3 const&(w_t::*)() const) &w_t::num, ccr())
        .def("den", (int const&(w_t::*)() const) &w_t::den, ccr())
        .def("__eq__", &w_t::operator==)
        .def("__ne__", &w_t::operator!=)
        .def("is_valid", &w_t::is_valid)
        .def("is_zero", &w_t::is_zero)
        .def("new_denominator", &w_t::new_denominator, (arg_("new_den")))
        .def("scale", &w_t::scale, (arg_("factor")))
        .def("mod_positive", &w_t::mod_positive)
        .def("mod_short", &w_t::mod_short)
        .def("cancel", &w_t::cancel)
        .def("as_double", &w_t::as_double)
        .def("plus", &w_t::plus, (arg_("rhs")))
        .def("minus", &w_t::minus, (arg_("rhs")))
        .def("as_string", &w_t::as_string, as_string_overloads((
          arg_("decimal")=false,
          arg_("separator")=",")))
        .def("__str__", str)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_tr_vec()
  {
    tr_vec_wrappers::wrap();
  }

}}} // namespace cctbx::sgtbx::boost_python
