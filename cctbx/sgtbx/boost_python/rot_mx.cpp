/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Created (rwgk)
 */

#include <boost/python/class.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <cctbx/sgtbx/rot_mx_info.h>

namespace cctbx { namespace sgtbx { namespace boost_python {

namespace {

  struct rot_mx_wrappers
  {
    typedef rot_mx w_t;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      order_overloads, order, 0, 1)
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      accumulate_overloads, accumulate, 0, 1)

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("rot_mx", no_init)
        .def(init<optional<int, int> >())
        .def(init<sg_mat3 const&, optional<int> >())
        .def("num", (sg_mat3 const&(w_t::*)() const) &w_t::num, ccr())
        .def("den", (int const&(w_t::*)() const) &w_t::den, ccr())
        .def("__eq__", &w_t::operator==)
        .def("__ne__", &w_t::operator!=)
        .def("is_valid", &w_t::is_valid)
        .def("is_unit_mx", &w_t::is_unit_mx)
        .def("minus_unit_mx", &w_t::minus_unit_mx)
        .def("new_denominator", &w_t::new_denominator)
        .def("scale", &w_t::scale)
        .def("determinant", &w_t::determinant)
        .def("inverse", &w_t::inverse)
        .def("cancel", &w_t::cancel)
        .def("inverse_cancel", &w_t::inverse_cancel)
        .def("multiply", (rot_mx (w_t::*)(rot_mx const&) const) &w_t::multiply)
        .def("multiply", (tr_vec (w_t::*)(tr_vec const&) const) &w_t::multiply)
        .def("divide", &w_t::divide)
        .def("type", &w_t::type)
        .def("order", &w_t::order, order_overloads())
        .def("accumulate", &w_t::accumulate, accumulate_overloads())
        .def("info", &w_t::info)
        .def("as_double", &w_t::as_double)
        .def("__float__", &w_t::as_double)
      ;
    }
  };

  struct rot_mx_info_wrappers
  {
    typedef rot_mx_info w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("rot_mx_info", no_init)
        .def(init<rot_mx const&>())
        .def("type", &w_t::type)
        .def("ev", &w_t::ev, ccr())
        .def("sense", &w_t::sense)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_rot_mx()
  {
    rot_mx_wrappers::wrap();
    rot_mx_info_wrappers::wrap();
  }

}}} // namespace cctbx::sgtbx::boost_python
