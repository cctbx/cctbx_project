/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Created (rwgk)
 */

#include <cctbx/sgtbx/tr_vec.h>
#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>

namespace cctbx { namespace sgtbx { namespace boost_python {

namespace {

  struct tr_vec_wrappers
  {
    typedef tr_vec w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("tr_vec")
        .def(init<int>())
        .def(init<sg_vec3 const&, optional<int> >())
        .def("num", (sg_vec3 const&(w_t::*)() const) &w_t::num, ccr())
        .def("den", (int const&(w_t::*)() const) &w_t::den, ccr())
        .def("__eq__", &w_t::operator==)
        .def("__ne__", &w_t::operator!=)
        .def("is_valid", &w_t::is_valid)
        .def("is_zero", &w_t::is_zero)
        .def("new_denominator", &w_t::new_denominator)
        .def("scale", &w_t::scale)
        .def("mod_positive", &w_t::mod_positive)
        .def("mod_short", &w_t::mod_short)
        .def("cancel", &w_t::cancel)
        .def("as_double", &w_t::as_double)
        .def("__float__", &w_t::as_double)
        .def("plus", &w_t::plus)
        .def("minus", &w_t::minus)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_tr_vec()
  {
    tr_vec_wrappers::wrap();
  }

}}} // namespace cctbx::sgtbx::boost_python
