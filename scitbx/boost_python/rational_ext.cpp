/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2003 Jan: Created (R.W. Grosse-Kunstleve)
 */

#include <boost/rational.hpp>
#include <scitbx/boost_python/utils.h>
#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/operators.hpp>

namespace scitbx { namespace boost_python { namespace {

  struct rational_int_wrappers
  {
    typedef boost::rational<int> w_t;

    static boost::python::tuple
    as_tuple(w_t const& o)
    {
      return boost::python::make_tuple(o.numerator(), o.denominator());
    }

    static bool eq_rr(w_t const& lhs, w_t const& rhs) { return lhs == rhs; }
    static bool ne_rr(w_t const& lhs, w_t const& rhs) { return lhs != rhs; }
    static bool lt_rr(w_t const& lhs, w_t const& rhs) { return lhs < rhs; }
    static bool gt_rr(w_t const& lhs, w_t const& rhs) { return lhs > rhs; }
    static bool le_rr(w_t const& lhs, w_t const& rhs) { return lhs <= rhs; }
    static bool ge_rr(w_t const& lhs, w_t const& rhs) { return lhs >= rhs; }

    static bool eq_ri(w_t const& lhs, int rhs) { return lhs == rhs; }
    static bool ne_ri(w_t const& lhs, int rhs) { return lhs != rhs; }
    static bool lt_ri(w_t const& lhs, int rhs) { return lhs < rhs; }
    static bool gt_ri(w_t const& lhs, int rhs) { return lhs > rhs; }
    static bool le_ri(w_t const& lhs, int rhs) { return lhs <= rhs; }
    static bool ge_ri(w_t const& lhs, int rhs) { return lhs >= rhs; }

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("int")
        .def(init<int, optional<int> >())
        .def("numerator", &w_t::numerator)
        .def("denominator", &w_t::denominator)
        .def("as_tuple", as_tuple)
        .def(self + self)
        .def(self - self)
        .def(self * self)
        .def(self / self)
        .def(self + int())
        .def(self - int())
        .def(self * int())
        .def(self / int())
        .def(int() + self)
        .def(int() - self)
        .def(int() * self)
        .def(int() / self)
        .def("__eq__", eq_rr)
        .def("__ne__", ne_rr)
        .def("__lt__", lt_rr)
        .def("__gt__", gt_rr)
        .def("__le__", le_rr)
        .def("__ge__", ge_rr)
        .def("__eq__", eq_ri)
        .def("__ne__", ne_ri)
        .def("__lt__", lt_ri)
        .def("__gt__", gt_ri)
        .def("__le__", le_ri)
        .def("__ge__", ge_ri)
      ;
    }
  };

  void init_module()
  {
    using namespace boost::python;

    scope().attr("__version__") = scitbx::boost_python::cvs_revision(
      "$Revision$");

    rational_int_wrappers::wrap();
    def("gcd", (int(*)(int,int))boost::gcd);
    def("lcm", (int(*)(int,int))boost::lcm);
  }

}}} // namespace scitbx::boost_python::<anonymous>

BOOST_PYTHON_MODULE(rational_ext)
{
  scitbx::boost_python::init_module();
}
