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
#include <boost/python/return_internal_reference.hpp>
#include <cctbx/sgtbx/rt_mx.h>

namespace cctbx { namespace sgtbx { namespace boost_python {

namespace {

  struct rt_mx_wrappers
  {
    typedef rt_mx w_t;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      as_xyz_overloads, as_xyz, 0, 4)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      new_denominators_overloads, new_denominators, 1, 2)

    static std::string
    str(w_t const& o) { return o.as_xyz(); }

    static scitbx::vec3<double>
    mul(w_t const& o, scitbx::vec3<double> const& rhs) { return o * rhs; }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_internal_reference<> rir;
      class_<w_t>("rt_mx", no_init)
        .def(init<optional<int, int> >())
        .def(init<rot_mx const&, tr_vec const&>())
        .def(init<rot_mx const&, optional<int> >())
        .def(init<tr_vec const&, optional<int> >())
        .def(init<parse_string&, optional<const char*, int, int> >())
        .def(init<std::string const&, optional<const char*, int, int> >())
        .def("r", (rot_mx const&(w_t::*)() const) &w_t::r, rir())
        .def("t", (tr_vec const&(w_t::*)() const) &w_t::t, rir())
        .def("__eq__", &w_t::operator==)
        .def("__ne__", &w_t::operator!=)
        .def("is_valid", &w_t::is_valid)
        .def("unit_mx", &w_t::unit_mx)
        .def("is_unit_mx", &w_t::is_unit_mx)
        .def("as_xyz", &w_t::as_xyz, as_xyz_overloads())
        .def("__str__", str)
        .def("as_int_array", &w_t::as_int_array)
        .def("as_double_array", &w_t::as_double_array)
        .def("new_denominators",
          (rt_mx(w_t::*)(int, int) const) 0,
            new_denominators_overloads())
        .def("new_denominators",
          (rt_mx(w_t::*)(rt_mx const&) const)
            &w_t::new_denominators)
        .def("mod_positive", &w_t::mod_positive)
        .def("mod_short", &w_t::mod_short)
        .def("inverse", &w_t::inverse)
        .def("refine_gridding", (sg_vec3(w_t::*)(sg_vec3 const&) const)
          &w_t::refine_gridding)
        .def("cancel", &w_t::cancel)
        .def("inverse_cancel", &w_t::inverse_cancel)
        .def("multiply", &w_t::multiply)
        .def("__mul__", mul)
      ;
    }
  };

  struct translation_part_info_wrappers
  {
    typedef translation_part_info w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("translation_part_info", no_init)
        .def(init<rt_mx const&>())
        .def("intrinsic_part", &w_t::intrinsic_part, ccr())
        .def("location_part", &w_t::location_part, ccr())
        .def("origin_shift", &w_t::origin_shift, ccr())
      ;
    }
  };

} // namespace <anoymous>

  void wrap_rt_mx()
  {
    rt_mx_wrappers::wrap();
    translation_part_info_wrappers::wrap();
  }

}}} // namespace cctbx::sgtbx::boost_python
