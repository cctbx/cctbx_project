/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Created (rwgk)
 */

#include <cctbx/sgtbx/space_group.h>
#include <boost/python/class.hpp>
#include <boost/python/overloads.hpp>

namespace cctbx { namespace sgtbx { namespace boost_python {

namespace {

  struct phase_info_wrappers
  {
    typedef phase_info w_t;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      ht_angle_overloads, ht_angle, 0, 1)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      is_valid_phase_overloads, is_valid_phase, 1, 3)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      nearest_valid_phase_overloads, nearest_valid_phase, 1, 2)

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("phase_info", no_init)
        .def(init<space_group const&,
                  miller::index<> const&,
                  optional<bool> >())
        .def("sys_abs_was_tested", &w_t::sys_abs_was_tested)
        .def("is_sys_absent", &w_t::is_sys_absent)
        .def("is_centric", &w_t::is_centric)
        .def("ht", &w_t::ht)
        .def("t_den", &w_t::t_den)
        .def("ht_angle", &w_t::ht_angle, ht_angle_overloads())
        .def("is_valid_phase",
          &w_t::is_valid_phase,
          is_valid_phase_overloads())
        .def("nearest_valid_phase",
          &w_t::nearest_valid_phase,
          nearest_valid_phase_overloads())
      ;
    }
  };

} // namespace <anoymous>

  void wrap_phase_info()
  {
    phase_info_wrappers::wrap();
  }

}}} // namespace cctbx::sgtbx::boost_python
