/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Created (rwgk)
 */

#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/sgtbx/change_of_basis_op.h>
#include <boost/python/tuple.hpp>
#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>

namespace cctbx { namespace sgtbx { namespace boost_python {

namespace {

  struct change_of_basis_op_wrappers : boost::python::pickle_suite
  {
    typedef change_of_basis_op w_t;

    static w_t
    new_denominators_int_int(w_t const& o, int r_den, int t_den)
    {
      return o.new_denominators(r_den, t_den);
    }

    static w_t
    new_denominators_w_t(w_t const& o, w_t const& other)
    {
      return o.new_denominators(other);
    }

    static void
    update_w_t(w_t& o, w_t const& other)
    {
      o.update(other);
    }

    static boost::python::tuple
    getinitargs(w_t const& o)
    {
      return boost::python::make_tuple(o.c().as_xyz());
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("change_of_basis_op", no_init)
        .def(init<rt_mx const&, rt_mx const&>())
        .def(init<rt_mx const&>())
        .def(init<std::string const&, optional<const char*, int, int> >())
        .def(init<optional<int, int> >())
        .def("is_valid", &w_t::is_valid)
        .def("identity_op", &w_t::identity_op)
        .def("is_identity_op", &w_t::is_identity_op)
        .def("new_denominators", new_denominators_int_int)
        .def("new_denominators", new_denominators_w_t)
        .def("c", &w_t::c, ccr())
        .def("c_inv", &w_t::c_inv, ccr())
        .def("select", &w_t::select, ccr())
        .def("inverse", &w_t::inverse)
        .def("mod_positive_in_place", &w_t::mod_positive_in_place)
        .def("mod_short_in_place", &w_t::mod_short_in_place)
        .def("apply",
          (rt_mx(w_t::*)(rt_mx const&) const)
          &w_t::apply)
        .def("apply",
          (uctbx::unit_cell(w_t::*)(uctbx::unit_cell const&) const)
          &w_t::apply)
        .def("apply",
          (miller::index<>(w_t::*)(miller::index<> const&) const)
          &w_t::apply)
        .def("apply",
          (af::shared<miller::index<> >(w_t::*)
            (af::const_ref<miller::index<> > const&) const)
              &w_t::apply)
        .def("__call__",
          (fractional<>(w_t::*)(fractional<> const&) const)
          &w_t::operator())
        .def("update", update_w_t)
        .def("__mul__", &w_t::operator*)
        .def_pickle(change_of_basis_op_wrappers())
      ;
    }
  };

} // namespace <anoymous>

  void wrap_change_of_basis_op()
  {
    change_of_basis_op_wrappers::wrap();
  }

}}} // namespace cctbx::sgtbx::boost_python
