/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Created (rwgk)
 */

#include <cctbx/sgtbx/space_group_type.h>
#include <boost/python/tuple.hpp>
#include <boost/python/class.hpp>
#include <boost/python/overloads.hpp>
#include <scitbx/boost_python/utils.h>

namespace cctbx { namespace sgtbx { namespace boost_python {

  // Split for Visual C++ 7
  void wrap_space_group_2(boost::python::class_<space_group>& cl);

namespace {

  struct space_group_wrappers : boost::python::pickle_suite
  {
    typedef space_group w_t;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      reset_overloads, reset, 0, 1)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      parse_hall_symbol_overloads, parse_hall_symbol, 1, 3)

    static rt_mx
    getitem(w_t const& o, std::size_t i_op)
    {
      if (i_op >= o.order_z()) scitbx::boost_python::raise_index_error();
      return o(i_op);
    }

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      z2p_op_overloads, z2p_op, 0, 2)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      construct_z2p_op_overloads, construct_z2p_op, 0, 2)

    static boost::python::tuple
    getinitargs(w_t const& o)
    {
      return boost::python::make_tuple(o.type().hall_symbol());
    }

    static void
    wrap()
    {
      using namespace boost::python;
      wrap_space_group_2(
      class_<w_t>("space_group")
        .def(init<parse_string&, optional<bool, bool, bool> >())
        .def(init<std::string const&, optional<bool, bool, bool> >())
        .def(init<space_group_symbols const&>())
        .def("reset", &w_t::reset, reset_overloads())
        .def("expand_ltr", &w_t::expand_ltr)
        .def("expand_inv", &w_t::expand_inv)
        .def("expand_smx", (void(w_t::*)(rt_mx const&)) &w_t::expand_smx)
        .def("expand_conventional_centring_type",
          &w_t::expand_conventional_centring_type)
        .def("parse_hall_symbol",
          &w_t::parse_hall_symbol,
          parse_hall_symbol_overloads())
        .def("change_basis", &w_t::change_basis)
        .def("r_den", &w_t::r_den)
        .def("t_den", &w_t::t_den)
        .def("order_p", &w_t::order_p)
        .def("order_z", &w_t::order_z)
        .def("__len__", &w_t::order_z)
        .def("n_ltr", &w_t::n_ltr)
        .def("is_centric", (bool(w_t::*)() const) &w_t::is_centric)
        .def("is_origin_centric", &w_t::is_origin_centric)
        .def("f_inv", &w_t::f_inv)
        .def("n_smx", &w_t::n_smx)
        .def("__call__",
          (rt_mx(w_t::*)(std::size_t, std::size_t, std::size_t) const)
          &w_t::operator())
        .def("__call__",
          (rt_mx(w_t::*)(std::size_t) const)
          &w_t::operator())
        .def("__getitem__", getitem)
        .def("make_tidy", &w_t::make_tidy)
        .def("__eq__", &w_t::operator==)
        .def("__ne__", &w_t::operator!=)
        .def("conventional_centring_type_symbol",
          &w_t::conventional_centring_type_symbol)
        .def("z2p_op", &w_t::z2p_op, z2p_op_overloads())
        .def("construct_z2p_op",
          &w_t::construct_z2p_op,
          construct_z2p_op_overloads())
        .def("is_chiral", &w_t::is_chiral)
        .def_pickle(space_group_wrappers())
      );
    }
  };

} // namespace <anoymous>

  void wrap_space_group()
  {
    space_group_wrappers::wrap();
  }

}}} // namespace cctbx::sgtbx::boost_python
