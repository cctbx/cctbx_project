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
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>

namespace cctbx { namespace sgtbx { namespace boost_python {

namespace {

  struct space_group_type_wrappers : boost::python::pickle_suite
  {
    typedef space_group_type w_t;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      hall_symbol_overloads, hall_symbol, 0, 1)

    static boost::python::tuple
    getinitargs(w_t const& o)
    {
      return boost::python::make_tuple("Hall: " + o.hall_symbol());
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_internal_reference<> rir;
      class_<w_t>("space_group_type")
        .def(init<std::string const&, optional<std::string const&> >())
        .def(init<space_group const&, optional<bool, int, int> >())
        .def("group", &w_t::group, rir())
        .def("number", &w_t::number)
        .def("cb_op", &w_t::cb_op, ccr())
        .def("addl_generators_of_euclidean_normalizer",
          &w_t::addl_generators_of_euclidean_normalizer)
        .def("expand_addl_generators_of_euclidean_normalizer",
          &w_t::expand_addl_generators_of_euclidean_normalizer)
        .def("is_enantiomorphic", &w_t::is_enantiomorphic)
        .def("change_of_hand_op", &w_t::change_of_hand_op)
        .def("hall_symbol", &w_t::hall_symbol, hall_symbol_overloads())
        .def("lookup_symbol", &w_t::lookup_symbol)
        .def_pickle(space_group_type_wrappers())
      ;
    }
  };

} // namespace <anoymous>

  void wrap_space_group_type()
  {
    space_group_type_wrappers::wrap();
  }

}}} // namespace cctbx::sgtbx::boost_python
