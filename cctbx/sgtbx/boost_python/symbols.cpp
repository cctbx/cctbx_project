/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Created (rwgk)
 */

#include <cctbx/sgtbx/symbols.h>
#include <scitbx/boost_python/iterator_wrappers.h>
#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>

namespace cctbx { namespace sgtbx { namespace boost_python {

namespace {

  struct space_group_symbols_wrappers
  {
    typedef space_group_symbols w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("space_group_symbols", no_init)
        .def(init<std::string const&, optional<std::string const&> >())
        .def(init<int, optional<std::string const&, std::string const&> >())
        .def("number", &w_t::number)
        .def("schoenflies", &w_t::schoenflies, ccr())
        .def("qualifier", &w_t::qualifier, ccr())
        .def("hermann_mauguin", &w_t::hermann_mauguin, ccr())
        .def("extension", &w_t::extension)
        .def("extended_hermann_mauguin", &w_t::extended_hermann_mauguin, ccr())
        .def("hall", &w_t::hall, ccr())
      ;
    }
  };

} // namespace <anoymous>

  void wrap_symbols()
  {
    space_group_symbols_wrappers::wrap();
    scitbx::boost_python::iterator_wrappers<
      space_group_symbols,
      space_group_symbol_iterator>::wrap(
        "space_group_symbol_iterator");
  }

}}} // namespace cctbx::sgtbx::boost_python
