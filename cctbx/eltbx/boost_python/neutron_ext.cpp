/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Created (rwgk)
 */

#include <cctbx/eltbx/neutron.h>
#include <scitbx/boost_python/utils.h>
#include <scitbx/boost_python/iterator_wrappers.h>
#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/class.hpp>

namespace cctbx { namespace eltbx { namespace neutron {
namespace boost_python {

namespace {

  struct neutron_news_1992_table_wrappers
  {
    typedef neutron_news_1992_table w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("neutron_news_1992_table", no_init)
        .def(init<std::string const&, optional<bool> >())
        .def("label", &w_t::label)
        .def("bound_coh_scatt_length", &w_t::bound_coh_scatt_length)
        .def("abs_cross_sect", &w_t::abs_cross_sect)
      ;
    }
  };

  void init_module()
  {
    using namespace boost::python;

    scope().attr("__version__") = scitbx::boost_python::cvs_revision(
      "$Revision$");

    neutron_news_1992_table_wrappers::wrap();
    scitbx::boost_python::iterator_wrappers<
      neutron_news_1992_table,
      neutron_news_1992_table_iterator>::wrap(
        "neutron_news_1992_table_iterator");
  }

} // namespace <anonymous>
}}}} // namespace cctbx::eltbx::neutron::boost_python

BOOST_PYTHON_MODULE(neutron_ext)
{
  cctbx::eltbx::neutron::boost_python::init_module();
}
