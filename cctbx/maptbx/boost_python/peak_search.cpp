/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Created (rwgk)
 */

#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/maptbx/peak_search.h>
#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <scitbx/boost_python/container_conversions.h>

namespace cctbx { namespace maptbx { namespace boost_python {

namespace {

  struct peak_list_wrappers
  {
    typedef peak_list<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("peak_list", no_init)
        .def(init<af::const_ref<float, af::c_grid_padded<3> > const&,
                  af::ref<long, af::c_grid<3> > const&,
                  optional<int, std::size_t> >())
        .def(init<af::const_ref<double, af::c_grid_padded<3> > const&,
                  af::ref<long, af::c_grid<3> > const&,
                  optional<int, std::size_t> >())
        .def("entries", &w_t::entries, ccr())
      ;
    }
  };

  struct peak_list_indexed_value_wrappers
  {
    typedef peak_list<>::indexed_value_type w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("peak_list_indexed_value", no_init)
        .def_readonly("index", &w_t::index)
        .def_readonly("value", &w_t::value)
      ;

      using namespace scitbx::boost_python::container_conversions;
      to_python_converter<af::shared<w_t>,
                 to_tuple<af::shared<w_t> > >();
    }
  };

} // namespace <anoymous>

  void wrap_peak_list()
  {
    peak_list_wrappers::wrap();
    peak_list_indexed_value_wrappers::wrap();
  }

}}} // namespace cctbx::maptbx::boost_python
