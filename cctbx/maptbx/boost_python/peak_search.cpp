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
#include <scitbx/boost_python/utils.h>

namespace cctbx { namespace maptbx { namespace boost_python {

namespace {

  struct peak_list_wrappers
  {
    typedef peak_list<> w_t;

    static af::tiny<long, 3>
    grid_indices(w_t const& self, long i)
    {
      using scitbx::boost_python::positive_getitem_index;
      std::size_t j = positive_getitem_index(i, self.grid_indices().size());
      return self.grid_indices()[j];
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("peak_list", no_init)
        .def(init<af::const_ref<float, af::c_grid_padded<3> > const&,
                  af::ref<long, af::c_grid<3> > const&,
                  int, std::size_t, bool>())
        .def(init<af::const_ref<float, af::c_grid_padded<3> > const&,
                  af::ref<long, af::c_grid<3> > const&,
                  int, float, std::size_t, bool>())
        .def(init<af::const_ref<double, af::c_grid_padded<3> > const&,
                  af::ref<long, af::c_grid<3> > const&,
                  int, std::size_t, bool>())
        .def(init<af::const_ref<double, af::c_grid_padded<3> > const&,
                  af::ref<long, af::c_grid<3> > const&,
                  int, double, std::size_t, bool>())
        .def("gridding", &w_t::gridding, ccr())
        .def("size", &w_t::size)
        .def("grid_indices", grid_indices)
        .def("grid_heights", &w_t::grid_heights)
        .def("sites", &w_t::sites)
        .def("heights", &w_t::heights)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_peak_list()
  {
    peak_list_wrappers::wrap();
  }

}}} // namespace cctbx::maptbx::boost_python
