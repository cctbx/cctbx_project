/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Created (rwgk)
 */

#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/maptbx/structure_factors.h>
#include <boost/python/class.hpp>

namespace cctbx { namespace maptbx { namespace structure_factors {

namespace {

  struct to_map_wrappers
  {
    typedef to_map<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("structure_factors_to_map", no_init)
        .def(init<sgtbx::space_group const&,
                  bool,
                  af::const_ref<miller::index<> > const&,
                  af::const_ref<std::complex<double> > const&,
                  af::flex_grid<> const&,
                  bool>())
        .def("complex_map", &w_t::complex_map)
      ;
    }
  };

  struct from_map_wrappers
  {
    typedef from_map<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("structure_factors_from_map", no_init)
        .def(init<uctbx::unit_cell const&,
                  sgtbx::space_group_type const&,
                  bool,
                  double,
                  af::const_ref<std::complex<double>,
                                af::c_grid_padded<3> > const&,
                  bool>())
        .def(init<bool,
                  af::const_ref<miller::index<> > const&,
                  af::const_ref<std::complex<double>,
                                af::c_grid_padded<3> > const&,
                  bool>())
        .def("miller_indices", &w_t::miller_indices)
        .def("data", &w_t::data)
      ;
    }
  };

}} // namespace structure_factors::<anoymous>

namespace boost_python {

  void wrap_structure_factors()
  {
    structure_factors::to_map_wrappers::wrap();
    structure_factors::from_map_wrappers::wrap();
  }

}}} // namespace cctbx::maptbx::boost_python
