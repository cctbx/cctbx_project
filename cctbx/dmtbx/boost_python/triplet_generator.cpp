/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2003 Jun: Created based on triplet.cpp (Ralf W. Grosse-Kunstleve)
 */

#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/dmtbx/triplet_generator.h>
#include <boost/python/class.hpp>

namespace cctbx { namespace dmtbx { namespace boost_python {
namespace {

  struct triplet_generator_wrappers
  {
    typedef triplet_generator<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("triplet_generator", no_init)
        .def(init<sgtbx::space_group const&,
                  af::const_ref<miller::index<> > const&,
                  optional<bool, bool> >())
        .def("t_den", &w_t::t_den)
        .def("sigma_2_only", &w_t::sigma_2_only)
        .def("discard_weights", &w_t::discard_weights)
        .def("n_relations", &w_t::n_relations)
        .def("relations_for", &w_t::relations_for)
        .def("sums_of_amplitude_products", &w_t::sums_of_amplitude_products)
        .def("raw_apply_tangent_formula", &w_t::apply_tangent_formula)
      ;
    }
  };

} // namespace <anonymous>

  void wrap_triplet_generator()
  {
    triplet_generator_wrappers::wrap();
  }

}}} // namespace cctbx::dmtbx::boost_python
