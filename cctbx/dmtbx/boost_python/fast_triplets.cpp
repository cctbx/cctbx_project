/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2003 Jun: Created (Ralf W. Grosse-Kunstleve)
 */

#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/dmtbx/fast_triplets.h>
#include <boost/python/class.hpp>

namespace cctbx { namespace dmtbx { namespace boost_python {
namespace {

  struct fast_triplets_wrappers
  {
    typedef fast_triplets<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("fast_triplets", no_init)
        .def(init<sgtbx::space_group_type const&,
                  af::const_ref<miller::index<> > const&,
                  optional<bool> >())
        .def("sigma_2_only", &w_t::sigma_2_only)
        .def("number_of_weighted_triplets", &w_t::number_of_weighted_triplets)
        .def("total_number_of_triplets", &w_t::total_number_of_triplets)
        .def("average_number_of_triplets_per_reflection",
          &w_t::average_number_of_triplets_per_reflection)
        .def("dump_triplets", &w_t::dump_triplets)
        .def("pack_triplets", &w_t::pack_triplets)
        .def("n_relations", &w_t::n_relations)
        .def("sum_of_e_products", &w_t::sum_of_e_products)
        .def("apply_tangent_formula", &w_t::apply_tangent_formula)
      ;
    }
  };

} // namespace <anonymous>

  void wrap_fast_triplets()
  {
    fast_triplets_wrappers::wrap();
  }

}}} // namespace cctbx::dmtbx::boost_python
