/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2003 Jun: Created based on triplet.cpp (Ralf W. Grosse-Kunstleve)
 */

#include <boost/python/class.hpp>
#include <cctbx/dmtbx/triplet_phase_relation.h>
#include <scitbx/boost_python/container_conversions.h>
#include <scitbx/array_family/shared.h>

namespace cctbx { namespace dmtbx { namespace boost_python {
namespace {

  struct weighted_triplet_phase_relation_wrappers
  {
    typedef weighted_triplet_phase_relation w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("weighted_triplet_phase_relation", no_init)
        .def("ik", &w_t::ik)
        .def("friedel_flag_k", &w_t::friedel_flag_k)
        .def("ihmk", &w_t::ihmk)
        .def("friedel_flag_hmk", &w_t::friedel_flag_hmk)
        .def("ht_sum", &w_t::ht_sum)
        .def("weight", &w_t::weight)
      ;
    }
  };

  void register_tuple_mappings()
  {
    using namespace scitbx::boost_python::container_conversions;
    tuple_mapping<
      scitbx::af::shared<weighted_triplet_phase_relation>,
      variable_capacity_policy>();
  }

} // namespace <anonymous>

  void wrap_triplet_phase_relation()
  {
    weighted_triplet_phase_relation_wrappers::wrap();
    register_tuple_mappings();
  }

}}} // namespace cctbx::dmtbx::boost_python
