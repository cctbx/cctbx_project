// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Jul: Created (R.W. Grosse-Kunstleve)
 */

#include <boost/python/cross_module.hpp>
#include <cctbx/miller_bpl.h>
#include <cctbx/dmtbx/triplet.h>

namespace {

  using namespace cctbx;

  void
  inplace_sort(af::shared<Miller::Index> miller_indices,
               af::shared<double> data,
               bool reverse)
  {
    if (reverse) {
      Miller::inplace_sort(miller_indices, data, std::greater<double>());
    }
    else {
      Miller::inplace_sort(miller_indices, data, std::less<double>());
    }
  }

  af::shared<double>
  refine_phases_3(dmtbx::triplet_invariants<double> const& ti,
                  af::shared<Miller::Index> miller_indices,
                  af::shared<double> e_values,
                  af::shared<double> phases)
  {
    return ti.refine_phases(miller_indices, e_values, phases);
  }

# include <cctbx/basic/from_bpl_import.h>

  void init_module(python::module_builder& this_module)
  {
    const std::string Revision = "$Revision$";
    this_module.add(ref(to_python(
        Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

    python::import_converters<sgtbx::SpaceGroupInfo>
    py_SpaceGroup("cctbx_boost.sgtbx", "SpaceGroupInfo");

    python::import_converters<af::shared<double> >
    py_shared_double("cctbx_boost.arraytbx.shared", "double");

    python::import_converters<af::shared<Miller::Index> >
    py_shared_Miller_Index(
      "cctbx_boost.arraytbx.shared", "Miller_Index");

    class_builder<dmtbx::triplet_invariants<double> >
    py_triplet_invariants(this_module, "triplet_invariants");

    py_triplet_invariants.def(constructor<>());
    py_triplet_invariants.def(constructor<
      sgtbx::SpaceGroupInfo const&,
      af::shared<Miller::Index>,
      af::shared<double> >());
    py_triplet_invariants.def(
      &dmtbx::triplet_invariants<double>::number_of_weighted_triplets,
                                         "number_of_weighted_triplets");
    py_triplet_invariants.def(
      &dmtbx::triplet_invariants<double>::total_number_of_triplets,
                                         "total_number_of_triplets");
    py_triplet_invariants.def(
      &dmtbx::triplet_invariants<double
        >::average_number_of_triplets_per_reflection,
          "average_number_of_triplets_per_reflection");
    py_triplet_invariants.def(
      &dmtbx::triplet_invariants<double>::dump_triplets,
                                         "dump_triplets");
    py_triplet_invariants.def(
      &dmtbx::triplet_invariants<double>::unique_triplets,
                                         "unique_triplets");
    py_triplet_invariants.def(
      &dmtbx::triplet_invariants<double>::refine_phases,
                                         "refine_phases");
    py_triplet_invariants.def(
                                         &refine_phases_3,
                                         "refine_phases");

    this_module.def(inplace_sort, "inplace_sort");
  }

}

BOOST_PYTHON_MODULE_INIT(dmtbx)
{
  boost::python::module_builder this_module("dmtbx");
  init_module(this_module);
}
