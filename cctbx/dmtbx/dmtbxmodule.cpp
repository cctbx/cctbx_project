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
#include <cctbx/miller/sort.h>
#include <cctbx/dmtbx/triplet.h>

#include <cctbx/array_family/shared_bpl_.h>

namespace cctbx { namespace af { namespace bpl {

  void import_flex()
  {
    CCTBX_ARRAY_FAMILY_FLEX_IMPORT(double, "double")
    CCTBX_ARRAY_FAMILY_FLEX_IMPORT(cctbx::miller::Index, "miller_Index")
  }

}}} // namespace cctbx::af::bpl

CCTBX_ARRAY_FAMILY_IMPLICIT_SHARED_CONVERTERS(double)
CCTBX_ARRAY_FAMILY_IMPLICIT_SHARED_CONVERTERS(cctbx::miller::Index)

namespace {

  using namespace cctbx;

  void
  inplace_sort(af::shared<miller::Index> miller_indices,
               af::shared<double> data,
               bool reverse)
  {
    if (reverse) {
      miller::inplace_sort(miller_indices, data, std::greater<double>());
    }
    else {
      miller::inplace_sort(miller_indices, data, std::less<double>());
    }
  }

# include <cctbx/basic/from_bpl_import.h>

  void init_module(python::module_builder& this_module)
  {
    const std::string Revision = "$Revision$";
    this_module.add(ref(to_python(
        Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

    python::import_converters<sgtbx::SpaceGroupInfo>
    py_SpaceGroup("cctbx_boost.sgtbx", "SpaceGroupInfo");

    af::bpl::import_flex();

    class_builder<dmtbx::triplet_invariants<double> >
    py_triplet_invariants(this_module, "triplet_invariants");

    py_triplet_invariants.def(constructor<>());
    py_triplet_invariants.def(constructor<
      sgtbx::SpaceGroupInfo const&,
      af::shared<miller::Index>,
      bool,
      bool>());
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
      &dmtbx::triplet_invariants<double>::n_relations,
                                         "n_relations");
    py_triplet_invariants.def(
      &dmtbx::triplet_invariants<double>::apply_tangent_formula,
                                         "apply_tangent_formula");
    py_triplet_invariants.def(
      &dmtbx::triplet_invariants<double>::weights_and_epsilon,
                                         "weights_and_epsilon");

    this_module.def(inplace_sort, "inplace_sort");
  }

}

BOOST_PYTHON_MODULE_INIT(dmtbx)
{
  boost::python::module_builder this_module("dmtbx");
  init_module(this_module);
}
