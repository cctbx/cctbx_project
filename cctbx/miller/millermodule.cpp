// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jul 2002: Created, based on fragments from sharedmodule.cpp (rwgk)
 */

#include <boost/python/cross_module.hpp>
#include <cctbx/bpl_utils.h>
#include <cctbx/miller_bpl.h>
#include <cctbx/hendrickson_lattman_bpl.h>
#include <cctbx/array_family/shared_bpl.h>
#include <cctbx/miller/sym_equiv.h>
#include <cctbx/miller/asu.h>
#include <cctbx/miller/build.h>
#include <cctbx/miller/span.h>
#include <cctbx/miller/join.h>
#include <cctbx/miller/expand.h>
#include <cctbx/miller/bins.h>
#include <cctbx/miller/math.h>

#include <cctbx/array_family/flex_bpl.h>

namespace cctbx { namespace af { namespace bpl { namespace {

  typedef tiny<std::size_t, 2> tiny_size_t_2;

  void import_flex()
  {
    CCTBX_ARRAY_FAMILY_FLEX_IMPORT(bool, "bool")
    CCTBX_ARRAY_FAMILY_FLEX_IMPORT(std::size_t, "size_t")
    CCTBX_ARRAY_FAMILY_FLEX_IMPORT(double, "double")
    CCTBX_ARRAY_FAMILY_FLEX_IMPORT(cctbx::miller::Index, "miller_Index")
    CCTBX_ARRAY_FAMILY_FLEX_IMPORT(tiny_size_t_2, "tiny_size_t_2")
    CCTBX_ARRAY_FAMILY_FLEX_IMPORT(std::complex<double>, "complex_double");
    CCTBX_ARRAY_FAMILY_FLEX_IMPORT(
      cctbx::hendrickson_lattman<double>,
      "hendrickson_lattman");
  }

}}}} // namespace cctbx::af::bpl<anonymous>

CCTBX_ARRAY_FAMILY_IMPLICIT_SHARED_CONVERTERS(bool)
CCTBX_ARRAY_FAMILY_IMPLICIT_SHARED_CONVERTERS(std::size_t)
CCTBX_ARRAY_FAMILY_IMPLICIT_SHARED_CONVERTERS(double)
CCTBX_ARRAY_FAMILY_IMPLICIT_SHARED_CONVERTERS(cctbx::miller::Index)
CCTBX_ARRAY_FAMILY_IMPLICIT_SHARED_CONVERTERS(cctbx::af::bpl::tiny_size_t_2)
CCTBX_ARRAY_FAMILY_IMPLICIT_SHARED_CONVERTERS(std::complex<double>)
CCTBX_ARRAY_FAMILY_IMPLICIT_SHARED_CONVERTERS(
  cctbx::hendrickson_lattman<double>)

namespace {

  using namespace cctbx;
  using namespace cctbx::miller;

  double SymEquivIndex_HT_angle_0(
    SymEquivIndex const& SEI)
  {
    return SEI.HT_angle();
  }

  SymEquivIndex
  SymEquivIndex_Mate_0(
    const SymEquivIndex& SEI)
  {
    return SEI.Mate();
  }

  double SymEquivIndex_phase_eq_1(
    SymEquivIndex const& SEI,
    double phi_in)
  {
    return SEI.phase_eq(phi_in);
  }
  double SymEquivIndex_phase_eq_2(
    SymEquivIndex const& SEI,
    double phi_in,
    bool deg)
  {
    return SEI.phase_eq(phi_in, deg);
  }
  double SymEquivIndex_phase_in_1(
    const SymEquivIndex& SEI,
    double phi_eq)
  {
    return SEI.phase_in(phi_eq);
  }
  double SymEquivIndex_phase_in_2(
    const SymEquivIndex& SEI,
    double phi_eq,
    bool deg)
  {
    return SEI.phase_in(phi_eq, deg);
  }
  std::complex<double>
  SymEquivIndex_complex_eq(
    const SymEquivIndex& SEI,
    const std::complex<double>& f_in)
  {
    return SEI.complex_eq(f_in);
  }
  std::complex<double>
  SymEquivIndex_complex_in(
    const SymEquivIndex& SEI,
    const std::complex<double>& f_eq)
  {
    return SEI.complex_in(f_eq);
  }
  hendrickson_lattman<double>
  SymEquivIndex_hl_eq(
    SymEquivIndex const& SEI,
    hendrickson_lattman<double> const& coeff_in)
  {
    return SEI.hl_eq(coeff_in);
  }
  hendrickson_lattman<double>
  SymEquivIndex_hl_in(
    SymEquivIndex const& SEI,
    hendrickson_lattman<double> const& coeff_eq)
  {
    return SEI.hl_in(coeff_eq);
  }

  SymEquivIndex
  SymEquivIndices_getitem(const SymEquivIndices& SEMI,
                                std::size_t key) {
    if (key >= SEMI.N()) bpl_utils::raise_index_error();
    return SEMI[key];
  }
  SymEquivIndex
  SymEquivIndices_call_1(const SymEquivIndices& SEMI,
                               int iIL) {
    return SEMI(iIL);
  }
  SymEquivIndex
  SymEquivIndices_call_2(const SymEquivIndices& SEMI,
                               int iInv, int iList) {
    return SEMI(iInv, iList);
  }
  bool
  SymEquivIndices_isValidPhase_2(const SymEquivIndices& SEMI,
                                       double phi, bool deg) {
    return SEMI.isValidPhase(phi, deg);
  }
  bool
  SymEquivIndices_isValidPhase_1(const SymEquivIndices& SEMI,
                                       double phi) {
    return SEMI.isValidPhase(phi);
  }

  Index
  IndexGenerator_getitem(IndexGenerator& MIG, std::size_t)
  {
    Index result = MIG.next();
    if (result.is000()) {
      PyErr_SetString(PyExc_IndexError, "End of list.");
      boost::python::throw_error_already_set();
    }
    return result;
  }

  af::shared<Index>
  py_BuildIndices_Resolution_d_min(
    const uctbx::UnitCell& UC,
    const sgtbx::SpaceGroupInfo& SgInfo,
    bool FriedelFlag,
    double Resolution_d_min)
  {
    af::shared<Index> result;
    IndexGenerator(
      UC, SgInfo, FriedelFlag, Resolution_d_min).AddToArray(result);
    return result;
  }
  af::shared<Index>
  py_BuildIndices_MaxIndex(
    const sgtbx::SpaceGroupInfo& SgInfo,
    bool FriedelFlag,
    const Index& MaxIndex)
  {
    af::shared<Index> result;
    IndexGenerator(
      SgInfo, FriedelFlag, MaxIndex).AddToArray(result);
    return result;
  }

  template <typename ValueType>
  struct map_to_asu_wrappers
  {
    static
    void
    map_to_asu_no_bool(
      sgtbx::SpaceGroupInfo const& sginfo,
      bool friedel_flag,
      af::shared<Index> miller_indices,
      af::shared<ValueType> data_array)
    {
      map_to_asu(sginfo, friedel_flag, miller_indices, data_array);
    }

    static
    void
    map_to_asu_with_bool(
      sgtbx::SpaceGroupInfo const& sginfo,
      bool friedel_flag,
      af::shared<Index> miller_indices,
      af::shared<ValueType> data_array,
      bool deg)
    {
      map_to_asu(sginfo, friedel_flag, miller_indices, data_array, deg);
    }
  };

  af::shared<double>
  join_sets_plus(
    join_sets const& js, af::shared<double> data0, af::shared<double> data1)
  {
    return js.plus(data0, data1);
  }
  af::shared<double>
  join_sets_minus(
    join_sets const& js, af::shared<double> data0, af::shared<double> data1)
  {
    return js.minus(data0, data1);
  }
  af::shared<double>
  join_sets_multiplies(
    join_sets const& js, af::shared<double> data0, af::shared<double> data1)
  {
    return js.multiplies(data0, data1);
  }
  af::shared<double>
  join_sets_divides(
    join_sets const& js, af::shared<double> data0, af::shared<double> data1)
  {
    return js.divides(data0, data1);
  }

  af::shared<double>
  join_sets_additive_sigmas(
    join_sets const& js,
    af::shared<double> sigmas0,
    af::shared<double> sigmas1)
  {
    return js.additive_sigmas(sigmas0, sigmas1);
  }

  af::shared<double>
  join_bijvoet_mates_minus(
    join_bijvoet_mates const& jbm,
    af::shared<double> data)
  {
    return jbm.minus(data);
  }

  af::shared<double>
  join_bijvoet_mates_additive_sigmas(
    join_bijvoet_mates const& jbm,
    af::shared<double> sigmas)
  {
    return jbm.additive_sigmas(sigmas);
  }

  af::shared<double>
  join_bijvoet_mates_average(
    join_bijvoet_mates const& jbm,
    af::shared<double> sigmas)
  {
    return jbm.average(sigmas);
  }

  struct binning_
  {
    static
    boost::python::ref
    range_used(binning const& bng)
    {
      return boost::python::ref(PyRange_New(1, bng.n_bins_used(), 1, 1));
    }

    static
    boost::python::ref
    range_all(binning const& bng)
    {
      return boost::python::ref(PyRange_New(0, bng.n_bins_all(), 1, 1));
    }

    static std::size_t
    get_i_bin_d(binning const& bng, double d_star_sq)
    {
      return bng.get_i_bin(d_star_sq);
    }

    static std::size_t
    get_i_bin_i(binning const& bng, Index const& h)
    {
      return bng.get_i_bin(h);
    }
  };

  void
  py_expand_to_p1_4(
    sgtbx::SpaceGroup const& SgOps,
    bool friedel_flag,
    const af::shared<Index>& in,
    af::versa<Index, af::flex_grid<> >& out)
  {
    af::shared<Index> out_ = bpl_utils::as_base_array(out);
    expand_to_p1(SgOps, friedel_flag, in, out_);
    out.resize(af::flex_grid<>(out_.size()));
  }

  void
  py_expand_to_p1_9(
    sgtbx::SpaceGroup const& SgOps,
    bool friedel_flag,
    af::shared<Index> const& h_in,
    af::shared<double> const& ampl_in,
    af::shared<double> const& phase_in,
    af::versa<Index, af::flex_grid<> >& h_out,
    af::versa<double, af::flex_grid<> >& ampl_out,
    af::versa<double, af::flex_grid<> >& phase_out,
    bool phase_degrees)
  {
    af::shared<Index> h_out_ = bpl_utils::as_base_array(h_out);
    af::shared<double> ampl_out_ = bpl_utils::as_base_array(ampl_out);
    af::shared<double> phase_out_ = bpl_utils::as_base_array(phase_out);
    expand_to_p1(
      SgOps, friedel_flag,
      h_in, ampl_in, phase_in,
      h_out_, ampl_out_, phase_out_,
      phase_degrees);
    h_out.resize(af::flex_grid<>(h_out_.size()));
    ampl_out.resize(af::flex_grid<>(ampl_out_.size()));
    phase_out.resize(af::flex_grid<>(phase_out_.size()));
  }

  void
  py_expand_to_p1_8(
    sgtbx::SpaceGroup const& SgOps,
    bool friedel_flag,
    af::shared<Index> const& h_in,
    af::shared<double> const& ampl_in,
    af::shared<double> const& phase_in,
    af::versa<Index, af::flex_grid<> >& h_out,
    af::versa<double, af::flex_grid<> >& ampl_out,
    af::versa<double, af::flex_grid<> >& phase_out)
  {
    py_expand_to_p1_9(
      SgOps, friedel_flag,
      h_in, ampl_in, phase_in,
      h_out, ampl_out, phase_out, false);
  }

  double
  statistical_mean_double(
    sgtbx::SpaceGroup const& SgOps,
    bool friedel_flag,
    af::shared<Index> miller_indices,
    af::shared<double> data)
  {
    return statistical_mean(SgOps, friedel_flag, miller_indices, data);
  }

# include <cctbx/basic/from_bpl_import.h>

  void init_module(python::module_builder& this_module)
  {
    const std::string Revision = "$Revision$";
    this_module.add(ref(to_python(
        Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

    python::import_converters<uctbx::UnitCell>
    py_UnitCell("cctbx_boost.uctbx", "UnitCell");

    python::import_converters<sgtbx::PhaseInfo>
    py_PhaseInfo("cctbx_boost.sgtbx", "PhaseInfo");

    python::import_converters<sgtbx::SpaceGroup>
    py_SpaceGroup("cctbx_boost.sgtbx", "SpaceGroup");

    python::import_converters<sgtbx::SpaceGroupInfo>
    py_SpaceGroupInfo("cctbx_boost.sgtbx", "SpaceGroupInfo");

    python::import_converters<sgtbx::ReciprocalSpaceASU>
    py_ReciprocalSpaceASU("cctbx_boost.sgtbx", "ReciprocalSpaceASU");

    af::bpl::import_flex();

    class_builder<SymEquivIndex>
    py_SymEquivIndex(this_module, "SymEquivIndex");

    class_builder<SymEquivIndices>
    py_SymEquivIndices(this_module, "SymEquivIndices");

    class_builder<IndexTableLayoutAdaptor>
    py_IndexTableLayoutAdaptor(this_module, "IndexTableLayoutAdaptor");

    class_builder<AsymIndex>
    py_AsymIndex(this_module, "AsymIndex");

    class_builder<IndexGenerator>
    py_IndexGenerator(this_module, "IndexGenerator");

    class_builder<index_span>
    py_index_span(this_module, "index_span");

    class_builder<join_sets>
    py_join_sets(this_module, "join_sets");

    class_builder<join_bijvoet_mates>
    py_join_bijvoet_mates(this_module, "join_bijvoet_mates");

    class_builder<binning>
    py_binning(this_module, "binning");

    class_builder<binner>
    py_binner(this_module, "binner");

    py_IndexTableLayoutAdaptor.declare_base(
      py_SymEquivIndex, python::without_downcast);

    py_AsymIndex.declare_base(
      py_SymEquivIndex, python::without_downcast);

    py_binner.declare_base(
      py_binning, python::without_downcast);

    py_SymEquivIndex.def(constructor<>());
    py_SymEquivIndex.def(
      constructor<const Index&, int, int, bool>());
    py_SymEquivIndex.def(&SymEquivIndex::H, "H");
    py_SymEquivIndex.def(&SymEquivIndex::HR, "HR");
    py_SymEquivIndex.def(&SymEquivIndex::HT, "HT");
    py_SymEquivIndex.def(&SymEquivIndex::TBF, "TBF");
    py_SymEquivIndex.def(&SymEquivIndex::HT_angle, "HT_angle");
    py_SymEquivIndex.def(SymEquivIndex_HT_angle_0, "HT_angle");
    py_SymEquivIndex.def(
      &SymEquivIndex::FriedelFlag, "FriedelFlag");
    py_SymEquivIndex.def(SymEquivIndex_Mate_0, "Mate");
    py_SymEquivIndex.def(&SymEquivIndex::Mate, "Mate");
    py_SymEquivIndex.def(
      SymEquivIndex_phase_eq_1, "phase_eq");
    py_SymEquivIndex.def(
      SymEquivIndex_phase_eq_2, "phase_eq");
    py_SymEquivIndex.def(
      SymEquivIndex_phase_in_1, "phase_in");
    py_SymEquivIndex.def(
      SymEquivIndex_phase_in_2, "phase_in");
    py_SymEquivIndex.def(
      SymEquivIndex_complex_eq, "complex_eq");
    py_SymEquivIndex.def(
      SymEquivIndex_complex_in, "complex_in");
    py_SymEquivIndex.def(
      SymEquivIndex_hl_eq, "hl_eq");
    py_SymEquivIndex.def(
      SymEquivIndex_hl_in, "hl_in");

    py_SymEquivIndices.def(constructor<>());
    py_SymEquivIndices.def(constructor<
      sgtbx::SpaceGroup const&,
      Index const&>());
    py_SymEquivIndices.def(
      &SymEquivIndices::getPhaseRestriction, "getPhaseRestriction");
    py_SymEquivIndices.def(
      &SymEquivIndices::isCentric, "isCentric");
    py_SymEquivIndices.def(&SymEquivIndices::N, "N");
    py_SymEquivIndices.def(&SymEquivIndices::M, "M");
    py_SymEquivIndices.def(&SymEquivIndices::fMates, "fMates");
    py_SymEquivIndices.def(&SymEquivIndices::epsilon, "epsilon");
    py_SymEquivIndices.def(&SymEquivIndices::N, "__len__");
    py_SymEquivIndices.def(SymEquivIndices_getitem, "__getitem__");
    py_SymEquivIndices.def(SymEquivIndices_call_1, "__call__");
    py_SymEquivIndices.def(SymEquivIndices_call_2, "__call__");
    py_SymEquivIndices.def(
      &SymEquivIndices::isValidPhase, "isValidPhase");
    py_SymEquivIndices.def(
      SymEquivIndices_isValidPhase_2, "isValidPhase");
    py_SymEquivIndices.def(
      SymEquivIndices_isValidPhase_1, "isValidPhase");

    py_IndexTableLayoutAdaptor.def(constructor<>());
    py_IndexTableLayoutAdaptor.def(
      &IndexTableLayoutAdaptor::H, "H");
    py_IndexTableLayoutAdaptor.def(
      &IndexTableLayoutAdaptor::iColumn, "iColumn");

    py_AsymIndex.def(constructor<>());
    py_AsymIndex.def(constructor<
      const sgtbx::SpaceGroup&,
      const sgtbx::ReciprocalSpaceASU&,
      const Index&>());
    py_AsymIndex.def(constructor<
      const sgtbx::SpaceGroup&,
      const Index&>());
    py_AsymIndex.def(constructor<
      const SymEquivIndices&>());
    py_AsymIndex.def(
      &AsymIndex::one_column, "one_column");
    py_AsymIndex.def(
      &AsymIndex::two_column, "two_column");

    py_IndexGenerator.def(constructor<>());
    py_IndexGenerator.def(constructor<const uctbx::UnitCell&,
                                      const sgtbx::SpaceGroupInfo&,
                                      bool,
                                      double>());
    py_IndexGenerator.def(constructor<const sgtbx::SpaceGroupInfo&,
                                      bool,
                                      const Index&>());
    py_IndexGenerator.def(&IndexGenerator::next, "next");
    py_IndexGenerator.def(IndexGenerator_getitem, "__getitem__");
    py_IndexGenerator.def(&IndexGenerator::ASU, "ASU");

    this_module.def(py_BuildIndices_Resolution_d_min, "BuildIndices");
    this_module.def(py_BuildIndices_MaxIndex, "BuildIndices");

    this_module.def(
       map_to_asu_wrappers<double>::map_to_asu_no_bool,
      "map_to_asu");
    this_module.def(
      map_to_asu_wrappers<double>::map_to_asu_with_bool,
      "map_to_asu");
    this_module.def(
       map_to_asu_wrappers<std::complex<double> >::map_to_asu_no_bool,
      "map_to_asu");
    this_module.def(
       map_to_asu_wrappers<hendrickson_lattman<double> >::map_to_asu_no_bool,
      "map_to_asu");

    py_index_span.def(constructor<>());
    py_index_span.def(constructor<af::shared<Index> >());
    py_index_span.def(&index_span::min, "min");
    py_index_span.def(&index_span::max, "max");
    py_index_span.def(&index_span::abs_range, "abs_range");
    py_index_span.def(&index_span::map_grid, "map_grid");
    py_index_span.def(&index_span::is_in_domain, "is_in_domain");

    py_join_sets.def(constructor<>());
    py_join_sets.def(constructor<
      af::shared<Index>,
      af::shared<Index> >());
    py_join_sets.def(&join_sets::pairs, "pairs");
    py_join_sets.def(&join_sets::singles, "singles");
    py_join_sets.def(&join_sets::have_singles, "have_singles");
    py_join_sets.def(&join_sets::size_processed, "size_processed");
    py_join_sets.def(&join_sets::pair_selection,
                                "pair_selection");
    py_join_sets.def(&join_sets::single_selection,
                                "single_selection");
    py_join_sets.def(&join_sets::paired_miller_indices,
                                "paired_miller_indices");
    py_join_sets.def(join_sets_plus, "plus");
    py_join_sets.def(join_sets_minus, "minus");
    py_join_sets.def(join_sets_multiplies, "multiplies");
    py_join_sets.def(join_sets_divides, "divides");
    py_join_sets.def(join_sets_additive_sigmas, "additive_sigmas");

    py_join_bijvoet_mates.def(constructor<>());
    py_join_bijvoet_mates.def(constructor<
      sgtbx::SpaceGroupInfo const&, af::shared<Index> >());
    py_join_bijvoet_mates.def(constructor<
      cctbx::sgtbx::ReciprocalSpaceASU const&, af::shared<Index> >());
    py_join_bijvoet_mates.def(constructor<
      af::shared<Index> >());
    py_join_bijvoet_mates.def(&join_bijvoet_mates::pairs, "pairs");
    py_join_bijvoet_mates.def(&join_bijvoet_mates::singles, "singles");
    py_join_bijvoet_mates.def(&join_bijvoet_mates::have_singles,
                                                  "have_singles");
    py_join_bijvoet_mates.def(&join_bijvoet_mates::size_processed,
                                                  "size_processed");
    py_join_bijvoet_mates.def(
      &join_bijvoet_mates::miller_indices_in_hemisphere,
                          "miller_indices_in_hemisphere");
    py_join_bijvoet_mates.def(join_bijvoet_mates_minus, "minus");
    py_join_bijvoet_mates.def(join_bijvoet_mates_additive_sigmas,
                                                "additive_sigmas");
    py_join_bijvoet_mates.def(join_bijvoet_mates_average, "average");

    py_binning.def(constructor<>());
    py_binning.def(constructor<
      uctbx::UnitCell const&,
      std::size_t,
      af::shared<Index>,
      double,
      double,
      double>());
    py_binning.def(constructor<
      uctbx::UnitCell const&,
      std::size_t,
      af::shared<Index>,
      double,
      double>());
    py_binning.def(constructor<
      uctbx::UnitCell const&,
      std::size_t,
      af::shared<Index>,
      double>());
    py_binning.def(constructor<
      uctbx::UnitCell const&,
      std::size_t,
      af::shared<Index> >());
    py_binning.def(constructor<
      uctbx::UnitCell const&,
      std::size_t,
      double,
      double>());
    py_binning.def(&binning::unit_cell, "unit_cell");
    py_binning.def(&binning::n_bins_used, "n_bins_used");
    py_binning.def(&binning::n_bins_all, "n_bins_all");
    py_binning.def(binning_::range_used, "range_used");
    py_binning.def(binning_::range_all, "range_all");
    py_binning.def(&binning::i_bin_d_too_large, "i_bin_d_too_large");
    py_binning.def(&binning::i_bin_d_too_small, "i_bin_d_too_small");
    py_binning.def(&binning::d_max, "d_max");
    py_binning.def(&binning::d_min, "d_min");
    py_binning.def(&binning::bin_d_range, "bin_d_range");
    py_binning.def(&binning::bin_d_min, "bin_d_min");
    py_binning.def(&binning::limits, "limits");
    py_binning.def(binning_::get_i_bin_d, "get_i_bin");
    py_binning.def(binning_::get_i_bin_i, "get_i_bin");

    py_binner.def(constructor<>());
    py_binner.def(constructor<binning const&, af::shared<Index> >());
    py_binner.def(&binner::bin_indices, "bin_indices");
    py_binner.def(&binner::count, "count");
    py_binner.def(&binner::counts, "counts");
    py_binner.def(&binner::operator(), "__call__");
    py_binner.def(&binner::array_indices, "array_indices");

    this_module.def(py_expand_to_p1_4, "expand_to_p1");
    this_module.def(py_expand_to_p1_9, "expand_to_p1");
    this_module.def(py_expand_to_p1_8, "expand_to_p1");

    this_module.def(statistical_mean_double, "statistical_mean");
  }

}

BOOST_PYTHON_MODULE_INIT(miller)
{
  boost::python::module_builder this_module("miller");
  init_module(this_module);
}
