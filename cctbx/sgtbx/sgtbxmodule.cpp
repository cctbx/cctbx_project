// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Oct 12: SpecialPosition -> SiteSymmetry (R.W. Grosse-Kunstleve)
     2001 Sep 13: SpaceGroupType -> SpaceGroupInfo (R.W. Grosse-Kunstleve)
     2001 Jul 02: Merged from CVS branch sgtbx_special_pos (rwgk)
     2001 May 31: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <cctbx/sgtbx/symbols.h>
#include <cctbx/sgtbx/groups.h>
#include <cctbx/sgtbx/coordinates.h>
#include <cctbx/sgtbx/asym_units.h>
#include <cctbx/sgtbx/miller_asu.h>
#include <cctbx/sgtbx/seminvariant.h>
#include <boost/python/cross_module.hpp>
#include <cctbx/coordinates_bpl.h>
#include <cctbx/miller_bpl.h>
#include <cctbx/hendrickson_lattman_bpl.h>

using namespace cctbx;
using namespace cctbx::sgtbx;

namespace {

# include <cctbx/basic/from_bpl_import.h>

  int rational_numerator(const rational& r) {
    return r.numerator();
  }
  int rational_denominator(const rational& r) {
    return r.denominator();
  }

  std::string rational_format_0(const rational& r) {
    return r.format();
  }
  std::string rational_format_1(const rational& r, bool Decimal) {
    return r.format(Decimal);
  }

  SpaceGroupSymbols
  SpaceGroupSymbolIterator_getitem(SpaceGroupSymbolIterator& iter,
                                   std::size_t) {
    SpaceGroupSymbols result = iter.next();
    if (result.SgNumber() == 0) {
      PyErr_SetString(PyExc_IndexError, "At end of table.");
      throw python::error_already_set();
    }
    return result;
  }

  RTMx TranslationComponents_IntrinsicPart(const TranslationComponents& tc) {
    return RTMx(RotMx(1, 0), tc.IntrinsicPart());
  }
  RTMx TranslationComponents_LocationPart(const TranslationComponents& tc) {
    return RTMx(RotMx(1, 0), tc.LocationPart());
  }
  RTMx TranslationComponents_OriginShift(const TranslationComponents& tc) {
    return RTMx(RotMx(1, 0), tc.OriginShift());
  }

  RTMx RTMx_newBaseFactors(const RTMx& M, int RBF, int TBF) {
    return M.newBaseFactors(RBF, TBF);
  }

  std::string RTMx_as_xyz_0(const RTMx& M) {
    return M.as_xyz();
  }

  tuple RTMx_as_tuple_0(const RTMx& M) {
    tuple result(2);
    result.set_item(0,
      ref(to_python(M.Rpart().as_array(double()))));
    result.set_item(1,
      ref(to_python(M.Tpart().as_array(double()))));
    return result;
  }

  tuple RTMx_as_tuple_2(const RTMx& M, int RBF, int TBF) {
    RTMx N = M.newBaseFactors(RBF, TBF);
    tuple result(2);
    result.set_item(0, ref(to_python(N.Rpart().vec())));
    result.set_item(1, ref(to_python(N.Tpart().vec())));
    return result;
  }

  RTMx RTMx_from_tuple(const tuple& t, int RBF, int TBF) {
    if (t.size() != 2) {
      PyErr_SetString(PyExc_ValueError,
        "First argument must be a tuple that contains two tuples"
        " (rotation part and translation part)");
      throw boost::python::error_already_set();
    }
    af::int9 Rpart = from_python(
      t[0].get(), boost::python::type<const af::int9&>());
    af::int3 Tpart = from_python(
      t[1].get(), boost::python::type<const af::int3&>());
    return RTMx(RotMx(Rpart, RBF), TrVec(Tpart, TBF));
  }

  RTMx RTMx_mul_RTMx(const RTMx& lhs, const RTMx& rhs) {
    return lhs * rhs;
  }
  RTMx RTMx_multiply_RTMx(const RTMx& lhs, const RTMx& rhs) {
    return lhs.multiply(rhs);
  }
  af::double3
  RTMx_multiply_int3(const RTMx& lhs, const af::double3& rhs) {
    return lhs * rhs;
  }

  fractional<double>
  ChOfBasisOp_call(const ChOfBasisOp& CBOp, const fractional<double>& X) {
    return CBOp(X);
  }

  RTMx
  ChOfBasisOp_apply_RTMx(const ChOfBasisOp& CBOp, const RTMx& RT) {
    return CBOp.apply(RT);
  }
  uctbx::UnitCell
  ChOfBasisOp_apply_UnitCell(const ChOfBasisOp& CBOp,
                             const uctbx::UnitCell& uc) {
    return CBOp.apply(uc);
  }

  int SpaceGroup_ParseHallSymbol_parse_string(SpaceGroup& SgOps,
                                              parse_string& HSym) {
    return SgOps.ParseHallSymbol(HSym);
  }
  int SpaceGroup_ParseHallSymbol_std_string(SpaceGroup& SgOps,
                                            const std::string& HSym) {
    parse_string HSymPS(HSym);
    return SgOps.ParseHallSymbol(HSymPS);
  }

  void SpaceGroup_expandLTr(SpaceGroup& SgOps,
                            const af::int3& vector, int base_factor)
  {
    SgOps.expandLTr(TrVec(vector, base_factor).newBaseFactor(SgOps.TBF()));
  }

  RTMx SpaceGroup_call_1(const SpaceGroup& SgOps, int iLIS) {
    return SgOps(iLIS);
  }
  RTMx
  SpaceGroup_call_3(const SpaceGroup& SgOps, int iLTr, int iInv, int iSMx) {
    return SgOps(iLTr, iInv, iSMx);
  }
  RTMx SpaceGroup_getitem(const SpaceGroup& SgOps, std::size_t key) {
    try {
      return SgOps(key);
    }
    catch (const error_index& e) {
      PyErr_SetString(PyExc_IndexError, e.what());
      throw boost::python::error_already_set();
    }
  }

  bool
  SpaceGroup_isSysAbsent_s(const SpaceGroup& SgOps,
                         const Miller::Index& H) {
    return SgOps.isSysAbsent(H);
  }
  af::shared<bool>
  SpaceGroup_isSysAbsent_a(const SpaceGroup& SgOps,
                         const af::shared<Miller::Index>& H) {
    return SgOps.isSysAbsent(H);
  }

  bool SpaceGroup_isCentric_0(const SpaceGroup& SgOps) {
    return SgOps.isCentric();
  }
  bool
  SpaceGroup_isCentric_s(const SpaceGroup& SgOps,
                         const Miller::Index& H) {
    return SgOps.isCentric(H);
  }
  af::shared<bool>
  SpaceGroup_isCentric_a(const SpaceGroup& SgOps,
                         const af::shared<Miller::Index>& H) {
    return SgOps.isCentric(H);
  }

  bool SpaceGroup_isValidPhase_3(const SpaceGroup& SgOps,
                                 const Miller::Index& H, double phi,
                                 bool deg) {
    return SgOps.isValidPhase(H, phi, deg);
  }
  bool SpaceGroup_isValidPhase_2(const SpaceGroup& SgOps,
                                 const Miller::Index& H, double phi) {
    return SgOps.isValidPhase(H, phi);
  }

  int
  SpaceGroup_epsilon_s(const SpaceGroup& SgOps,
                       const Miller::Index& H) {
    return SgOps.epsilon(H);
  }
  af::shared<int>
  SpaceGroup_epsilon_a(const SpaceGroup& SgOps,
                       const af::shared<Miller::Index>& H) {
    return SgOps.epsilon(H);
  }

  int
  SpaceGroup_multiplicity_s(const SpaceGroup& SgOps,
                            const Miller::Index& H,
                            bool FriedelFlag) {
    return SgOps.multiplicity(H, FriedelFlag);
  }
  af::shared<int>
  SpaceGroup_multiplicity_a(const SpaceGroup& SgOps,
                            const af::shared<Miller::Index>& H,
                            bool FriedelFlag) {
    return SgOps.multiplicity(H, FriedelFlag);
  }

  void SpaceGroup_CheckMetricalMatrix_1(const SpaceGroup& SgOps,
                                        const af::double9& G) {
    SgOps.CheckMetricalMatrix(G);
  }
  void SpaceGroup_CheckUnitCell_1(const SpaceGroup& SgOps,
                                  const uctbx::UnitCell& uc) {
    SgOps.CheckUnitCell(uc);
  }

  int SpaceGroup_cmp_equal(const SpaceGroup& lhs, const SpaceGroup& rhs) {
    if (lhs == rhs) return 0;
    return 1;
  }

  ChOfBasisOp SpaceGroup_getZ2POp_0(const SpaceGroup& SgOps) {
    return SgOps.getZ2POp();
  }

  af::int3 SpaceGroup_refine_gridding_0(
    const SpaceGroup& SgOps) {
    return SgOps.refine_gridding();
  }
  af::int3 SpaceGroup_refine_gridding_1(
    const SpaceGroup& SgOps,
    const af::int3& grid) {
    return SgOps.refine_gridding(grid);
  }

  python::string SpaceGroup_repr(const SpaceGroup& SgOps)
  {
    std::string result;
    char buf[256];
    int i;
    sprintf(buf, "nLTr=%d\n", SgOps.nLTr());
    result += buf;
    for (i = 0; i < SgOps.nLTr(); i++)
      result += "  " + SgOps(i, 0, 0).as_xyz() + "\n";
    sprintf(buf, "fInv=%d\n", SgOps.fInv());
    result += buf;
    if (SgOps.isCentric())
      result += "  " + SgOps(0, 1, 0).as_xyz() + "\n";
    sprintf(buf, "nSMx=%d\n", SgOps.nSMx());
    result += buf;
    for (i = 0; i < SgOps.nSMx(); i++)
      result += "  " + SgOps(0, 0, i).as_xyz() + "\n";
    return python::string(result.c_str());
  }

  // Support for pickle.
  tuple SpaceGroup_getinitargs(const SpaceGroup& SgOps) {
    tuple initargs(1);
    initargs.set_item(0,
      ref(to_python(SpaceGroupInfo(SgOps).BuildHallSymbol())));
    return initargs;
  }

  std::string SpaceGroupInfo_BuildHallSymbol_0(const SpaceGroupInfo& SgInfo) {
    return SgInfo.BuildHallSymbol();
  }
  std::string SpaceGroupInfo_BuildHallSymbol_1(const SpaceGroupInfo& SgInfo,
                                               bool TidyCBOp) {
    return SgInfo.BuildHallSymbol(TidyCBOp);
  }

  double PhaseInfo_HT_angle_0(const PhaseInfo& phinfo)
  {
    return phinfo.HT_angle();
  }

  bool PhaseInfo_isValidPhase_2(const PhaseInfo& phinfo, double phi, bool deg)
  {
    return phinfo.isValidPhase(phi, deg);
  }
  bool PhaseInfo_isValidPhase_1(const PhaseInfo& phinfo, double phi)
  {
    return phinfo.isValidPhase(phi);
  }
  double PhaseInfo_nearest_valid_phase_1(const PhaseInfo& phinfo, double phi)
  {
    return phinfo.nearest_valid_phase(phi);
  }

  double Miller_SymEquivIndex_HT_angle_0(
    Miller::SymEquivIndex const& SEI)
  {
    return SEI.HT_angle();
  }

  Miller::SymEquivIndex Miller_SymEquivIndex_Mate_0(
    const Miller::SymEquivIndex& SEI)
  {
    return SEI.Mate();
  }

  double Miller_SymEquivIndex_phase_eq_1(
    const Miller::SymEquivIndex& SEI,
    double phi_in)
  {
    return SEI.phase_eq(phi_in);
  }
  double Miller_SymEquivIndex_phase_eq_2(
    const Miller::SymEquivIndex& SEI,
    double phi_in,
    bool deg)
  {
    return SEI.phase_eq(phi_in, deg);
  }
  double Miller_SymEquivIndex_phase_in_1(
    const Miller::SymEquivIndex& SEI,
    double phi_eq)
  {
    return SEI.phase_in(phi_eq);
  }
  double Miller_SymEquivIndex_phase_in_2(
    const Miller::SymEquivIndex& SEI,
    double phi_eq,
    bool deg)
  {
    return SEI.phase_in(phi_eq, deg);
  }
  std::complex<double>
  Miller_SymEquivIndex_complex_eq(
    const Miller::SymEquivIndex& SEI,
    const std::complex<double>& f_in)
  {
    return SEI.complex_eq(f_in);
  }
  std::complex<double>
  Miller_SymEquivIndex_complex_in(
    const Miller::SymEquivIndex& SEI,
    const std::complex<double>& f_eq)
  {
    return SEI.complex_in(f_eq);
  }
  hendrickson_lattman<double>
  Miller_SymEquivIndex_hl_eq(
    Miller::SymEquivIndex const& SEI,
    hendrickson_lattman<double> const& coeff_in)
  {
    return SEI.hl_eq(coeff_in);
  }
  hendrickson_lattman<double>
  Miller_SymEquivIndex_hl_in(
    Miller::SymEquivIndex const& SEI,
    hendrickson_lattman<double> const& coeff_eq)
  {
    return SEI.hl_in(coeff_eq);
  }

  Miller::SymEquivIndex
  SymEquivMillerIndices_getitem(const SymEquivMillerIndices& SEMI,
                                std::size_t key) {
    if (key >= SEMI.N()) {
      PyErr_SetString(PyExc_IndexError,
                      "SymEquivMillerIndices: index out of range");
      throw boost::python::error_already_set();
    }
    return SEMI[key];
  }
  Miller::SymEquivIndex
  SymEquivMillerIndices_call_1(const SymEquivMillerIndices& SEMI,
                               int iIL) {
    return SEMI(iIL);
  }
  Miller::SymEquivIndex
  SymEquivMillerIndices_call_2(const SymEquivMillerIndices& SEMI,
                               int iInv, int iList) {
    return SEMI(iInv, iList);
  }
  bool
  SymEquivMillerIndices_isValidPhase_2(const SymEquivMillerIndices& SEMI,
                                       double phi, bool deg) {
    return SEMI.isValidPhase(phi, deg);
  }
  bool
  SymEquivMillerIndices_isValidPhase_1(const SymEquivMillerIndices& SEMI,
                                       double phi) {
    return SEMI.isValidPhase(phi);
  }

  RTMx WyckoffPosition_getitem(const WyckoffPosition& WP,
                               std::size_t key) {
    WP.CheckExpanded();
    try {
      return WP(key);
    }
    catch (const error_index& e) {
      PyErr_SetString(PyExc_IndexError, e.what());
      throw boost::python::error_already_set();
    }
  }

  WyckoffPosition WyckoffTable_call_size_t(const WyckoffTable& WTab,
                                           std::size_t i) {
    return WTab(i);
  }
  WyckoffPosition WyckoffTable_call_char(const WyckoffTable& WTab,
                                         char Letter) {
    return WTab(Letter);
  }

  WyckoffPosition WyckoffTable_getitem(const WyckoffTable& WTab,
                                       std::size_t key) {
    try {
      return WTab(key);
    }
    catch (const error_index& e) {
      PyErr_SetString(PyExc_IndexError, e.what());
      throw boost::python::error_already_set();
    }
  }

  WyckoffMapping
  WyckoffTable_getWyckoffMapping_SS(const WyckoffTable& WTab,
                                    const SiteSymmetry& SS) {
    return WTab.getWyckoffMapping(SS);
  }
  WyckoffMapping
  WyckoffTable_getWyckoffMapping_3(const WyckoffTable& WTab,
                                   const uctbx::UnitCell& uc,
                                   const SpaceGroup& SgOps,
                                   const fractional<double>& X) {
    return WTab.getWyckoffMapping(uc, SgOps, X);
  }
  WyckoffMapping
  WyckoffTable_getWyckoffMapping_4(const WyckoffTable& WTab,
                                   const uctbx::UnitCell& uc,
                                   const SpaceGroup& SgOps,
                                   const fractional<double>& X,
                                   double SnapRadius) {
    return WTab.getWyckoffMapping(uc, SgOps, X, SnapRadius);
  }

  const char* SiteSymmetry_PointGroupType(const SiteSymmetry& SS) {
    return SS.PointGroupType().Label();
  }

  fractional<double>
  SiteSymmetry_ApplySpecialOp(const SiteSymmetry& SS,
                              const fractional<double>& X) {
    return SS.ApplySpecialOp(X);
  }

  bool SiteSymmetry_isCompatibleUstar_1(const SiteSymmetry& SS,
                                        const af::double6& Ustar) {
    return SS.isCompatibleUstar(Ustar);
  }
  bool SiteSymmetry_isCompatibleUstar_2(const SiteSymmetry& SS,
                                        const af::double6& Ustar,
                                        double tolerance) {
    return SS.isCompatibleUstar(Ustar, tolerance);
  }
  void SiteSymmetry_CheckUstar_1(const SiteSymmetry& SS,
                                 const af::double6& Ustar) {
    SS.CheckUstar(Ustar);
  }
  void SiteSymmetry_CheckUstar_2(const SiteSymmetry& SS,
                                 const af::double6& Ustar,
                                 double tolerance) {
    SS.CheckUstar(Ustar, tolerance);
  }
  af::double6
  SiteSymmetry_AverageUstar(const SiteSymmetry& SS,
                            const af::double6& Ustar) {
    return SS.AverageUstar(Ustar);
  }

  RTMx SiteSymmetry_getitem(const SiteSymmetry& SS,
                            std::size_t key) {
    try {
      return SS(key);
    }
    catch (const error_index& e) {
      PyErr_SetString(PyExc_IndexError, e.what());
      throw boost::python::error_already_set();
    }
  }

  fractional<double>
  SymEquivCoordinates_getitem(const SymEquivCoordinates<double>& SymEqCoor,
                              std::size_t key) {
    try {
      return SymEqCoor(key);
    }
    catch (const error_index& e) {
      PyErr_SetString(PyExc_IndexError, e.what());
      throw boost::python::error_already_set();
    }
  }

  const char*
  ReferenceReciprocalSpaceASU_LaueGroupCode(
                                   const ReferenceReciprocalSpaceASU& StdASU) {
    return StdASU.LaueGroupCode().Label();
  }

  Miller::Index
  MillerIndexGenerator_getitem(MillerIndexGenerator& MIG, std::size_t) {
    Miller::Index result = MIG.next();
    if (result.is000()) {
      PyErr_SetString(PyExc_IndexError, "End of list.");
      throw python::error_already_set();
    }
    return result;
  }

  af::int3 StructureSeminvariant_refine_gridding_0(
    const StructureSeminvariant& ssVM) {
    return ssVM.refine_gridding();
  }
  af::int3 StructureSeminvariant_refine_gridding_1(
    const StructureSeminvariant& ssVM,
    const af::int3& grid) {
    return ssVM.refine_gridding(grid);
  }

  void
  py_expand_to_p1_4(
    const SpaceGroup& SgOps,
    bool friedel_flag,
    const af::shared<Miller::Index>& in,
    af::shared<Miller::Index>& out)
  {
    expand_to_p1(SgOps, friedel_flag, in, out);
  }

  void
  py_expand_to_p1_9(
    SpaceGroup const& SgOps,
    bool friedel_flag,
    af::shared<Miller::Index> const& h_in,
    af::shared<double> const& ampl_in,
    af::shared<double> const& phase_in,
    af::shared<Miller::Index> h_out,
    af::shared<double> ampl_out,
    af::shared<double> phase_out,
    bool phase_degrees)
  {
    expand_to_p1(
      SgOps, friedel_flag,
      h_in, ampl_in, phase_in,
      h_out, ampl_out, phase_out,
      phase_degrees);
  }
  void
  py_expand_to_p1_8(
    SpaceGroup const& SgOps,
    bool friedel_flag,
    af::shared<Miller::Index> const& h_in,
    af::shared<double> const& ampl_in,
    af::shared<double> const& phase_in,
    af::shared<Miller::Index> h_out,
    af::shared<double> ampl_out,
    af::shared<double> phase_out)
  {
    expand_to_p1(
      SgOps, friedel_flag,
      h_in, ampl_in, phase_in,
      h_out, ampl_out, phase_out);
  }

} // namespace <anonymous>

namespace boost { namespace python {
  template class
  class_builder<ReferenceReciprocalSpaceASU>; // explicitly instantiate
}} // namespace boost::python

BOOST_PYTHON_BEGIN_CONVERSION_NAMESPACE
  PyObject* to_python(const ReferenceReciprocalSpaceASU* p)
  {
    return boost::python::python_extension_class_converters<
      ReferenceReciprocalSpaceASU>::smart_ptr_to_python(
      const_cast<ReferenceReciprocalSpaceASU*>(p));
  }
BOOST_PYTHON_END_CONVERSION_NAMESPACE


BOOST_PYTHON_MODULE_INIT(sgtbx)
{
# include <cctbx/basic/from_bpl_import.h>

  python::module_builder this_module("sgtbx");

  const std::string Revision = "$Revision$";
  this_module.add(ref(to_python(
      Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

  python::import_converters<uctbx::UnitCell>
  py_UnitCell("cctbx_boost.uctbx", "UnitCell");

  python::import_converters<af::shared<bool> >
  py_shared_bool("cctbx_boost.arraytbx.shared", "bool");

  python::import_converters<af::shared<int> >
  py_shared_int("cctbx_boost.arraytbx.shared", "int");

  python::import_converters<af::shared<double> >
  py_shared_double("cctbx_boost.arraytbx.shared", "double");

  python::import_converters<af::shared<Miller::Index> >
  py_shared_Miller_Index("cctbx_boost.arraytbx.shared", "Miller_Index");

  python::import_converters<af::shared<RTMx> >
  py_shared_RTMx("cctbx_boost.arraytbx.shared", "RTMx");

  this_module.add(ref(to_python(STBF)), "STBF");
  this_module.add(ref(to_python(CRBF)), "CRBF");
  this_module.add(ref(to_python(CTBF)), "CTBF");

  class_builder<rational>
  py_rational(this_module, "rational");
  class_builder<SpaceGroupSymbols>
  py_SpaceGroupSymbols(this_module, "SpaceGroupSymbols");
  class_builder<SpaceGroupSymbolIterator>
  py_SpaceGroupSymbolIterator(this_module, "SpaceGroupSymbolIterator");
  class_builder<parse_string>
  py_parse_string(this_module, "parse_string");
  class_builder<RotMxInfo>
  py_RotMxInfo(this_module, "RotMxInfo");
  class_builder<TranslationComponents>
  py_TranslationComponents(this_module, "TranslationComponents");
  class_builder<RTMx>
  py_RTMx(this_module, "RTMx");
  python::export_converters(py_RTMx);
  class_builder<ChOfBasisOp>
  py_ChOfBasisOp(this_module, "ChOfBasisOp");
  class_builder<Miller::SymEquivIndex>
  py_Miller_SymEquivIndex(this_module, "Miller_SymEquivIndex");
  class_builder<PhaseInfo>
  py_PhaseInfo(this_module, "PhaseInfo");
  class_builder<SymEquivMillerIndices>
  py_SymEquivMillerIndices(this_module, "SymEquivMillerIndices");
  class_builder<SpaceGroup>
  py_SpaceGroup(this_module, "SpaceGroup");
  python::export_converters(py_SpaceGroup);
  class_builder<SpaceGroupInfo>
  py_SpaceGroupInfo(this_module, "SpaceGroupInfo");
  python::export_converters(py_SpaceGroupInfo);
  class_builder<WyckoffPosition>
  py_WyckoffPosition(this_module, "WyckoffPosition");
  class_builder<WyckoffMapping>
  py_WyckoffMapping(this_module, "WyckoffMapping");
  class_builder<WyckoffTable>
  py_WyckoffTable(this_module, "WyckoffTable");
  class_builder<SpecialPositionSnapParameters>
  py_SpecialPositionSnapParameters(this_module,
    "SpecialPositionSnapParameters");
  class_builder<SpecialPositionTolerances>
  py_SpecialPositionTolerances(this_module, "SpecialPositionTolerances");
  class_builder<SiteSymmetry>
  py_SiteSymmetry(this_module, "SiteSymmetry");
  python::export_converters(py_SiteSymmetry);
  class_builder<SymEquivCoordinates<double> >
  py_SymEquivCoordinates(this_module, "SymEquivCoordinates");
  class_builder<BrickPoint>
  py_BrickPoint(this_module, "BrickPoint");
  class_builder<Brick>
  py_Brick(this_module, "Brick");
  class_builder<ReferenceReciprocalSpaceASU>
  py_ReferenceReciprocalSpaceASU(this_module, "ReferenceReciprocalSpaceASU");
  class_builder<ReciprocalSpaceASU>
  py_ReciprocalSpaceASU(this_module, "ReciprocalSpaceASU");
  class_builder<MillerIndexGenerator>
  py_MillerIndexGenerator(this_module, "MillerIndexGenerator");
  class_builder<Miller::IndexTableLayoutAdaptor>
  py_Miller_IndexTableLayoutAdaptor(this_module,
    "Miller_IndexTableLayoutAdaptor");
  class_builder<Miller::AsymIndex>
  py_Miller_AsymIndex(this_module, "Miller_AsymIndex");
  class_builder<StructureSeminvariant>
  py_StructureSeminvariant(this_module, "StructureSeminvariant");

  py_Miller_IndexTableLayoutAdaptor.declare_base(
    py_Miller_SymEquivIndex, python::without_downcast);
  py_Miller_AsymIndex.declare_base(
    py_Miller_SymEquivIndex, python::without_downcast);

  py_rational.def(constructor<>());
  py_rational.def(constructor<int>());
  py_rational.def(constructor<int, int>());
  py_rational.def(rational_numerator, "numerator");
  py_rational.def(rational_denominator, "denominator");
  py_rational.def(rational_format_0, "format");
  py_rational.def(rational_format_0, "__repr__");
  py_rational.def(rational_format_1, "format");

  py_SpaceGroupSymbols.def(constructor<>());
  py_SpaceGroupSymbols.def(constructor<const std::string&>());
  py_SpaceGroupSymbols.def(constructor<const std::string&,
                                       const std::string&>());
  py_SpaceGroupSymbols.def(constructor<int>());
  py_SpaceGroupSymbols.def(constructor<int,
                                       const std::string&>());
  py_SpaceGroupSymbols.def(constructor<int,
                                       const std::string&,
                                       const std::string&>());
  py_SpaceGroupSymbols.def(&SpaceGroupSymbols::SgNumber, "SgNumber");
  py_SpaceGroupSymbols.def(&SpaceGroupSymbols::Schoenflies, "Schoenflies");
  py_SpaceGroupSymbols.def(&SpaceGroupSymbols::Qualifier, "Qualifier");
  py_SpaceGroupSymbols.def(&SpaceGroupSymbols::Hermann_Mauguin,
                                              "Hermann_Mauguin");
  py_SpaceGroupSymbols.def(&SpaceGroupSymbols::Extension, "Extension");
  py_SpaceGroupSymbols.def(&SpaceGroupSymbols::ExtendedHermann_Mauguin,
                                              "ExtendedHermann_Mauguin");
  py_SpaceGroupSymbols.def(&SpaceGroupSymbols::Hall, "Hall");

  py_SpaceGroupSymbolIterator.def(constructor<>());
  py_SpaceGroupSymbolIterator.def(&SpaceGroupSymbolIterator::next, "next");
  py_SpaceGroupSymbolIterator.def(SpaceGroupSymbolIterator_getitem,
    "__getitem__");

  py_parse_string.def(constructor<>());
  py_parse_string.def(constructor<const std::string&>());
  py_parse_string.def(&parse_string::string, "string");
  py_parse_string.def(&parse_string::where, "where");

  py_RotMxInfo.def(constructor<>());
  py_RotMxInfo.def(&RotMxInfo::Rtype, "Rtype");
  py_RotMxInfo.def(&RotMxInfo::EV, "EV");
  py_RotMxInfo.def(&RotMxInfo::SenseOfRotation, "SenseOfRotation");

  py_TranslationComponents.def(constructor<>());
  py_TranslationComponents.def(
     TranslationComponents_IntrinsicPart, "IntrinsicPart");
  py_TranslationComponents.def(
     TranslationComponents_LocationPart, "LocationPart");
  py_TranslationComponents.def(
     TranslationComponents_OriginShift, "OriginShift");

  this_module.def(RTMx_from_tuple, "RTMx_from_tuple");
  py_RTMx.def(constructor<>());
  py_RTMx.def(constructor<parse_string&>());
  py_RTMx.def(constructor<parse_string&, const char*>());
  py_RTMx.def(constructor<parse_string&, const char*, int>());
  py_RTMx.def(constructor<parse_string&, const char*, int, int>());
  py_RTMx.def(constructor<const std::string&>());
  py_RTMx.def(constructor<const std::string&, const char*>());
  py_RTMx.def(constructor<const std::string&, const char*, int>());
  py_RTMx.def(constructor<const std::string&, const char*, int, int>());
  py_RTMx.def(&RTMx::isValid, "isValid");
  py_RTMx.def(&RTMx::RBF, "RBF");
  py_RTMx.def(&RTMx::TBF, "TBF");
  py_RTMx.def(&RTMx::Unit, "Unit");
  py_RTMx.def(&RTMx::isUnit, "isUnit");
  py_RTMx.def(RTMx_newBaseFactors, "newBaseFactors");
  py_RTMx.def(&RTMx::cancel, "cancel");
  py_RTMx.def(&RTMx::modPositiveInPlace, "modPositiveInPlace");
  py_RTMx.def(&RTMx::modPositive, "modPositive");
  py_RTMx.def(&RTMx::modShortInPlace, "modShortInPlace");
  py_RTMx.def(&RTMx::modShort, "modShort");
  py_RTMx.def(&RTMx::inverse, "inverse");
  py_RTMx.def(&RTMx::inverse_with_cancel, "inverse_with_cancel");
  py_RTMx.def(RTMx_as_xyz_0, "__repr__");
  py_RTMx.def(RTMx_as_xyz_0, "as_xyz");
  py_RTMx.def(&RTMx::as_xyz, "as_xyz");
  py_RTMx.def(RTMx_as_tuple_0, "as_tuple");
  py_RTMx.def(RTMx_as_tuple_2, "as_tuple");
  py_RTMx.def(&RTMx::getRotMxInfo, "getRotMxInfo");
  py_RTMx.def(&RTMx::isPerpendicular, "isPerpendicular");
  py_RTMx.def(&RTMx::analyzeTpart, "analyzeTpart");
  py_RTMx.def(RTMx_mul_RTMx, "__mul__");
  py_RTMx.def(RTMx_multiply_RTMx, "multiply");
  py_RTMx.def(RTMx_multiply_int3, "multiply");

  py_ChOfBasisOp.def(constructor<const RTMx&, const RTMx&>());
  py_ChOfBasisOp.def(constructor<const RTMx&>());
  py_ChOfBasisOp.def(constructor<>());
  py_ChOfBasisOp.def(constructor<int>());
  py_ChOfBasisOp.def(constructor<int, int>());
  py_ChOfBasisOp.def(constructor<const std::string&>());
  py_ChOfBasisOp.def(constructor<const std::string&, const char*>());
  py_ChOfBasisOp.def(constructor<const std::string&, const char*,int>());
  py_ChOfBasisOp.def(constructor<const std::string&, const char*,int,int>());
  py_ChOfBasisOp.def(&ChOfBasisOp::isValid, "isValid");
  py_ChOfBasisOp.def(&ChOfBasisOp::M, "M");
  py_ChOfBasisOp.def(&ChOfBasisOp::InvM, "InvM");
  py_ChOfBasisOp.def(&ChOfBasisOp::swap, "swap");
  py_ChOfBasisOp.def(ChOfBasisOp_call, "__call__");
  py_ChOfBasisOp.def(ChOfBasisOp_apply_RTMx, "apply");
  py_ChOfBasisOp.def(ChOfBasisOp_apply_UnitCell, "apply");

  py_Miller_SymEquivIndex.def(constructor<>());
  py_Miller_SymEquivIndex.def(
    constructor<const Miller::Index&, int, int, bool>());
  py_Miller_SymEquivIndex.def(&Miller::SymEquivIndex::H, "H");
  py_Miller_SymEquivIndex.def(&Miller::SymEquivIndex::HR, "HR");
  py_Miller_SymEquivIndex.def(&Miller::SymEquivIndex::HT, "HT");
  py_Miller_SymEquivIndex.def(&Miller::SymEquivIndex::TBF, "TBF");
  py_Miller_SymEquivIndex.def(&Miller::SymEquivIndex::HT_angle, "HT_angle");
  py_Miller_SymEquivIndex.def(Miller_SymEquivIndex_HT_angle_0, "HT_angle");
  py_Miller_SymEquivIndex.def(
    &Miller::SymEquivIndex::FriedelFlag, "FriedelFlag");
  py_Miller_SymEquivIndex.def(Miller_SymEquivIndex_Mate_0, "Mate");
  py_Miller_SymEquivIndex.def(&Miller::SymEquivIndex::Mate, "Mate");
  py_Miller_SymEquivIndex.def(
    Miller_SymEquivIndex_phase_eq_1, "phase_eq");
  py_Miller_SymEquivIndex.def(
    Miller_SymEquivIndex_phase_eq_2, "phase_eq");
  py_Miller_SymEquivIndex.def(
    Miller_SymEquivIndex_phase_in_1, "phase_in");
  py_Miller_SymEquivIndex.def(
    Miller_SymEquivIndex_phase_in_2, "phase_in");
  py_Miller_SymEquivIndex.def(
    Miller_SymEquivIndex_complex_eq, "complex_eq");
  py_Miller_SymEquivIndex.def(
    Miller_SymEquivIndex_complex_in, "complex_in");
  py_Miller_SymEquivIndex.def(
    Miller_SymEquivIndex_hl_eq, "hl_eq");
  py_Miller_SymEquivIndex.def(
    Miller_SymEquivIndex_hl_in, "hl_in");

  py_PhaseInfo.def(constructor<>());
  py_PhaseInfo.def(
    constructor<SpaceGroup const&, Miller::Index const&>());
  py_PhaseInfo.def(
    constructor<SpaceGroup const&, Miller::Index const&, bool>());
  py_PhaseInfo.def(&PhaseInfo::SysAbsWasTested, "SysAbsWasTested");
  py_PhaseInfo.def(&PhaseInfo::isSysAbsent, "isSysAbsent");
  py_PhaseInfo.def(&PhaseInfo::isCentric, "isCentric");
  py_PhaseInfo.def(&PhaseInfo::HT, "HT");
  py_PhaseInfo.def(&PhaseInfo::TBF, "TBF");
  py_PhaseInfo.def(&PhaseInfo::HT_angle, "HT_angle");
  py_PhaseInfo.def(PhaseInfo_HT_angle_0, "HT_angle");
  py_PhaseInfo.def(&PhaseInfo::isValidPhase, "isValidPhase");
  py_PhaseInfo.def(PhaseInfo_isValidPhase_2, "isValidPhase");
  py_PhaseInfo.def(PhaseInfo_isValidPhase_1, "isValidPhase");
  py_PhaseInfo.def(&PhaseInfo::nearest_valid_phase, "nearest_valid_phase");
  py_PhaseInfo.def(PhaseInfo_nearest_valid_phase_1, "nearest_valid_phase");

  py_SymEquivMillerIndices.def(constructor<>());
  py_SymEquivMillerIndices.def(
    &SymEquivMillerIndices::getPhaseRestriction, "getPhaseRestriction");
  py_SymEquivMillerIndices.def(
    &SymEquivMillerIndices::isCentric, "isCentric");
  py_SymEquivMillerIndices.def(&SymEquivMillerIndices::N, "N");
  py_SymEquivMillerIndices.def(&SymEquivMillerIndices::M, "M");
  py_SymEquivMillerIndices.def(&SymEquivMillerIndices::fMates, "fMates");
  py_SymEquivMillerIndices.def(&SymEquivMillerIndices::epsilon, "epsilon");
  py_SymEquivMillerIndices.def(&SymEquivMillerIndices::N, "__len__");
  py_SymEquivMillerIndices.def(SymEquivMillerIndices_getitem, "__getitem__");
  py_SymEquivMillerIndices.def(SymEquivMillerIndices_call_1, "__call__");
  py_SymEquivMillerIndices.def(SymEquivMillerIndices_call_2, "__call__");
  py_SymEquivMillerIndices.def(
    &SymEquivMillerIndices::isValidPhase, "isValidPhase");
  py_SymEquivMillerIndices.def(
    SymEquivMillerIndices_isValidPhase_2, "isValidPhase");
  py_SymEquivMillerIndices.def(
    SymEquivMillerIndices_isValidPhase_1, "isValidPhase");

  py_SpaceGroup.def(constructor<>());
  py_SpaceGroup.def(constructor<parse_string&>());
  py_SpaceGroup.def(constructor<parse_string&, bool>());
  py_SpaceGroup.def(constructor<parse_string&, bool, bool>());
  py_SpaceGroup.def(constructor<const std::string&>());
  py_SpaceGroup.def(constructor<const std::string&, bool>());
  py_SpaceGroup.def(constructor<const std::string&, bool, bool>());
  py_SpaceGroup.def(constructor<const SpaceGroupSymbols&>());
  py_SpaceGroup.def(SpaceGroup_ParseHallSymbol_parse_string,
                              "ParseHallSymbol");
  py_SpaceGroup.def(SpaceGroup_ParseHallSymbol_std_string,
                              "ParseHallSymbol");
  py_SpaceGroup.def(&SpaceGroup::ParseHallSymbol, "ParseHallSymbol");
  py_SpaceGroup.def(&SpaceGroup::ChangeBasis, "ChangeBasis");
  py_SpaceGroup.def(&SpaceGroup::expandSMx, "expandSMx");
  py_SpaceGroup.def(SpaceGroup_expandLTr, "expandLTr");
  py_SpaceGroup.def(&SpaceGroup::expandConventionalCentringType,
                                "expandConventionalCentringType");
  py_SpaceGroup.def(&SpaceGroup::RBF, "RBF");
  py_SpaceGroup.def(&SpaceGroup::TBF, "TBF");
  py_SpaceGroup.def(&SpaceGroup::nLTr, "nLTr");
  py_SpaceGroup.def(&SpaceGroup::fInv, "fInv");
  py_SpaceGroup.def(&SpaceGroup::nSMx, "nSMx");
  py_SpaceGroup.def(&SpaceGroup::OrderP, "OrderP");
  py_SpaceGroup.def(&SpaceGroup::OrderZ, "OrderZ");
  py_SpaceGroup.def(SpaceGroup_call_1, "__call__");
  py_SpaceGroup.def(SpaceGroup_call_3, "__call__");
  py_SpaceGroup.def(&SpaceGroup::OrderZ, "__len__");
  py_SpaceGroup.def(SpaceGroup_getitem, "__getitem__");
  py_SpaceGroup.def(&SpaceGroup::getConventionalCentringTypeSymbol,
                                "getConventionalCentringTypeSymbol");
  py_SpaceGroup.def(&SpaceGroup::isChiral, "isChiral");
  py_SpaceGroup.def(SpaceGroup_isSysAbsent_a, "isSysAbsent");
  py_SpaceGroup.def(SpaceGroup_isSysAbsent_s, "isSysAbsent");
  py_SpaceGroup.def(SpaceGroup_isCentric_0, "isCentric");
  py_SpaceGroup.def(SpaceGroup_isCentric_a, "isCentric");
  py_SpaceGroup.def(SpaceGroup_isCentric_s, "isCentric");
  py_SpaceGroup.def(&SpaceGroup::isOriginCentric, "isOriginCentric");
  py_SpaceGroup.def(&SpaceGroup::getPhaseRestriction, "getPhaseRestriction");
  py_SpaceGroup.def(&SpaceGroup::isValidPhase, "isValidPhase");
  py_SpaceGroup.def(SpaceGroup_isValidPhase_3, "isValidPhase");
  py_SpaceGroup.def(SpaceGroup_isValidPhase_2, "isValidPhase");
  py_SpaceGroup.def(SpaceGroup_epsilon_a, "epsilon");
  py_SpaceGroup.def(SpaceGroup_epsilon_s, "epsilon");
  py_SpaceGroup.def(SpaceGroup_multiplicity_a, "multiplicity");
  py_SpaceGroup.def(SpaceGroup_multiplicity_s, "multiplicity");
  py_SpaceGroup.def(&SpaceGroup::getEquivMillerIndices,
                                "getEquivMillerIndices");
  py_SpaceGroup.def(SpaceGroup_CheckMetricalMatrix_1, "CheckMetricalMatrix");
  py_SpaceGroup.def(&SpaceGroup::CheckMetricalMatrix, "CheckMetricalMatrix");
  py_SpaceGroup.def(SpaceGroup_CheckUnitCell_1, "CheckUnitCell");
  py_SpaceGroup.def(&SpaceGroup::CheckUnitCell, "CheckUnitCell");
  py_SpaceGroup.def(SpaceGroup_getZ2POp_0, "getZ2POp");
  py_SpaceGroup.def(&SpaceGroup::getZ2POp, "getZ2POp");
  py_SpaceGroup.def(&SpaceGroup::makeTidy, "makeTidy");
  py_SpaceGroup.def(SpaceGroup_cmp_equal, "__cmp__");
  py_SpaceGroup.def(&SpaceGroup::MatchTabulatedSettings,
                                "MatchTabulatedSettings");
  py_SpaceGroup.def(&SpaceGroup::Info, "Info");
  py_SpaceGroup.def(SpaceGroup_refine_gridding_0, "refine_gridding");
  py_SpaceGroup.def(SpaceGroup_refine_gridding_1, "refine_gridding");
  py_SpaceGroup.def(SpaceGroup_repr, "__repr__");
  py_SpaceGroup.def(SpaceGroup_getinitargs, "__getinitargs__");

  py_SpaceGroupInfo.def(constructor<>());
  py_SpaceGroupInfo.def(constructor<const SpaceGroup&>());
  py_SpaceGroupInfo.def(constructor<const SpaceGroup&, bool>());
  py_SpaceGroupInfo.def(constructor<const SpaceGroup&, bool, int>());
  py_SpaceGroupInfo.def(constructor<const SpaceGroup&, bool, int, int>());
  py_SpaceGroupInfo.def(&SpaceGroupInfo::SgOps, "SgOps");
  py_SpaceGroupInfo.def(&SpaceGroupInfo::SgNumber, "SgNumber");
  py_SpaceGroupInfo.def(&SpaceGroupInfo::CBOp, "CBOp");
  py_SpaceGroupInfo.def(
    &SpaceGroupInfo::getAddlGeneratorsOfEuclideanNormalizer,
                    "getAddlGeneratorsOfEuclideanNormalizer");
  py_SpaceGroupInfo.def(
    &SpaceGroupInfo::expandAddlGeneratorsOfEuclideanNormalizer,
                    "expandAddlGeneratorsOfEuclideanNormalizer");
  py_SpaceGroupInfo.def(&SpaceGroupInfo::isEnantiomorphic,
                                        "isEnantiomorphic");
  py_SpaceGroupInfo.def(&SpaceGroupInfo::getChangeOfHandOp,
                                        "getChangeOfHandOp");
  py_SpaceGroupInfo.def(SpaceGroupInfo_BuildHallSymbol_0, "BuildHallSymbol");
  py_SpaceGroupInfo.def(SpaceGroupInfo_BuildHallSymbol_1, "BuildHallSymbol");
  py_SpaceGroupInfo.def(&SpaceGroupInfo::BuildLookupSymbol,
                                        "BuildLookupSymbol");

  py_WyckoffPosition.def(constructor<>());
  py_WyckoffPosition.def(&WyckoffPosition::M, "M");
  py_WyckoffPosition.def(&WyckoffPosition::Letter, "Letter");
  py_WyckoffPosition.def(&WyckoffPosition::SpecialOp, "SpecialOp");
  py_WyckoffPosition.def(&WyckoffPosition::isExpanded, "isExpanded");
  py_WyckoffPosition.def(&WyckoffPosition::CheckExpanded,"CheckExpanded");
  py_WyckoffPosition.def(&WyckoffPosition::operator(), "__call__");
  py_WyckoffPosition.def(&WyckoffPosition::M, "__len__");
  py_WyckoffPosition.def(WyckoffPosition_getitem, "__getitem__");

  py_WyckoffMapping.def(constructor<>());
  py_WyckoffMapping.def(&WyckoffMapping::WP, "WP");
  py_WyckoffMapping.def(&WyckoffMapping::Mapping, "Mapping");
  py_WyckoffMapping.def(&WyckoffMapping::snap_to_representative,
                                        "snap_to_representative");
  py_WyckoffMapping.def(&WyckoffMapping::snap, "snap");

  py_WyckoffTable.def(constructor<>());
  py_WyckoffTable.def(constructor<const SpaceGroupInfo&>());
  py_WyckoffTable.def(constructor<const SpaceGroupInfo&, bool>());
  py_WyckoffTable.def(&WyckoffTable::expand, "expand");
  py_WyckoffTable.def(&WyckoffTable::N, "N");
  py_WyckoffTable.def(WyckoffTable_call_size_t, "__call__");
  py_WyckoffTable.def(WyckoffTable_call_char, "__call__");
  py_WyckoffTable.def(&WyckoffTable::N, "__len__");
  py_WyckoffTable.def(WyckoffTable_getitem, "__getitem__");
  py_WyckoffTable.def(&WyckoffTable::LookupIndex, "LookupIndex");
  py_WyckoffTable.def(WyckoffTable_getWyckoffMapping_SS,"getWyckoffMapping");
  py_WyckoffTable.def(WyckoffTable_getWyckoffMapping_3, "getWyckoffMapping");
  py_WyckoffTable.def(WyckoffTable_getWyckoffMapping_4, "getWyckoffMapping");

  py_SpecialPositionSnapParameters.def(constructor<const uctbx::UnitCell&,
                                                   const SpaceGroup&>());
  py_SpecialPositionSnapParameters.def(constructor<const uctbx::UnitCell&,
                                                   const SpaceGroup&,
                                                   bool>());
  py_SpecialPositionSnapParameters.def(constructor<const uctbx::UnitCell&,
                                                   const SpaceGroup&,
                                                   bool,
                                                   double>());

  py_SpecialPositionTolerances.def(constructor<const uctbx::UnitCell&,
                                               const SpaceGroup&>());
  py_SpecialPositionTolerances.def(constructor<const uctbx::UnitCell&,
                                               const SpaceGroup&,
                                               double>());
  py_SpecialPositionTolerances.def(constructor<const uctbx::UnitCell&,
                                               const SpaceGroup&,
                                               double,
                                               double>());

  py_SiteSymmetry.def(constructor<>());
  py_SiteSymmetry.def(
    constructor<const SpecialPositionSnapParameters&,
                const fractional<double>&>());
  py_SiteSymmetry.def(
    constructor<const SpecialPositionSnapParameters&,
                const fractional<double>&,
                bool>());
  py_SiteSymmetry.def(&SiteSymmetry::OriginalPosition, "OriginalPosition");
  py_SiteSymmetry.def(&SiteSymmetry::SnapPosition, "SnapPosition");
  py_SiteSymmetry.def(&SiteSymmetry::DistanceMoved2, "DistanceMoved2");
  py_SiteSymmetry.def(&SiteSymmetry::DistanceMoved, "DistanceMoved");
  py_SiteSymmetry.def(&SiteSymmetry::ShortestDistance2, "ShortestDistance2");
  py_SiteSymmetry.def(&SiteSymmetry::ShortestDistance, "ShortestDistance");
  py_SiteSymmetry.def(&SiteSymmetry::isWellBehaved, "isWellBehaved");
  py_SiteSymmetry.def(&SiteSymmetry::M, "M");
  py_SiteSymmetry.def(&SiteSymmetry::SpecialOp, "SpecialOp");
  py_SiteSymmetry.def(SiteSymmetry_PointGroupType, "PointGroupType");
  py_SiteSymmetry.def(SiteSymmetry_ApplySpecialOp, "ApplySpecialOp");
  py_SiteSymmetry.def(SiteSymmetry_isCompatibleUstar_1, "isCompatibleUstar");
  py_SiteSymmetry.def(SiteSymmetry_isCompatibleUstar_2, "isCompatibleUstar");
  py_SiteSymmetry.def(SiteSymmetry_CheckUstar_1, "CheckUstar");
  py_SiteSymmetry.def(SiteSymmetry_CheckUstar_2, "CheckUstar");
  py_SiteSymmetry.def(SiteSymmetry_AverageUstar, "AverageUstar");
  py_SiteSymmetry.def(&SiteSymmetry::expand, "expand");
  py_SiteSymmetry.def(&SiteSymmetry::isExpanded, "isExpanded");
  py_SiteSymmetry.def(&SiteSymmetry::CheckExpanded, "CheckExpanded");
  py_SiteSymmetry.def(&SiteSymmetry::operator(), "__call__");
  py_SiteSymmetry.def(&SiteSymmetry::M, "__len__");
  py_SiteSymmetry.def(SiteSymmetry_getitem, "__getitem__");

  py_SymEquivCoordinates.def(constructor<>());
  py_SymEquivCoordinates.def(
    constructor<const SpecialPositionSnapParameters&,
                const fractional<double>&>());
  py_SymEquivCoordinates.def(constructor<const SiteSymmetry&>());
  py_SymEquivCoordinates.def(
    constructor<const WyckoffMapping&,
    const fractional<double>&>());
  py_SymEquivCoordinates.def(
    constructor<const WyckoffPosition&,
    const fractional<double>&>());
  py_SymEquivCoordinates.def(
    constructor<const SpecialPositionTolerances&,
    const fractional<double>&>());
  py_SymEquivCoordinates.def(
    constructor<const SpaceGroup&,
    const fractional<double>&>());
  py_SymEquivCoordinates.def(&SymEquivCoordinates<double>::M, "M");
  py_SymEquivCoordinates.def(
    &SymEquivCoordinates<double>::operator(), "__call__");
  py_SymEquivCoordinates.def(&SymEquivCoordinates<double>::M, "__len__");
  py_SymEquivCoordinates.def(SymEquivCoordinates_getitem, "__getitem__");
  py_SymEquivCoordinates.def(
    &SymEquivCoordinates<double>::getShortestDifference,
                                 "getShortestDifference");
  py_SymEquivCoordinates.def(&SymEquivCoordinates<double>
     ::getShortestDifferenceUnderAllowedOriginShifts,
      "getShortestDifferenceUnderAllowedOriginShifts");
  py_SymEquivCoordinates.def(
    &SymEquivCoordinates<double>::getShortestDistance2,
                                 "getShortestDistance2");
  py_SymEquivCoordinates.def(
    &SymEquivCoordinates<double>::getShortestDistance,
                                 "getShortestDistance");
  py_SymEquivCoordinates.def(
    &SymEquivCoordinates<double>::StructureFactor,
                                 "StructureFactor");

  py_BrickPoint.def(constructor<>());
  py_BrickPoint.def(&BrickPoint::Point, "Point");
  py_BrickPoint.def(&BrickPoint::Off, "Off");

  py_Brick.def(constructor<>());
  py_Brick.def(constructor<const SpaceGroupInfo&>());
  py_Brick.def(&Brick::operator(), "__call__");
  py_Brick.def(&Brick::format, "format");
  py_Brick.def(&Brick::format, "__repr__");

  py_ReferenceReciprocalSpaceASU.def(constructor<>());
  py_ReferenceReciprocalSpaceASU.def(
     ReferenceReciprocalSpaceASU_LaueGroupCode, "LaueGroupCode");
  py_ReferenceReciprocalSpaceASU.def(
    &ReferenceReciprocalSpaceASU::isInASU, "isInASU");
  py_ReferenceReciprocalSpaceASU.def(
    &ReferenceReciprocalSpaceASU::representation, "representation");
  py_ReferenceReciprocalSpaceASU.def(
    &ReferenceReciprocalSpaceASU::getCutParameters, "getCutParameters");

  py_ReciprocalSpaceASU.def(constructor<>());
  py_ReciprocalSpaceASU.def(constructor<const SpaceGroupInfo&>());
  py_ReciprocalSpaceASU.def(
    &ReciprocalSpaceASU::ReferenceASU, "ReferenceASU");
  py_ReciprocalSpaceASU.def(&ReciprocalSpaceASU::CBOp, "CBOp");
  py_ReciprocalSpaceASU.def(
    &ReciprocalSpaceASU::isReferenceASU, "isReferenceASU");
  py_ReciprocalSpaceASU.def(&ReciprocalSpaceASU::isInASU, "isInASU");

  py_MillerIndexGenerator.def(constructor<>());
  py_MillerIndexGenerator.def(constructor<const uctbx::UnitCell&,
                                          const SpaceGroupInfo&,
                                          bool,
                                          double>());
  py_MillerIndexGenerator.def(constructor<const SpaceGroupInfo&,
                                          bool,
                                          const Miller::Index&>());
  py_MillerIndexGenerator.def(&MillerIndexGenerator::next, "next");
  py_MillerIndexGenerator.def(MillerIndexGenerator_getitem, "__getitem__");
  py_MillerIndexGenerator.def(&MillerIndexGenerator::ASU, "ASU");

  py_Miller_IndexTableLayoutAdaptor.def(constructor<>());
  py_Miller_IndexTableLayoutAdaptor.def(
    &Miller::IndexTableLayoutAdaptor::H, "H");
  py_Miller_IndexTableLayoutAdaptor.def(
    &Miller::IndexTableLayoutAdaptor::iColumn, "iColumn");

  py_Miller_AsymIndex.def(constructor<>());
  py_Miller_AsymIndex.def(constructor<
    const SpaceGroup&,
    const ReciprocalSpaceASU&,
    const Miller::Index&>());
  py_Miller_AsymIndex.def(constructor<
    const SpaceGroup&,
    const Miller::Index&>());
  py_Miller_AsymIndex.def(constructor<
    const SymEquivMillerIndices&>());
  py_Miller_AsymIndex.def(
    &Miller::AsymIndex::one_column, "one_column");
  py_Miller_AsymIndex.def(
    &Miller::AsymIndex::two_column, "two_column");

  py_StructureSeminvariant.def(constructor<>());
  py_StructureSeminvariant.def(constructor<const SpaceGroup&>());
  py_StructureSeminvariant.def(&StructureSeminvariant::size, "size");
  py_StructureSeminvariant.def(&StructureSeminvariant::size, "__len__");
  py_StructureSeminvariant.def(&StructureSeminvariant::V, "V");
  py_StructureSeminvariant.def(&StructureSeminvariant::M, "M");
  py_StructureSeminvariant.def(&StructureSeminvariant::is_ss, "is_ss");
  //XXX: expose
  //py_StructureSeminvariant.def(&StructureSeminvariant::apply_mod,
  //                                                    "apply_mod");
  py_StructureSeminvariant.def(
    StructureSeminvariant_refine_gridding_0, "refine_gridding");
  py_StructureSeminvariant.def(
    StructureSeminvariant_refine_gridding_1, "refine_gridding");

  this_module.def(py_expand_to_p1_4, "expand_to_p1");
  this_module.def(py_expand_to_p1_9, "expand_to_p1");
  this_module.def(py_expand_to_p1_8, "expand_to_p1");

  sgtbx::sanity_check();
}
