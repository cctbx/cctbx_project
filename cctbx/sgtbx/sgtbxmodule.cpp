// $Id$

#include <cctbx/sgtbx/symbols.h>
#include <cctbx/sgtbx/groups.h>
#include <boost/python/cross_module.hpp>
#include <cctbx/basic/boost_array_bpl.h>
#include <cctbx/miller_bpl.h>

using namespace sgtbx;

namespace {

# include <cctbx/basic/from_bpl_import.h>

  std::string RTMx_as_xyz_0(const RTMx& M) {
    return M.as_xyz();
  }

  tuple RTMx_as_tuple(const RTMx& M) {
    tuple result(2);
    result.set_item(0,
      ref(to_python(M.Rpart().as_array(static_cast<double>(0)))));
    result.set_item(1,
      ref(to_python(M.Tpart().as_array(static_cast<double>(0)))));
    return result;
  }

  int SgOps_ParseHallSymbol_parse_string(SgOps& sgo, parse_string& HSym) {
    return sgo.ParseHallSymbol(HSym);
  }
  int SgOps_ParseHallSymbol_std_string(SgOps& sgo, const std::string& HSym) {
    parse_string HSymPS(HSym);
    return sgo.ParseHallSymbol(HSymPS);
  }

  RTMx SgOps_call_1(const SgOps& sgo, int iLIS) {
    return sgo(iLIS);
  }
  RTMx SgOps_call_3(const SgOps& sgo, int iLTr, int iInv, int iSMx) {
    return sgo(iLTr, iInv, iSMx);
  }

  bool SgOps_isCentric_0(const SgOps& sgo) {
    return sgo.isCentric();
  }
  bool SgOps_isCentric_1(const SgOps& sgo, const Miller::Index& H) {
    return sgo.isCentric(H);
  }

  bool SgOps_isValidPhase_rad_2(const SgOps& sgo,
                                const Miller::Index& H, double phi) {
    return sgo.isValidPhase_rad(H, phi);
  }
  bool SgOps_isValidPhase_rad_3(const SgOps& sgo,
                                const Miller::Index& H, double phi,
                                double tolerance) {
    return sgo.isValidPhase_rad(H, phi, tolerance);
  }
  bool SgOps_isValidPhase_deg_2(const SgOps& sgo,
                                const Miller::Index& H, double phi) {
    return sgo.isValidPhase_deg(H, phi);
  }
  bool SgOps_isValidPhase_deg_3(const SgOps& sgo,
                                const Miller::Index& H, double phi,
                                double tolerance) {
    return sgo.isValidPhase_deg(H, phi, tolerance);
  }

  void SgOps_CheckMetricalMatrix_1(const SgOps& sgo, const uctbx::Mx33& G) {
    sgo.CheckMetricalMatrix(G);
  }
  void SgOps_CheckUnitCell_1(const SgOps& sgo, const uctbx::UnitCell& uc) {
    sgo.CheckUnitCell(uc);
  }

  Miller::MasterIndex
  SgOps_getMasterIndex_2(const SgOps& sgo,
                         const Miller::Index& H,
                         bool Pretty) {
    return sgo.getMasterIndex(H, Pretty);
  }
  Miller::MasterIndex
  SgOps_getMasterIndex_3(const SgOps& sgo,
                         const Miller::Index& H,
                         const Miller::Vec3& CutP,
                         bool Pretty) {
    return sgo.getMasterIndex(H, CutP, Pretty);
  }

  int SgOps_cmp_equal(const SgOps& lhs, const SgOps& rhs) {
    if (lhs == rhs) return 0;
    return 1;
  }

  ChOfBasisOp SgOps_getZ2POp_0(const SgOps& sgo) {
    return sgo.getZ2POp();
  }

  int PhaseRestriction_HT_0(const PhaseRestriction& PR) {
    return PR.HT();
  }
  double PhaseRestriction_HT_1(const PhaseRestriction& PR, double Period) {
    return PR.HT(Period);
  }

  bool PhaseRestriction_isValidPhase_rad_1(const PhaseRestriction& PR,
                                           double phi) {
    return PR.isValidPhase_rad(phi);
  }
  bool PhaseRestriction_isValidPhase_rad_2(const PhaseRestriction& PR,
                                           double phi, double tolerance) {
    return PR.isValidPhase_rad(phi, tolerance);
  }
  bool PhaseRestriction_isValidPhase_deg_1(const PhaseRestriction& PR,
                                           double phi) {
    return PR.isValidPhase_deg(phi);
  }
  bool PhaseRestriction_isValidPhase_deg_2(const PhaseRestriction& PR,
                                           double phi, double tolerance) {
    return PR.isValidPhase_deg(phi, tolerance);
  }

  Miller::SymEquivIndex
  SymEquivMillerIndices_getitem(const SymEquivMillerIndices& SEMI,
                                int iList) {
    if (iList < 0 || iList >= SEMI.N()) {
      throw error("Index out of range.");
    }
    return SEMI[iList];
  }
  Miller::Index
  SymEquivMillerIndices_call_1(const SymEquivMillerIndices& SEMI,
                               int iIL) {
    return SEMI(iIL);
  }
  Miller::Index
  SymEquivMillerIndices_call_2(const SymEquivMillerIndices& SEMI,
                               int iInv, int iList) {
    return SEMI(iInv, iList);
  }
  bool
  SymEquivMillerIndices_isValidPhase_rad_1(const SymEquivMillerIndices& SEMI,
                                           double phi) {
    return SEMI.isValidPhase_rad(phi);
  }
  bool
  SymEquivMillerIndices_isValidPhase_rad_2(const SymEquivMillerIndices& SEMI,
                                           double phi, double tolerance) {
    return SEMI.isValidPhase_rad(phi, tolerance);
  }
  bool
  SymEquivMillerIndices_isValidPhase_deg_1(const SymEquivMillerIndices& SEMI,
                                           double phi) {
    return SEMI.isValidPhase_deg(phi);
  }
  bool
  SymEquivMillerIndices_isValidPhase_deg_2(const SymEquivMillerIndices& SEMI,
                                           double phi, double tolerance) {
    return SEMI.isValidPhase_deg(phi, tolerance);
  }
  Miller::MasterIndex
  SymEquivMillerIndices_getMasterIndex_1(const SymEquivMillerIndices& SEMI,
                                         bool Pretty) {
    return SEMI.getMasterIndex(Pretty);
  }
  Miller::MasterIndex
  SymEquivMillerIndices_getMasterIndex_2(const SymEquivMillerIndices& SEMI,
                                         const Miller::Vec3& CutP,
                                         bool Pretty) {
    return SEMI.getMasterIndex(CutP, Pretty);
  }

  double Miller_MasterIndex_Phase_rad(const Miller::MasterIndex& Master,
                                      double phi, bool FriedelSym) {
    return Master.Phase_rad(phi, FriedelSym);
  }
  double Miller_MasterIndex_Phase_deg(const Miller::MasterIndex& Master,
                                      double phi, bool FriedelSym) {
    return Master.Phase_deg(phi, FriedelSym);
  }
  std::complex<double>
  Miller_MasterIndex_ShiftPhase(const Miller::MasterIndex& Master,
                                const std::complex<double>& F,
                                bool FriedelSym) {
    return Master.ShiftPhase(F, FriedelSym);
  }
}

BOOST_PYTHON_MODULE_INIT(sgtbx)
{
  try
  {
#   include <cctbx/basic/from_bpl_import.h>

    python::module_builder this_module("sgtbx");

    const std::string Revision = "$Revision$";
    this_module.add(ref(to_python(
        Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

    this_module.add(ref(to_python(STBF)), "STBF");
    this_module.add(ref(to_python(CRBF)), "CRBF");
    this_module.add(ref(to_python(CTBF)), "CTBF");

    class_builder<SpaceGroupSymbols>
    py_SpaceGroupSymbols(this_module, "SpaceGroupSymbols");
    class_builder<parse_string>
    py_parse_string(this_module, "parse_string");
    class_builder<RotMxInfo>
    py_RotMxInfo(this_module, "RotMxInfo");
    class_builder<RTMx>
    py_RTMx(this_module, "RTMx");
    class_builder<ChOfBasisOp>
    py_ChOfBasisOp(this_module, "ChOfBasisOp");
    class_builder<Miller::SymEquivIndex>
    py_Miller_SymEquivIndex(this_module, "Miller_SymEquivIndex");
    class_builder<Miller::MasterIndex>
    py_Miller_MasterIndex(this_module, "Miller_MasterIndex");
    class_builder<PhaseRestriction>
    py_PhaseRestriction(this_module, "PhaseRestriction");
    class_builder<SymEquivMillerIndices>
    py_SymEquivMillerIndices(this_module, "SymEquivMillerIndices");
    class_builder<SgOps>
    py_SgOps(this_module, "SgOps");
    python::export_converters(py_SgOps);

    python::import_converters<uctbx::UnitCell>
    UnitCell_converters("uctbx", "UnitCell");

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
    py_SpaceGroupSymbols.def(&SpaceGroupSymbols::Hall, "Hall");

    py_parse_string.def(constructor<const std::string&>());
    py_parse_string.def(&parse_string::string, "string");
    py_parse_string.def(&parse_string::where, "where");

    py_RotMxInfo.def(&RotMxInfo::Rtype, "Rtype");
    py_RotMxInfo.def(&RotMxInfo::EV, "EV");
    py_RotMxInfo.def(&RotMxInfo::SenseOfRotation, "SenseOfRotation");

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
    py_RTMx.def(&RTMx::normalize, "normalize");
    py_RTMx.def(&RTMx::inverse, "inverse");
    py_RTMx.def(RTMx_as_xyz_0, "as_xyz");
    py_RTMx.def(&RTMx::as_xyz, "as_xyz");
    py_RTMx.def(RTMx_as_tuple, "as_tuple");
    py_RTMx.def(&RTMx::getRotMxInfo, "getRotMxInfo");
    py_RTMx.def(operators<python::op_mul>());

    py_ChOfBasisOp.def(&ChOfBasisOp::M, "M");
    py_ChOfBasisOp.def(&ChOfBasisOp::InvM, "InvM");
    py_ChOfBasisOp.def(&ChOfBasisOp::swap, "swap");

    py_Miller_SymEquivIndex.def(&Miller::SymEquivIndex::HR, "HR");
    py_Miller_SymEquivIndex.def(&Miller::SymEquivIndex::HT, "HT");

    py_Miller_MasterIndex.def(&Miller::MasterIndex::TBF, "TBF");
    py_Miller_MasterIndex.def(&Miller::MasterIndex::haveCutP, "haveCutP");
    py_Miller_MasterIndex.def(&Miller::MasterIndex::CutP, "CutP");
    py_Miller_MasterIndex.def(&Miller::MasterIndex::Pretty, "Pretty");
    py_Miller_MasterIndex.def(&Miller::MasterIndex::H, "H");
    py_Miller_MasterIndex.def(&Miller::MasterIndex::iMate, "iMate");
    py_Miller_MasterIndex.def(&Miller::MasterIndex::HR, "HR");
    py_Miller_MasterIndex.def(&Miller::MasterIndex::HT, "HT");
    py_Miller_MasterIndex.def(Miller_MasterIndex_Phase_rad, "Phase_rad");
    py_Miller_MasterIndex.def(Miller_MasterIndex_Phase_deg, "Phase_deg");
    py_Miller_MasterIndex.def(Miller_MasterIndex_ShiftPhase, "ShiftPhase");

    py_PhaseRestriction.def(&PhaseRestriction::isCentric, "isCentric");
    py_PhaseRestriction.def(PhaseRestriction_HT_0, "HT");
    py_PhaseRestriction.def(&PhaseRestriction::TBF, "TBF");
    py_PhaseRestriction.def(PhaseRestriction_HT_1, "HT");
    py_PhaseRestriction.def(&PhaseRestriction::HT_rad, "HT_rad");
    py_PhaseRestriction.def(&PhaseRestriction::HT_deg, "HT_deg");
    py_PhaseRestriction.def(
      PhaseRestriction_isValidPhase_rad_1, "isValidPhase_rad");
    py_PhaseRestriction.def(
      PhaseRestriction_isValidPhase_rad_2, "isValidPhase_rad");
    py_PhaseRestriction.def(
      PhaseRestriction_isValidPhase_deg_1, "isValidPhase_deg");
    py_PhaseRestriction.def(
      PhaseRestriction_isValidPhase_deg_2, "isValidPhase_deg");

    py_SymEquivMillerIndices.def(
      &SymEquivMillerIndices::getPhaseRestriction, "getPhaseRestriction");
    py_SymEquivMillerIndices.def(
      &SymEquivMillerIndices::isCentric, "isCentric");
    py_SymEquivMillerIndices.def(&SymEquivMillerIndices::N, "N");
    py_SymEquivMillerIndices.def(&SymEquivMillerIndices::M, "M");
    py_SymEquivMillerIndices.def(&SymEquivMillerIndices::fMates, "fMates");
    py_SymEquivMillerIndices.def(&SymEquivMillerIndices::epsilon, "epsilon");
    py_SymEquivMillerIndices.def(SymEquivMillerIndices_getitem, "__getitem__");
    py_SymEquivMillerIndices.def(SymEquivMillerIndices_call_1, "__call__");
    py_SymEquivMillerIndices.def(SymEquivMillerIndices_call_2, "__call__");
    py_SymEquivMillerIndices.def(
      SymEquivMillerIndices_isValidPhase_rad_1, "isValidPhase_rad");
    py_SymEquivMillerIndices.def(
      SymEquivMillerIndices_isValidPhase_rad_2, "isValidPhase_rad");
    py_SymEquivMillerIndices.def(
      SymEquivMillerIndices_isValidPhase_deg_1, "isValidPhase_deg");
    py_SymEquivMillerIndices.def(
      SymEquivMillerIndices_isValidPhase_deg_2, "isValidPhase_deg");
    py_SymEquivMillerIndices.def(
      SymEquivMillerIndices_getMasterIndex_1, "getMasterIndex");
    py_SymEquivMillerIndices.def(
      SymEquivMillerIndices_getMasterIndex_2, "getMasterIndex");

    py_SgOps.def(constructor<>());
    py_SgOps.def(constructor<parse_string&>());
    py_SgOps.def(constructor<parse_string&, bool>());
    py_SgOps.def(constructor<parse_string&, bool, bool>());
    py_SgOps.def(constructor<const std::string&>());
    py_SgOps.def(constructor<const std::string&, bool>());
    py_SgOps.def(constructor<const std::string&, bool, bool>());
    py_SgOps.def(SgOps_ParseHallSymbol_parse_string, "ParseHallSymbol");
    py_SgOps.def(SgOps_ParseHallSymbol_std_string, "ParseHallSymbol");
    py_SgOps.def(&SgOps::ParseHallSymbol, "ParseHallSymbol");
    py_SgOps.def(&SgOps::ChangeBasis, "ChangeBasis");
    py_SgOps.def(&SgOps::expandSMx, "expandSMx");
    py_SgOps.def(&SgOps::RBF, "RBF");
    py_SgOps.def(&SgOps::TBF, "TBF");
    py_SgOps.def(&SgOps::nLTr, "nLTr");
    py_SgOps.def(&SgOps::fInv, "fInv");
    py_SgOps.def(&SgOps::nSMx, "nSMx");
    py_SgOps.def(&SgOps::OrderP, "OrderP");
    py_SgOps.def(&SgOps::OrderZ, "OrderZ");
    py_SgOps.def(SgOps_call_1, "__call__");
    py_SgOps.def(SgOps_call_3, "__call__");
    py_SgOps.def(&SgOps::isChiral, "isChiral");
    py_SgOps.def(&SgOps::isSysAbsent, "isSysAbsent");
    py_SgOps.def(SgOps_isCentric_0, "isCentric");
    py_SgOps.def(SgOps_isCentric_1, "isCentric");
    py_SgOps.def(&SgOps::getPhaseRestriction, "getPhaseRestriction");
    py_SgOps.def(SgOps_isValidPhase_rad_2, "isValidPhase_rad");
    py_SgOps.def(SgOps_isValidPhase_rad_3, "isValidPhase_rad");
    py_SgOps.def(SgOps_isValidPhase_deg_2, "isValidPhase_deg");
    py_SgOps.def(SgOps_isValidPhase_deg_3, "isValidPhase_deg");
    py_SgOps.def(&SgOps::epsilon, "epsilon");
    py_SgOps.def(&SgOps::multiplicity, "multiplicity");
    py_SgOps.def(&SgOps::getEquivMillerIndices, "getEquivMillerIndices");
    py_SgOps.def(&SgOps::getCutParameters, "getCutParameters");
    py_SgOps.def(SgOps_getMasterIndex_2, "getMasterIndex");
    py_SgOps.def(SgOps_getMasterIndex_3, "getMasterIndex");
    py_SgOps.def(SgOps_CheckMetricalMatrix_1, "CheckMetricalMatrix");
    py_SgOps.def(&SgOps::CheckMetricalMatrix, "CheckMetricalMatrix");
    py_SgOps.def(SgOps_CheckUnitCell_1, "CheckUnitCell");
    py_SgOps.def(&SgOps::CheckUnitCell, "CheckUnitCell");
    py_SgOps.def(SgOps_getZ2POp_0, "getZ2POp");
    py_SgOps.def(&SgOps::getZ2POp, "getZ2POp");
    py_SgOps.def(&SgOps::makeTidy, "makeTidy");
    py_SgOps.def(SgOps_cmp_equal, "__cmp__");

    sgtbx::sanity_check();
  }
  catch(...)
  {
    boost::python::handle_exception(); // Deal with the exception for Python
  }
}
