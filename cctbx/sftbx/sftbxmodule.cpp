// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Oct 22: Created (R.W. Grosse-Kunstleve)
 */

#include <boost/python/cross_module.hpp>
#include <cctbx/miller_bpl.h>
#include <cctbx/coordinates_bpl.h>
#include <cctbx/sftbx/xray_scatterer.h>
#include <cctbx/sftbx/sfmap.h>

namespace {

  using namespace cctbx;

  typedef sftbx::XrayScatterer<double, eltbx::CAASF_WK1995> ex_xray_scatterer;

  void
  xray_scatterer_set_Occ_1(ex_xray_scatterer& site,
                           double Occ) {
    site.set_Occ(Occ);
  }
  void
  xray_scatterer_set_Occ_2(ex_xray_scatterer& site,
                           double Occ,
                           const sgtbx::SpaceGroup& SgOps) {
    site.set_Occ(Occ, SgOps);
  }

  ex_xray_scatterer
  xray_scatterer_copy(const ex_xray_scatterer& site) {
    return site;
  }

  fractional<double>
  py_least_squares_shift(
    const uctbx::UnitCell& ucell,
    const af::shared<ex_xray_scatterer>& sites1,
    const af::shared<ex_xray_scatterer>& sites2)
  {
    return sftbx::least_squares_shift(
      ucell, sites1.const_ref(), sites2.const_ref());
  }

  double
  py_rms_coordinates_plain(
    const uctbx::UnitCell& ucell,
    const af::shared<ex_xray_scatterer>& sites1,
    const af::shared<ex_xray_scatterer>& sites2)
  {
    return sftbx::rms_coordinates(
      ucell, sites1.const_ref(), sites2.const_ref());
  }
  double
  py_rms_coordinates_shift(
    const uctbx::UnitCell& ucell,
    const af::shared<ex_xray_scatterer>& sites1,
    const af::shared<ex_xray_scatterer>& sites2,
    const fractional<double>& shift)
  {
    return sftbx::rms_coordinates(
      ucell, sites1.const_ref(), sites2.const_ref(), shift);
  }

  std::complex<double>
  py_StructureFactor_plain(const sgtbx::SpaceGroup& SgOps,
                           const Miller::Index& H,
                           const fractional<double>& X) {
    return sftbx::StructureFactor(SgOps, H, X);
  }
  std::complex<double>
  py_StructureFactor_iso(const sgtbx::SpaceGroup& SgOps,
                         const uctbx::UnitCell& UCell,
                         const Miller::Index& H,
                         const fractional<double>& X,
                         double Uiso) {
    return sftbx::StructureFactor(SgOps, UCell, H, X, Uiso);
  }
  std::complex<double>
  py_StructureFactor_aniso(const sgtbx::SpaceGroup& SgOps,
                           const Miller::Index& H,
                           const fractional<double>& X,
                           const af::double6& Ustar) {
    return sftbx::StructureFactor(SgOps, H, X, Ustar);
  }

  void
  py_StructureFactorArray_add(
    const uctbx::UnitCell& UC,
    const sgtbx::SpaceGroup& SgOps,
    const af::shared<Miller::Index>& H,
    const af::shared<ex_xray_scatterer>& Sites,
    af::shared<std::complex<double> > Fcalc)
  {
    sftbx::StructureFactorArray(
      UC, SgOps, H.const_ref(), Sites.const_ref(),
      Fcalc.ref());
  }

  af::shared<std::complex<double> >
  py_StructureFactorArray_new(
    const uctbx::UnitCell& UC,
    const sgtbx::SpaceGroup& SgOps,
    const af::shared<Miller::Index>& H,
    const af::shared<ex_xray_scatterer>& Sites)
  {
    af::shared<std::complex<double> > fcalc(H.size());
    sftbx::StructureFactorArray(
      UC, SgOps, H.const_ref(), Sites.const_ref(), fcalc.ref());
    return fcalc;
  }

  void
  py_StructureFactor_dT_dX_Array_add(
    const uctbx::UnitCell& UC,
    const sgtbx::SpaceGroup& SgOps,
    const af::shared<Miller::Index>& H,
    const af::shared<std::complex<double> >& dTarget_dFcalc,
    const af::shared<ex_xray_scatterer>& Sites,
    af::shared<af::double3> dT_dX)
  {
    sftbx::StructureFactor_dT_dX_Array(
      UC, SgOps, H.const_ref(), dTarget_dFcalc.const_ref(), Sites.const_ref(),
      dT_dX.ref());
  }

  af::shared<af::double3>
  py_StructureFactor_dT_dX_Array_new(
    const uctbx::UnitCell& UC,
    const sgtbx::SpaceGroup& SgOps,
    const af::shared<Miller::Index>& H,
    const af::shared<std::complex<double> >& dTarget_dFcalc,
    const af::shared<ex_xray_scatterer>& Sites)
  {
    af::shared<af::double3> dT_dX(Sites.size());
    dT_dX.fill(af::double3(0.,0.,0.));
    sftbx::StructureFactor_dT_dX_Array(
      UC, SgOps, H.const_ref(), dTarget_dFcalc.const_ref(), Sites.const_ref(),
      dT_dX.ref());
    return dT_dX;
  }

  void
  py_StructureFactor_dT_dOcc_Array_add(
    const uctbx::UnitCell& UC,
    const sgtbx::SpaceGroup& SgOps,
    const af::shared<Miller::Index>& H,
    const af::shared<std::complex<double> >& dTarget_dFcalc,
    const af::shared<ex_xray_scatterer>& Sites,
    af::shared<double> dT_dOcc)
  {
    sftbx::StructureFactor_dT_dOcc_Array(
      UC, SgOps, H.const_ref(), dTarget_dFcalc.const_ref(), Sites.const_ref(),
      dT_dOcc.ref());
  }

  af::shared<double>
  py_StructureFactor_dT_dOcc_Array_new(
    const uctbx::UnitCell& UC,
    const sgtbx::SpaceGroup& SgOps,
    const af::shared<Miller::Index>& H,
    const af::shared<std::complex<double> >& dTarget_dFcalc,
    const af::shared<ex_xray_scatterer>& Sites)
  {
    af::shared<double> dT_dOcc(Sites.size());
    sftbx::StructureFactor_dT_dOcc_Array(
      UC, SgOps, H.const_ref(), dTarget_dFcalc.const_ref(), Sites.const_ref(),
      dT_dOcc.ref());
    return dT_dOcc;
  }

  void
  py_StructureFactor_dT_dUiso_Array_add(
    const uctbx::UnitCell& UC,
    const sgtbx::SpaceGroup& SgOps,
    const af::shared<Miller::Index>& H,
    const af::shared<std::complex<double> >& dTarget_dFcalc,
    const af::shared<ex_xray_scatterer>& Sites,
    af::shared<double> dT_dUiso)
  {
    sftbx::StructureFactor_dT_dUiso_Array(
      UC, SgOps, H.const_ref(), dTarget_dFcalc.const_ref(), Sites.const_ref(),
      dT_dUiso.ref());
  }

  af::shared<double>
  py_StructureFactor_dT_dUiso_Array_new(
    const uctbx::UnitCell& UC,
    const sgtbx::SpaceGroup& SgOps,
    const af::shared<Miller::Index>& H,
    const af::shared<std::complex<double> >& dTarget_dFcalc,
    const af::shared<ex_xray_scatterer>& Sites)
  {
    af::shared<double> dT_dUiso(Sites.size());
    sftbx::StructureFactor_dT_dUiso_Array(
      UC, SgOps, H.const_ref(), dTarget_dFcalc.const_ref(), Sites.const_ref(),
      dT_dUiso.ref());
    return dT_dUiso;
  }

  void py_apply_special_position_ops(
    af::shared<af::double3> dT_dX,
    af::shared<sgtbx::RTMx> special_position_ops)
  {
    cctbx_assert(dT_dX.size() == special_position_ops.size());
    for(std::size_t i=0;i<dT_dX.size();i++) {
      dT_dX[i] = MatrixLite::vector_mul_matrix(
        dT_dX[i], special_position_ops[i].Rpart().as_array(double()));
    }
  }

  af::shared<Miller::Index>
  py_BuildMillerIndices_Resolution_d_min(
    const uctbx::UnitCell& UC,
    const sgtbx::SpaceGroupInfo& SgInfo,
    double Resolution_d_min)
  {
    af::shared<Miller::Index> result;
    sgtbx::MillerIndexGenerator(
      UC, SgInfo, Resolution_d_min).AddToArray(result);
    return result;
  }
  af::shared<Miller::Index>
  py_BuildMillerIndices_MaxIndex(
    const sgtbx::SpaceGroupInfo& SgInfo,
    const Miller::Index& MaxIndex)
  {
    af::shared<Miller::Index> result;
    sgtbx::MillerIndexGenerator(
      SgInfo, MaxIndex).AddToArray(result);
    return result;
  }

  std::size_t
  pack_parameters_generic(
    const uctbx::UnitCell* UCell,
    const af::shared<ex_xray_scatterer>& sites,
    af::shared<double> x,
    bool coordinates, bool occ, bool uiso)
  {
    if (coordinates) {
      if (!UCell) {
        for(std::size_t i=0;i<sites.size();i++) {
          const af::double3& c = sites[i].Coordinates();
          x.insert(x.end(), c.begin(), c.end());
        }
      }
      else {
        for(std::size_t i=0;i<sites.size();i++) {
          const af::double3 c = UCell->orthogonalize(sites[i].Coordinates());
          x.insert(x.end(), c.begin(), c.end());
        }
      }
    }
    if (occ) {
      for(std::size_t i=0;i<sites.size();i++) {
        x.push_back(sites[i].Occ());
      }
    }
    if (uiso) {
      for(std::size_t i=0;i<sites.size();i++) {
        x.push_back(sites[i].Uiso());
      }
    }
    return x.size();
  }

  std::size_t
  unpack_parameters_generic(
    const uctbx::UnitCell* UCell,
    af::shared<double> x,
    std::size_t start,
    af::shared<ex_xray_scatterer> sites,
    bool coordinates, bool occ, bool uiso)
  {
    if (coordinates) {
      cctbx_assert(x.size() >= start + sites.size() * 3);
      if (!UCell) {
        for(std::size_t i=0;i<sites.size();i++) {
          sites[i].set_Coordinates(fractional<double>(&x[i * 3]));
        }
      }
      else {
        for(std::size_t i=0;i<sites.size();i++) {
          sites[i].set_Coordinates(
            UCell->fractionalize(cartesian<double>(&x[i * 3])));
        }
      }
      start += sites.size() * 3;
    }
    if (occ) {
      cctbx_assert(x.size() >= start + sites.size());
      for(std::size_t i=0;i<sites.size();i++) {
        sites[i].set_Occ(x[i]);
      }
      start += sites.size();
    }
    if (uiso) {
      cctbx_assert(x.size() >= start + sites.size());
      for(std::size_t i=0;i<sites.size();i++) {
        sites[i].set_Uiso(x[i]);
      }
      start += sites.size();
    }
    return start;
  }

  std::size_t
  pack_parameters_frac(
    const af::shared<ex_xray_scatterer>& sites,
    af::shared<double> x,
    bool coordinates, bool occ, bool uiso)
  {
    return pack_parameters_generic(
      0, sites, x, coordinates, occ, uiso);
  }

  std::size_t
  unpack_parameters_frac(
    af::shared<double> x,
    std::size_t start,
    af::shared<ex_xray_scatterer> sites,
    bool coordinates, bool occ, bool uiso)
  {
    return unpack_parameters_generic(
      0, x, start, sites, coordinates, occ, uiso);
  }

  std::size_t
  pack_parameters_cart(
    const uctbx::UnitCell& UCell,
    const af::shared<ex_xray_scatterer>& sites,
    af::shared<double> x,
    bool coordinates, bool occ, bool uiso)
  {
    return pack_parameters_generic(
      &UCell, sites, x, coordinates, occ, uiso);
  }

  std::size_t
  unpack_parameters_cart(
    const uctbx::UnitCell& UCell,
    af::shared<double> x,
    std::size_t start,
    af::shared<ex_xray_scatterer> sites,
    bool coordinates, bool occ, bool uiso)
  {
    return unpack_parameters_generic(
      &UCell, x, start, sites, coordinates, occ, uiso);
  }

  af::shared<double>
  flatten(af::shared<af::double3> a) {
    return af::shared<double>(a.handle());
  }

  void
  dT_dX_inplace_frac_as_cart(
    const uctbx::UnitCell& UCell,
    af::shared<af::double3> dT_dX)
  {
    for(std::size_t i=0;i<dT_dX.size();i++) {
      dT_dX[i] = MatrixLite::vector_mul_matrix(
        dT_dX[i], UCell.getFractionalizationMatrix());
    }
  }

  sftbx::sampled_model_density<double>
  py_sample_model_density(
    const uctbx::UnitCell& ucell,
    const af::shared<ex_xray_scatterer>& sites,
    double max_q,
    double resolution_factor,
    int max_prime,
    const af::int3& mandatory_factors)
  {
    return sftbx::sample_model_density<double>(
      ucell, sites.const_ref(),
      max_q, resolution_factor, max_prime, mandatory_factors);
  }

#   include <cctbx/basic/from_bpl_import.h>

  tuple
  py_StructureFactor_dT_dX_dUiso_Array_new(
    const uctbx::UnitCell& UC,
    const sgtbx::SpaceGroup& SgOps,
    const af::shared<Miller::Index>& H,
    const af::shared<std::complex<double> >& dTarget_dFcalc,
    const af::shared<ex_xray_scatterer>& Sites)
  {
    af::shared<af::double3> dT_dX(Sites.size());
    dT_dX.fill(af::double3(0.,0.,0.));
    af::shared<double> dT_dUiso(Sites.size());
    sftbx::StructureFactor_dT_dX_dUiso_Array(
      UC, SgOps, H.const_ref(), dTarget_dFcalc.const_ref(), Sites.const_ref(),
      dT_dX.ref(), dT_dUiso.ref());
    tuple result(2);
    result.set_item(0, ref(to_python(dT_dX)));
    result.set_item(1, ref(to_python(dT_dUiso)));
    return result;
  }

  void init_module(python::module_builder& this_module)
  {
    const std::string Revision = "$Revision$";
    this_module.add(ref(to_python(
        Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

    python::import_converters<uctbx::UnitCell>
    py_UnitCell("cctbx_boost.uctbx", "UnitCell");
    python::import_converters<sgtbx::SpaceGroup>
    py_SpaceGroup("cctbx_boost.sgtbx", "SpaceGroup");
    python::import_converters<sgtbx::SpaceGroupInfo>
    py_SpaceGroupInfo("cctbx_boost.sgtbx", "SpaceGroupInfo");
    python::import_converters<sgtbx::SiteSymmetry>
    py_SiteSymmetry("cctbx_boost.sgtbx", "SiteSymmetry");
    python::import_converters<eltbx::CAASF_WK1995>
    py_CAASF_WK1995("cctbx_boost.eltbx.caasf_wk1995", "CAASF_WK1995");

    python::import_converters<af::shared<double> >
    py_shared_double("cctbx_boost.arraytbx.shared", "double");

    python::import_converters<af::shared<std::complex<double> > >
    py_shared_complex_double(
      "cctbx_boost.arraytbx.shared", "complex_double");

    python::import_converters<af::shared<Miller::Index> >
    py_shared_Miller_Index(
      "cctbx_boost.arraytbx.shared", "Miller_Index");

    python::import_converters<af::shared<ex_xray_scatterer> >
    py_shared_XrayScatterer("cctbx_boost.arraytbx.shared", "XrayScatterer");

    python::import_converters<af::shared<af::double3> >
    py_shared_double3("cctbx_boost.arraytbx.shared", "double3");

    python::import_converters<af::shared<sgtbx::RTMx> >
    py_shared_RTMx(
      "cctbx_boost.arraytbx.shared", "RTMx");

    class_builder<ex_xray_scatterer>
    py_XrayScatterer(this_module, "XrayScatterer");
    python::export_converters(py_XrayScatterer);

    class_builder<sftbx::sampled_model_density<double> >
    py_sampled_model_density(this_module, "sampled_model_density");

    this_module.def(py_sample_model_density, "sample_model_density");

    this_module.def(py_StructureFactor_plain, "StructureFactor");
    this_module.def(py_StructureFactor_iso,   "StructureFactor");
    this_module.def(py_StructureFactor_aniso, "StructureFactor");

    py_XrayScatterer.def(constructor<>());
    py_XrayScatterer.def(constructor<
      const std::string&,
      const eltbx::CAASF_WK1995&,
      const std::complex<double>&,
      const fractional<double>&,
      const double&,
      const double&>());
    py_XrayScatterer.def(constructor<
      const std::string&,
      const eltbx::CAASF_WK1995&,
      const std::complex<double>&,
      const fractional<double>&,
      const double&,
      const af::double6&>());
    py_XrayScatterer.def(
      &ex_xray_scatterer::Label, "Label");
    py_XrayScatterer.def(
      &ex_xray_scatterer::CAASF, "CAASF");
    py_XrayScatterer.def(
      &ex_xray_scatterer::fpfdp, "fpfdp");
    py_XrayScatterer.def(
      &ex_xray_scatterer::Coordinates, "Coordinates");
    py_XrayScatterer.def(
      &ex_xray_scatterer::Occ, "Occ");
    py_XrayScatterer.def(
      &ex_xray_scatterer::isAnisotropic, "isAnisotropic");
    py_XrayScatterer.def(
      &ex_xray_scatterer::Uiso, "Uiso");
    py_XrayScatterer.def(
      &ex_xray_scatterer::Uaniso, "Uaniso");
    py_XrayScatterer.def(
      &ex_xray_scatterer::M, "M");
    py_XrayScatterer.def(
      &ex_xray_scatterer::w, "w");
    py_XrayScatterer.def(
      &ex_xray_scatterer::CheckMultiplicity, "CheckMultiplicity");
    py_XrayScatterer.def(
      &ex_xray_scatterer::difference, "difference");
    py_XrayScatterer.def(
      &ex_xray_scatterer::distance2, "distance2");
    py_XrayScatterer.def(
      &ex_xray_scatterer::distance, "distance");
    py_XrayScatterer.def(
      &ex_xray_scatterer::ApplySymmetry, "ApplySymmetry");
    py_XrayScatterer.def(
      &ex_xray_scatterer::set_Coordinates, "set_Coordinates");
    py_XrayScatterer.def(
      xray_scatterer_set_Occ_1, "set_Occ");
    py_XrayScatterer.def(
      xray_scatterer_set_Occ_2, "set_Occ");
    py_XrayScatterer.def(
      &ex_xray_scatterer::set_Uiso, "set_Uiso");
    py_XrayScatterer.def(
      &ex_xray_scatterer::StructureFactor, "StructureFactor");
    py_XrayScatterer.def(
      xray_scatterer_copy, "copy");

    py_sampled_model_density.def(constructor<>());
    py_sampled_model_density.def(
      &sftbx::sampled_model_density<double>::max_q,
                                            "max_q");
    py_sampled_model_density.def(
      &sftbx::sampled_model_density<double>::resolution_factor,
                                            "resolution_factor");
    py_sampled_model_density.def(
      &sftbx::sampled_model_density<double>::quality_factor,
                                            "quality_factor");
    py_sampled_model_density.def(
      &sftbx::sampled_model_density<double>::wing_cutoff,
                                            "wing_cutoff");
    py_sampled_model_density.def(
      &sftbx::sampled_model_density<double>::u_extra,
                                            "u_extra");

    this_module.def(py_least_squares_shift, "least_squares_shift");
    this_module.def(py_rms_coordinates_plain, "rms_coordinates");
    this_module.def(py_rms_coordinates_shift, "rms_coordinates");

    this_module.def(py_StructureFactorArray_add, "StructureFactorArray");
    this_module.def(py_StructureFactorArray_new, "StructureFactorArray");
    this_module.def(
      py_StructureFactor_dT_dX_Array_add, "StructureFactor_dT_dX_Array");
    this_module.def(
      py_StructureFactor_dT_dX_Array_new, "StructureFactor_dT_dX_Array");
    this_module.def(
      py_StructureFactor_dT_dOcc_Array_add, "StructureFactor_dT_dOcc_Array");
    this_module.def(
      py_StructureFactor_dT_dOcc_Array_new, "StructureFactor_dT_dOcc_Array");
    this_module.def(
      py_StructureFactor_dT_dUiso_Array_add, "StructureFactor_dT_dUiso_Array");
    this_module.def(
      py_StructureFactor_dT_dUiso_Array_new, "StructureFactor_dT_dUiso_Array");
    this_module.def(
      py_StructureFactor_dT_dX_dUiso_Array_new,
        "StructureFactor_dT_dX_dUiso_Array");

    this_module.def(py_apply_special_position_ops, "apply");

    this_module.def(py_BuildMillerIndices_Resolution_d_min,
                      "BuildMillerIndices");
    this_module.def(py_BuildMillerIndices_MaxIndex,
                      "BuildMillerIndices");

    this_module.def(pack_parameters_frac, "pack_parameters");
    this_module.def(unpack_parameters_frac, "unpack_parameters");
    this_module.def(pack_parameters_cart, "pack_parameters");
    this_module.def(unpack_parameters_cart, "unpack_parameters");
    this_module.def(flatten, "flatten");
    this_module.def(dT_dX_inplace_frac_as_cart, "dT_dX_inplace_frac_as_cart");
  }

}

BOOST_PYTHON_MODULE_INIT(sftbx)
{
  boost::python::module_builder this_module("sftbx");
  init_module(this_module);
}
