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
#include <cctbx/array_family/grid_accessor_bpl.h>
#include <cctbx/sftbx/xray_scatterer.h>
#include <cctbx/sftbx/sfmap.h>

#include <cctbx/maps/peak_search.h>

namespace {

  using namespace cctbx;

  typedef sftbx::XrayScatterer<double, eltbx::CAASF_WK1995> ex_xray_scatterer;

  void
  xray_scatterer_set_fpfdp(ex_xray_scatterer& site,
                           const std::complex<double>& fpfdp) {
    site.set_fpfdp(fpfdp);
  }

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
    bool FriedelFlag,
    double Resolution_d_min)
  {
    af::shared<Miller::Index> result;
    sgtbx::MillerIndexGenerator(
      UC, SgInfo, FriedelFlag, Resolution_d_min).AddToArray(result);
    return result;
  }
  af::shared<Miller::Index>
  py_BuildMillerIndices_MaxIndex(
    const sgtbx::SpaceGroupInfo& SgInfo,
    bool FriedelFlag,
    const Miller::Index& MaxIndex)
  {
    af::shared<Miller::Index> result;
    sgtbx::MillerIndexGenerator(
      SgInfo, FriedelFlag, MaxIndex).AddToArray(result);
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

  double
  py_calc_u_extra_4(
    double max_q, double grid_resolution_factor,
    double quality_factor,
    double max_u_extra)
  {
    return sftbx::calc_u_extra(
      max_q, grid_resolution_factor, quality_factor, max_u_extra);
  }

  double
  py_calc_u_extra_3(
    double max_q, double grid_resolution_factor,
    double quality_factor)
  {
    return sftbx::calc_u_extra(
      max_q, grid_resolution_factor, quality_factor);
  }

  double
  py_calc_u_extra_2(
    double max_q, double grid_resolution_factor)
  {
    return sftbx::calc_u_extra(
      max_q, grid_resolution_factor);
  }

  af::int3
  py_determine_grid(
    const uctbx::UnitCell& ucell,
    double max_q,
    double resolution_factor,
    int max_prime,
    const af::int3& mandatory_factors)
  {
    return maps::determine_grid<af::int3>(
      ucell, max_q, resolution_factor, max_prime, mandatory_factors);
  }

  af::shared<double>
  sampled_model_density_map_real_as_shared(
    sftbx::sampled_model_density<double>& smd)
  {
    return smd.map_real().as_base_array();
  }

  af::shared<std::complex<double> >
  sampled_model_density_map_complex_as_shared(
    sftbx::sampled_model_density<double>& smd)
  {
    return smd.map_complex().as_base_array();
  }

  af::shared<double>
  as_shared(const af::versa<double, af::grid<3> >& a) {
    return af::shared<double>(a.as_base_array());
  }

  bool
  grid_tags_verify_1(
    const maps::grid_tags<long>& tags,
    af::versa<double, af::grid<3> >& data) {
    return tags.verify(data);
  }
  bool
  grid_tags_verify_2(
    const maps::grid_tags<long>& tags,
    af::versa<double, af::grid<3> >& data,
    double min_correlation) {
    return tags.verify(data, min_correlation);
  }

  void
  py_eliminate_u_extra_5(
    const uctbx::UnitCell& ucell,
    double u_extra,
    const af::shared<Miller::Index>& miller_indices,
    af::shared<std::complex<double> > structure_factors,
    double norm)
  {
    sftbx::eliminate_u_extra(
      ucell, u_extra, miller_indices.const_ref(), structure_factors.ref(),
      norm);
  }

  void
  py_eliminate_u_extra_4(
    const uctbx::UnitCell& ucell,
    double u_extra,
    const af::shared<Miller::Index>& miller_indices,
    af::shared<std::complex<double> > structure_factors)
  {
    sftbx::eliminate_u_extra(
      ucell, u_extra, miller_indices.const_ref(), structure_factors.ref());
  }

  void sampled_model_density_apply_symmetry(
    sftbx::sampled_model_density<double>& dens,
    const maps::grid_tags<long>& tags)
  {
    dens.apply_symmetry(tags);
  }

  void sampled_model_density_eliminate_u_extra_and_normalize(
    const sftbx::sampled_model_density<double>& dens,
    const af::shared<Miller::Index>& miller_indices,
    af::shared<std::complex<double> > structure_factors)
  {
    dens.eliminate_u_extra_and_normalize(
      miller_indices.const_ref(), structure_factors.ref());
  }

  af::shared<std::complex<double> >
  py_structure_factor_map(
    const sgtbx::SpaceGroup& sgops,
    bool friedel_flag,
    const af::shared<Miller::Index>& h,
    const af::shared<std::complex<double> >& f,
    const af::tiny<long, 3>& n_complex,
    bool conjugate)
  {
    return sftbx::structure_factor_map(
      sgops, friedel_flag, h.const_ref(), f.const_ref(), n_complex,
      conjugate).as_base_array();
  }

  af::shared<std::complex<double> >
  py_collect_structure_factors_miller_indices_real(
    bool friedel_flag,
    const af::shared<Miller::Index>& miller_indices,
    const af::shared<double>& transformed_real_map,
    const af::tiny<long, 3>& n_complex,
    bool conjugate)
  {
    cctbx_assert(
         transformed_real_map.size()
      >= 2 * af::compile_time_product<3>::get(n_complex));
    af::const_ref<std::complex<double> > transformed_complex_map(
      reinterpret_cast<const std::complex<double>*>(
        transformed_real_map.begin()), transformed_real_map.size()/2);
    return sftbx::collect_structure_factors(
      friedel_flag,
      miller_indices.const_ref(),
      transformed_complex_map,
      n_complex,
      conjugate);
  }

  af::shared<std::complex<double> >
  py_collect_structure_factors_miller_indices_complex(
    bool friedel_flag,
    const af::shared<Miller::Index>& miller_indices,
    const af::shared<std::complex<double> >& transformed_complex_map,
    const af::tiny<long, 3>& n_complex,
    bool conjugate)
  {
    return sftbx::collect_structure_factors(
      friedel_flag,
      miller_indices.const_ref(),
      transformed_complex_map.const_ref(),
      n_complex,
      conjugate);
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

  tuple
  py_collect_structure_factors_max_q_real(
    const uctbx::UnitCell& ucell,
    const sgtbx::SpaceGroupInfo& sginfo,
    bool friedel_flag,
    double max_q,
    const af::shared<double>& transformed_real_map,
    const af::tiny<long, 3>& n_complex,
    bool conjugate)
  {
    cctbx_assert(
         transformed_real_map.size()
      >= 2 * af::compile_time_product<3>::get(n_complex));
    af::const_ref<std::complex<double> > transformed_complex_map(
      reinterpret_cast<const std::complex<double>*>(
        transformed_real_map.begin()), transformed_real_map.size()/2);
    std::pair<af::shared<Miller::Index>,
              af::shared<std::complex<double> > >
    indexed_structure_factors = sftbx::collect_structure_factors(
      ucell, sginfo, friedel_flag,
      max_q, transformed_complex_map, n_complex, conjugate);
    tuple result(2);
    result.set_item(0, indexed_structure_factors.first);
    result.set_item(1, indexed_structure_factors.second);
    return result;
  }

  tuple
  py_collect_structure_factors_max_q_complex(
    const uctbx::UnitCell& ucell,
    const sgtbx::SpaceGroupInfo& sginfo,
    bool friedel_flag,
    double max_q,
    const af::shared<std::complex<double> >& transformed_complex_map,
    const af::tiny<long, 3>& n_complex,
    bool conjugate)
  {
    std::pair<af::shared<Miller::Index>,
              af::shared<std::complex<double> > >
    indexed_structure_factors = sftbx::collect_structure_factors(
      ucell, sginfo, friedel_flag,
      max_q, transformed_complex_map.const_ref(), n_complex, conjugate);
    tuple result(2);
    result.set_item(0, indexed_structure_factors.first);
    result.set_item(1, indexed_structure_factors.second);
    return result;
  }

  std::size_t peak_list_size(const maps::peak_list<double>& pl) {
    return pl.size();
  }
  boost::python::dictionary
  peak_list_getitem(const maps::peak_list<double>& pl, std::size_t i) {
    boost::python::dictionary res;
    res.set_item("index", boost::python::ref(to_python(pl[i].index)));
    res.set_item("value", boost::python::ref(to_python(pl[i].value)));
    return res;
  }

  maps::peak_list<double>
  get_peak_list(
    af::int3 const& n_real,
    af::int3 const& m_real,
    af::shared<double> data,
    maps::grid_tags<long>& tags,
    int peak_search_level,
    std::size_t max_peaks)
  {
    cctbx_assert(af::product(n_real) <= af::product(m_real));
    cctbx_assert(af::product(n_real) == tags.size());
    cctbx_assert(af::product(m_real) == data.size());
    typedef af::grid<3> grid_type;
    typedef grid_type::index_type grid_point_type;
    af::ref<double, grid_type> data_3d(data.begin(), grid_type(m_real));
    af::versa<double, grid_type> data_compact_3d; // XXX why do we need the
    data_compact_3d.resize(grid_type(n_real));    // XXX separate resize()?
    af::nested_loop<grid_point_type> loop(data_compact_3d.accessor());
    for (const grid_point_type& pt = loop(); !loop.over(); loop.incr()) {
      data_compact_3d(pt) = data_3d(pt);
    }
    return maps::peak_list<double>(
      data_compact_3d, tags, peak_search_level, max_peaks);
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

    // XXX move to sgtbx
    class_builder<Miller::map_to_asym_index<std::complex<double> > >
    py_map_to_asym_index(this_module, "map_to_asym_index");

    class_builder<ex_xray_scatterer>
    py_XrayScatterer(this_module, "XrayScatterer");
    python::export_converters(py_XrayScatterer);

    class_builder<sftbx::sampled_model_density<double> >
    py_sampled_model_density(this_module, "sampled_model_density");

    class_builder<af::versa<double, af::grid<3> > >
    py_versa_double_grid_3(this_module, "versa_double_grid_3");

    class_builder<maps::map_symmetry_flags>
    py_map_symmetry_flags(this_module, "map_symmetry_flags");

    class_builder<maps::grid_tags<long> >
    py_grid_tags(this_module, "grid_tags");

    class_builder<maps::peak_list<double> >
    py_peak_list(this_module, "peak_list");

    class_builder<Miller::index_span>
    py_index_span(this_module, "index_span");

    py_map_to_asym_index.def(constructor<>());
    py_map_to_asym_index.def(constructor<
      const sgtbx::SpaceGroupInfo&,
      bool,
      af::shared<Miller::Index>,
      af::shared<std::complex<double> >,
      bool>());
    py_map_to_asym_index.def(constructor<
      const sgtbx::SpaceGroupInfo&,
      bool,
      af::shared<Miller::Index>,
      af::shared<std::complex<double> > >());
    py_map_to_asym_index.def(
      &Miller::map_to_asym_index<std::complex<double> >::friedel_flag,
                                                        "friedel_flag");
    //py_map_to_asym_index.def(
    //  &Miller::map_to_asym_index<std::complex<double> >::asu,
    //                                                    "asu");
    py_map_to_asym_index.def(
      &Miller::map_to_asym_index<std::complex<double> >::asym_miller_indices,
                                                        "asym_miller_indices");
    py_map_to_asym_index.def(
      &Miller::map_to_asym_index<std::complex<double> >::asym_data_array,
                                                        "asym_data_array");

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
      xray_scatterer_set_fpfdp, "set_fpfdp");
    py_XrayScatterer.def(
      &ex_xray_scatterer::set_Coordinates, "set_Coordinates");
    py_XrayScatterer.def(
      xray_scatterer_set_Occ_1, "set_Occ");
    py_XrayScatterer.def(
      xray_scatterer_set_Occ_2, "set_Occ");
    py_XrayScatterer.def(
      &ex_xray_scatterer::set_Uiso, "set_Uiso");
    py_XrayScatterer.def(
      &ex_xray_scatterer::set_Uaniso, "set_Uaniso");
    py_XrayScatterer.def(
      &ex_xray_scatterer::StructureFactor, "StructureFactor");
    py_XrayScatterer.def(
      xray_scatterer_copy, "copy");

    py_sampled_model_density.def(constructor<>());
    py_sampled_model_density.def(constructor<
      const uctbx::UnitCell&,
      const af::shared<ex_xray_scatterer>&,
      const af::tiny<long, 3>&,
      const af::tiny<long, 3>&>());
    py_sampled_model_density.def(constructor<
      const uctbx::UnitCell&,
      const af::shared<ex_xray_scatterer>&,
      const af::tiny<long, 3>&,
      const af::tiny<long, 3>&,
      double>());
    py_sampled_model_density.def(constructor<
      const uctbx::UnitCell&,
      const af::shared<ex_xray_scatterer>&,
      const af::tiny<long, 3>&,
      const af::tiny<long, 3>&,
      double,
      double,
      double>());
    py_sampled_model_density.def(constructor<
      const uctbx::UnitCell&,
      const af::shared<ex_xray_scatterer>&,
      const af::tiny<long, 3>&,
      const af::tiny<long, 3>&,
      double,
      double,
      double,
      bool>());
    py_sampled_model_density.def(constructor<
      const uctbx::UnitCell&,
      const af::shared<ex_xray_scatterer>&,
      const af::tiny<long, 3>&,
      const af::tiny<long, 3>&,
      double,
      double,
      double,
      bool,
      bool>());
    py_sampled_model_density.def(
      &sftbx::sampled_model_density<double>::ucell,
                                            "ucell");
    py_sampled_model_density.def(
      &sftbx::sampled_model_density<double>::u_extra,
                                            "u_extra");
    py_sampled_model_density.def(
      &sftbx::sampled_model_density<double>::wing_cutoff,
                                            "wing_cutoff");
    py_sampled_model_density.def(
      &sftbx::sampled_model_density<double>::exp_table_one_over_step_size,
                                            "exp_table_one_over_step_size");
    py_sampled_model_density.def(
      &sftbx::sampled_model_density<double>::exp_table_size,
                                            "exp_table_size");
    py_sampled_model_density.def(
      &sftbx::sampled_model_density<double>::n_passed_scatterers,
                                            "n_passed_scatterers");
    py_sampled_model_density.def(
      &sftbx::sampled_model_density<double>::n_contributing_scatterers,
                                            "n_contributing_scatterers");
    py_sampled_model_density.def(
      &sftbx::sampled_model_density<double>::n_anomalous_scatterers,
                                            "n_anomalous_scatterers");
    py_sampled_model_density.def(
      &sftbx::sampled_model_density<double>::max_shell_radii,
                                            "max_shell_radii");
    py_sampled_model_density.def(
      &sftbx::sampled_model_density<double>::max_shell_radii_frac,
                                            "max_shell_radii_frac");
    py_sampled_model_density.def(
      &sftbx::sampled_model_density<double>::friedel_flag,
                                            "friedel_flag");
    py_sampled_model_density.def(
      sampled_model_density_map_real_as_shared,
                           "map_real_as_shared");
    py_sampled_model_density.def(
      sampled_model_density_map_complex_as_shared,
                           "map_complex_as_shared");
    py_sampled_model_density.def(
      sampled_model_density_apply_symmetry,
                           "apply_symmetry");
    py_sampled_model_density.def(
      sampled_model_density_eliminate_u_extra_and_normalize,
                           "eliminate_u_extra_and_normalize");

    this_module.def(py_StructureFactor_plain, "StructureFactor");
    this_module.def(py_StructureFactor_iso,   "StructureFactor");
    this_module.def(py_StructureFactor_aniso, "StructureFactor");

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

    // XXX move to sgtbx
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

    this_module.def(py_calc_u_extra_4, "calc_u_extra");
    this_module.def(py_calc_u_extra_3, "calc_u_extra");
    this_module.def(py_calc_u_extra_2, "calc_u_extra");
    this_module.def(py_eliminate_u_extra_5, "eliminate_u_extra");
    this_module.def(py_eliminate_u_extra_4, "eliminate_u_extra");
    this_module.def(py_collect_structure_factors_max_q_real,
      "collect_structure_factors");
    this_module.def(py_collect_structure_factors_max_q_complex,
      "collect_structure_factors");
    this_module.def(py_collect_structure_factors_miller_indices_real,
      "collect_structure_factors");
    this_module.def(py_collect_structure_factors_miller_indices_complex,
      "collect_structure_factors");
    this_module.def(py_structure_factor_map, "structure_factor_map");

    py_versa_double_grid_3.def(constructor<>());
    py_versa_double_grid_3.def(as_shared, "as_shared");

    py_map_symmetry_flags.def(constructor<>());
    py_map_symmetry_flags.def(constructor<bool>());
    py_map_symmetry_flags.def(constructor<bool, bool>());
    py_map_symmetry_flags.def(constructor<bool, bool, bool>());
    py_map_symmetry_flags.def(
      &maps::map_symmetry_flags::use_space_group_symmetry,
                                "use_space_group_symmetry");
    py_map_symmetry_flags.def(
      &maps::map_symmetry_flags::use_normalizer_K2L,
                                "use_normalizer_K2L");
    py_map_symmetry_flags.def(
      &maps::map_symmetry_flags::use_structure_seminvariants,
                                "use_structure_seminvariants");
    py_map_symmetry_flags.def(
      &maps::map_symmetry_flags::select_sub_space_group,
                                "select_sub_space_group");
    py_map_symmetry_flags.def(
      &maps::map_symmetry_flags::get_grid_factors,
                                "get_grid_factors");

    py_grid_tags.def(constructor<>());
    py_grid_tags.def(constructor<const cctbx::af::grid<3>&>());
    py_grid_tags.def(&maps::grid_tags<long>::build, "build");
    py_grid_tags.def(&maps::grid_tags<long>::is_valid, "is_valid");
    py_grid_tags.def(&maps::grid_tags<long>::SgInfo, "SgInfo");
    py_grid_tags.def(&maps::grid_tags<long>::sym_flags, "sym_flags");
    py_grid_tags.def(&maps::grid_tags<long>::n_grid_misses, "n_grid_misses");
    py_grid_tags.def(&maps::grid_tags<long>::n_independent, "n_independent");
    py_grid_tags.def(grid_tags_verify_1, "verify");
    py_grid_tags.def(grid_tags_verify_2, "verify");

    py_peak_list.def(constructor<>());
    py_peak_list.def(peak_list_size, "size");
    py_peak_list.def(peak_list_size, "__len__");
    py_peak_list.def(peak_list_getitem, "__getitem__");

    this_module.def(get_peak_list, "get_peak_list");

    this_module.def(py_determine_grid, "determine_grid");

    py_index_span.def(constructor<>());
    py_index_span.def(constructor<af::shared<Miller::Index> >());
    py_index_span.def(&Miller::index_span::min, "min");
    py_index_span.def(&Miller::index_span::max, "max");
    py_index_span.def(&Miller::index_span::abs_range, "abs_range");
    py_index_span.def(&Miller::index_span::map_grid, "map_grid");
    py_index_span.def(&Miller::index_span::is_in_domain, "is_in_domain");
  }

}

BOOST_PYTHON_MODULE_INIT(sftbx)
{
  boost::python::module_builder this_module("sftbx");
  init_module(this_module);
}
