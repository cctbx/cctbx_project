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
#include <cctbx/xray_scatterer.h>

namespace {

  using namespace cctbx;

  boost::shared_ptr<std::vector<std::complex<double> > >
  py_StructureFactorVector(
    const uctbx::UnitCell& UC,
    const sgtbx::SpaceGroup& SgOps,
    const std::vector<Miller::Index>& H,
    const std::vector<sftbx::XrayScatterer<double, eltbx::CAASF_WK1995> >&
      Sites)
  {
    boost::shared_ptr<std::vector<std::complex<double> > >
    Fcalc(new std::vector<std::complex<double> >(H.size()));
    sftbx::StructureFactorVector(UC, SgOps, H, Sites, *Fcalc);
    return Fcalc;
  }

  boost::shared_ptr<std::vector<Miller::Index> >
  py_BuildMillerIndices_Resolution_d_min(
    const uctbx::UnitCell& UC,
    const sgtbx::SpaceGroupInfo& SgInfo,
    double Resolution_d_min)
  {
    boost::shared_ptr<std::vector<Miller::Index> >
    VectorOfH(new std::vector<Miller::Index>);
    sgtbx::MillerIndexGenerator(
      UC, SgInfo, Resolution_d_min).AddToVector(*VectorOfH);
    return VectorOfH;
  }
  boost::shared_ptr<std::vector<Miller::Index> >
  py_BuildMillerIndices_MaxIndex(
    const sgtbx::SpaceGroupInfo& SgInfo,
    const Miller::Index& MaxIndex)
  {
    boost::shared_ptr<std::vector<Miller::Index> >
    VectorOfH(new std::vector<Miller::Index>);
    sgtbx::MillerIndexGenerator(
      SgInfo, MaxIndex).AddToVector(*VectorOfH);
    return VectorOfH;
  }

#   include <cctbx/basic/from_bpl_import.h>

  void init_module(python::module_builder& this_module)
  {
    const std::string Revision = "$Revision$";
    this_module.add(ref(to_python(
        Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

    python::import_converters<uctbx::UnitCell>
    py_UnitCell("cctbx.uctbx", "UnitCell");
    python::import_converters<sgtbx::SpaceGroup>
    py_SpaceGroup("cctbx.sgtbx", "SpaceGroup");
    python::import_converters<sgtbx::SpaceGroupInfo>
    py_SpaceGroupInfo("cctbx.sgtbx", "SpaceGroupInfo");
    python::import_converters<eltbx::CAASF_WK1995>
    py_CAASF_WK1995("cctbx.eltbx.caasf_wk1995", "CAASF_WK1995");

    python::import_converters<std::vector<std::complex<double> > >
    py_std_vector_complex_double(
      "cctbx.arraytbx.std_vector", "complex_double");

    python::import_converters<std::vector<Miller::Index> >
    py_std_vector_Miller_Index(
      "cctbx.arraytbx.std_vector", "Miller_Index");

    python::import_converters<
      std::vector<sftbx::XrayScatterer<double, eltbx::CAASF_WK1995> > >
    py_vector_XrayScatterer("cctbx.arraytbx.std_vector", "XrayScatterer");

    class_builder<sftbx::XrayScatterer<double, eltbx::CAASF_WK1995> >
    py_XrayScatterer(this_module, "XrayScatterer");
    python::export_converters(py_XrayScatterer);

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
      const double6&>());
    py_XrayScatterer.def(
      &sftbx::XrayScatterer<double, eltbx::CAASF_WK1995
      >::Label, "Label");
    py_XrayScatterer.def(
      &sftbx::XrayScatterer<double, eltbx::CAASF_WK1995
      >::CAASF, "CAASF");
    py_XrayScatterer.def(
      &sftbx::XrayScatterer<double, eltbx::CAASF_WK1995
      >::fpfdp, "fpfdp");
    py_XrayScatterer.def(
      &sftbx::XrayScatterer<double, eltbx::CAASF_WK1995
      >::Coordinates, "Coordinates");
    py_XrayScatterer.def(
      &sftbx::XrayScatterer<double, eltbx::CAASF_WK1995
      >::Occ, "Occ");
    py_XrayScatterer.def(
      &sftbx::XrayScatterer<double, eltbx::CAASF_WK1995
      >::isAnisotropic, "isAnisotropic");
    py_XrayScatterer.def(
      &sftbx::XrayScatterer<double, eltbx::CAASF_WK1995
      >::Uiso, "Uiso");
    py_XrayScatterer.def(
      &sftbx::XrayScatterer<double, eltbx::CAASF_WK1995
      >::Uaniso, "Uaniso");
    py_XrayScatterer.def(
      &sftbx::XrayScatterer<double, eltbx::CAASF_WK1995
      >::M, "M");
    py_XrayScatterer.def(
      &sftbx::XrayScatterer<double, eltbx::CAASF_WK1995
      >::w, "w");
    py_XrayScatterer.def(
      &sftbx::XrayScatterer<double, eltbx::CAASF_WK1995
      >::ApplySymmetry, "ApplySymmetry");
    py_XrayScatterer.def(
      &sftbx::XrayScatterer<double, eltbx::CAASF_WK1995
      >::StructureFactor, "StructureFactor");

    this_module.def(py_StructureFactorVector, "StructureFactorVector");
    this_module.def(py_BuildMillerIndices_Resolution_d_min,
                      "BuildMillerIndices");
    this_module.def(py_BuildMillerIndices_MaxIndex,
                      "BuildMillerIndices");
  }

}

BOOST_PYTHON_MODULE_INIT(sftbx)
{
  boost::python::module_builder this_module("sftbx");
  init_module(this_module);
}
