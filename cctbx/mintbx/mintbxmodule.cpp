// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jul 2002: Created (R.W. Grosse-Kunstleve)
 */

#include <boost/python/cross_module.hpp>
#include <cctbx/miller_bpl.h>
#include <cctbx/mintbx/k_b_scaling.h>

using namespace cctbx;
using namespace cctbx::mintbx;

namespace {

# include <cctbx/basic/from_bpl_import.h>

  void init_module(python::module_builder& this_module)
  {
    const std::string Revision = "$Revision$";
    this_module.add(ref(to_python(
        Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

    python::import_converters<uctbx::UnitCell>
    py_UnitCell("cctbx_boost.uctbx", "UnitCell");

    python::import_converters<sgtbx::SpaceGroup>
    py_SpaceGroup("cctbx_boost.sgtbx", "SpaceGroup");

    python::import_converters<af::shared<int> >
    py_shared_int("cctbx_boost.arraytbx.shared", "int");

    python::import_converters<af::shared<double> >
    py_shared_double("cctbx_boost.arraytbx.shared", "double");

    python::import_converters<af::shared<miller::Index> >
    py_shared_miller_Index(
      "cctbx_boost.arraytbx.shared", "miller_Index");

    typedef k_b_scaling_target_and_gradients<double> kbstg;
    class_builder<kbstg>
    py_kbstg(this_module, "k_b_scaling_target_and_gradients");

    py_kbstg.def(constructor<>());
    py_kbstg.def(constructor<
      uctbx::UnitCell const&,
      af::shared<miller::Index>,
      af::shared<int>,
      af::shared<double>,
      af::shared<double>,
      double,
      double,
      bool,
      bool
      >());
    py_kbstg.def(constructor<
      uctbx::UnitCell const&,
      af::shared<miller::Index>,
      af::shared<int>,
      af::shared<double>,
      af::shared<double>,
      double,
      af::tiny<double, 6> const&,
      bool,
      bool
      >());
    py_kbstg.def(&kbstg::target, "target");
    py_kbstg.def(&kbstg::gradient_k, "gradient_k");
    py_kbstg.def(&kbstg::anisotropic, "anisotropic");
    py_kbstg.def(&kbstg::gradient_b_iso, "gradient_b_iso");
    py_kbstg.def(&kbstg::gradients_b_cif, "gradients_b_cif");
  }

}

BOOST_PYTHON_MODULE_INIT(mintbx)
{
  boost::python::module_builder this_module("mintbx");
  init_module(this_module);
}
