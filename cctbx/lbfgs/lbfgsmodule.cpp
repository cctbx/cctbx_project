// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Mar 2002: Created (R.W. Grosse-Kunstleve)
 */

#include <boost/python/cross_module.hpp>

#include <cctbx/error.h>
#include <cctbx/lbfgs.h>
#include <cctbx/array_family/shared.h>

namespace {

  using namespace cctbx;

  double py_minimizer_gnorm_0(
    const lbfgs::minimizer<double>& minimizer)
  {
    return minimizer.gnorm();
  }
  double py_minimizer_gnorm_1(
    const lbfgs::minimizer<double>& minimizer,
    const af::shared<double>& g)
  {
    cctbx_assert(g.size() == minimizer.n());
    return minimizer.gnorm(g.begin());
  }

  void py_minimizer_run_3(
    lbfgs::minimizer<double>& minimizer,
    af::shared<double> x,
    double f,
    const af::shared<double>& g)
  {
    cctbx_assert(x.size() == minimizer.n());
    cctbx_assert(g.size() == minimizer.n());
    minimizer.run(x.begin(), f, g.begin());
  }

  void py_minimizer_run_4(
    lbfgs::minimizer<double>& minimizer,
    af::shared<double> x,
    double f,
    const af::shared<double>& g,
    af::shared<double> diag)
  {
    cctbx_assert(x.size() == minimizer.n());
    cctbx_assert(g.size() == minimizer.n());
    cctbx_assert(diag.size() == minimizer.n());
    minimizer.run(x.begin(), f, g.begin(), diag.begin());
  }

# include <cctbx/basic/from_bpl_import.h>

  void init_module(python::module_builder& this_module)
  {
    const std::string Revision = "$Revision$";
    this_module.add(ref(to_python(
        Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

    python::import_converters<af::shared<double> >
    py_shared_double(
      "cctbx_boost.arraytbx.shared", "double");

    python::class_builder<lbfgs::minimizer<double> >
    py_minimizer(this_module, "minimizer");

    py_minimizer.def(constructor<>());
    py_minimizer.def(constructor<std::size_t>());
    py_minimizer.def(constructor<std::size_t, std::size_t>());
    py_minimizer.def(constructor<std::size_t, std::size_t, std::size_t,
      double>());
    py_minimizer.def(constructor<std::size_t, std::size_t, std::size_t,
      double, double>());
    py_minimizer.def(constructor<std::size_t, std::size_t, std::size_t,
      double, double, double>());
    py_minimizer.def(constructor<std::size_t, std::size_t, std::size_t,
      double, double, double, double>());
    py_minimizer.def(constructor<std::size_t, std::size_t, std::size_t,
      double, double, double, double, double>());
    py_minimizer.def(py_minimizer_run_3, "run");
    py_minimizer.def(py_minimizer_run_4, "run");
    py_minimizer.def(&lbfgs::minimizer<double>::n, "n");
    py_minimizer.def(&lbfgs::minimizer<double>::m, "m");
    py_minimizer.def(&lbfgs::minimizer<double>::maxfev, "maxfev");
    py_minimizer.def(&lbfgs::minimizer<double>::eps, "eps");
    py_minimizer.def(&lbfgs::minimizer<double>::gtol, "gtol");
    py_minimizer.def(&lbfgs::minimizer<double>::xtol, "xtol");
    py_minimizer.def(&lbfgs::minimizer<double>::stpmin, "stpmin");
    py_minimizer.def(&lbfgs::minimizer<double>::stpmax, "stpmax");
    py_minimizer.def(
      &lbfgs::minimizer<double>::is_converged, "is_converged");
    py_minimizer.def(
      &lbfgs::minimizer<double>::requests_f_and_g, "requests_f_and_g");
    py_minimizer.def(
      &lbfgs::minimizer<double>::requests_diag, "requests_diag");
    py_minimizer.def(&lbfgs::minimizer<double>::iter, "iter");
    py_minimizer.def(&lbfgs::minimizer<double>::nfun, "nfun");
    py_minimizer.def(py_minimizer_gnorm_0, "gnorm");
    py_minimizer.def(py_minimizer_gnorm_1, "gnorm");
    py_minimizer.def(&lbfgs::minimizer<double>::stp, "stp");
  }

}

BOOST_PYTHON_MODULE_INIT(lbfgs)
{
  boost::python::module_builder this_module("lbfgs");
  init_module(this_module);
}
