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

#include <cctbx/array_family/flex_bpl.h>

namespace cctbx { namespace af { namespace bpl { namespace {

  void import_flex()
  {
    CCTBX_ARRAY_FAMILY_FLEX_IMPORT(double, "double")
  }

}}}} // namespace cctbx::af::bpl<anonymous>

CCTBX_ARRAY_FAMILY_IMPLICIT_SHARED_CONVERTERS(double)

namespace {

  using namespace cctbx;

  double py_minimizer_euclidean_norm(
    const lbfgs::minimizer<double>& minimizer,
    const af::shared<double>& a)
  {
    cctbx_assert(a.size() == minimizer.n());
    return minimizer.euclidean_norm(a.begin());
  }

  bool py_minimizer_run_3(
    lbfgs::minimizer<double>& minimizer,
    af::shared<double> x,
    double f,
    const af::shared<double>& g)
  {
    cctbx_assert(x.size() == minimizer.n());
    cctbx_assert(g.size() == minimizer.n());
    return minimizer.run(x.begin(), f, g.begin());
  }

  bool py_minimizer_run_4(
    lbfgs::minimizer<double>& minimizer,
    af::shared<double> x,
    double f,
    const af::shared<double>& g,
    af::shared<double> diag)
  {
    cctbx_assert(x.size() == minimizer.n());
    cctbx_assert(g.size() == minimizer.n());
    cctbx_assert(diag.size() == minimizer.n());
    return minimizer.run(x.begin(), f, g.begin(), diag.begin());
  }

  bool py_traditional_convergence_test_call(
    const lbfgs::traditional_convergence_test<double>& is_converged,
    af::shared<double> x,
    const af::shared<double>& g)
  {
    cctbx_assert(x.size() == is_converged.n());
    cctbx_assert(g.size() == is_converged.n());
    return is_converged(x.begin(), g.begin());
  }

# include <cctbx/basic/from_bpl_import.h>

  void init_module(python::module_builder& this_module)
  {
    const std::string Revision = "$Revision$";
    this_module.add(ref(to_python(
        Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

    af::bpl::import_flex();

    python::class_builder<lbfgs::minimizer<double> >
    py_minimizer(this_module, "minimizer");

    python::class_builder<lbfgs::traditional_convergence_test<double> >
    py_traditional_convergence_test(
      this_module, "traditional_convergence_test");

    python::class_builder<lbfgs::drop_convergence_test<double> >
    py_drop_convergence_test(
      this_module, "drop_convergence_test");

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
    py_minimizer.def(py_minimizer_run_3, "run");
    py_minimizer.def(py_minimizer_run_4, "run");
    py_minimizer.def(&lbfgs::minimizer<double>::n, "n");
    py_minimizer.def(&lbfgs::minimizer<double>::m, "m");
    py_minimizer.def(&lbfgs::minimizer<double>::maxfev, "maxfev");
    py_minimizer.def(&lbfgs::minimizer<double>::gtol, "gtol");
    py_minimizer.def(&lbfgs::minimizer<double>::xtol, "xtol");
    py_minimizer.def(&lbfgs::minimizer<double>::stpmin, "stpmin");
    py_minimizer.def(&lbfgs::minimizer<double>::stpmax, "stpmax");
    py_minimizer.def(
      &lbfgs::minimizer<double>::requests_f_and_g, "requests_f_and_g");
    py_minimizer.def(
      &lbfgs::minimizer<double>::requests_diag, "requests_diag");
    py_minimizer.def(&lbfgs::minimizer<double>::iter, "iter");
    py_minimizer.def(&lbfgs::minimizer<double>::nfun, "nfun");
    py_minimizer.def(py_minimizer_euclidean_norm, "euclidean_norm");
    py_minimizer.def(&lbfgs::minimizer<double>::stp, "stp");

    py_traditional_convergence_test.def(constructor<>());
    py_traditional_convergence_test.def(constructor<std::size_t>());
    py_traditional_convergence_test.def(constructor<std::size_t, double>());
    py_traditional_convergence_test.def(
      &lbfgs::traditional_convergence_test<double>::n, "n");
    py_traditional_convergence_test.def(
      &lbfgs::traditional_convergence_test<double>::eps, "eps");
    py_traditional_convergence_test.def(
      py_traditional_convergence_test_call, "__call__");

    py_drop_convergence_test.def(constructor<>());
    py_drop_convergence_test.def(constructor<std::size_t>());
    py_drop_convergence_test.def(constructor<std::size_t, double>());
    py_drop_convergence_test.def(constructor<std::size_t, double, double>());
    py_drop_convergence_test.def(
      &lbfgs::drop_convergence_test<double>::p, "p");
    py_drop_convergence_test.def(
      &lbfgs::drop_convergence_test<double>::max_drop_eps, "max_drop_eps");
    py_drop_convergence_test.def(
      &lbfgs::drop_convergence_test<double>::iteration_coefficient,
                                            "iteration_coefficient");
    py_drop_convergence_test.def(
      &lbfgs::drop_convergence_test<double>::operator(), "__call__");
    py_drop_convergence_test.def(
      &lbfgs::drop_convergence_test<double>::objective_function_values,
                                            "objective_function_values");
    py_drop_convergence_test.def(
      &lbfgs::drop_convergence_test<double>::max_drop, "max_drop");
  }

}

BOOST_PYTHON_MODULE_INIT(lbfgs)
{
  boost::python::module_builder this_module("lbfgs");
  init_module(this_module);
}
