/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Ported from cctbx/lbfgs/lbfgsmodule.cpp (rwgk)
     2002 Mar: Created (R.W. Grosse-Kunstleve)
 */

#include <scitbx/error.h>
#include <scitbx/lbfgs.h>
#include <scitbx/lbfgs/drop_convergence_test.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/boost_python/utils.h>
#include <boost/python/module.hpp>
#include <boost/python/class.hpp>

namespace {

  using namespace scitbx;

  struct lbfgs_minimizer_wrappers
  {
    typedef lbfgs::minimizer<double> w_t;

    static bool
    run_4(
      w_t& minimizer,
      af::flex_double& x,
      double f,
      af::flex_double const& g,
      af::flex_double& diag)
    {
      SCITBX_ASSERT(x.size() == minimizer.n());
      SCITBX_ASSERT(g.size() == minimizer.n());
      SCITBX_ASSERT(diag.size() == minimizer.n());
      return minimizer.run(x.begin(), f, g.begin(), diag.begin());
    }

    static bool
    run_3(
      w_t& minimizer,
      af::flex_double& x,
      double f,
      af::flex_double const& g)
    {
      SCITBX_ASSERT(x.size() == minimizer.n());
      SCITBX_ASSERT(g.size() == minimizer.n());
      return minimizer.run(x.begin(), f, g.begin());
    }

    static double
    euclidean_norm(
      w_t const& minimizer,
      af::flex_double const& a)
    {
      SCITBX_ASSERT(a.size() == minimizer.n());
      return minimizer.euclidean_norm(a.begin());
    }

    static void
    def_all(boost::python::module& this_module)
    {
      using namespace boost::python;
      this_module.add(
        class_<w_t>("minimizer",
          args<>())
          .def_init(args<std::size_t>())
          .def_init(args<std::size_t, std::size_t>())
          .def_init(args<std::size_t, std::size_t, std::size_t>())
          .def_init(args<std::size_t, std::size_t, std::size_t,
            double>())
          .def_init(args<std::size_t, std::size_t, std::size_t,
            double, double>())
          .def_init(args<std::size_t, std::size_t, std::size_t,
            double, double, double>())
          .def_init(args<std::size_t, std::size_t, std::size_t,
            double, double, double, double>())
          .def("run", run_4)
          .def("run", run_3)
          .def("n", &w_t::n)
          .def("m", &w_t::m)
          .def("maxfev", &w_t::maxfev)
          .def("gtol", &w_t::gtol)
          .def("xtol", &w_t::xtol)
          .def("stpmin", &w_t::stpmin)
          .def("stpmax", &w_t::stpmax)
          .def("requests_f_and_g", &w_t::requests_f_and_g)
          .def("requests_diag", &w_t::requests_diag)
          .def("iter", &w_t::iter)
          .def("nfun", &w_t::nfun)
          .def("euclidean_norm", euclidean_norm)
          .def("stp", &w_t::stp)
      );
    }
  };

  struct traditional_convergence_test_wrappers
  {
    typedef lbfgs::traditional_convergence_test<double> w_t;

    static bool
    call(
      w_t const& is_converged,
      af::flex_double const& x,
      af::flex_double const& g)
    {
      SCITBX_ASSERT(x.size() == is_converged.n());
      SCITBX_ASSERT(g.size() == is_converged.n());
      return is_converged(x.begin(), g.begin());
    }

    static void
    def_all(boost::python::module& this_module)
    {
      using namespace boost::python;
      this_module.add(
        class_<w_t>("traditional_convergence_test",
          args<>())
          .def_init(args<std::size_t>())
          .def_init(args<std::size_t, double>())
          .def("n", &w_t::n)
          .def("eps", &w_t::eps)
          .def("__call__", call)
      );
    }
  };

  struct drop_convergence_test_wrappers
  {
    typedef lbfgs::drop_convergence_test<double> w_t;

    static void
    def_all(boost::python::module& this_module)
    {
      using namespace boost::python;
      this_module.add(
        class_<w_t>("drop_convergence_test",
          args<>())
          .def_init(args<std::size_t>())
          .def_init(args<std::size_t, double>())
          .def_init(args<std::size_t, double, double>())
          .def("p", &w_t::p)
          .def("max_drop_eps", &w_t::max_drop_eps)
          .def("iteration_coefficient", &w_t::iteration_coefficient)
          .def("__call__", &w_t::operator())
          .def("objective_function_values", &w_t::objective_function_values)
          .def("max_drop", &w_t::max_drop)
      );
    }
  };

  void init_module(boost::python::module& this_module)
  {
    this_module
      .setattr("__version__",
        scitbx::boost_python::cvs_revision("$Revision$"))
    ;

    scitbx::boost_python::import_module("scitbx_boost.array_family.flex");

    lbfgs_minimizer_wrappers::def_all(this_module);
    traditional_convergence_test_wrappers::def_all(this_module);
    drop_convergence_test_wrappers::def_all(this_module);
  }

}

BOOST_PYTHON_MODULE_INIT(lbfgs)
{
  boost::python::module this_module("lbfgs");
  init_module(this_module);
}
