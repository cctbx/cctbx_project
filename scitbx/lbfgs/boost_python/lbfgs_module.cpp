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
#include <scitbx/array_family/boost_python/utils.h>
#include <scitbx/boost_python/utils.h>
#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>

namespace scitbx { namespace lbfgs { namespace {

  struct minimizer_wrappers
  {
    typedef minimizer<double> w_t;

    static bool
    run_4(
      w_t& minimizer,
      af::flex_double& x,
      double f,
      af::flex_double const& g,
      af::flex_double& diag)
    {
      using namespace scitbx::af::boost_python;
      SCITBX_ASSERT(flex_as_base_array(x).size() == minimizer.n());
      SCITBX_ASSERT(flex_as_base_array(g).size() == minimizer.n());
      SCITBX_ASSERT(flex_as_base_array(diag).size() == minimizer.n());
      return minimizer.run(x.begin(), f, g.begin(), diag.begin());
    }

    static bool
    run_3(
      w_t& minimizer,
      af::flex_double& x,
      double f,
      af::flex_double const& g)
    {
      using namespace scitbx::af::boost_python;
      SCITBX_ASSERT(flex_as_base_array(x).size() == minimizer.n());
      SCITBX_ASSERT(flex_as_base_array(g).size() == minimizer.n());
      return minimizer.run(x.begin(), f, g.begin());
    }

    static double
    euclidean_norm(
      w_t const& minimizer,
      af::flex_double const& a)
    {
      using namespace scitbx::af::boost_python;
      SCITBX_ASSERT(flex_as_base_array(a).size() == minimizer.n());
      return minimizer.euclidean_norm(a.begin());
    }

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("minimizer", no_init)
        .def(init<optional<std::size_t, std::size_t, std::size_t,
          double, double, double, double> >())
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
      ;
    }
  };

  struct traditional_convergence_test_wrappers
  {
    typedef traditional_convergence_test<double> w_t;

    static bool
    call(
      w_t const& is_converged,
      af::flex_double const& x,
      af::flex_double const& g)
    {
      using namespace scitbx::af::boost_python;
      SCITBX_ASSERT(flex_as_base_array(x).size() == is_converged.n());
      SCITBX_ASSERT(flex_as_base_array(g).size() == is_converged.n());
      return is_converged(x.begin(), g.begin());
    }

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("traditional_convergence_test")
        .def(init<std::size_t, optional<double> >())
        .def("n", &w_t::n)
        .def("eps", &w_t::eps)
        .def("__call__", call)
      ;
    }
  };

  struct drop_convergence_test_wrappers
  {
    typedef drop_convergence_test<double> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("drop_convergence_test", no_init)
        .def(init<optional<std::size_t, double, double> >())
        .def("p", &w_t::p)
        .def("max_drop_eps", &w_t::max_drop_eps)
        .def("iteration_coefficient", &w_t::iteration_coefficient)
        .def("__call__", &w_t::operator())
        .def("objective_function_values", &w_t::objective_function_values)
        .def("max_drop", &w_t::max_drop)
      ;
    }
  };

  void init_module()
  {
    using namespace boost::python;

    scope().attr("__version__") = scitbx::boost_python::cvs_revision(
      "$Revision$");

    scitbx::boost_python::import_module(
      "scitbx_boost.array_family.flex_scitbx_ext");

    minimizer_wrappers::wrap();
    traditional_convergence_test_wrappers::wrap();
    drop_convergence_test_wrappers::wrap();
  }

}}} // namespace scitbx::lbfgs::<anonymous>

BOOST_PYTHON_MODULE(lbfgs_ext)
{
  scitbx::lbfgs::init_module();
}
