#include <scitbx/array_family/boost_python/flex_fwd.h>

#if defined(SCITBX_LBFGS_HAVE_LBFGS_FEM)
#include <scitbx/lbfgs_fem.hpp>
#endif

#include <scitbx/error.h>
#include <scitbx/lbfgs.h>
#include <scitbx/lbfgs/drop_convergence_test.h>
#include <scitbx/lbfgs/raw.h>
#include <scitbx/lbfgs/raw_reference.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/boost_python/utils.h>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>

namespace scitbx { namespace lbfgs { namespace ext {

  int
  fortran(
    int n,
    int m,
    af::ref<double> const& x,
    double f,
    af::const_ref<double> const& g,
    int diagco,
    af::ref<double> const& diag,
    af::tiny<int, 2> const& iprint,
    double eps,
    double xtol,
    af::ref<double> const& w,
    int iflag)
  {
    SCITBX_ASSERT(n > 0);
    SCITBX_ASSERT(m > 0);
    std::size_t n_ = static_cast<std::size_t>(n);
    std::size_t m_ = static_cast<std::size_t>(m);
    SCITBX_ASSERT(x.size() == n_);
    SCITBX_ASSERT(g.size() == n_);
    SCITBX_ASSERT(diagco == 0 || diagco == 1);
    SCITBX_ASSERT(diag.size() == n_);
    SCITBX_ASSERT(w.size() == n_*(2*m_+1)+2*m_);
#if defined(SCITBX_LBFGS_HAVE_LBFGS_FEM)
    static lbfgs_fem::common cmn;
    lbfgs_fem::blockdata_lb2(cmn);
    lbfgs_fem::lbfgs(
      cmn,
      n,
      m,
      x[0],
      f,
      g[0],
      diagco,
      diag[0],
      iprint[0],
      eps,
      xtol,
      w[0],
      iflag);
#else
    throw std::runtime_error("lbfgs_fem is not available.");
#endif
    return iflag;
  }

  int
  raw_reference(
    int n,
    int m,
    af::ref<double> const& x,
    double f,
    af::const_ref<double> const& g,
    int diagco,
    af::ref<double> const& diag,
    af::tiny<int, 2> const& iprint,
    double eps,
    double xtol,
    af::ref<double> const& w,
    int iflag)
  {
    SCITBX_ASSERT(n > 0);
    SCITBX_ASSERT(m > 0);
    std::size_t n_ = static_cast<std::size_t>(n);
    std::size_t m_ = static_cast<std::size_t>(m);
    SCITBX_ASSERT(x.size() == n_);
    SCITBX_ASSERT(g.size() == n_);
    SCITBX_ASSERT(diagco == 0 || diagco == 1);
    SCITBX_ASSERT(diag.size() == n_);
    SCITBX_ASSERT(w.size() == n_*(2*m_+1)+2*m_);
    using scitbx::lbfgs::raw_reference::const_ref1; // fully-qualified
    using scitbx::lbfgs::raw_reference::ref1;       // to work around
    scitbx::lbfgs::raw_reference::lbfgs(            // gcc 3.2 bug
      n,
      m,
      ref1<double>(x.begin(), n),
      f,
      const_ref1<double>(g.begin(), n),
      diagco,
      ref1<double>(diag.begin(), n),
      const_ref1<int>(iprint.begin(), 2),
      eps,
      xtol,
      ref1<double>(w.begin(), static_cast<int>(w.size())),
      iflag);
    return iflag;
  }

struct raw_lbfgs : boost::noncopyable {

  scitbx::lbfgs::raw::lbfgs lbfgs_obj;

  int
  operator()(
    int n,
    int m,
    af::ref<double> const& x,
    double f,
    af::const_ref<double> const& g,
    int diagco,
    af::ref<double> const& diag,
    af::tiny<int, 2> const& iprint,
    double eps,
    double xtol,
    af::ref<double> const& w,
    int iflag)
  {
    SCITBX_ASSERT(n > 0);
    SCITBX_ASSERT(m > 0);
    std::size_t n_ = static_cast<std::size_t>(n);
    std::size_t m_ = static_cast<std::size_t>(m);
    SCITBX_ASSERT(x.size() == n_);
    SCITBX_ASSERT(g.size() == n_);
    SCITBX_ASSERT(diagco >= 0);
    SCITBX_ASSERT(diagco <= 3);
    SCITBX_ASSERT(diag.size() == n_);
    SCITBX_ASSERT(w.size() == n_*(2*m_+1)+2*m_);
    using scitbx::lbfgs::raw::const_ref1; // fully-qualified
    using scitbx::lbfgs::raw::ref1;       // to work around
    lbfgs_obj(                            // gcc 3.2 bug
      n,
      m,
      ref1<double>(x.begin(), n),
      f,
      const_ref1<double>(g.begin(), n),
      diagco,
      ref1<double>(diag.begin(), n),
      const_ref1<int>(iprint.begin(), 2),
      eps,
      xtol,
      ref1<double>(w.begin(), static_cast<int>(w.size())),
      iflag);
    return iflag;
  }

  int
  nfun() const { return lbfgs_obj.nfun; }

  int
  iter() const { return lbfgs_obj.iter; }

  double
  stp() const { return lbfgs_obj.stp; }

  void
  set_stp(double value) { lbfgs_obj.stp = value; }

  af::shared<double>
  current_search_direction() const
  {
    return af::shared<double>(
      lbfgs_obj.current_search_direction.begin(),
      lbfgs_obj.current_search_direction.end());
  }
};

struct raw_lbfgs_wrappers
{
  typedef raw_lbfgs w_t;

  static
  void
  wrap()
  {
    using namespace boost::python;
    class_<w_t, boost::noncopyable>("raw_lbfgs", no_init)
      .def(init<>())
      .def("__call__", &w_t::operator(), (
        arg("n"),
        arg("m"),
        arg("x"),
        arg("f"),
        arg("g"),
        arg("diagco"),
        arg("diag"),
        arg("iprint"),
        arg("eps"),
        arg("xtol"),
        arg("w"),
        arg("iflag")))
      .def("nfun", &w_t::nfun)
      .def("iter", &w_t::iter)
      .def("stp", &w_t::stp)
      .def("set_stp", &w_t::set_stp, (arg("value")))
      .def("current_search_direction", &w_t::current_search_direction)
    ;
  }
};

  struct minimizer_wrappers
  {
    typedef minimizer<double> w_t;

    static bool
    run_4(
      w_t& minimizer,
      af::flex_double& x,
      double f,
      af::flex_double const& g,
      af::flex_double const& diag)
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
        .def(init<optional<std::size_t, double, double> >((
          arg("n_test_points")=5,
          arg("max_drop_eps")=1.e-5,
          arg("iteration_coefficient")=2)))
        .def("n_test_points", &w_t::n_test_points)
        .def("max_drop_eps", &w_t::max_drop_eps)
        .def("iteration_coefficient", &w_t::iteration_coefficient)
        .def("__call__", &w_t::operator())
        .def("objective_function_values", &w_t::objective_function_values, (
          arg("f")))
        .def("max_drop", &w_t::max_drop)
      ;
    }
  };

  void init_module()
  {
    using namespace boost::python;

    scope().attr("have_lbfgs_fem") =
#if defined(SCITBX_LBFGS_HAVE_LBFGS_FEM)
      true;
#else
      false;
#endif
    def("fortran", fortran, (
      arg("n"),
      arg("m"),
      arg("x"),
      arg("f"),
      arg("g"),
      arg("diagco"),
      arg("diag"),
      arg("iprint"),
      arg("eps"),
      arg("xtol"),
      arg("w"),
      arg("iflag")));
    def("raw_reference", raw_reference, (
      arg("n"),
      arg("m"),
      arg("x"),
      arg("f"),
      arg("g"),
      arg("diagco"),
      arg("diag"),
      arg("iprint"),
      arg("eps"),
      arg("xtol"),
      arg("w"),
      arg("iflag")));
    raw_lbfgs_wrappers::wrap();

    minimizer_wrappers::wrap();
    traditional_convergence_test_wrappers::wrap();
    drop_convergence_test_wrappers::wrap();
  }

}}} // namespace scitbx::lbfgs::ext

BOOST_PYTHON_MODULE(scitbx_lbfgs_ext)
{
  scitbx::lbfgs::ext::init_module();
}
