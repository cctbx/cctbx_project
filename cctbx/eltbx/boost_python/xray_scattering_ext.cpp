#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <cctbx/eltbx/gaussian_fit.h>
#include <scitbx/boost_python/iterator_wrappers.h>

namespace cctbx { namespace eltbx { namespace xray_scattering {
namespace boost_python {

namespace {

  struct gaussian_wrappers : boost::python::pickle_suite
  {
    typedef gaussian w_t;

    static
    boost::python::tuple
    getinitargs(w_t const& w)
    {
      return boost::python::make_tuple(w.a(), w.b(), w.c());
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("gaussian", no_init)
        .def(init<float>())
        .def(init<af::small<float, w_t::max_n_ab> const&,
                  af::small<float, w_t::max_n_ab> const&,
                  optional<float> >())
        .def("n_ab", &w_t::n_ab)
        .def("a", (af::small<float, w_t::max_n_ab> const&(w_t::*)()const)
          &w_t::a, ccr())
        .def("b", (af::small<float, w_t::max_n_ab> const&(w_t::*)()const)
          &w_t::b, ccr())
        .def("c", &w_t::c)
        .def("all_zero", &w_t::all_zero)
        .def("at_stol_sq", &w_t::at_stol_sq)
        .def("at_stol", &w_t::at_stol)
        .def("at_d_star_sq", (double(w_t::*)(double) const) &w_t::at_d_star_sq)
        .def("at_d_star_sq",
          (af::shared<double>(w_t::*)(af::const_ref<double> const&) const)
          &w_t::at_d_star_sq)
        .def("gradient_at_d_star", &w_t::gradient_at_d_star)
        .def("integral_at_d_star", &w_t::integral_at_d_star)
        .def("gradients_at_stol", &w_t::gradients_at_stol)
        .def_pickle(gaussian_wrappers())
      ;
    }
  };

  struct gaussian_fit_wrappers
  {
    typedef gaussian_fit w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t, bases<gaussian> >("gaussian_fit", no_init)
        .def(init<af::shared<double> const&,
                  af::shared<double> const&,
                  af::shared<double> const&,
                  gaussian const&>())
        .def(init<af::shared<double> const&,
                  gaussian const&,
                  af::shared<double> const&,
                  gaussian const&>())
        .def("stols", &w_t::stols)
        .def("target_values", &w_t::target_values)
        .def("sigmas", &w_t::sigmas)
        .def("fitted_values", &w_t::fitted_values)
        .def("differences", &w_t::differences)
        .def("apply_shifts", &w_t::apply_shifts)
        .def("target_function", &w_t::target_function)
        .def("gradients", &w_t::gradients)
      ;
    }
  };

  template <std::size_t N>
  struct base_wrappers
  {
    typedef base<N> w_t;

    static void
    wrap(const char* python_name)
    {
      using namespace boost::python;
      class_<w_t>(python_name, no_init)
        .def("table", &w_t::table)
        .def("label", &w_t::label)
        .def("a", (af::small<float, gaussian::max_n_ab>(w_t::*)()const)
          &w_t::a)
        .def("b", (af::small<float, gaussian::max_n_ab>(w_t::*)()const)
          &w_t::b)
        .def("c", &w_t::c)
        .def("fetch", &w_t::fetch)
        .def("at_stol_sq", &w_t::at_stol_sq)
        .def("at_stol", &w_t::at_stol)
        .def("at_d_star_sq", (double(w_t::*)(double) const) &w_t::at_d_star_sq)
        .def("at_d_star_sq",
          (af::shared<double>(w_t::*)(af::const_ref<double> const&) const)
          &w_t::at_d_star_sq)
        .def("gradient_at_d_star", &w_t::gradient_at_d_star)
        .def("integral_at_d_star", &w_t::integral_at_d_star)
        .def("gradients_at_stol", &w_t::gradients_at_stol)
      ;
    }
  };

  struct it1992_wrappers
  {
    typedef it1992 w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      base_wrappers<4>::wrap("base_4");
      class_<w_t, bases<base<4> > >("it1992", no_init)
        .def(init<std::string const&, optional<bool> >())
      ;
    }
  };

  struct wk1995_wrappers
  {
    typedef wk1995 w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      base_wrappers<5>::wrap("base_5");
      class_<w_t, bases<base<5> > >("wk1995", no_init)
        .def(init<std::string const&, optional<bool> >())
      ;
    }
  };

  void init_module()
  {
    using namespace boost::python;

    gaussian_wrappers::wrap();
    gaussian_fit_wrappers::wrap();

    it1992_wrappers::wrap();
    scitbx::boost_python::iterator_wrappers<
      it1992, it1992_iterator>::wrap("it1992_iterator");

    wk1995_wrappers::wrap();
    scitbx::boost_python::iterator_wrappers<
      wk1995, wk1995_iterator>::wrap("wk1995_iterator");
  }

} // namespace <anonymous>
}}}} // namespace cctbx::eltbx::xray_scattering::boost_python

BOOST_PYTHON_MODULE(xray_scattering_ext)
{
  cctbx::eltbx::xray_scattering::boost_python::init_module();
}
