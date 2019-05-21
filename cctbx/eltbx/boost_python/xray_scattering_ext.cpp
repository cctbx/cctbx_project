#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <cctbx/eltbx/xray_scattering.h>
#include <cctbx/eltbx/xray_scattering/n_gaussian.h>
#include <scitbx/boost_python/iterator_wrappers.h>
#include <scitbx/boost_python/array_as_list.h>
#include <boost_adaptbx/optional_conversions.h>
#include <boost/optional.hpp>

namespace cctbx { namespace eltbx { namespace xray_scattering {
namespace boost_python {

namespace {

  boost::python::object
  standard_labels_list()
  {
    boost::python::object result = scitbx::boost_python::array_as_list(
      standard_labels, sizeof(standard_labels) / sizeof(char*) - 1);
    return result;
  }

  template<class Derived>
  struct isotropic_form_factor_mixin_wrapper
  {
    typedef isotropic_form_factor_mixin<Derived> w_t;

    static void wrap(char const *name)
    {
      using namespace boost::python;
      class_<w_t>(name, no_init)
        .def("at_stol_sq", &w_t::at_stol_sq)
        .def("at_stol", &w_t::at_stol)
        .def("at_d_star_sq",
          (double(w_t::*)(double) const)
            &w_t::at_d_star_sq)
        .def("at_d_star_sq",
          (af::shared<double>(w_t::*)(af::const_ref<double> const&) const)
            &w_t::at_d_star_sq)
        .def("at_d_star", &w_t::at_d_star)
        ;
    }
  };

  struct gaussian_wrappers
  {
    typedef gaussian w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t,
            bases<w_t::base_t, isotropic_form_factor_mixin<w_t> >
            >("gaussian", no_init)
        .def(init<w_t::base_t const&>())
        .def(init<double, optional<bool> >())
        .def(init<
          af::small<double, gaussian::max_n_terms> const&,
          af::small<double, gaussian::max_n_terms> const&,
          optional<double, bool> >())
        .def(init<
          af::const_ref<double> const&,
          optional<double, bool> >())
      ;
      boost_adaptbx::optional_conversions::to_and_from_python<
        boost::optional<gaussian> >();
    }
  };

  struct n_gaussian_table_entry_wrappers
  {
    typedef n_gaussian::table_entry w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_internal_reference<> rir;
      class_<w_t>("n_gaussian_table_entry", no_init)
        .def(init<std::size_t, std::size_t>())
        .def(init<std::string const&, std::size_t>())
        .def(init<std::size_t, double, double>())
        .def(init<std::string const&, double, double>())
        .def("label", &w_t::label, ccr())
        .def("gaussian", &w_t::gaussian, rir())
        .def("max_stol", &w_t::max_stol)
        .def("d_min", &w_t::d_min)
        .def("max_relative_error", &w_t::max_relative_error)
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
        .def("fetch", &w_t::fetch)
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

    def("standard_labels_list", standard_labels_list);
    def("get_standard_label", get_standard_label, (
      arg("label"),
      arg("exact")=false,
      arg("optional")=false));

    isotropic_form_factor_mixin_wrapper<gaussian>::wrap(
      "gaussian_mixin");
    gaussian_wrappers::wrap();

    def("n_gaussian_table_size", n_gaussian::table_size);
    def("n_gaussian_table_index", n_gaussian::table_index);
    n_gaussian_table_entry_wrappers::wrap();

    it1992_wrappers::wrap();
    scitbx::boost_python::iterator_wrappers<
      it1992, it1992_iterator>::wrap("it1992_iterator");

    wk1995_wrappers::wrap();
    scitbx::boost_python::iterator_wrappers<
      wk1995, wk1995_iterator>::wrap("wk1995_iterator");
  }

} // namespace <anonymous>
}}}} // namespace cctbx::eltbx::xray_scattering::boost_python

BOOST_PYTHON_MODULE(cctbx_eltbx_xray_scattering_ext)
{
  cctbx::eltbx::xray_scattering::boost_python::init_module();
}
