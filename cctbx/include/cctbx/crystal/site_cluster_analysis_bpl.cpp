#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/overloads.hpp>
#include <cctbx/crystal/site_cluster_analysis.h>

namespace cctbx { namespace crystal {
namespace {

  struct site_cluster_analysis_wrappers
  {
    typedef site_cluster_analysis<> w_t;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      process_sites_frac_3_overloads, process_sites_frac, 2, 3)
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      process_sites_frac_2_overloads, process_sites_frac, 1, 2)
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      process_sites_cart_3_overloads, process_sites_cart, 2, 3)
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      process_sites_cart_2_overloads, process_sites_cart, 1, 2)

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("site_cluster_analysis", no_init)
        .def(init<
          sgtbx::space_group const&,
          direct_space_asu::float_asu<> const&,
          double const&,
          optional<
            unsigned,
            double const&,
            double const&> >((
              arg_("space_group"),
              arg_("asu"),
              arg_("distance_cutoff"),
              arg_("estimated_reduction_factor")=4,
              arg_("asu_mappings_buffer_thickness")=-1,
              arg_("cubicle_epsilon")=-1)))
        .def_readwrite("estimated_reduction_factor",
                  &w_t::estimated_reduction_factor)
        .def_readwrite("min_distance_sym_equiv",
                  &w_t::min_distance_sym_equiv)
        .def_readwrite("assert_min_distance_sym_equiv",
                  &w_t::assert_min_distance_sym_equiv)
        .def("asu_mappings", &w_t::asu_mappings)
        .def("process_site_frac",
          (bool(w_t::*)(fractional<> const&, sgtbx::site_symmetry_ops const&))
            &w_t::process_site_frac, (
              arg_("original_site"), arg_("site_symmetry_ops")))
        .def("process_site_frac",
          (bool(w_t::*)(fractional<> const&))
            &w_t::process_site_frac, (
              arg_("original_site")))
        .def("process_sites_frac",
          (af::shared<std::size_t>(w_t::*)(
            af::const_ref<scitbx::vec3<double> > const&,
            sgtbx::site_symmetry_table const&,
            std::size_t))
              &w_t::process_sites_frac, process_sites_frac_3_overloads((
                arg_("original_sites"),
                arg_("site_symmetry_table"),
                arg_("max_clusters")=0)))
        .def("process_sites_frac",
          (af::shared<std::size_t>(w_t::*)(
            af::const_ref<scitbx::vec3<double> > const&,
            std::size_t))
              &w_t::process_sites_frac, process_sites_frac_2_overloads((
                arg_("original_sites"),
                arg_("max_clusters")=0)))
        .def("process_sites_cart",
          (af::shared<std::size_t>(w_t::*)(
            af::const_ref<scitbx::vec3<double> > const&,
            sgtbx::site_symmetry_table const&,
            std::size_t))
              &w_t::process_sites_cart, process_sites_cart_3_overloads((
                arg_("original_sites"),
                arg_("site_symmetry_table"),
                arg_("max_clusters")=0)))
        .def("process_sites_cart",
          (af::shared<std::size_t>(w_t::*)(
            af::const_ref<scitbx::vec3<double> > const&,
            std::size_t))
              &w_t::process_sites_cart, process_sites_cart_2_overloads((
                arg_("original_sites"),
                arg_("max_clusters")=0)))
      ;
    }
  };

  void
  wrap_all()
  {
    site_cluster_analysis_wrappers::wrap();
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_site_cluster_analysis() { wrap_all(); }

}}} // namespace cctbx::crystal::boost_python
