#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <scitbx/boost_python/is_polymorphic_workaround.h>
#include <cctbx/crystal/incremental_pairs.h>
#include <cctbx/crystal/workarounds_bpl.h>

SCITBX_BOOST_IS_POLYMORPHIC_WORKAROUND(
  cctbx::crystal::incremental_pairs<>)

namespace cctbx { namespace crystal {
namespace {

  struct incremental_pairs_wrappers
  {
    typedef incremental_pairs<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("incremental_pairs", no_init)
        .def(init<
          sgtbx::space_group const&,
          direct_space_asu::float_asu<> const&,
          double const&,
          double const&,
          double const&>((
            arg("space_group"),
            arg("asu"),
            arg("distance_cutoff"),
            arg("asu_mappings_buffer_thickness")=-1,
            arg("cubicle_epsilon")=-1)))
        .def_readwrite("min_distance_sym_equiv",
                  &w_t::min_distance_sym_equiv)
        .def_readwrite("assert_min_distance_sym_equiv",
                  &w_t::assert_min_distance_sym_equiv)
        .def("asu_mappings", &w_t::asu_mappings)
        .def("pair_asu_table", &w_t::pair_asu_table)
        .def("process_site_frac",
          (void(w_t::*)(fractional<> const&, sgtbx::site_symmetry_ops const&))
            &w_t::process_site_frac, (
              arg("original_site"), arg("site_symmetry_ops")))
        .def("process_site_frac",
          (void(w_t::*)(fractional<> const&))
            &w_t::process_site_frac, (
              arg("original_site")))
        .def("process_sites_frac",
          (void(w_t::*)(
            af::const_ref<scitbx::vec3<double> > const&,
            sgtbx::site_symmetry_table const&))
              &w_t::process_sites_frac, (
                arg("original_sites"), arg("site_symmetry_table")))
        .def("process_sites_frac",
          (void(w_t::*)(af::const_ref<scitbx::vec3<double> > const&))
            &w_t::process_sites_frac, (
              arg("original_sites")))
        .def("process_sites_cart",
          (void(w_t::*)(
            af::const_ref<scitbx::vec3<double> > const&,
            sgtbx::site_symmetry_table const&))
              &w_t::process_sites_cart, (
                arg("original_sites"), arg("site_symmetry_table")))
        .def("process_sites_cart",
          (void(w_t::*)(af::const_ref<scitbx::vec3<double> > const&))
            &w_t::process_sites_cart, (
              arg("original_sites")))
        .def("cubicle_size_counts", &w_t::cubicle_size_counts)
      ;
    }
  };

  void
  wrap_all()
  {
    incremental_pairs_wrappers::wrap();
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_incremental_pairs() { wrap_all(); }

}}} // namespace cctbx::crystal::boost_python
