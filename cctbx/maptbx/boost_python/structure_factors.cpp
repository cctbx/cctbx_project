#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/maptbx/structure_factors.h>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>

namespace cctbx { namespace maptbx { namespace structure_factors {

namespace {

  struct to_map_wrappers
  {
    typedef to_map<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("structure_factors_to_map", no_init)
        .def(init<
          sgtbx::space_group const&,
          bool,
          af::const_ref<miller::index<> > const&,
          af::const_ref<std::complex<double> > const&,
          af::int3 const&,
          af::flex_grid<> const&,
          bool,
          optional<bool> >((
            arg_("space_group"),
            arg_("anomalous_flag"),
            arg_("miller_indices"),
            arg_("structure_factors"),
            arg_("n_real"),
            arg_("map_grid"),
            arg_("conjugate_flag"),
            arg_("treat_restricted")=true)))
        .def("complex_map", &w_t::complex_map)
      ;
    }
  };

  struct from_map_wrappers
  {
    typedef from_map<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("structure_factors_from_map", no_init)
        .def(init<
          uctbx::unit_cell const&,
          sgtbx::space_group_type const&,
          bool,
          double,
          af::const_ref<std::complex<double>,
            af::c_grid_padded<3> > const&,
          bool,
          optional<bool> >((
            arg_("unit_cell"),
            arg_("space_group_type"),
            arg_("anomalous_flag"),
            arg_("d_min"),
            arg_("complex_map"),
            arg_("conjugate_flag"),
            arg_("discard_indices_affected_by_aliasing")=false)))
        .def(init<
          bool,
          af::const_ref<miller::index<> > const&,
          af::const_ref<std::complex<double>,
            af::c_grid_padded<3> > const&,
          bool,
          optional<bool> >((
            arg_("anomalous_flag"),
            arg_("miller_indices"),
            arg_("complex_map"),
            arg_("conjugate_flag"),
            arg_("allow_miller_indices_outside_map")=false)))
        .def(init<
          sgtbx::space_group const&,
          bool,
          af::const_ref<miller::index<> > const&,
          af::const_ref<std::complex<double>,
          af::c_grid_padded<3> > const&,
          bool>((
            arg_("space_group"),
            arg_("anomalous_flag"),
            arg_("miller_indices"),
            arg_("complex_map"),
            arg_("conjugate_flag"))))
        .def("miller_indices", &w_t::miller_indices)
        .def("data", &w_t::data)
        .def("n_indices_affected_by_aliasing",
             &w_t::n_indices_affected_by_aliasing)
        .def("outside_map", &w_t::outside_map)
      ;
    }
  };

}} // namespace structure_factors::<anoymous>

namespace boost_python {

  void wrap_structure_factors()
  {
    structure_factors::to_map_wrappers::wrap();
    structure_factors::from_map_wrappers::wrap();
  }

}}} // namespace cctbx::maptbx::boost_python
