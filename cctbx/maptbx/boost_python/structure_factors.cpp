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
            arg("space_group"),
            arg("anomalous_flag"),
            arg("miller_indices"),
            arg("structure_factors"),
            arg("n_real"),
            arg("map_grid"),
            arg("conjugate_flag"),
            arg("treat_restricted")=true)))
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
            arg("unit_cell"),
            arg("space_group_type"),
            arg("anomalous_flag"),
            arg("d_min"),
            arg("complex_map"),
            arg("conjugate_flag"),
            arg("discard_indices_affected_by_aliasing")=false)))
        .def(init<
          bool,
          af::const_ref<miller::index<> > const&,
          af::const_ref<std::complex<double>,
            af::c_grid_padded<3> > const&,
          bool,
          optional<bool> >((
            arg("anomalous_flag"),
            arg("miller_indices"),
            arg("complex_map"),
            arg("conjugate_flag"),
            arg("allow_miller_indices_outside_map")=false)))
        .def(init<
          sgtbx::space_group const&,
          bool,
          af::const_ref<miller::index<> > const&,
          af::const_ref<std::complex<double>,
          af::c_grid_padded<3> > const&,
          bool>((
            arg("space_group"),
            arg("anomalous_flag"),
            arg("miller_indices"),
            arg("complex_map"),
            arg("conjugate_flag"))))
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
