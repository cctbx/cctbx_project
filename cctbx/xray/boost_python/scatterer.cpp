#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/scatterer.h>
#include <boost/python/class.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

namespace cctbx { namespace xray { namespace boost_python {

namespace {

  struct scatterer_wrappers
  {
    typedef scatterer<> w_t;
    typedef w_t::float_type flt_t;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      apply_symmetry_overloads, apply_symmetry, 2, 6)

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      typedef default_call_policies dcp;
      class_<w_t>("scatterer", no_init)
        .def(init<std::string const&,
                  fractional<flt_t> const&,
                  flt_t const&,
                  flt_t const&,
                  std::string const&,
                  flt_t const&,
                  flt_t const&>())
        .def(init<std::string const&,
                  fractional<flt_t> const&,
                  scitbx::sym_mat3<flt_t> const&,
                  flt_t const&,
                  std::string const&,
                  flt_t const&,
                  flt_t const&>())
        .add_property("label", make_getter(&w_t::label, rbv()),
                               make_setter(&w_t::label, dcp()))
        .add_property("scattering_type",
          make_getter(&w_t::scattering_type, rbv()),
          make_setter(&w_t::scattering_type, dcp()))
        .add_property("fp", make_getter(&w_t::fp, rbv()),
                            make_setter(&w_t::fp, dcp()))
        .add_property("fdp", make_getter(&w_t::fdp, rbv()),
                             make_setter(&w_t::fdp, dcp()))
        .add_property("site", make_getter(&w_t::site, rbv()),
                              make_setter(&w_t::site, dcp()))
        .add_property("occupancy", make_getter(&w_t::occupancy, rbv()),
                                   make_setter(&w_t::occupancy, dcp()))
        .add_property("anisotropic_flag",
          make_getter(&w_t::anisotropic_flag, rbv()),
          make_setter(&w_t::anisotropic_flag, dcp()))
        .add_property("u_iso", make_getter(&w_t::u_iso, rbv()),
                               make_setter(&w_t::u_iso, dcp()))
        .add_property("u_star", make_getter(&w_t::u_star, rbv()),
                                make_setter(&w_t::u_star, dcp()))
        .def("apply_symmetry", &w_t::apply_symmetry,apply_symmetry_overloads())
        .def("update_weight", &w_t::update_weight)
        .def("multiplicity", &w_t::multiplicity)
        .def("weight_without_occupancy", &w_t::weight_without_occupancy)
        .def("weight", &w_t::weight)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_scatterer()
  {
    scatterer_wrappers::wrap();
  }

}}} // namespace cctbx::xray::boost_python
