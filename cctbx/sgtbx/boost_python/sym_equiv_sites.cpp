#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/sgtbx/sym_equiv_sites.h>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>

namespace cctbx { namespace sgtbx { namespace boost_python {

namespace {

  struct sym_equiv_sites_wrappers
  {
    typedef sym_equiv_sites<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_internal_reference<> rir;
      class_<w_t>("sym_equiv_sites", no_init)
        .def(init<site_symmetry const&>())
        .def(init<wyckoff::mapping const&>(
          (arg_("wyckoff_mapping"))))
        .def(init<sgtbx::space_group const&,
                  fractional<> const&,
                  optional<uctbx::unit_cell const&> >(
          (arg_("space_group"), arg_("original_site"), arg_("unit_cell"))))
        .def(init<uctbx::unit_cell const&,
                  sgtbx::space_group const&,
                  fractional<> const&,
                  rt_mx const&>(
          (arg_("unit_cell"), arg_("space_group"), arg_("original_site"),
           arg_("special_op"))))
        .def(init<uctbx::unit_cell const&,
                  sgtbx::space_group const&,
                  fractional<> const&,
                  optional<double, double> >(
          (arg_("unit_cell"), arg_("space_group"), arg_("original_site"),
           arg_("minimum_distance"), arg_("tolerance"))))
        .def("unit_cell", &w_t::unit_cell, rir())
        .def("space_group", &w_t::space_group, rir())
        .def("original_site", &w_t::original_site, ccr())
        .def("special_op", &w_t::special_op, ccr())
        .def("max_accepted_tolerance", &w_t::max_accepted_tolerance)
        .def("coordinates", &w_t::coordinates)
        .def("sym_op", &w_t::sym_op, (arg_("i_coor")))
        .def("sym_op_indices", &w_t::sym_op_indices)
        .def("is_special_position", &w_t::is_special_position)
      ;
    }
  };

  struct min_sym_equiv_distance_info_wrappers
  {
    typedef min_sym_equiv_distance_info<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("min_sym_equiv_distance_info", no_init)
        .def(init<sym_equiv_sites<> const&,
                  fractional<> const&,
                  optional<af::tiny<bool, 3> const&> >(
          (arg_("reference_sites"),
           arg_("other"),
           arg_("principal_continuous_allowed_origin_shift_flags")
             =no_continuous_allowed_shifts)))
        .def(init<sym_equiv_sites<> const&,
                  af::const_ref<scitbx::vec3<double> > const&,
                  optional<af::tiny<bool, 3> const&> >(
          (arg_("reference_sites"),
           arg_("others"),
           arg_("principal_continuous_allowed_origin_shift_flags")
             =no_continuous_allowed_shifts)))
        .def("i_other", &w_t::i_other)
        .def("sym_op", &w_t::sym_op, ccr())
        .def("continuous_shifts", &w_t::continuous_shifts, ccr())
        .def("diff", &w_t::diff, ccr())
        .def("dist", &w_t::dist)
        .def("apply", &w_t::apply, (arg_("sites_frac")))
      ;
    }
  };

} // namespace <anoymous>

  void wrap_sym_equiv_sites()
  {
    sym_equiv_sites_wrappers::wrap();
    min_sym_equiv_distance_info_wrappers::wrap();
  }

}}} // namespace cctbx::sgtbx::boost_python
