#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/optional.hpp>
#include <cctbx/sgtbx/rt_mx.h>
#include <cctbx/sgtbx/rt_mx_hash.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/stl/vector_wrapper.h>
#include <scitbx/boost_python/container_conversions.h>
#include <boost_adaptbx/optional_conversions.h>
#include <boost_adaptbx/hash.h>
#include <tbxx/optional_copy.hpp>

namespace cctbx { namespace sgtbx { namespace boost_python {

namespace {

  struct rt_mx_wrappers : boost_adaptbx::py_hashable<rt_mx>
  {
    typedef rt_mx w_t;

    static rt_mx*
    unpickle_init(af::const_ref<int> const& values)
    {
      CCTBX_ASSERT(values.size() == 14);
      rt_mx result(values[12], values[13]);
      for(std::size_t i=0;i<9;i++) {
        result.r().num()[i] = values[i];
      }
      for(std::size_t i=0;i<3;i++) {
        result.t().num()[i] = values[i+9];
      }
      return new rt_mx(result);
    }

    static std::string
    str(w_t const& o) { return o.as_xyz(); }

    static scitbx::vec3<double>
    mul_double(
      w_t const& o, scitbx::vec3<double> const& rhs) { return o * rhs; }

    static vec3_rat
    mul_rat(
      w_t const& o, vec3_rat const& rhs) { return o * rhs; }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_internal_reference<> rir;
      class_<w_t>("rt_mx", no_init)
        .enable_pickling()
        .def("__init__", boost::python::make_constructor(unpickle_init))
        .def(init<optional<int, int> >((
          arg("r_den")=1,
          arg("t_den")=sg_t_den)))
        .def(init<rot_mx const&, tr_vec const&>((arg("r"), arg("t"))))
        .def(init<rot_mx const&, optional<int> >((
          arg("r"),
          arg("t_den")=sg_t_den)))
        .def(init<tr_vec const&, optional<int> >((
          arg("t"),
          arg("r_den")=1)))
        .def(init<parse_string&, optional<const char*, int, int> >((
          arg("symbol"),
          arg("stop_chars")="",
          arg("r_den")=1,
          arg("t_den")=sg_t_den)))
        .def(init<std::string const&, optional<const char*, int, int> >((
          arg("symbol"),
          arg("stop_chars")="",
          arg("r_den")=1,
          arg("t_den")=sg_t_den)))
        .def(init<
          scitbx::mat3<double> const&,
          scitbx::vec3<double> const&,
          optional<int, int> >((
            arg("r"),
            arg("t"),
            arg("r_den")=1,
            arg("t_den")=sg_t_den)))
        .def("r", (rot_mx const&(w_t::*)() const) &w_t::r, rir())
        .def("t", (tr_vec const&(w_t::*)() const) &w_t::t, rir())
        .def("__eq__", &w_t::operator==)
        .def("__ne__", &w_t::operator!=)
        .def("__hash__", py_hash)
        .def("is_valid", &w_t::is_valid)
        .def("unit_mx", &w_t::unit_mx)
        .def("is_unit_mx", &w_t::is_unit_mx)
        .def("as_xyz", &w_t::as_xyz, (
          arg("decimal")=false,
          arg("t_first")=false,
          arg("symbol_letters")="xyz",
          arg("separator")=","))
        .def("__str__", str)
        .def("as_int_array", &w_t::as_int_array)
        .def("as_double_array", &w_t::as_double_array)
        .def("new_denominators",
          (rt_mx(w_t::*)(int, int) const) &w_t::new_denominators, (
            arg("r_den"),
            arg("t_den")=0))
        .def("new_denominators",
          (rt_mx(w_t::*)(rt_mx const&) const) &w_t::new_denominators, (
            arg("other")))
        .def("mod_positive", &w_t::mod_positive)
        .def("mod_short", &w_t::mod_short)
        .def("inverse", &w_t::inverse)
        .def("refine_gridding", (sg_vec3(w_t::*)(sg_vec3 const&) const)
          &w_t::refine_gridding, (arg("grid")))
        .def("cancel", &w_t::cancel)
        .def("inverse_cancel", &w_t::inverse_cancel)
        .def("multiply", &w_t::multiply, (arg("rhs")))
        .def("__mul__", mul_double)
        .def("__mul__", mul_rat)
        .def("__call__",
             (scitbx::vec3<double>
              (w_t::*)(scitbx::vec3<double> const &) const)
             &w_t::operator())
        .def("__call__",
             (scitbx::sym_mat3<double>
              (w_t::*)(scitbx::sym_mat3<double> const &) const)
             &w_t::operator())
        .def("__add__", (rt_mx(w_t::*)(sg_vec3 const&) const)&w_t::operator+)
        .def("__add__", (rt_mx(w_t::*)(tr_vec const&) const)&w_t::operator+)
        .def("unit_shifts_minimum_distance", (
          scitbx::vec3<int>(w_t::*)(
            fractional<double> const&,
            fractional<double> const&) const)
              &w_t::unit_shifts_minimum_distance, (
                arg("site_frac_1"), arg("site_frac_2")))
        .def("add_unit_shifts_minimum_distance", (
          rt_mx(w_t::*)(
            fractional<double> const&,
            fractional<double> const&) const)
              &w_t::add_unit_shifts_minimum_distance, (
                arg("site_frac_1"), arg("site_frac_2")))
      ;

      scitbx::stl::boost_python::vector_wrapper<rt_mx>::wrap(
        "stl_vector_rt_mx");
      {
        using namespace scitbx::boost_python::container_conversions;
        tuple_mapping_variable_capacity<af::shared<rt_mx> >();
      }
      {
        using boost_adaptbx::optional_conversions::to_and_from_python;
        // used in cctbx/geometry_restraints
        to_and_from_python<
          tbxx::optional_container<
            af::shared<cctbx::sgtbx::rt_mx> > >();
        to_and_from_python<
          tbxx::optional_copy<cctbx::sgtbx::rt_mx> >();
      }
    }
  };

  struct translation_part_info_wrappers
  {
    typedef translation_part_info w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("translation_part_info", no_init)
        .def(init<rt_mx const&>((arg("s"))))
        .def("intrinsic_part", &w_t::intrinsic_part, ccr())
        .def("location_part", &w_t::location_part, ccr())
        .def("origin_shift", &w_t::origin_shift, ccr())
      ;
    }
  };

} // namespace <anoymous>

  void wrap_rt_mx()
  {
    rt_mx_wrappers::wrap();
    translation_part_info_wrappers::wrap();
  }

}}} // namespace cctbx::sgtbx::boost_python
