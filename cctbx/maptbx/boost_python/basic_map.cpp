// copyright (c) Jacob N. Smith & Erik McKee; leave this here; use at your whim
#include <cctbx/boost_python/flex_fwd.h>
#include <cctbx/maptbx/basic_map.h>
#include <cctbx/maptbx/generic_grid.h>
#include <cctbx/maptbx/map_out_of_bounds.h>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <chiltbx/handle.h>

namespace cctbx { namespace maptbx { namespace boost_python {

namespace {

struct basic_map_wrapper {

  typedef double                                        FloatType;
  typedef signed long                                   IntType;

  typedef out_of_bounds<void,FloatType,IntType>         out_of_intf;
  typedef chiltbx::handle::handle<out_of_intf>          out_of_handle;
  typedef out_of_bounds<raise,FloatType,IntType>        out_of_raise;
  typedef out_of_bounds<substitute,FloatType,IntType>   out_of_substitute;
  typedef out_of_bounds<clamp,FloatType,IntType>        out_of_clamp;
  typedef out_of_bounds<interpolate,FloatType,IntType>  out_of_interpolate;

  typedef basic_map<FloatType,IntType>                  basic_map_type;

  typedef af::flex_grid<>                               flex_grid;
  typedef af::c_grid_padded<dimension_3>                c_grid_padded;
  typedef af::versa<FloatType,flex_grid>                af_versa;
  typedef af::versa<FloatType,c_grid_padded>            af_c_grid;
  typedef af::tiny<IntType,dimension_3>                 af_extents;
  typedef cdsa::float_asu<FloatType>                    cdsa_float_asu;
  typedef sgtbx::space_group                            tbx_space_group;
  typedef cctbx::uctbx::unit_cell                       tbx_unit_cell;

  typedef scitbx::mat3<FloatType>                       mat3;

  typedef fractional<FloatType>                         frac_type;
  typedef cartesian<FloatType>                          cart_type;
  typedef grid_point<IntType>                           grid_type;
  typedef transform<frac_type,grid_type>                f2g_type;
  typedef transform<frac_type,cart_type>                f2c_type;
  typedef transform<grid_type,frac_type>                g2f_type;
  typedef transform<cart_type,frac_type>                c2f_type;
  typedef transform<grid_type,cart_type>                g2c_type;
  typedef transform<cart_type,grid_type>                c2g_type;

  static void wrap () {

    using namespace boost::python;

    typedef return_value_policy<return_by_value> rbv;

    class_<asu>("basic_map_asu_flag",no_init)
      .def(init<>());
    class_<unit_cell>("basic_map_unit_cell_flag",no_init)
      .def(init<>());
    class_<non_symmetric>("basic_map_non_symmetric_flag",no_init)
      .def(init<>());

    class_<out_of_handle>("out_of_bounds_handle",no_init);
    class_<out_of_raise>("out_of_bounds_raise",init<>())
      .def("as_handle",&out_of_raise::as_handle);
    class_<out_of_substitute>("out_of_bounds_substitute",no_init)
      .def(init<FloatType>())
      .def("as_handle",&out_of_substitute::as_handle);
    class_<out_of_clamp>("out_of_bounds_clamp",no_init)
      .def(init<FloatType>())
      .def("as_handle",&out_of_clamp::as_handle);
    class_<out_of_interpolate>("out_of_bounds_interpolate",init<>())
      .def("as_handle",&out_of_interpolate::as_handle);

    class_<basic_map_type>("basic_map",no_init)
      .def(init<asu const&,
          af_versa const&,
          tbx_space_group const&,
          cdsa_float_asu const&,
          af_extents const&,
          mat3 const&,
          out_of_handle const&,
          tbx_unit_cell const&,
          FloatType const&,
          bool>())
      .def(init<unit_cell const&,
          af_versa const&,
          af_extents const&,
          mat3 const&,
          out_of_handle const&,
          tbx_unit_cell const&>())
      .def(init<unit_cell const&,
          af_c_grid const&,
          af_extents const&,
          mat3 const&,
          out_of_handle const&,
          tbx_unit_cell const&>())
      .def(init<non_symmetric const&,
          af_versa const&,
          af_extents const&,
          mat3 const&,
          out_of_handle const&,
          tbx_unit_cell const&>())
      .def(init<basic_map_type>())
      .def("as_asu",
        (void(basic_map_type::*)
          (af_versa const&
          ,tbx_space_group const&
          ,cdsa_float_asu const&
          ,af_extents const&
          ,FloatType const&
          ,bool))
        &basic_map_type::set_grid_handle,
        (arg("data")
        ,arg("space_group")
        ,arg("float_asu")
        ,arg("grid_length")
        ,arg("min_distance_sym_equiv")
        ,arg("assert_min_distance_sym_equiv")))
      .def("as_unit_cell",
        (void(basic_map_type::*)
          (af_versa const&))
        &basic_map_type::set_grid_handle,
        arg("data"))
      .def("as_non_symmetric",
        (void(basic_map_type::*)
          (af_versa const&
          ,af_extents const&))
        &basic_map_type::set_grid_handle,
        (arg("data")
        ,arg("grid_length")))
      .def("set_out_of_bounds_handle",
        &basic_map_type::set_out_of_bounds_handle,arg("handle"))
      .def("get_cart_value",
        &basic_map_type::get_cart_value,arg("coordinate"))
      .def("get_frac_value",
        &basic_map_type::get_frac_value,arg("coordinate"))
      .def("get_grid_value",
        &basic_map_type::get_grid_value,arg("coordinate"))
      .def("get_cart_values",
        &basic_map_type::get_cart_values,arg("coordinates"))
      .def("get_frac_values",
        &basic_map_type::get_frac_values,arg("coordinates"))
      .def("get_grid_values",
        &basic_map_type::get_grid_values,arg("coordinates"))
      .def("set_grid_value",
        &basic_map_type::set_grid_value,
        (arg("coordinate")
        ,arg("value")))
      .def("cart2frac",&basic_map_type::cart_to_frac)
      .def("cart2grid",&basic_map_type::cart_to_grid)
      .def("grid2frac",&basic_map_type::grid_to_frac)
      .def("cart2frac",&basic_map_type::cart_to_frac)
      .def("frac2grid",&basic_map_type::frac_to_grid)
      .def("frac2cart",&basic_map_type::frac_to_cart)
      .def("nearest_grid_point_cartesian",
        &basic_map_type::cart_nearest_grid_point,arg("coordinate"))
      .def("nearest_grid_point_fractional",
        &basic_map_type::frac_nearest_grid_point,arg("coordinate"))
      .def("value_at_nearest_grid_point_fractional",
        &basic_map_type::frac_value_at_nearest_grid_point,arg("coordinate"))
      .def("value_at_nearest_grid_point_cartesian",
        &basic_map_type::cart_value_at_nearest_grid_point,arg("coordinate"))
      .def("remap_cartesian",
        &basic_map_type::remap_cart_coordinate,arg("coordinate"))
      .def("remap_fractional",
        &basic_map_type::remap_frac_coordinate,arg("coordinate"))
      .def("remap_grid",
        &basic_map_type::remap_grid_coordinate,arg("coordinate"))
      .def("remap",&basic_map_type::remap,arg("coordinate"))
      .def("is_inside_grid",
        (bool(basic_map_type::*)(grid_type const&)const)
        &basic_map_type::is_inside,
        arg("coordinate"))
      .def("is_inside_frac",
        (bool(basic_map_type::*)(frac_type const&)const)
        &basic_map_type::is_inside,
        arg("coordinate"))
      .def("is_inside_cart",
        (bool(basic_map_type::*)(cart_type const&)const)
        &basic_map_type::is_inside,
        arg("coordinate"))
      .def("rebuild_transformers",&basic_map_type::rebuild_transformers,
        (arg("extents"),arg("matrix")))
      .def("extents",&basic_map_type::extents)
      .def("unit_cell",&basic_map_type::unit_cell)
      ;

  }

};

} // namespace <anoymous>

void wrap_basic_map() {
  basic_map_wrapper::wrap();
}

}}} // namespace cctbx::maptbx::boost_python
