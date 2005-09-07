// copyright (c) Jacob N. Smith & Erik McKee; leave this here; use at your whim
#include <cctbx/boost_python/flex_fwd.h>
#include <cctbx/maptbx/mapper.h>
#include <chiltbx/handle.h>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

namespace cctbx { namespace maptbx { namespace boost_python {

namespace {

struct mapper_wrappers {

  typedef double                      FloatType;
  typedef signed long                    IntType;

  typedef fractional<FloatType>              fractional;
  typedef grid_point<IntType>                unit_shift;
  typedef IntType                      symmetry_operation;
  typedef scitbx::mat3<FloatType>              matrix3;
  typedef sgtbx::space_group                space_group;
  typedef crystal::direct_space_asu::float_asu<FloatType>  float_asu;

  typedef basic_mapper<void,FloatType,IntType>      void_basic_mapper;
  typedef basic_mapper<non_symmetric,FloatType,IntType>  ns_basic_mapper;
  typedef basic_mapper<unit_cell,FloatType,IntType>    uc_basic_mapper;
  typedef basic_mapper<asu,FloatType,IntType>        as_basic_mapper;

  typedef mapper_factory<void,FloatType,IntType>      vv_factory;
  typedef mapper_factory<non_symmetric,FloatType,IntType>  ns_factory;
  typedef mapper_factory<unit_cell,FloatType,IntType>    uc_factory;
  typedef mapper_factory<asu,FloatType,IntType>      as_factory;

  typedef chiltbx::handle::handle<vv_factory>        factory_handle;

    static void_basic_mapper
    get_non_symmetric_mapper ( fractional const& f ) {
        return ns_basic_mapper(f);
    }
    static void_basic_mapper
    get_unit_cell_mapper ( fractional const& f ) {
        return uc_basic_mapper(f);
    }
    static void_basic_mapper
    get_asu_mapper ( fractional const& f,
        space_group const& sg, float_asu const& fa ) {
        return as_basic_mapper(f,sg,fa);
    }

  static void wrap() {

    using namespace boost::python;
    typedef return_value_policy<return_by_value> rbv;

    class_<void_basic_mapper>("basic_mapper",no_init)
      .add_property("mapped_coordinate",
        make_getter(&void_basic_mapper::mapped_coordinate, rbv()))
      .add_property("unit_shift",
        make_getter(&void_basic_mapper::unit_shift, rbv()))
      .add_property("symmetry_operation",
        make_getter(&void_basic_mapper::symmetry_operation, rbv()));

        def("get_non_symmetric_mapper",
            &mapper_wrappers::get_non_symmetric_mapper,
            arg_("coordinate"));

        def("get_unit_cell_mapper",
            &mapper_wrappers::get_unit_cell_mapper,
            arg_("coordinate"));

        def("get_asu_mapper",
            &mapper_wrappers::get_asu_mapper,
            (arg_("coordinate")
            ,arg_("space_group")
            ,arg_("float_asu")));

    class_<factory_handle>("factory_handle",no_init);

    class_<ns_factory>("non_symmetric_factory",init<>())
      .def("map",&ns_factory::map,arg_("coordinate"))
            .def("as_handle",&ns_factory::as_handle);

    class_<uc_factory>("unit_cell_factory",init<>())
      .def("map",&uc_factory::map,arg_("coordinate"))
            .def("as_handle",&uc_factory::as_handle);

    class_<as_factory>("asu_factory",no_init)
      .def(init<space_group,float_asu>())
      .def("map",&as_factory::map,arg_("coordinate"))
            .def("as_handle",&as_factory::as_handle);

  }
};

} // namespace <anoymous>

void wrap_mappers() {
  mapper_wrappers::wrap();
}

}}} // namespace cctbx::maptbx::boost_python
