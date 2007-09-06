#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>

#include <scitbx/iso_surface.h>

namespace scitbx { namespace iso_surface { namespace boost_python {

  template <class CoordinatesType, class ValueType>
  struct triangulation_wrapper
  {
    typedef triangulation<CoordinatesType, ValueType> wt;

    static void wrap(const char *name) {
      using namespace boost::python;

      class_<wt>(name, no_init)
        .def(init<typename wt::map_const_ref_type,
                  ValueType, scitbx::vec3<CoordinatesType> const&>
                  ((arg("map"), arg("iso_level"), arg("grid_size"))))
        .add_property("vertices", &wt::vertices)
        .add_property("triangles", &wt::triangles)
        .add_property("normals", &wt::normals)
      ;
    }

  };

  void init_module() {
    triangulation_wrapper<double, double>::wrap("triangulation");
  }

}}}

BOOST_PYTHON_MODULE(scitbx_iso_surface_ext)
{
  scitbx::iso_surface::boost_python::init_module();
}

