#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/tuple.hpp>

#include <gltbx/quadrics.h>

namespace gltbx { namespace quadrics { namespace boost_python {

  #define COMMON_ARGS arg("draw_style")=GLU_FILL,\
                      arg("orientation")=GLU_OUTSIDE,\
                      arg("normals")=GLU_SMOOTH

  struct proto_cylinder_wrapper
  {
    typedef proto_cylinder wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt>("proto_cylinder", no_init)
        .def(init<GLdouble, GLint, GLint, GLenum, GLenum, GLenum>(
             (arg("top_to_base_radius_ratio"), arg("slices"), arg("stacks"),
              COMMON_ARGS)))
        .def(init<GLint, GLenum, GLenum, GLenum>(
                  (arg("slices"), COMMON_ARGS)))
        .def("draw", &wt::draw, (arg("start"), arg("end"), arg("base_radius")))
      ;
    }
  };


  struct ellipsoid_to_sphere_transform_wrapper
  {
    typedef ellipsoid_to_sphere_transform wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt>("ellipsoid_to_sphere_transform", no_init)
        .def(init<scitbx::vec3<GLdouble> const &,
                  scitbx::sym_mat3<GLdouble> const &>
                ((arg("centre"), arg("metrics"))))
        .def("non_positive_definite", &wt::non_positive_definite)
        .def("linear_part", &wt::linear_part)
        .def("translation_part", &wt::translation_part)
      ;
    }

  };

  struct proto_ellipsoid_wrapper
  {
    typedef proto_ellipsoid wt;

    static void wrap() {
      using namespace boost::python;
      void (wt::*draw_1)(wt::vec3_t const &, wt::metrics_t const &) = &wt::draw;
      void (wt::*draw_2)(ellipsoid_to_sphere_transform const &) = &wt::draw     ;
      class_<wt>("proto_ellipsoid", no_init)
        .def(init<GLint, GLint, GLenum, GLenum, GLenum>(
             (arg("slices"), arg("stacks"), COMMON_ARGS)))
        .def("draw", draw_1, (arg("centre"), arg("metrics")))
        .def("draw", draw_2, arg("ellipsoid_to_sphere_transform"))
      ;
    }
  };

  namespace {
    void init_module() {
      proto_cylinder_wrapper::wrap();
      proto_ellipsoid_wrapper::wrap();
      ellipsoid_to_sphere_transform_wrapper::wrap();
    }
  }

}}} // gltbx::quadrics::boost_python

BOOST_PYTHON_MODULE(gltbx_quadrics_ext)
{
  gltbx::quadrics::boost_python::init_module();
}
