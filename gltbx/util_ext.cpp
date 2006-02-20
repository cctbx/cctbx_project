#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/list.hpp>

#include <scitbx/mat3.h>

#include <gltbx/error.h>

namespace gltbx { namespace util {

  inline
  void
  translate(scitbx::vec3<double> const& t)
  {
    glTranslated(t[0], t[1], t[2]);
  }

  template <typename ElementType>
  boost::python::list
  as_python_list(ElementType* elements, unsigned size)
  {
    boost::python::list result;
    for(unsigned i=0;i<size;i++) result.append(elements[i]);
    return result;
  }

  template <unsigned Size>
  struct gl_vector_as_python_list
  {
    static
    boost::python::list
    int_(GLenum pname)
    {
      GLint vector[Size];
      glGetIntegerv(pname, vector);
      return as_python_list(vector, Size);
    }

    static
    boost::python::list
    double_(GLenum pname)
    {
      GLdouble vector[Size];
      glGetDoublev(pname, vector);
      return as_python_list(vector, Size);
    }
  };

  inline
  boost::python::list
  get_gl_modelview_matrix()
  {
    return gl_vector_as_python_list<16>::double_(GL_MODELVIEW_MATRIX);
  }

  inline
  boost::python::list
  get_gl_projection_matrix()
  {
    return gl_vector_as_python_list<16>::double_(GL_PROJECTION_MATRIX);
  }

  inline
  boost::python::list
  get_gl_viewport()
  {
    return gl_vector_as_python_list<4>::int_(GL_VIEWPORT);
  }

  inline
  scitbx::mat3<double>
  extract_rotation_from_gl_matrix(GLdouble* m)
  {
    return scitbx::mat3<double>(
      m[0], m[4], m[8],
      m[1], m[5], m[9],
      m[2], m[6], m[10]);
  }

  inline
  scitbx::vec3<double>
  extract_translation_from_gl_matrix(GLdouble* m)
  {
    return scitbx::vec3<double>(m[12], m[13], m[14]);
  }

  inline
  scitbx::mat3<double>
  extract_rotation_from_gl_modelview_matrix()
  {
    GLdouble mvm[16];
    glGetDoublev(GL_MODELVIEW_MATRIX, mvm);
    return extract_rotation_from_gl_matrix(mvm);
  }

  inline
  scitbx::vec3<double>
  extract_translation_from_gl_modelview_matrix()
  {
    GLdouble mvm[16];
    glGetDoublev(GL_MODELVIEW_MATRIX, mvm);
    return extract_translation_from_gl_matrix(mvm);
  }

  inline
  scitbx::vec3<double>
  object_as_eye_coordinates(
    scitbx::vec3<double> const& object_coordinates)
  {
    GLdouble mvm[16];
    glGetDoublev(GL_MODELVIEW_MATRIX, mvm);
    return extract_rotation_from_gl_matrix(mvm) * object_coordinates
         + extract_translation_from_gl_matrix(mvm);
  }

  inline
  void
  translate_object(double eye_x, double eye_y, double eye_z)
  {
    GLdouble mvm[16];
    glGetDoublev(GL_MODELVIEW_MATRIX, mvm);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslated(eye_x, eye_y, eye_z);
    glMultMatrixd(mvm);
  }

  inline
  void
  translate_object(scitbx::vec3<double> const& eye_vector)
  {
    translate_object(eye_vector[0], eye_vector[1], eye_vector[2]);
  }

  inline
  void
  translate_object(
    double s,
    double x,
    double y,
    double mousex,
    double mousey)
  {
    translate_object(s * (x - mousex), s * (mousey - y), 0.0);
  }

  inline
  void
  rotate_object_about_eye_x_and_y(
    double s,
    double xcenter,
    double ycenter,
    double zcenter,
    double x,
    double y,
    double mousex,
    double mousey)
  {
    GLdouble mvm[16];
    glGetDoublev(GL_MODELVIEW_MATRIX, mvm);
    scitbx::vec3<double> eye_center =
        extract_rotation_from_gl_matrix(mvm)
          * scitbx::vec3<double>(xcenter, ycenter, zcenter)
      + extract_translation_from_gl_matrix(mvm);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    translate(eye_center);
    glRotated((s * (y - mousey)), 1.0, 0.0, 0.0);
    glRotated((s * (x - mousex)), 0.0, 1.0, 0.0);
    translate(-eye_center);
    glMultMatrixd(mvm);
  }

  inline
  void
  rotate_object_about_eye_vector(
    double xcenter,
    double ycenter,
    double zcenter,
    double xvector,
    double yvector,
    double zvector,
    double angle)
  {
    GLdouble mvm[16];
    glGetDoublev(GL_MODELVIEW_MATRIX, mvm);
    scitbx::vec3<double> eye_center =
        extract_rotation_from_gl_matrix(mvm)
          * scitbx::vec3<double>(xcenter, ycenter, zcenter)
      + extract_translation_from_gl_matrix(mvm);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    translate(eye_center);
    glRotated(angle, xvector, yvector, zvector);
    translate(-eye_center);
    glMultMatrixd(mvm);
  }

  inline
  void
  TranslateScene(
    double s,
    double x,
    double y,
    double mousex,
    double mousey)
  {
    GLdouble mat[16];
    glGetDoublev(GL_MODELVIEW_MATRIX, mat);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslated((s * (x - mousex)), (s * (mousey - y)), 0.0);
    glMultMatrixd(mat);
  }

  inline
  void
  RotateScene(
    double s,
    double xcenter,
    double ycenter,
    double zcenter,
    double x,
    double y,
    double mousex,
    double mousey)
  {
    GLdouble mat[16];
    glGetDoublev(GL_MODELVIEW_MATRIX, mat);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslated(xcenter, ycenter, zcenter);
    glRotated((s * (y - mousey)), 1.0, 0.0, 0.0);
    glRotated((s * (x - mousex)), 0.0, 1.0, 0.0);
    glTranslated(-xcenter, -ycenter, -zcenter);
    glMultMatrixd(mat);
  }

  inline
  void
  RotateAboutVector(
    double xcenter,
    double ycenter,
    double zcenter,
    double xvector,
    double yvector,
    double zvector,
    double angle)
  {
    glMatrixMode(GL_MODELVIEW);
    GLdouble mat[16];
    glGetDoublev(GL_MODELVIEW_MATRIX, mat);
    glLoadIdentity();
    glTranslated(xcenter, ycenter, zcenter);
    glRotated(angle, xvector, yvector, zvector);
    glTranslated(-xcenter, -ycenter, -zcenter);
    glMultMatrixd(mat);
  }

  void
  init_module()
  {
    using namespace boost::python;
    def("handle_error", handle_error);
    def("get_gl_modelview_matrix", get_gl_modelview_matrix);
    def("get_gl_projection_matrix", get_gl_projection_matrix);
    def("get_gl_viewport", get_gl_viewport);
    def("extract_rotation_from_gl_modelview_matrix",
      extract_rotation_from_gl_modelview_matrix);
    def("object_as_eye_coordinates", object_as_eye_coordinates, (
      arg_("object_coordinates")));
    def("translate_object",
      (void(*)(double, double, double)) translate_object, (
        arg_("eye_x"),
        arg_("eye_y"),
        arg_("eye_z")));
    def("translate_object",
      (void(*)(scitbx::vec3<double> const&)) translate_object, (
        arg_("eye_vector")));
    def("translate_object",
      (void(*)(double, double, double, double, double)) translate_object, (
        arg_("s"),
        arg_("x"),
        arg_("y"),
        arg_("mousex"),
        arg_("mousey")));
    def("rotate_object_about_eye_x_and_y", rotate_object_about_eye_x_and_y, (
      arg_("s"),
      arg_("xcenter"),
      arg_("ycenter"),
      arg_("zcenter"),
      arg_("x"),
      arg_("y"),
      arg_("mousex"),
      arg_("mousey")));
    def("rotate_object_about_eye_vector", rotate_object_about_eye_vector, (
      arg_("xcenter"),
      arg_("ycenter"),
      arg_("zcenter"),
      arg_("xvector"),
      arg_("yvector"),
      arg_("zvector"),
      arg_("angle")));
    def("TranslateScene", TranslateScene, (
      arg_("s"),
      arg_("x"),
      arg_("y"),
      arg_("mousex"),
      arg_("mousey")));
    def("RotateScene", RotateScene, (
      arg_("s"),
      arg_("xcenter"),
      arg_("ycenter"),
      arg_("zcenter"),
      arg_("x"),
      arg_("y"),
      arg_("mousex"),
      arg_("mousey")));
    def("RotateAboutVector", RotateAboutVector, (
      arg_("xcenter"),
      arg_("ycenter"),
      arg_("zcenter"),
      arg_("xvector"),
      arg_("yvector"),
      arg_("zvector"),
      arg_("angle")));
  }

}} // namespace gltbx::util

BOOST_PYTHON_MODULE(gltbx_util_ext)
{
  gltbx::util::init_module();
}
