#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/list.hpp>

#include <gltbx/error.h>

namespace gltbx { namespace util {

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
      glGetIntegerv(pname, vector); handle_error();
      return as_python_list(vector, Size);
    }

    static
    boost::python::list
    double_(GLenum pname)
    {
      GLdouble vector[Size];
      glGetDoublev(pname, vector); handle_error();
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
  void
  TranslateScene(
    double s,
    double x,
    double y,
    double mousex,
    double mousey)
  {
    GLdouble mat[16];
    glGetDoublev(GL_MODELVIEW_MATRIX, mat); handle_error();
    glMatrixMode(GL_MODELVIEW); handle_error();
    glLoadIdentity(); handle_error();
    glTranslated((s * (x - mousex)), (s * (mousey - y)), 0.0); handle_error();
    glMultMatrixd(mat); handle_error();
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
    glGetDoublev(GL_MODELVIEW_MATRIX, mat); handle_error();
    glMatrixMode(GL_MODELVIEW); handle_error();
    glLoadIdentity(); handle_error();
    glTranslated(xcenter, ycenter, zcenter); handle_error();
    glRotated((s * (y - mousey)), 1.0, 0.0, 0.0); handle_error();
    glRotated((s * (x - mousex)), 0.0, 1.0, 0.0); handle_error();
    glTranslated(-xcenter, -ycenter, -zcenter); handle_error();
    glMultMatrixd(mat); handle_error();
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
    glMatrixMode(GL_MODELVIEW); handle_error();
    GLdouble mat[16];
    glGetDoublev(GL_MODELVIEW_MATRIX, mat); handle_error();
    glLoadIdentity(); handle_error();
    glTranslated(xcenter, ycenter, zcenter); handle_error();
    glRotated(angle, xvector, yvector, zvector); handle_error();
    glTranslated(-xcenter, -ycenter, -zcenter); handle_error();
    glMultMatrixd(mat); handle_error();
  }

  void
  init_module()
  {
    using namespace boost::python;
    def("get_gl_modelview_matrix", get_gl_modelview_matrix);
    def("get_gl_projection_matrix", get_gl_projection_matrix);
    def("get_gl_viewport", get_gl_viewport);
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
