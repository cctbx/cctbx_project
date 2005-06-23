#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include <gltbx/error.h>

namespace gltbx { namespace util {

  void
  TranslateScene(
    double s,
    double x,
    double y,
    double mousex,
    double mousey)
  {
    glMatrixMode(GL_MODELVIEW); handle_error();
    GLdouble mat[16];
    glGetDoublev(GL_MODELVIEW_MATRIX, mat); handle_error();
    glLoadIdentity(); handle_error();
    glTranslated((s * (x - mousex)), (s * (mousey - y)), 0.0); handle_error();
    glMultMatrixd(mat); handle_error();
  }

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
    glMatrixMode(GL_MODELVIEW); handle_error();
    GLdouble mat[16];
    glGetDoublev(GL_MODELVIEW_MATRIX, mat); handle_error();
    glLoadIdentity(); handle_error();
    glTranslated(xcenter, ycenter, zcenter); handle_error();
    glRotated((s * (y - mousey)), 1.0, 0.0, 0.0); handle_error();
    glRotated((s * (x - mousex)), 0.0, 1.0, 0.0); handle_error();
    glTranslated(-xcenter, -ycenter, -zcenter); handle_error();
    glMultMatrixd(mat); handle_error();
  }

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
