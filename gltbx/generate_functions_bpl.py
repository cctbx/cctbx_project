from __future__ import absolute_import, division, print_function
from libtbx.utils import write_this_is_auto_generated
from libtbx.str_utils import line_breaker
import libtbx.load_env
import libtbx.path
import os
import sys
from six.moves import range

this = "gltbx.generate_functions_bpl"

return_types = {
  "GLenum": 0,
  "GLboolean": 0,
  "GLint": 0,
  "GLuint": 0,
  "const GLubyte*": 0,
  "GLUnurbs*": 0,
  "GLUquadric*": 0,
  "GLUtesselator*": 0,
}

arg_types = {
  "const void*": 0,
  "GLbitfield": 0,
  "GLboolean": 0,
  "GLboolean*": 0,
  "GLbyte": 0,
  "GLclampd": 0,
  "GLclampf": 0,
  "GLdouble": 0,
  "GLdouble*": 0,
  "GLenum": 0,
  "GLfloat": 0,
  "GLfloat*": 0,
  "GLint": 0,
  "GLint*": 0,
  "GLshort": 0,
  "GLsizei": 0,
  "GLubyte": 0,
  "GLubyte*": 0,
  "GLuint": 0,
  "GLuint*": 0,
  "GLushort": 0,
  "GLushort*": 0,
  "GLvoid*": 0,
  "GLvoid**": 0,
  "const GLboolean*": 0,
  "const GLbyte*": 0,
  "const GLclampf*": 0,
  "const GLdouble*": 0,
  "const GLfloat*": 0,
  "const GLint*": 0,
  "const GLshort*": 0,
  "const GLubyte*": 0,
  "const GLuint*": 0,
  "const GLushort*": 0,
  "const GLvoid*": 0,
  "GLUnurbs*": 0,
  "GLUquadric*": 0,
  "GLUtesselator*": 0,
  "glu_function_pointer": 0,
}

opaque_pointers = [
  "GLUnurbs*",
  "GLUquadric*",
  "GLUtesselator*",
]

pointee_sizes = {
  "glAreTexturesResident textures": 0,
  "glAreTexturesResident residences": 0,
  "glBitmap bitmap": 0,
  "glCallLists lists": "?n*sizeof(type)",
  "glClipPlane equation": 4,
  "glColor3bv v": 3,
  "glColor3dv v": 3,
  "glColor3fv v": 3,
  "glColor3iv v": 3,
  "glColor3sv v": 3,
  "glColor3ubv v": 3,
  "glColor3uiv v": 3,
  "glColor3usv v": 3,
  "glColor4bv v": 4,
  "glColor4dv v": 4,
  "glColor4fv v": 4,
  "glColor4iv v": 4,
  "glColor4sv v": 4,
  "glColor4ubv v": 4,
  "glColor4uiv v": 4,
  "glColor4usv v": 4,
  "glColorPointer pointer": 0,
  "glDeleteTextures textures": 0,
  "glDrawElements indices": 0,
  "glDrawPixels pixels": 0,
  "glEdgeFlagv flag": 1,
  "glEdgeFlagPointer pointer": 0,
  "glEvalCoord1dv u": 1,
  "glEvalCoord1fv u": 1,
  "glEvalCoord2dv u": 2,
  "glEvalCoord2fv u": 2,
  "glFeedbackBuffer buffer": "size",
  "glFogfv params": "?pname=GL_FOG_COLOR: 4, default: 1",
  "glFogiv params": "?pname=GL_FOG_COLOR: 4, default: 1",
  "glGenTextures textures": "n",
  "glGetClipPlane equation": 4,
  "glGetBooleanv params": "?1..16 depending on pname",
  "glGetDoublev params": "?1..16 depending on pname",
  "glGetFloatv params": "?1..16 depending on pname",
  "glGetIntegerv params": "?1..16 depending on pname",
  "glGetLightfv params": "?1..4 depending on pname",
  "glGetLightiv params": "?1..4 depending on pname",
  "glGetMapdv v": 0,
  "glGetMapfv v": 0,
  "glGetMapiv v": 0,
  "glGetMaterialfv params": "?1..4 depending on pname",
  "glGetMaterialiv params": 0,
  "glGetPixelMapfv values": "?glGet(map)",
  "glGetPixelMapuiv values": 0,
  "glGetPixelMapusv values": 0,
  "glGetPointerv params": 0,
  "glGetPolygonStipple mask": 0,
  "glGetTexEnvfv params": 0,
  "glGetTexEnviv params": 0,
  "glGetTexGendv params": 0,
  "glGetTexGenfv params": 0,
  "glGetTexGeniv params": 0,
  "glGetTexImage pixels": 0,
  "glGetTexLevelParameterfv params": 0,
  "glGetTexLevelParameteriv params": 0,
  "glGetTexParameterfv params": 0,
  "glGetTexParameteriv params": 0,
  "glIndexdv c": 0,
  "glIndexfv c": 0,
  "glIndexiv c": 0,
  "glIndexsv c": 0,
  "glIndexubv c": 0,
  "glIndexPointer pointer": 0,
  "glInterleavedArrays pointer": 0,
  "glLightfv params": 0,
  "glLightiv params": 0,
  "glLightModelfv params": 0,
  "glLightModeliv params": 0,
  "glLoadMatrixd m": 0,
  "glLoadMatrixf m": 0,
  "glMap1d points": 0,
  "glMap1f points": 0,
  "glMap2d points": 0,
  "glMap2f points": 0,
  "glMaterialfv params": 0,
  "glMaterialiv params": 0,
  "glMultMatrixd m": 0,
  "glMultMatrixf m": 0,
  "glNormal3bv v": 3,
  "glNormal3dv v": 3,
  "glNormal3fv v": 3,
  "glNormal3iv v": 3,
  "glNormal3sv v": 3,
  "glNormalPointer pointer": 0,
  "glPixelMapfv values": 0,
  "glPixelMapuiv values": 0,
  "glPixelMapusv values": 0,
  "glPolygonStipple mask": 0,
  "glPrioritizeTextures textures": 0,
  "glPrioritizeTextures priorities": 0,
  "glRasterPos2dv v": 2,
  "glRasterPos2fv v": 2,
  "glRasterPos2iv v": 2,
  "glRasterPos2sv v": 2,
  "glRasterPos3dv v": 3,
  "glRasterPos3fv v": 3,
  "glRasterPos3iv v": 3,
  "glRasterPos3sv v": 3,
  "glRasterPos4dv v": 4,
  "glRasterPos4fv v": 4,
  "glRasterPos4iv v": 4,
  "glRasterPos4sv v": 4,
  "glReadPixels pixels": 0,
  "glRectdv v1": 2,
  "glRectdv v2": 2,
  "glRectfv v1": 2,
  "glRectfv v2": 2,
  "glRectiv v1": 2,
  "glRectiv v2": 2,
  "glRectsv v1": 2,
  "glRectsv v2": 2,
  "glSelectBuffer buffer": 0,
  "glTexCoord1dv v": 1,
  "glTexCoord1fv v": 1,
  "glTexCoord1iv v": 1,
  "glTexCoord1sv v": 1,
  "glTexCoord2dv v": 2,
  "glTexCoord2fv v": 2,
  "glTexCoord2iv v": 2,
  "glTexCoord2sv v": 2,
  "glTexCoord3dv v": 3,
  "glTexCoord3fv v": 3,
  "glTexCoord3iv v": 3,
  "glTexCoord3sv v": 3,
  "glTexCoord4dv v": 4,
  "glTexCoord4fv v": 4,
  "glTexCoord4iv v": 4,
  "glTexCoord4sv v": 4,
  "glTexCoordPointer pointer": 0,
  "glTexEnvfv params": 0,
  "glTexEnviv params": 0,
  "glTexGendv params": 0,
  "glTexGenfv params": 0,
  "glTexGeniv params": 0,
  "glTexImage1D pixels": 0,
  "glTexImage2D pixels": 0,
  "glTexParameterfv params": 0,
  "glTexParameteriv params": 0,
  "glTexSubImage1D pixels": 0,
  "glTexSubImage2D pixels": 0,
  "gluBeginCurve nurb": 0,
  "gluEndCurve nurb": 0,
  "gluBeginPolygon tess": 0,
  "gluEndPolygon tess": 0,
  "gluBeginSurface nurb": 0,
  "gluEndSurface nurb": 0,
  "gluBeginTrim nurb": 0,
  "gluEndTrim nurb": 0,
  "gluBuild1DMipmaps data": 0,
  "gluBuild2DMipmaps data": 0,
  "gluCylinder quad": 0,
  "gluDeleteNurbsRenderer nurb": 0,
  "gluDeleteQuadric quad": 0,
  "gluDeleteTess tess": 0,
  "gluDisk quad": 0,
  "gluGetNurbsProperty nurb": 0,
  "gluGetNurbsProperty data": 0,
  "gluGetTessProperty tess": 0,
  "gluGetTessProperty data": 0,
  "gluLoadSamplingMatrices nurb": 0,
  "gluLoadSamplingMatrices model": 16,
  "gluLoadSamplingMatrices perspective": 16,
  "gluLoadSamplingMatrices view": 4,
  "gluNextContour tess": 0,
  "gluNurbsCallbackDataEXT nurb": 0,
  "gluNurbsCallbackDataEXT userData": 0,
  "gluNurbsCallback nurb": 0,
  "gluNurbsCurve nurb": 0,
  "gluNurbsCurve knots": 0,
  "gluNurbsCurve control": 0,
  "gluNurbsProperty nurb": 0,
  "gluNurbsSurface nurb": 0,
  "gluNurbsSurface sKnots": 0,
  "gluNurbsSurface tKnots": 0,
  "gluNurbsSurface control": 0,
  "gluPartialDisk quad": 0,
  "gluPickMatrix viewport": 4,
  "gluProject model": 16,
  "gluProject proj": 16,
  "gluProject view": 4,
  "gluProject winX": 1,
  "gluProject winY": 1,
  "gluProject winZ": 1,
  "gluPwlCurve nurb": 0,
  "gluPwlCurve data": 0,
  "gluQuadricCallback quad": 0,
  "gluQuadricDrawStyle quad": 0,
  "gluQuadricNormals quad": 0,
  "gluQuadricOrientation quad": 0,
  "gluQuadricTexture quad": 0,
  "gluScaleImage dataIn": 0,
  "gluScaleImage dataOut": 0,
  "gluSphere quad": 0,
  "gluTessBeginContour tess": 0,
  "gluTessEndContour tess": 0,
  "gluTessBeginPolygon tess": 0,
  "gluTessBeginPolygon data": 0,
  "gluTessCallback tess": 0,
  "gluTessEndPolygon tess": 0,
  "gluTessNormal tess": 0,
  "gluTessProperty tess": 0,
  "gluTessVertex tess": 0,
  "gluTessVertex location": 0,
  "gluTessVertex data": 0,
  "gluUnProject model": 16,
  "gluUnProject proj": 16,
  "gluUnProject view": 4,
  "gluUnProject objX": 1,
  "gluUnProject objY": 1,
  "gluUnProject objZ": 1,
  "glVertex2dv v": 2,
  "glVertex2fv v": 2,
  "glVertex2iv v": 2,
  "glVertex2sv v": 2,
  "glVertex3dv v": 3,
  "glVertex3fv v": 3,
  "glVertex3iv v": 3,
  "glVertex3sv v": 3,
  "glVertex4dv v": 4,
  "glVertex4fv v": 4,
  "glVertex4iv v": 4,
  "glVertex4sv v": 4,
  "glVertexPointer pointer": 0,
}

version_guards = {
  "glBlendColorEXT": "GL_XXX",
  "glEdgeFlagPointer": "GLTBX_XXX",
  "gluNurbsCallbackDataEXT": "GL_XXX",
}

special_wrappers = {

"glGetString": [
"""\
  boost::python::str
  gl_GetString(boost::python::object const& py_name)
  {
    boost::python::extract<GLenum> name_proxy(py_name);
    GLenum name = name_proxy();
    boost::python::str result(
      reinterpret_cast<const char*>(glGetString(name)));
    return result;
  }
""",
None
],

"gluGetString": [
"""\
  boost::python::str
  glu_GetString(boost::python::object const& py_name)
  {
    boost::python::extract<GLenum> name_proxy(py_name);
    GLenum name = name_proxy();
    boost::python::str result(
      reinterpret_cast<const char*>(gluGetString(name)));
    return result;
  }
""",
None
],

"gluErrorString": [
"""\
  boost::python::str
  glu_ErrorString(boost::python::object const& py_error)
  {
    boost::python::extract<GLenum> error_proxy(py_error);
    GLenum error = error_proxy();
    return boost::python::str(
      reinterpret_cast<const char*>(gluErrorString(error)));
  }
""",
None
],

}

def bytes_converters(signature, expected_size="0", post_extract=""):
  assert signature.return_type == "void"
  function_name = signature.function_name
  arg_type = signature.args[-1].type
  arg_name = signature.args[-1].name
  arg_type_name = arg_type+" "+arg_name
  is_const = arg_type.startswith("const ")
  call = "\n".join(signature.format_call(
    return_directly=is_const,
    prefix="    "))
  if (not is_const):
    call += "\n    %s_proxy.write_back();" % arg_name
    is_const = "false"
  else:
    is_const = "true"
  extracts = [""]
  for arg in signature.args[:-1]:
    assert not arg.type.startswith("const ")
    extracts.append("boost::python::extract<%s> %s_proxy(py_%s);" % (
      arg.type, arg.name, arg.name))
    extracts.append("%s %s = %s_proxy();" % (
      arg.type, arg.name, arg.name))
  extracts = "\n  ".join(extracts)
  return """\
%(extracts)s%(post_extract)s
  if      (type == GL_BYTE) {
    boost_python::converter_str<GLubyte> %(arg_name)s_proxy(
      "%(arg_name)s", py_%(arg_name)s, %(expected_size)s, %(is_const)s);
    %(arg_type_name)s = reinterpret_cast<%(arg_type)s>(
      %(arg_name)s_proxy.get());
%(call)s
  }
  else if (type == GL_UNSIGNED_BYTE) {
    boost_python::converter_str<GLbyte> %(arg_name)s_proxy(
      "%(arg_name)s", py_%(arg_name)s, %(expected_size)s, %(is_const)s);
    %(arg_type_name)s = reinterpret_cast<%(arg_type)s>(
      %(arg_name)s_proxy.get());
%(call)s
  }
  else {
    throw std::runtime_error(
      "Conversion not implemented for given GLenum type:"
      " %(function_name)s(): %(arg_type_name)s");
  }""" % vars()

def glReadPixels_wrapper_body(signature):
  return bytes_converters(
    signature=signature,
    expected_size="expected_size",
    post_extract="""
  boost::python::ssize_t expected_size = glReadPixels_pixels_expected_size(
    width, height, format, type);""")

special_wrapper_bodies = {

"glCallLists": bytes_converters,
"glDrawPixels": bytes_converters,
"glGetTexImage": bytes_converters,
"glReadPixels": glReadPixels_wrapper_body,
"glTexImage1D": bytes_converters,
"glTexImage2D": bytes_converters,
"glTexSubImage1D": bytes_converters,
"glTexSubImage2D": bytes_converters,

}

class argument:

  def __init__(self, function_name, string):
    fields = string.split()
    self.type = " ".join(fields[:-1])
    self.name = fields[-1]
    assert self.type in arg_types
    arg_types[self.type] += 1
    if (self.type[-1] != "*"):
      self.pointee_size = None
    else:
      self.pointee_size = pointee_sizes[function_name + " " + self.name]

class signature:

  def __init__(self, string):
    assert string.endswith(" )")
    fields = string[:-2].split("(")
    assert len(fields) == 2
    arg_strings = []
    for arg in fields[1].split(","):
      arg_strings.append(
        " ".join(arg.replace("*", " * ").split()).replace(" *", "*"))
    fields = fields[0].split()
    self.return_type = " ".join(" ".join(fields[:-1])
      .replace("*", " * ").split()).replace(" *", "*")
    self.function_name = fields[-1]
    if (self.return_type != "void"):
      assert self.return_type in return_types
      return_types[self.return_type] += 1
    self.args = []
    if (arg_strings != ["void"]):
      for arg in arg_strings:
        self.args.append(argument(self.function_name, arg))
    self.version_guard = version_guards.get(self.function_name, None)
    self.have_opaque_pointer = self.return_type in opaque_pointers
    if (not self.have_opaque_pointer):
      for arg in self.args:
        if (arg.type in opaque_pointers):
          self.have_opaque_pointer = True

  def show(self, f=None):
    if (f is None): f = sys.stdout
    print("function name:", self.function_name, file=f)
    print("  return type:", self.return_type, file=f)
    for arg in self.args:
      print("  arg type:", arg.type, "name:", arg.name, file=f)

  def wrapper_function_name(self):
    i = 2
    if (self.function_name.startswith("glu")): i = 3
    return self.function_name[:i]+"_"+self.function_name[i:]

  def write_no_opaque_pointers_guard_if(self, f):
    if (self.have_opaque_pointer):
      print("#if !defined(GLTBX_NO_OPAQUE_POINTERS)", file=f)

  def write_no_opaque_pointers_guard_endif(self, f):
    if (self.have_opaque_pointer):
      print("#endif", file=f)

  def write_version_guard_if(self, f):
    if (self.version_guard is not None):
      print("#if defined(%s)" % self.version_guard, file=f)

  def write_version_guard_endif(self, f):
    if (self.version_guard is not None):
      print("#endif", file=f)

  def format_call(self, return_directly, prefix):
    s = ""
    if (self.return_type != "void"):
      if (return_directly):
        s += "return "
      else:
        s += self.return_type + " result = "
    s += self.function_name+"("
    s += ", ".join([arg.name for arg in self.args])
    s += ");"
    result = []
    indent = ""
    for line in line_breaker(s, 70):
      result.append(prefix+indent+line)
      indent = "  "
    return result

  def write_wrapper(self, f):
    special = special_wrappers.get(self.function_name, None)
    if (special is not None and special[0] is not None):
      print(special[0], file=f)
      return
    lines = [
      self.return_type,
      self.wrapper_function_name()+"("
    ]
    for arg in self.args:
      lines.append("  %s %s," % (
        "boost::python::object const&", "py_"+arg.name))
    if (lines[-1][-1] == ","):
      lines[-1] = lines[-1][:-1]
    lines[-1] += ")"
    lines.append("{")
    special_body = special_wrapper_bodies.get(self.function_name, None)
    if (special_body is not None):
      lines.extend(special_body(self).splitlines())
    else:
      not_implemented = [
        "const void*",
        "GLvoid*",
        "GLvoid**",
        "const GLvoid*",
        "glu_function_pointer"]
      to_write_back = []
      ss = ""
      for arg in self.args:
        if ((arg.pointee_size is not None
             and arg.type not in opaque_pointers)
            or arg.type == "glu_function_pointer"):
          if (arg.type in not_implemented):
            lines.append(ss+"  throw std::runtime_error(")
            lines.append(ss+'    "Conversion not implemented:"')
            lines.append(ss+'    " %s(): %s %s");' % (
                self.function_name, arg.type, arg.name))
            ss = "//"
            lines.append(ss+"  %s %s = 0;" % (arg.type, arg.name))
          else:
            expected_size = arg.pointee_size
            if (isinstance(expected_size, str)):
              if (expected_size[0] == "?"):
                expected_size = "0"
            else:
              assert isinstance(expected_size, int)
              expected_size = str(expected_size)
            if (arg.type.startswith("const ")):
              is_const = "true"
              converter_t = arg.type[6:-1]
            else:
              is_const = "false"
              converter_t = arg.type[:-1]
            if (   arg.type.endswith("GLbyte*")
                or arg.type.endswith("GLubyte*")):
              converter = "boost_python::converter_str"
            else:
              converter = "boost_python::converter"
            lines.append(ss+'  %s<%s> %s_proxy(' % (
              converter,
              converter_t,
              arg.name))
            lines.append(ss+'    "%s", py_%s, %s, %s);' % (
              arg.name,
              arg.name,
              expected_size,
              is_const))
            lines.append(ss+"  %s %s = %s_proxy.get();" % (
              arg.type, arg.name, arg.name))
            if (is_const == "false"): to_write_back.append(arg)
        else:
          assert not arg.type.startswith("const ")
          lines.append(ss+'  boost::python::extract<%s> %s_proxy(py_%s);' % (
            arg.type, arg.name, arg.name))
          lines.append(ss+"  %s %s = %s_proxy();" % (
            arg.type, arg.name, arg.name))
      return_directly = len(to_write_back) == 0
      lines.extend([ss+line for line in
        self.format_call(return_directly=return_directly, prefix="  ")])
      for arg in to_write_back:
        lines.append(ss+"  %s_proxy.write_back();" % arg.name)
      if (self.return_type != "void" and not return_directly):
        lines.append(ss+"  return result;")
    lines.append("}")
    self.write_no_opaque_pointers_guard_if(f=f)
    self.write_version_guard_if(f=f)
    for line in lines:
      print(" ", line, file=f)
    self.write_version_guard_endif(f=f)
    self.write_no_opaque_pointers_guard_endif(f=f)
    print(file=f)

  def write_def(self, f):
    special = special_wrappers.get(self.function_name, None)
    if (special is not None and special[1] is not None):
      print(special[1], file=f)
      return
    return_opaque = self.return_type in opaque_pointers
    def_args = (self.function_name, self.wrapper_function_name())
    self.write_no_opaque_pointers_guard_if(f=f)
    self.write_version_guard_if(f=f)
    if (len(self.args) == 0):
      if (not return_opaque):
        print('    def("%s", %s);' % def_args, file=f)
      else:
        print('    def("%s", %s,' % def_args, file=f)
        print("      return_value_policy<return_opaque_pointer>());", file=f)
    else:
      assert not return_opaque
      print('    def("%s", %s, (' % def_args, file=f)
      s = ""
      for arg in self.args:
        s += ', arg("%s")' % arg.name
      s = s[2:] + "));"
      for line in line_breaker(s, 73):
        print("      "+line, file=f)
    self.write_version_guard_endif(f=f)
    self.write_no_opaque_pointers_guard_endif(f=f)

def get_signatures():
  result = []
  specs_file = libtbx.env.under_dist("gltbx", "opengl_specs.txt")
  for line in open(specs_file).read().splitlines():
    if (not (line.startswith("GL_") or line.startswith("GLU_"))):
      result.append(signature(line))
  return result

def write_function_wrappers(f, namespace, signatures, i_fragment):
  write_this_is_auto_generated(f, this)
  print("""\
#include <gltbx/special_wrapper_support.h>
#include <gltbx/pointer_args_bpl.h>
#include <gltbx/error.h>
""", file=f)
  if (namespace == "glu"):
    print("#if defined(__GNUC__) && __GNUC__ == 2 \\", file=f)
    print("     && __GNUC_MINOR__ == 96 && __GNUC_PATCHLEVEL__ == 0", file=f)
    print("#define GLTBX_NO_OPAQUE_POINTERS", file=f)
    print("#else", file=f)
    print("#include <boost/python/return_value_policy.hpp>", file=f)
    print("#include <boost/python/return_opaque_pointer.hpp>", file=f)
    for opaque_pointer in opaque_pointers:
      print("BOOST_PYTHON_OPAQUE_SPECIALIZED_TYPE_ID(%s)" % (
        opaque_pointer[:-1]), file=f)
    print("#endif", file=f)
    print(file=f)
  print("""\
namespace gltbx { namespace %s { namespace {
""" % namespace, file=f)
  for signature in signatures:
    signature.write_wrapper(f=f)
  print("""\
} // namespace <anonymous>

namespace boost_python {

  void
  wrap_functions_%02d()
  {
    using namespace boost::python;""" % i_fragment, file=f)
  for signature in signatures:
    signature.write_def(f=f)
  print("""\
  }

}}} // namespace gltbx::%s::boost_python""" % namespace, file=f)

def run(target_dir):
  if (not os.path.isdir(target_dir)):
    os.makedirs(target_dir)
  gl_signatures = []
  glu_signatures = []
  for signature in get_signatures():
    if (signature.function_name.startswith("glu")):
      glu_signatures.append(signature)
    else:
      gl_signatures.append(signature)
  for namespace,signatures,n_fragments in [("gl", gl_signatures, 16),
                                           ("glu", glu_signatures, 4)]:
    block_size = len(signatures) // n_fragments
    if (block_size * n_fragments < len(signatures)):
      block_size += 1
    for i_fragment in range(n_fragments):
      file_name = libtbx.path.norm_join(
        target_dir, namespace+"_functions_%02d_bpl.cpp" % i_fragment)
      with open(file_name, "w") as f:
        write_function_wrappers(
          f=f,
          namespace=namespace,
          signatures=signatures[i_fragment*block_size:(i_fragment+1)*block_size],
          i_fragment=i_fragment)

if __name__ == "__main__":
  run(".")
