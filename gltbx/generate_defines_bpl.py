from __future__ import absolute_import, division, print_function
from libtbx.utils import write_this_is_auto_generated
import libtbx.load_env
import libtbx.path
import os
from six.moves import range

this = "gltbx.generate_defines_bpl"

additional_defines = """
GL_VERSION_1_1
GL_VERSION_1_2
GL_VERSION_1_3
GL_VERSION_1_4
GL_VERSION_1_5
GL_VERSION_2_0
GL_ARB_imaging
GLU_VERSION_1_1
GLU_VERSION_1_2
GLU_VERSION_1_3
GLU_VERSION_1_4
GLU_VERSION_1_5
GLU_VERSION_2_0
GL_MULTISAMPLE
GL_SAMPLE_BUFFERS
GL_SAMPLES
"""

def write_one(f, define):
  print('#if defined(%s)' % define, file=f)
  print('    scope.attr("%s") = %s;' % (define, define), file=f)
  print('#endif', file=f)

def write_define_wrappers(f, namespace, defines, i_fragment):
  write_this_is_auto_generated(f, this)
  print("""\
#include <boost/python/def.hpp>
#include <boost/python/scope.hpp>
#include <gltbx/include_opengl.h>

namespace gltbx { namespace %s { namespace boost_python {

  void
  wrap_defines_%02d(boost::python::scope scope)
  {""" % (namespace, i_fragment), file=f)
  for define in defines:
    write_one(f, define)
  if (namespace == "gl"):
    for i in range(4): write_one(f, "GL_AUX%d" % i)
    for i in range(8): write_one(f, "GL_LIGHT%d" % i)
    for i in range(32): write_one(f, "GL_TEXTURE%d" % i)
  print("""\
  }

}}} // gltbx::%s::boost_python""" % namespace, file=f)

def import_opengl_defines():
  result = []
  specs_file = libtbx.env.under_dist("gltbx", "opengl_specs.txt")
  for line in open(specs_file).read().splitlines():
    if (line.startswith("GL_") or line.startswith("GLU_")):
      result.append(line)
  return result

def run(target_dir):
  if not os.path.isdir(target_dir):
    os.makedirs(target_dir)
  gl_defines = []
  glu_defines = []
  for define in additional_defines.split() + import_opengl_defines():
    if (define.startswith("GLU_")):
      glu_defines.append(define)
    else:
      gl_defines.append(define)
  for namespace,defines,n_fragments in [("gl", gl_defines, 8),
                                        ("glu", glu_defines, 2)]:
    block_size = len(defines) // n_fragments
    if (block_size * n_fragments < len(defines)):
      block_size += 1
    for i_fragment in range(n_fragments):
      file_name = libtbx.path.norm_join(
        target_dir, namespace+"_defines_%02d_bpl.cpp" % i_fragment)
      with open(file_name, "w") as f:
        write_define_wrappers(
          f=f,
          namespace=namespace,
          defines=defines[i_fragment*block_size:(i_fragment+1)*block_size],
          i_fragment=i_fragment)

if __name__ == "__main__":
  run(".")
