def touch_init(env, target, source):
  if (str(target[0]) != "libtbx/phaser_boost/__init__.py"):
    import warnings
    warnings.warn("""\
The touch_init() function is obsolete and will be
removed in the future. The target for all Boost.Python
extension modules should be a file at the top level
of the #libtbx directory. E.g.:

*_ext.cpp file:

  BOOST_PYTHON_MODULE(scitbx_math_ext)

SConscript:

  env_scitbx_boost_python_ext.SharedLibrary(
    target="#libtbx/scitbx_math_ext",
    source=["math_ext.cpp"])

*.py file:

  import libtbx.boost_python
  ext = libtbx.boost_python.import_ext("scitbx_math_ext")
  from scitbx_math_ext import *

Note that it is important to use the fully-qualified long name
in order to avoid name clashes (e.g. with a math_ext in another
package.)
""",
    DeprecationWarning)
  assert len(target) == 1
  print "Creating:", str(target[0])
  open(str(target[0]), "w").close()
