# $Id$

import makefile_generator

class write_makefiles(makefile_generator.write_makefiles):

  def dependencies(self):

    self.prefix_macro = "BOOST"

    self.files = (
      "libs/python/src/classes.cpp",
      "libs/python/src/conversions.cpp",
      "libs/python/src/extension_class.cpp",
      "libs/python/src/functions.cpp",
      "libs/python/src/init_function.cpp",
      "libs/python/src/module_builder.cpp",
      "libs/python/src/objects.cpp",
      "libs/python/src/types.cpp",
      "libs/python/src/cross_module.cpp",
    )

    self.libraries = {
      "boost_python": (
        "classes",
        "conversions",
        "extension_class",
        "functions",
        "init_function",
        "module_builder",
        "objects",
        "types",
        "cross_module",
      )
    }
