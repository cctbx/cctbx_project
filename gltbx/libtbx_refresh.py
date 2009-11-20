from gltbx import generate_defines_bpl
from gltbx import generate_functions_bpl
from gltbx import generate_fonts_ucs_cpp
import sys, os

if (self.env.is_ready_for_build()):
  target_dir = self.env.under_build("gltbx")
  print '  Generating C++ files in:\n    "%s"' % target_dir
  generate_defines_bpl.run(target_dir=target_dir)
  generate_functions_bpl.run(target_dir=target_dir)
  generate_fonts_ucs_cpp.run(target_dir=target_dir)

  sources = ["lib/libpng.a", "include/png.h", "include/pngconf.h"]
  for file_name in sources :
    source = self.env.under_build("base/%s" % file_name)
    if (os.path.isfile(source)):
      target = self.env.under_build(file_name)
      print "  Copying: %s" % file_name
      open(target, "wb").write(open(source, "rb").read())
