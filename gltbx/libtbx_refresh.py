from __future__ import absolute_import, division, print_function
from gltbx import generate_defines_bpl
from gltbx import generate_functions_bpl
from gltbx import generate_fonts_ucs_cpp

if self.env.is_ready_for_build():
  target_dir = self.env.under_build("gltbx")
  print('  Generating C++ files in:\n    "%s"' % target_dir)
  generate_defines_bpl.run(target_dir=target_dir)
  generate_functions_bpl.run(target_dir=target_dir)
  generate_fonts_ucs_cpp.run(target_dir=target_dir)
