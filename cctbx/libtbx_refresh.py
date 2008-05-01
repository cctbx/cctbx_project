import os

if (self.env.is_ready_for_build()):
  message_template = '  Generating C++ files in:\n    "%s"'

  # eltbx
  from cctbx.source_generators.eltbx import generate_henke_cpp
  from cctbx.source_generators.eltbx import generate_sasaki_cpp
  target_dir = self.env.under_build("cctbx/eltbx")
  print message_template % target_dir
  for label,generator_module in [("Henke", generate_henke_cpp),
                                 ("Sasaki", generate_sasaki_cpp)]:
    if (os.path.isdir(generator_module.reference_tables_directory)):
      if (not os.path.isdir(target_dir)):
        os.makedirs(target_dir)
      generator_module.run(target_dir=target_dir)
    else:
      print "*"*79
      print "Warning: directory with %s tables is missing:" % label
      print " ", repr(generator_module.reference_tables_directory)
      print "*"*79

  # flex_fwd.h
  from cctbx.source_generators import flex_fwd_h
  target_dir = self.env.under_build("include/cctbx/boost_python")
  print message_template % target_dir
  if not os.path.isdir(target_dir):
    os.makedirs(target_dir)
  flex_fwd_h.run(target_dir)
