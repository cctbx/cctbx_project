from cctbx.source_generators.eltbx import generate_henke_cpp
from cctbx.source_generators.eltbx import generate_sasaki_cpp
import libtbx.load_env
import os

if (libtbx.env.is_ready_for_build()):
  target_dir = libtbx.env.under_build("cctbx/eltbx")
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
