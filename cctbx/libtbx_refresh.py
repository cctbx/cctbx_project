from cctbx.source_generators.eltbx import generate_henke_cpp
from cctbx.source_generators.eltbx import generate_sasaki_cpp
import libtbx.load_env
import os

target_dir = libtbx.env.under_build("cctbx/eltbx")
if (not os.path.isdir(target_dir)):
  os.makedirs(target_dir)
generate_henke_cpp.run(target_dir=target_dir)
generate_sasaki_cpp.run(target_dir=target_dir)
