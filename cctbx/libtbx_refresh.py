import os
from libtbx.utils import warn_if_unexpected_md5_hexdigest

self.remove_obsolete_pyc_if_possible(pyc_file_names=[
  "eltbx/xray_scattering.pyc", # XXX backward compatibility 2011-01-14
])

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

  # reference_table.cpp : checking that it is up-to-date
  for f,sig in [
      ("reference_table.py", "0e7cc97096ecb8f6dab2abce9fe98dba"),
      ("short_cuts.py", "92ad3a575f1ad230f7b9e89e947e7f5c"),
      ("proto/generate_cpp_asu_table.py", "78796aea988d52e0faa95aa0c40dba86")]:
    fn = "sgtbx/direct_space_asu/" + f
    warn_if_unexpected_md5_hexdigest(
      path=self.env.under_dist( module_name="cctbx", path=fn),
      expected_md5_hexdigests=[ sig ],
      hints=[
        "  Files to review:",
        "    "+fn,
        "    cctbx/libtbx_refresh.py"])
