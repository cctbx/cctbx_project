from scitbx.source_generators.array_family import generate_all
from scitbx.source_generators import flex_fwd_h
import os

if (self.env.is_ready_for_build()):
  message_template = '  Generating C++ files in:\n    "%s"'

  # array_family
  target_dir = self.env.under_build("include/scitbx/array_family/detail")
  if (not os.path.isdir(target_dir)):
    os.makedirs(target_dir)
  generate_all.refresh(array_family=os.path.dirname(target_dir))

  # flex_fwd.h
  target_dir = self.env.under_build("include/scitbx/array_family/boost_python")
  print message_template % target_dir
  if not os.path.isdir(target_dir):
    os.makedirs(target_dir)
  flex_fwd_h.run(target_dir)
