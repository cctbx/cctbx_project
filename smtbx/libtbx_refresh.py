from __future__ import absolute_import, division, print_function
import os

if self.env.is_ready_for_build():
  message_template = '  Generating C++ files in:\n    "%s"'
  #flex_fwd.h
  from smtbx.source_generators import flex_fwd_h
  target_dir = self.env.under_build("include/smtbx/boost_python")
  print(message_template % target_dir)
  if not os.path.isdir(target_dir):
    os.makedirs(target_dir)
  flex_fwd_h.run(target_dir)
