from scitbx.source_generators.array_family import generate_all
import os

if (self.env.is_ready_for_build()):
  target_dir = self.env.under_build("include/scitbx/array_family/detail")
  if (not os.path.isdir(target_dir)):
    os.makedirs(target_dir)
  generate_all.refresh(array_family=os.path.dirname(target_dir))
