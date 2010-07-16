import os
op = os.path

if (self.env.is_ready_for_build()):
  from scitbx.source_generators.array_family import generate_all
  target_dir = self.env.under_build("include/scitbx/array_family/detail")
  if (not op.isdir(target_dir)):
    os.makedirs(target_dir)
  generate_all.refresh(array_family=op.dirname(target_dir))

  from scitbx.source_generators import flex_fwd_h
  from libtbx.str_utils import show_string
  target_dir = self.env.under_build("include/scitbx/array_family/boost_python")
  print "  Generating C++ files in:\n    %s" % show_string(target_dir)
  if not op.isdir(target_dir):
    os.makedirs(target_dir)
  flex_fwd_h.run(target_dir)

  from scitbx.source_generators import lbfgs_fem
  print "  Using fable to convert", op.join("scitbx", "lbfgs.f")
  lbfgs_fem.run()
