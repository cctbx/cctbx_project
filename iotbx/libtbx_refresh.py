import os
op = os.path

self.remove_obsolete_pyc_if_possible(pyc_file_names=[
  "mtz/extract_from_symop_lib.pyc", # XXX backward compatibility 2010-10-29
])

if (self.env.is_ready_for_build()):
  f = self.env.under_dist("iotbx", path="pdb/hybrid_36_f.f")
  from libtbx.path import tail_levels
  print "  Using fable to convert", tail_levels(f, 3)
  import fable.cout
  cpp_lines = fable.cout.process(file_names=[f], fem_do_safe=False)
  d = self.env.under_build("iotbx/pdb")
  if (not op.isdir(d)):
    os.makedirs(d)
  t = op.join(d, "hybrid_36_fem.cpp")
  print "    Writing:", tail_levels(t, 3)
  print >> open(t, "w"), "\n".join(cpp_lines)
