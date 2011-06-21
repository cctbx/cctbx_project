from libtbx.str_utils import show_string
import os
op = os.path
pycbf_py = self.env.under_dist(module_name="cbflib", path="pycbf/pycbf.py")
if (not op.isfile(pycbf_py)):
  raise RuntimeError("Missing file: %s" % show_string(pycbf_py))
target_dir = self.env.under_build(path="lib")
if (not op.isdir(target_dir)):
  os.makedirs(target_dir)
print "  Copying to lib: %s" % show_string(pycbf_py)
open(op.join(target_dir, "pycbf.py"), "w").write(open(pycbf_py).read())
