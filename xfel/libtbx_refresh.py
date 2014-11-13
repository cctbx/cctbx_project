from __future__ import division
import libtbx.load_env, os
cxi_user = libtbx.env.find_in_repositories(
             relative_path="cxi_user",
             test=os.path.isdir)

if cxi_user is None or not os.path.exists(cxi_user):
  print "  Creating cxi_user directory"

  sources_root = libtbx.env.find_in_repositories(
                   relative_path=".",
                   test=os.path.isdir)

  cxi_user = os.path.join(sources_root, "cxi_user")
  os.mkdir(cxi_user)

init = os.path.join(cxi_user, "__init__.py")
if not os.path.exists(init):
  print "  Creating cxi_user/__init__.py"
  with open(init, "w") as f:
    f.write("from xfel.mono_simulation.mono_treatment import post_outlier_rejection\n")
    f.write("from xfel.mono_simulation.mono_treatment import pre_get_predictions\n")
