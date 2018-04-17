from __future__ import division
prios = ["boost", "scitbx", "cctbx", "rstbx", "labelit"]
prios = dict(zip(prios, range(len(prios))))
priority_modules = ['dials_algorithms_profile_model_modeller_ext']

def cmp_so(a, b):
  if a in priority_modules:
    return -1
  elif b in priority_modules:
    return 1

  def prio(s):
    k = s.split("_")[0]
    return prios.get(k, len(prios))
  pa, pb = prio(a), prio(b)
  result = cmp(pa, pb)
  if (result == 0):
    result = cmp(a, b)
  return result

def import_modules():
  # Search for __init__.py file and import those modules. This ensures all dependencies are loaded
  # prior to import individual so files
  import libtbx.load_env, os, traceback
  # Modules in the form of X/X/__init__.py
  doubled = ["elbow","phaser","phenix"]
  # These modules are not dependencies not maintained by cctbx or should be skipped for other reasons
  skip_modules = ["boost","cbflib","crys3d","PyQuante","phenix_html","tntbx","reel",
                  "amber", # external
                  ]
  # These modules fail to import
  failing_modules = ["dials.framework.dftbx","dials.nexus"]
  # These modules can only import if Eigen is available
  ok_to_fail = ["scitbx.examples.bevington","cctbx.examples.merging","cctbx.examples.merging.samosa"]

  for root_module in libtbx.env.module_list:
    if root_module.name in skip_modules:
      continue
    root_path = libtbx.env.find_in_repositories(root_module.name)
    if root_path is None:
      continue
    if root_module.name in doubled:
      root_path = os.path.join(root_path, root_module.name)
    for dirpath, dirnames, filenames in os.walk(root_path):
      if "__init__.py" in filenames:
        full_module = root_module.name + ".".join(dirpath.split(root_path)[-1].split(os.path.sep))
        if full_module in failing_modules:
          continue
        print full_module
        try:
          exec("import %s" % full_module)
        except ImportError as e:
          if full_module not in ok_to_fail:
            print traceback.format_exc()
            raise e

def run(args):
  assert len(args) == 0
  import_modules()
  import time
  t_start = time.time()
  from libtbx import introspection
  import os
  import sysconfig
  print "After script imports:"
  print "  wall clock time: %.2f" % (time.time() - t_start)
  print
  mb = 1024 * 1024
  lib = os.path.join(
    os.environ["LIBTBX_BUILD"], # intentionally not using libtbx.env for speed
    "lib")
  ext_so = []
  for node in os.listdir(lib):
    pylibext = sysconfig.get_config_vars("SO")[0]
    if (node.endswith("_ext" + pylibext)):
      ext_so.append(node.split(".")[0])
  ext_so.sort(cmp_so)
  print "Before importing extensions:"
  vmi = introspection.virtual_memory_info()
  vmi.show(prefix="  ")
  prev_vms = vmi.get_bytes('VmSize:')
  prev_rss = vmi.get_bytes('VmRSS:')
  print "  wall clock time: %.2f" % (time.time() - t_start)
  print
  for so in ext_so:
    t0 = time.time()
    exec("import %s" % so)
    vmi = introspection.virtual_memory_info()
    vms = vmi.get_bytes('VmSize:')
    rss = vmi.get_bytes('VmRSS:')
    if (vms is not None) : # won't work on Mac
      print "%.2f %3.0f %3.0f %s" % (
        time.time()-t0, (vms-prev_vms)/mb, (rss-prev_rss)/mb, so)
    else :
      assert (sys.platform in ["darwin", "win32"])
      print "%.2f %s" % (time.time()-t0, so)
    prev_vms = vms
    prev_rss = rss
  print
  print "After importing all extensions:"
  introspection.virtual_memory_info().show(prefix="  ")
  print "  wall clock time: %.2f" % (time.time() - t_start)
  print

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
