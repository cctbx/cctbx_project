prios = ["boost", "scitbx", "cctbx", "rstbx", "labelit"]
prios = dict(zip(prios, range(len(prios))))

def cmp_so(a, b):
  def prio(s):
    k = s.split("_")[0]
    return prios.get(k, len(prios))
  pa, pb = prio(a), prio(b)
  result = cmp(pa, pb)
  if (result == 0):
    result = cmp(a, b)
  return result

def run(args):
  assert len(args) == 0
  import time
  t_start = time.time()
  from libtbx import introspection
  import os
  print "After script imports:"
  print "  wall clock time: %.2f" % (time.time() - t_start)
  print
  mb = 1024 * 1024
  lib = os.path.join(
    os.environ["LIBTBX_BUILD"], # intentionally not using libtbx.env for speed
    "lib")
  ext_so = []
  for node in os.listdir(lib):
    if (node.endswith("_ext.so")):
      ext_so.append(node[:-3])
  ext_so.sort(cmp_so)
  print "Before importing extensios:"
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
    print "%.2f %3.0f %3.0f %s" % (
      time.time()-t0, (vms-prev_vms)/mb, (rss-prev_rss)/mb, so)
    prev_vms = vms
    prev_rss = rss
  print
  print "After importing all extensios:"
  introspection.virtual_memory_info().show(prefix="  ")
  print "  wall clock time: %.2f" % (time.time() - t_start)
  print

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
