from __future__ import division
import sys, os
op = os.path

def eval_logs(file_names, out=None):
  if (out is None): out = sys.stdout
  from scitbx.array_family import flex
  from libtbx.str_utils import format_value
  min_secs_epoch = None
  max_secs_epoch = None
  gaps = flex.double()
  infos = flex.std_string()
  n_stale = 0
  n_exception = 0
  n_traceback = 0
  n_abort = 0
  seconds = []
  space_groups_by_cod_code = {}
  for file_name in file_names:
    cod_code = None
    n_scatt = None
    iso = None
    file_str = open(file_name).read()
    if (file_str.find(chr(0)) >= 0):
      n_stale += 1
      continue
    for line in file_str.splitlines():
      if (line.startswith("cod_code: ")):
        cod_code = line[10:]
        iso = None
      elif (line.startswith("Space group: ")):
        assert cod_code is not None
        space_group = line.split(None, 2)[2]
        tabulated = space_groups_by_cod_code.setdefault(cod_code, space_group)
        assert tabulated == space_group
      elif (line.startswith("Number of scatterers: ")):
        assert cod_code is not None
        n_scatt = int(line.split(": ",1)[1])
      elif (line.startswith("iso          cc, r1: ")):
        assert cod_code is not None
        assert iso is None
        iso = line.split(": ",1)[1]
      elif (   line.startswith("dev          cc, r1: ")
            or line.startswith("ls_simple    cc, r1: ")
            or line.startswith("ls_lm        cc, r1: ")
            or line.startswith("shelxl_fm    cc, r1: ")
            or line.startswith("shelxl_cg    cc, r1: ")
            or line.startswith("shelx76      cc, r1: ")):
        assert iso is not None
        ref = line.split(": ",1)[1]
        gap = float(ref.split()[1]) - float(iso.split()[1])
        gaps.append(gap)
        infos.append(" : ".join([
          cod_code, iso, ref, "%.3f" % gap, str(n_scatt),
          space_groups_by_cod_code[cod_code]]))
        cod_code = None
        n_scatt = None
        iso = None
      else:
        def get_secs_epoch():
          return float(line.split()[-2][1:])
        if (line.find("EXCEPTION") >= 0):
          n_exception += 1
        if (line.startswith("Traceback")):
          n_traceback += 1
        if (line.find("Abort") >= 0):
          n_abort += 1
        if (line.startswith("wall clock time: ")):
          if (line.endswith(" seconds")):
            secs = float(line.split()[-2])
          else:
            _, fld = line.split("(", 1)
            assert fld.endswith(" seconds total)")
            secs = float(fld.split()[0])
          seconds.append(secs)
        elif (line.startswith("TIME BEGIN cod_refine: ")):
          s = get_secs_epoch()
          if (min_secs_epoch is None or s < min_secs_epoch): min_secs_epoch = s
        elif (line.startswith("TIME END cod_refine: ")):
          s = get_secs_epoch()
          if (max_secs_epoch is None or s > max_secs_epoch): max_secs_epoch = s
  perm = flex.sort_permutation(gaps)
  gaps = gaps.select(perm)
  print >> out, "Number of results:", gaps.size()
  print >> out, "Stale, Exceptions, Tracebacks, Abort:", \
    n_stale, n_exception, n_traceback, n_abort
  if (len(seconds) != 0):
    if (min_secs_epoch is not None and max_secs_epoch is not None):
      g = max_secs_epoch - min_secs_epoch
    else:
      g = None
    print >> out, "min, max, global seconds: %.2f %.2f %s" % (
      min(seconds), max(seconds), format_value("%.2f", g))
  print >> out
  def stats(f):
    n = f.count(True)
    return "%6d = %5.2f %%" % (n, 100 * n / max(1,gaps.size()))
  print >> out, "gaps below -0.05:", stats(gaps < -0.05)
  print >> out, "gaps below -0.01:", stats(gaps < -0.01)
  print >> out, "gaps below  0.01:", stats(gaps <  0.01)
  print >> out, "gaps above  0.01:", stats(gaps >  0.01)
  print >> out, "gaps above  0.05:", stats(gaps >  0.05)
  print >> out
  print >> out, "Histogram of gaps:"
  flex.histogram(gaps, n_slots=10).show(f=out)
  print >> out
  infos = infos.select(perm)
  for info in infos:
    print >> out, info
  return (min_secs_epoch, max_secs_epoch)

def run(args):
  file_names = []
  dir_names = []
  for arg in args:
    if (op.isfile(arg)):
      file_names.append(arg)
    elif (op.isdir(arg)):
      dir_names.append(arg)
  assert len(file_names) == 0 or len(dir_names) == 0
  min_max_secs_epoch = [None, None]
  def track_times(min_max):
    a, b = min_max_secs_epoch[0], min_max[0]
    if (b is not None):
      if (a is None or a > b): min_max_secs_epoch[0] = b
    a, b = min_max_secs_epoch[1], min_max[1]
    if (b is not None):
      if (a is None or a < b): min_max_secs_epoch[1] = b
  if (len(file_names) != 0):
    track_times(eval_logs(file_names))
  else:
    for dir_name in dir_names:
      file_names = []
      for node in sorted(os.listdir(dir_name)):
        if (node.startswith("log")):
          path = op.join(dir_name, node)
          if (op.isfile(path)):
            file_names.append(path)
      if (len(file_names) != 0):
        outfn = op.join(dir_name, "stats")
        print outfn
        sys.stdout.flush()
        track_times(eval_logs(file_names, out=open(outfn, "w")))
  if (min_max_secs_epoch.count(None) == 0):
    print "global seconds: %.2f" % (
      min_max_secs_epoch[1] - min_max_secs_epoch[0])

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
