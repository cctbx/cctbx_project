from __future__ import division
import sys, os
op = os.path

def eval_logs(file_names, out=None):
  if (out is None): out = sys.stdout
  from scitbx.array_family import flex
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
  perm = flex.sort_permutation(gaps)
  gaps = gaps.select(perm)
  print >> out, "Number of results:", gaps.size()
  print >> out, "Stale, Exceptions, Tracebacks, Abort:", \
    n_stale, n_exception, n_traceback, n_abort
  if (len(seconds) != 0):
    print >> out, "min, max seconds:", min(seconds), max(seconds)
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

def run(args):
  file_names = []
  dir_names = []
  for arg in args:
    if (op.isfile(arg)):
      file_names.append(arg)
    elif (op.isdir(arg)):
      dir_names.append(arg)
  assert len(file_names) == 0 or len(dir_names) == 0
  if (len(file_names) != 0):
    eval_logs(file_names)
    return
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
      eval_logs(file_names, out=open(outfn, "w"))

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
