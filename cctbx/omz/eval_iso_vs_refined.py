from __future__ import division

def run(args):
  from scitbx.array_family import flex
  gaps = flex.double()
  infos = flex.std_string()
  n_stale = 0
  n_exception = 0
  n_traceback = 0
  n_abort = 0
  seconds = []
  space_groups_by_cod_code = {}
  for file_name in args:
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
  print "Number of results:", gaps.size()
  print "Stale, Exceptions, Tracebacks, Abort:", \
    n_stale, n_exception, n_traceback, n_abort
  if (len(seconds) != 0):
    print "min, max seconds:", min(seconds), max(seconds)
  print
  def stats(f):
    n = f.count(True)
    return "%6d = %5.2f %%" % (n, 100 * n / max(1,gaps.size()))
  print "gaps below -0.05:", stats(gaps < -0.05)
  print "gaps below -0.01:", stats(gaps < -0.01)
  print "gaps below  0.01:", stats(gaps <  0.01)
  print "gaps above  0.01:", stats(gaps >  0.01)
  print "gaps above  0.05:", stats(gaps >  0.05)
  print
  print "Histogram of gaps:"
  flex.histogram(gaps, n_slots=10).show()
  print
  infos = infos.select(perm)
  for info in infos:
    print info

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
