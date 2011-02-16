from __future__ import division

def run(args):
  from scitbx.array_family import flex
  gaps = flex.double()
  infos = flex.std_string()
  for file_name in args:
    cod_code = None
    n_scatt = None
    iso = None
    for line in open(file_name).read().splitlines():
      if (line.startswith("cod_code: ")):
        cod_code = line[10:]
        iso = None
      elif (line.startswith("Number of scatterers: ")):
        assert cod_code is not None
        n_scatt = int(line.split(": ",1)[1])
      elif (line.startswith("iso          cc, r1: ")):
        assert cod_code is not None
        assert iso is None
        iso = line.split(": ",1)[1]
      elif (   line.startswith("ls12         cc, r1: ")
            or line.startswith("dev          cc, r1: ")):
        assert iso is not None
        ref = line.split(": ",1)[1]
        gap = float(ref.split()[1]) - float(iso.split()[1])
        gaps.append(gap)
        infos.append(" : ".join([
          cod_code, iso, ref, "%.3f" % gap, str(n_scatt)]))
        cod_code = None
        n_scatt = None
        iso = None
  perm = flex.sort_permutation(gaps)
  gaps = gaps.select(perm)
  print "Number of results:", gaps.size()
  print
  def stats(f):
    n = f.count(True)
    return "%6d = %5.2f %%" % (n, 100 * n / gaps.size())
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
