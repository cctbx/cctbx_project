from __future__ import absolute_import, division, print_function
def raise_usage():
  from libtbx.utils import Usage
  raise Usage("cctbx.wavelength_units 1A|1keV [...]")

def run(args):
  if (len(args) == 0): raise_usage()
  from cctbx import factor_kev_angstrom
  for arg in args:
    l = arg.lower()
    if (l.endswith("a")):
      try: a = eval(arg[:-1])
      except Exception: raise_usage()
      k = factor_kev_angstrom / a
    elif (l.endswith("kev")):
      try: k = eval(arg[:-3])
      except Exception: raise_usage()
      a = factor_kev_angstrom / k
    else:
      raise_usage()
    print("%.5f A = %.5f keV" % (a, k))

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
