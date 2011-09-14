def raise_usage():
  from libtbx.utils import Usage
  import libtbx.load_env
  raise Usage("%s 1A|1keV [...]" % libtbx.env.dispatcher_name)

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
    print "%.5f A = %.5f keV" % (a, k)

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
