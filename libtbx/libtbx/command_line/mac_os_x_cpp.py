#! /usr/bin/env python

import sys, os

def compile(argv):
  assert argv[-4] == "-c"
  assert argv[-3] == "-o"
  out = argv[-2]
  assert out.endswith(".os")
  cmd = " ".join(["g++"] + argv[1:])
  os.system(cmd)
  sys.stderr.flush()
  sys.stdout.flush()
  if (os.path.isfile(out)):
    print "One-phase success."
    return True
  ppout = out[:-3] + "_pp.cpp"
  cmd = " ".join(["g++"] + argv[1:-4] + ["-E", argv[-1], ">", ppout])
  os.system(cmd)
  sys.stderr.flush()
  sys.stdout.flush()
  args = ["g++"]
  for arg in argv[1:-4]:
    if (not arg.startswith("-D") and not arg.startswith("-I")):
      args.append(arg)
  args.extend(["-c", "-o", out, ppout])
  cmd = " ".join(args)
  os.system(cmd)
  sys.stderr.flush()
  sys.stdout.flush()
  if (os.path.isfile(out)):
    print "Two-phase success."
    return True
  return False

if (not compile(sys.argv)):
  try: i = sys.argv.index("-O3")
  except: sys.exit(1)
  print "Compiling without optimization."
  sys.stdout.flush()
  argv = sys.argv[:i] + ["-O0"] + sys.argv[i+1:]
  if (not compile(argv)):
    print "All attempts failed."
    sys.exit(1)
