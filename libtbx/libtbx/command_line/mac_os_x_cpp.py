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
    return True
  print "Trying compilation of preprocessor output in intermediate file."
  sys.stdout.flush()
  ppout = out[:-3] + "_pp.cpp"
  cmd = " ".join(["g++"] + argv[1:-4] + ["-E", argv[-1], ">", ppout])
  os.system(cmd)
  sys.stderr.flush()
  sys.stdout.flush()
  if (not os.path.isfile(ppout)):
    return False
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
    return True
  return False

def disable_optimization(argv):
  argv_mod = []
  for arg in argv:
    if (arg[:2] == "-O"):
      try: i = int(arg[2:])
      except: pass
      else:
        if (i > 0):
          arg = "-O0"
    argv_mod.append(arg)
  return argv_mod

def exit_fail():
  print "All attempts failed."
  sys.exit(1)

if (not compile(sys.argv)):
  argv_mod = disable_optimization(sys.argv)
  if (argv_mod == sys.argv):
    exit_fail()
  print "Trying without optimization."
  sys.stdout.flush()
  if (not compile(argv_mod)):
    exit_fail()
