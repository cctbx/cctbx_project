#! /usr/bin/env python
import sys, os
assert sys.argv[-4] == "-c"
assert sys.argv[-3] == "-o"
out = sys.argv[-2]
assert out.endswith(".os")
ppout = out[:-3] + "_pp.cpp"
cmd = " ".join(["g++"] + sys.argv[1:-4] + ["-E", sys.argv[-1], ">", ppout])
os.system(cmd)
args = ["g++"]
for arg in sys.argv[1:-4]:
  if (not arg.startswith("-D") and not arg.startswith("-I")):
    args.append(arg)
args.extend(["-c", "-o", out, ppout])
cmd = " ".join(args)
os.system(cmd)
