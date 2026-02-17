
from __future__ import absolute_import, division, print_function
import os

pdbmtz = os.getenv('PDBMTZ')

def process(f):
  with open(f,"r") as fo:
    rls = fo.readlines()
    full = ''.join(rls)
    result = None
    for l in rls:
      if l.startswith("  File "):
        result = l
  example = "******************* EXAMPLE: %s\n"%f
  return result, example+full

def run():
  epath = pdbmtz.replace("/mtz_files", "/errors/")
  d = {}
  cntr = 0
  for f in os.listdir(epath):
    if not f.endswith(".log"): continue
    cntr+=1
    #if cntr==20: break
    key, full = process(epath+f)
    d[key]=full
  #
  print("TOTAL UNIQUE ERRORS:", len(d.keys()))
  print()
  for k, v in d.items():
    print(v)

if (__name__ == "__main__"):
  run()
