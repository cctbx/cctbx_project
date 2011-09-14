def iobs_as_fobs(iobs, isigma):
  import math
  if (iobs >= 0):
    fobs = math.sqrt(iobs)
    if (isigma < iobs):
      fsigma = fobs - math.sqrt(iobs - isigma)
    else:
      fsigma = fobs
  else:
    fobs = 0
    fsigma = 0
  return fobs, fsigma

def dtrek_as_cns_hkl(file_object, file_name=None):
  info = []
  reflections = []
  mode = 0
  for line in file_object.xreadlines():
    if (line == ""): break
    if (mode == 0):
      if (line.startswith("CRYSTAL_")):
        info.append(line[:-1])
      if (not line.startswith(" ")): continue
      mode = 1
    flds = line.split()
    iobs = float(flds[3])
    isigma = float(flds[4])
    if (iobs < 0 or isigma < 0): continue
    fobs, sigma = iobs_as_fobs(iobs, isigma)
    reflections.append("INDEX %s %s %s FOBS %.6g SIGMA %.6g" % (
      flds[0], flds[1], flds[2], fobs, sigma))
  if (file_name): print "{ file:", file_name, "}"
  for line in info: print "{", line, "}"
  print "NREFlections=%d" % (len(reflections),)
  print "ANOMalous=FALSe"
  print "DECLare NAME=FOBS  DOMAin=RECIprocal TYPE=REAL END"
  print "DECLare NAME=SIGMA DOMAin=RECIprocal TYPE=REAL END"
  for line in reflections: print line

if (__name__ == "__main__"):
  import sys, os.path
  if (len(sys.argv) == 1):
    dtrek_as_cns_hkl(sys.stdin)
  elif (len(sys.argv) == 2):
    f = open(sys.argv[1], "r")
    dtrek_as_cns_hkl(f, os.path.abspath(sys.argv[1]))
    f.close()
  else:
    raise RuntimeError, (
      "usage: %s [d*trek_file_name]" % (os.path.basename(sys.argv[0]),))
