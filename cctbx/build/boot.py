import sys, os, os.path, string

def expand_cwd(cf):
  cwd = os.getcwd()
  for i in xrange(len(cf)):
    c = string.find(cf[i], "@(CWD)")
    if (c >= 0):
      cf[i] = cf[i][:c] + cwd + cf[i][c+6:]
      e = string.find(cf[i], "=")
      if (e >= 0):
        e = e + 1
        cf[i] = cf[i][:e] + os.path.normpath(cf[i][e:])

def read_configuration():
  f = open("configuration", "r")
  cf = f.readlines()
  f.close()
  for i in xrange(len(cf)): cf[i] = cf[i][:-1]
  expand_cwd(cf)
  return cf

def create_makefiles(path_cctbx):
  configuration = read_configuration()
  for subdir in ("eltbx", "sgtbx", "uctbx", "examples/cpp"):
    exe_g = {}
    exe_l = {}
    execfile(path_cctbx + "/" + subdir + "/MakeMakefiles.py", exe_g, exe_l)
    try:
      os.makedirs(subdir)
    except OSError:
      pass
    h = exe_l['write_makefiles'](configuration)
    f = open(subdir + "/Makefile", "w")
    h.write(f)
    f.close()

if (__name__ == "__main__"):
  path_cctbx = os.path.dirname(os.path.dirname(sys.argv[0]))
  sys.path.insert(0, path_cctbx + "/build")
  create_makefiles(path_cctbx)
