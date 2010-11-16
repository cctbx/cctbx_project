import os

def norm_join(path1, path2):
  return os.path.normpath(os.path.join(path1, path2))

def join_open(path1, path2, mode, verbose=0):
  path = norm_join(path1, path2)
  if (mode == "w" and verbose):
    print "     ", path2
  return open(path, mode)
