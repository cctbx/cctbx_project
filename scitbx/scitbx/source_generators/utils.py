import os.path

def norm_join(path1, path2):
  return os.path.normpath(os.path.join(path1, path2))

def join_open(path1, path2, mode):
  path = norm_join(path1, path2)
  if (mode == "w"):
    print "Generating:", path
  return open(path, mode)
