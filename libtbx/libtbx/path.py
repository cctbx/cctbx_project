import os

def norm_join(*args):
  return os.path.normpath(apply(os.path.join, args))
