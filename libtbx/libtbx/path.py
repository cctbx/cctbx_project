import os

def norm_join(*args):
  return os.path.normpath(apply(os.path.join, args))

def create_target_dir(target_file):
  target_dir = os.path.split(target_file)[0]
  if (not os.path.isdir(target_dir)):
    os.makedirs(target_dir)
