#! /usr/bin/env python

import sys, os

def norm(path):
  return os.path.normpath(os.path.normcase(path))

def run(args):
  assert len(args) == 3
  mode = args[0]
  assert mode in ("prepend", "append", "delete")
  env_key = args[1]
  arg_paths = args[2].split(os.pathsep)
  try: env_val = os.environ[env_key]
  except: env_val = ""
  env_paths = env_val.split(os.pathsep)
  remaining_env_paths = []
  if (mode == "delete"):
    arg_paths_norm = []
    for path in arg_paths:
      if (path == ""): continue
      arg_paths_norm.append(norm(path))
    for path in env_paths:
      if (path == ""): continue
      if (not norm(path) in arg_paths_norm):
        remaining_env_paths.append(path)
  else:
    remaining_env_paths_norm = []
    if (mode == "prepend"):
      all_paths = arg_paths + env_paths
    else:
      all_paths = env_paths + arg_paths
    for path in all_paths:
      if (path == ""): continue
      if (not norm(path) in remaining_env_paths_norm):
        remaining_env_paths.append(path)
        remaining_env_paths_norm.append(norm(path))
  return os.pathsep.join(remaining_env_paths)

if (__name__ == "__main__"):
  print run(sys.argv[1:])
