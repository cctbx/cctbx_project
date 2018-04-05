from __future__ import absolute_import, division, print_function

import os
import sys

def norm(path):
  return os.path.normcase(os.path.normpath(path))

def run(args):
  assert len(args) >= 2
  mode = args[0]
  env_key = args[1]
  if mode in ("prepend", "append", "delete"):
    assert len(args) == 3
    if args[2] == "__CWD__":
      arg_paths = [os.getcwd()]
    else:
      arg_paths = args[2].split(os.pathsep)
  elif mode == "tidy":
    assert len(args) == 2
    arg_paths = []
  else:
    raise RuntimeError('Unknown mode: "%s"' % mode)
  env_val = os.environ.get(env_key, "")
  env_paths = env_val.split(os.pathsep)
  if os.name == "nt":
    unquoted_paths = []
    for path in env_paths:
      if len(path) >= 2 and path[:1] == '"' and path[-1:] == '"':
        path = path[1:-1]
      unquoted_paths.append(path)
    env_paths = unquoted_paths
  remaining_env_paths = []
  if mode == "delete":
    arg_paths_norm = []
    for path in arg_paths:
      if path == "": continue
      arg_paths_norm.append(norm(path))
    for path in env_paths:
      if path == "": continue
      if norm(path) not in arg_paths_norm:
        remaining_env_paths.append(path)
  else:
    remaining_env_paths_norm = []
    if mode == "prepend":
      all_paths = arg_paths + env_paths
    else:
      all_paths = env_paths + arg_paths
    for path in all_paths:
      if path == "": continue
      if norm(path) not in remaining_env_paths_norm:
        remaining_env_paths.append(path)
        remaining_env_paths_norm.append(norm(path))
  if len(remaining_env_paths) == 0:
    return "L_I_B_T_B_X_E_M_P_T_Y"
  return os.pathsep.join(remaining_env_paths)

if __name__ == "__main__":
  print(run(sys.argv[1:]))
