from __future__ import absolute_import, division, print_function
import libtbx.load_env
import sys, os
from libtbx.utils import Usage

"""
Search libtbx_config files to show the complete dependency
chain of any cctbx modules. Optional modules are skipped,
as are modules tagged with 'exclude_from_binary_bundle'.
"""

def run():
  if len(sys.argv) != 2:
    raise Usage("%s modulename"%libtbx.env.dispatcher_name)

  initial_root = sys.argv[1]

  # save a list of all modules found so far so as to not enter loops
  all_found = []

  def show_dependencies(root, depth):
    """
    Recursivly search known modules for dependencies
    @param root module to start with
    @param depth recursion depth
    """
    indent = "  " * depth
    if root in all_found:
      # dependency already found
      return
    all_found.append(root)

    root_path = libtbx.env.find_in_repositories(root)
    if root_path is None or not os.path.exists(root_path):
      print(indent, "Dependency not found", root)
      return

    config_path = os.path.join(root_path, "libtbx_config")
    if not os.path.exists(config_path):
      print(indent, root, "has no libtbx_config")
      return

    f = open(config_path)
    config = eval(" ".join(f.readlines()), {}, {})
    f.close()

    all_mods = []
    for key, mods in config.items():
      if "optional" in key: continue
      elif key == "exclude_from_binary_bundle": continue
      else:
        for mod in mods:
          if mod not in all_mods:
            all_mods.append(mod)

    if len(all_mods) == 0:
      print(indent, root, "has no dependencies")
    else:
      print(indent, root, "depends on", " ".join(all_mods))
      for mod in all_mods:
        show_dependencies(mod, depth+1)

  show_dependencies(initial_root, 0)

if (__name__ == "__main__"):
  run()
