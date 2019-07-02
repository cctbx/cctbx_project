from __future__ import absolute_import, division, print_function
import libtbx.load_env
import os

def run():
  version = None
  for tag_file in ["TAG", "cctbx_bundle_TAG"]:
    tag_path = libtbx.env.under_dist("libtbx", os.path.join("..", tag_file))
    if (os.path.isfile(tag_path)):
      try: version = open(tag_path).read().strip()
      except KeyboardInterrupt: raise
      except Exception: pass
      else: break
  if (version is None):
    version = libtbx.env.command_version_suffix
  print(version)

if (__name__ == "__main__"):
  run()
