from __future__ import absolute_import, division, print_function
import libtbx.load_env
import libtbx.env_config
import libtbx.path
import libtbx.utils
import sys, os

def copy_tags(target_root):
  if (os.path.isdir(target_root)):
    flds = []
    for path in libtbx.env.repository_paths:
      src = libtbx.path.norm_join(abs(path), "TAG")
      if (os.path.isfile(src)):
        flds.extend(open(src).read().split())
    if (len(flds) > 0):
      dest = libtbx.path.norm_join(target_root, "TAG")
      print(" ".join(flds), file=open(dest, "w"))

def write_bundle_info(target_root, write_build_options=False):
  copy_tags(target_root)
  if (os.path.isdir(target_root)):
    dest = libtbx.path.norm_join(target_root, "bundle_info")
    f = open(dest, "w")
    print("date_and_time:", libtbx.utils.date_and_time(), file=f)
    print("hostname:", libtbx.env_config.get_hostname(), file=f)
    print("os.name:", os.name, file=f)
    print("sys.platform:", sys.platform, file=f)
    print("sys.executable:", sys.executable, file=f)
    print("sys.version:", " ".join(sys.version.splitlines()), file=f)
    print("repository_paths:", libtbx.env.repository_paths, file=f)
    if (write_build_options):
      libtbx.env.build_options.report(f=f)
    f.close()
