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
      print >> open(dest, "w"), " ".join(flds)

def write_bundle_info(target_root, write_build_options=False):
  copy_tags(target_root)
  if (os.path.isdir(target_root)):
    dest = libtbx.path.norm_join(target_root, "bundle_info")
    f = open(dest, "w")
    print >> f, "date_and_time:", libtbx.utils.date_and_time()
    print >> f, "hostname:", libtbx.env_config.get_hostname()
    print >> f, "os.name:", os.name
    print >> f, "sys.platform:", sys.platform
    print >> f, "sys.executable:", sys.executable
    print >> f, "sys.version:", " ".join(sys.version.splitlines())
    print >> f, "repository_paths:", libtbx.env.repository_paths
    if (write_build_options):
      libtbx.env.build_options.report(f=f)
    f.close()
