import libtbx.config
import libtbx.path
import shutil
import time
import sys, os

def copy_tag_file(src_root, target_root):
  if (os.path.isdir(target_root)):
    src = libtbx.path.norm_join(src_root, "TAG")
    if (os.path.isfile(src)):
      dest = libtbx.path.norm_join(target_root, "TAG")
      shutil.copy(src, dest)

def write_bundle_info(src_root, target_root):
  copy_tag_file(src_root, target_root)
  if (os.path.isdir(target_root)):
    dest = libtbx.path.norm_join(target_root, "bundle_info")
    f = open(dest, "w")
    print >> f, "%04d/%02d/%02d %02d:%02d:%02d" % time.localtime()[:6]
    print >> f, "time zone:", time.tzname
    print >> f, "hostname:", libtbx.config.get_hostname()
    print >> f, "os.name:", os.name
    print >> f, "src_root:", src_root
    print >> f, "sys.platform:", sys.platform
    print >> f, "sys.executable:", sys.executable
    print >> f, "sys.version:", sys.version
    f.close()
