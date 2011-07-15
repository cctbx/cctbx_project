# svn co https://cbflib.svn.sourceforge.net/svnroot/cbflib/trunk/CBFlib_bleeding_edge sourceforge_cbflib

class info_counters(object):
  def __init__(O):
    O.mkdir = 0
    O.copied = 0
    O.updated = 0
    O.already_up_to_date = 0
  def report(O):
    print "Directories created:", O.mkdir
    print "Files copied:", O.copied
    print "Files updated:", O.updated
    print "Files already up-to-date:", O.already_up_to_date

def run(args):
  assert len(args) == 1, "<path>/sourceforge_cbflib"
  sourceforge_cbflib = args[0]
  from shutil import copyfile
  import os
  op = os.path
  assert op.isdir(sourceforge_cbflib)
  counters = info_counters()
  def copy_from_directory(dname, fnames=None, h_c_only=False):
    if (not op.isdir(dname)):
      os.mkdir(dname)
      counters.mkdir += 1
    dpath = op.join(sourceforge_cbflib, dname)
    if (fnames is None): fnames = os.listdir(dpath)
    for fname in fnames:
      if (not h_c_only or fname.endswith(".h") or fname.endswith(".c")):
        src = op.join(dpath, fname)
        dst = op.join(dname, fname)
        src_bytes = open(src, "rb").read()
        if (not op.isfile(dst)):
          counters.copied += 1
        else:
          dst_bytes = open(dst, "rb").read()
          if (dst_bytes == src_bytes):
            counters.already_up_to_date += 1
            src = None
          else:
            counters.updated += 1
        if (src is not None):
          copyfile(src=src, dst=dst)
  for dname in ["include", "src"]:
    copy_from_directory(dname, h_c_only=True)
  copy_from_directory("src", ["cbf.stx.y"])
  copy_from_directory("pycbf", ["pycbf_wrap.c", "pycbf.py"])
  copy_from_directory("examples", ["img.h", "img.c", "fit2d_data.cbf"])
  copy_from_directory("doc", ["lgpl.txt"])
  fnames = ["README"]
  if (op.isfile(op.join(sourceforge_cbflib, "TAG"))):
    fnames.append("TAG")
  copy_from_directory(".", fnames)
  counters.report()
  print "Done."

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
