# svn co https://cbflib.svn.sourceforge.net/svnroot/cbflib/trunk/CBFlib_bleeding_edge sourceforge_cbflib

def run(args):
  assert len(args) == 1, "<path>/sourceforge_cbflib"
  sourceforge_cbflib = args[0]
  from shutil import copyfile
  import os
  op = os.path
  assert op.isdir(sourceforge_cbflib)
  def copy_from_directory(dname, fnames=None, h_c_only=False):
    if (not op.isdir(dname)):
      os.mkdir(dname)
    dpath = op.join(sourceforge_cbflib, dname)
    if (fnames is None): fnames = os.listdir(dpath)
    for fname in fnames:
      if (not h_c_only or fname.endswith(".h") or fname.endswith(".c")):
        copyfile(src=op.join(dpath, fname), dst=op.join(dname, fname))
  for dname in ["include", "src"]:
    copy_from_directory(dname, h_c_only=True)
  copy_from_directory("pycbf", ["pycbf_wrap.c", "pycbf.py"])
  copy_from_directory("examples", ["img.h", "img.c", "fit2d_data.cbf"])
  copy_from_directory("doc", ["lgpl.txt"])
  fnames = ["README"]
  if (op.isfile(op.join(sourceforge_cbflib, "TAG"))):
    fnames.append("TAG")
  copy_from_directory(".", fnames)

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
