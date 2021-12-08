# LIBTBX_SET_DISPATCHER_NAME cctbx.make_install_HKLviewer_ChimeraX

from __future__ import absolute_import, division, print_function

import subprocess, sys, shutil, os
from os.path import join as opj
import libtbx.load_env

# Bundle up a python wheel of the HKLviewer for chimeraX and subsequently install it into the local
# chimeraX installation
# Usage:
# cctbx.make_install_HKLviewer_ChimeraX "A:\Program Files\ChimeraX 1.2.5\bin\ChimeraX-console.exe"

if (__name__ == "__main__"):
  print("Bundling up and installing a ChimeraX wheel of cctbx.HKLviewer")
  if len(sys.argv) < 2:
    print("Specify the full path to the chimeraX commandline compiler as argument")
    exit()

  builddir = libtbx.env.under_root(os.path.join("build","ChimeraX_tools","HKLviewer","src"))
  wheeldir = libtbx.env.under_root(os.path.join("build","ChimeraX_tools","HKLviewer"))
  distdir = libtbx.env.under_dist("crys3d","hklviewer")

  shutil.copytree(opj(distdir,"chimeraX_wheel_src"), builddir, dirs_exist_ok=True)
  shutil.copy(opj(builddir,"bundle_info.xml"), wheeldir )
  os.remove(opj(builddir,"bundle_info.xml"))

  shutil.copy(opj(distdir,"helpers.py"), builddir)
  shutil.copy(opj(distdir,"HKLviewer.py"), builddir)
  shutil.copy(opj(distdir,"HKLviewerGui.py"), builddir)
  shutil.copy(opj(distdir,"qt.py"), builddir)
  shutil.copy(opj(distdir,"QtChromiumCheck.py"), builddir)

  cmdargs =  r'"%s" --nogui --cmd "devel install . ; exit"' %sys.argv[1]
  curdir = os.getcwd()
  os.chdir(wheeldir)
  buildproc = subprocess.Popen( cmdargs, shell=True, cwd=wheeldir,
                                    universal_newlines=True,
                                    stdin=subprocess.PIPE,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
  out, err = buildproc.communicate()
  print(out)
  print(err)
  os.chdir(curdir)
