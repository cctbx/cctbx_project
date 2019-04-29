from __future__ import absolute_import, division, print_function
from libtbx.bundle import copy_all
from libtbx.bundle import install_csh
from libtbx.bundle import install_bat
import sys, os

def run(args):
  if (len(os.listdir(".")) != 0):
    print("Please use this command only in an empty directory.")
    return
  if (len(args) != 2) and (len(args) != 3):
    print("usage: libtbx.start_binary_bundle bundle_name top_modules [--single_directory]")
    return
  single_dir = False
  if (len(args) == 2):
    bundle_name, top_modules = args
  else :
    if (args[2] != "--single_directory"):
      print("usage: libtbx.start_binary_bundle bundle_name top_modules [--single_directory]")
      return
    bundle_name, top_modules = args[:2]
    single_dir = True
  if (single_dir):
    os.mkdir(bundle_name)
    os.chdir(bundle_name)
  copy_all.run(bundle_name)
  if (os.name == "nt"):
    install_script = bundle_name+"_install_script.bat"
    open(install_script, "w").write(
      install_bat.create_script(
        bundle=bundle_name,
        top_modules=top_modules,
        single_dir=single_dir))
  else:
    install_script = bundle_name+"_install_script.csh"
    open(install_script, "w").write(
      install_csh.create_script(
        bundle=bundle_name,
        top_modules=top_modules))
    os.chmod(install_script, 0o755)

if (__name__ == "__main__"):
  run(sys.argv[1:])
