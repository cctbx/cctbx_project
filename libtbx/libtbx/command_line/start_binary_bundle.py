from libtbx.bundle import copy_all
from libtbx.bundle import install_csh
from libtbx.bundle import install_bat
import sys, os

def run(args):
  if (len(os.listdir(".")) != 0):
    print "Please use this command only in an empty directory."
    return
  if (len(args) != 2):
    print "usage: libtbx.start_binary_bundle bundle_name top_modules"
    return
  bundle_name, top_modules = args
  copy_all.run(bundle_name)
  if (os.name == "nt"):
    install_script = bundle_name+"_install_script.bat"
    open(install_script, "w").write(
      install_bat.create_script(
        bundle=bundle_name,
        top_modules=top_modules))
  else:
    install_script = bundle_name+"_install_script.csh"
    open(install_script, "w").write(
      install_csh.create_script(
        bundle=bundle_name,
        top_modules=top_modules))
    os.chmod(install_script, 0755)

if (__name__ == "__main__"):
  run(sys.argv[1:])
