from libtbx.bundle import copy_all
from libtbx.bundle import install_csh
from libtbx.bundle import install_bat
import sys, os

def run(args):
  if (len(os.listdir(".")) != 0):
    print "Please use this command only in an empty directory."
    return
  if (len(args) != 2 and (len(args) != 3 or args[0] != "--application")):
    print "usage: libtbx.start_binary_bundle [--application]",
    print "bundle_name top_modules"
    return
  if (len(args) == 2):
    bundle_name, top_modules = args
    prefer_usr_bin_python = 00000
  else:
    bundle_name, top_modules = args[1:]
    prefer_usr_bin_python = 0001
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
        top_modules=top_modules,
        prefer_usr_bin_python=prefer_usr_bin_python))
    os.chmod(install_script, 0755)

if (__name__ == "__main__"):
  run(sys.argv[1:])
