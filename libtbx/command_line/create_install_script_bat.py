import sys

def run(args):
  from libtbx.utils import Usage
  import libtbx.load_env
  if (len(args) != 2):
    raise Usage("%s bundle_name top_modules" % libtbx.env.dispatcher_name)
  bundle_name, top_modules = args
  install_script = bundle_name+"_install_script.bat"
  from libtbx.bundle import install_bat
  open(install_script, "w").write(
    install_bat.create_script(
      bundle=bundle_name,
      top_modules=top_modules))

if (__name__ == "__main__"):
  run(sys.argv[1:])
