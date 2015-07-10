import os, sys
import shutil

import libtbx.load_env
from libtbx.utils import Sorry
from libtbx import easy_run

def run():
  phenix_dir = libtbx.env.dist_path("phenix")
  print '\n  Phenix modules found :',phenix_dir
  amberhome = os.environ.get("AMBERHOME", None)
  if amberhome is None:
    raise Sorry("$AMBERHOME not set")
  print "\n  $AMBERHOME set :",amberhome
  print "\n  Linking AMBERHOME to Phenix modules"
  assert os.path.exists(amberhome)
  target = os.path.join(os.path.dirname(phenix_dir), "amber")
  if os.path.exists(target):
    print '\n  Not linking because it already exists : %s' % target
  else:
    os.symlink(amberhome, target)

  dispatcher_include = os.path.join(os.path.dirname(phenix_dir),
                                    "amber_adaptbx",
                                    "dispatcher_include_amber.sh",
                                    )
  build_dir = os.path.dirname(os.path.dirname(phenix_dir))
  build_dir = os.path.join(build_dir, "build")
  shutil.copyfile(dispatcher_include,
                  os.path.join(build_dir,
                               "dispatcher_include_amber.sh",
                               )
    )
  print "Building amber_adaptbx"
  cmd = "libtbx.configure amber amber_adaptbx"
  print "\n  ~> %s\n" % cmd
  easy_run.call(cmd)
  os.chdir(build_dir)
  cmd = "libtbx.scons -j 1"
  print "\n  ~> %s\n" % cmd
  easy_run.call(cmd)
  

if __name__=="__main__":
  run()#sys.argv[1])
  
