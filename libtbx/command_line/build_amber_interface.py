from __future__ import division
import os, sys
import shutil

import libtbx.load_env
from libtbx.utils import Sorry
from libtbx import easy_run

dispatcher_include_str = """#
# Environment additions for AMBER interface
#

# include before command
if [ ! -z "$AMBERHOME" ]; then
  if [ "$LIBTBX_OS_NAME" = "Darwin" ]; then
    if [ -z "$DYLD_LIBRARY_PATH" ]; then
      export DYLD_LIBRARY_PATH=${AMBERHOME}/lib
    else
      export DYLD_LIBRARY_PATH=${AMBERHOME}/lib:${DYLD_LIBRARY_PATH}
    fi
  else
    if [ -z "$LD_LIBRARY_PATH" ]; then
      export LD_LIBRARY_PATH=${AMBERHOME}/lib
    else
      export LD_LIBRARY_PATH=${AMBERHOME}/lib:${LD_LIBRARY_PATH}
    fi
    if [ -z "$PYTHONPATH" ]; then
      export PYTHONPATH="${AMBERHOME}/lib/%s/site-packages"
    else
      export PYTHONPATH="${AMBERHOME}/lib/%s/site-packages:${PYTHONPATH}"
    fi
  fi
else
  echo 'This Phenix build has been configured to use Amber'
  echo 'Therefore environment variable AMBERHOME needs set'
  exit
fi
"""

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

  build_dir = os.path.dirname(os.path.dirname(phenix_dir))
  build_dir = os.path.join(build_dir, "build")

  if 0: # don't do this as the python version can change
    dispatcher_include = os.path.join(os.path.dirname(phenix_dir),
                                      "amber_adaptbx",
                                      "dispatcher_include_amber.sh",
      )
    shutil.copyfile(dispatcher_include,
                    os.path.join(build_dir,
                                 "dispatcher_include_amber.sh",
                                 )
      )
  else:
    os.chdir(build_dir)
    for filename in os.listdir(os.path.join(amberhome, "lib")):
      if filename.startswith("python"):
        break
    f=file("dispatcher_include_amber.sh", "wb")
    f.write( dispatcher_include_str % (filename, filename) )
    f.close()

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
