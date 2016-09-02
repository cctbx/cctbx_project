# LIBTBX_SET_DISPATCHER_NAME phenix.build_divcon_interface
from __future__ import division
import os
import libtbx.load_env
from libtbx.utils import Sorry
from libtbx import easy_run

dispatcher_include_str = """#
# Environment additions for DivCon interface
# include before command
if [ ! -z "$QB_PYTHONPATH" ]; then
  export PYTHONPATH=$PYTHONPATH:$QB_PYTHONPATH
fi
if [ ! -z "$QB_LD_LIBRARY_PATH" ]; then
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$QB_LD_LIBRARY_PATH
fi
if [ ! -z "$QB_DYLD_LIBRARY_PATH" ]; then
  export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$QB_DYLD_LIBRARY_PATH
fi
"""

def run():
  phenix_dir = libtbx.env.dist_path("phenix")
  print '\n  Phenix modules found :',phenix_dir
  qbhome = os.environ.get("QBHOME", None)
  if qbhome is None:
    raise Sorry("$QBHOME not set")
  print "\n  $QBHOME set :",qbhome

  for env in ["QB_PYTHONPATH",
              "QB_LD_LIBRARY_PATH",
              "QB_DYLD_LIBRARY_PATH",
              ]:
    if not os.environ.get(env, False):
      print """
    Environment variable %s not set.
    Need to source file in $QBHOME/etc

    csh
      source %s

    bash
      source %s
    """ % (env,
           os.path.join(qbhome,
                        "etc",
                        "qbenv.csh"),
           os.path.join(qbhome,
                        "etc",
                        "qbenv.sh"),
    )
      raise Sorry('QB env. var. error')

  build_dir = os.path.dirname(os.path.dirname(phenix_dir))
  build_dir = os.path.join(build_dir, "build")

  os.chdir(build_dir)
  f=file("dispatcher_include_divcon.sh", "wb")
  f.write(dispatcher_include_str)
  f.close()

  print "Configure Phenix with DivCon"
  cmd = "libtbx.configure phenix"
  print "\n  ~> %s\n" % cmd
  easy_run.call(cmd)
  os.chdir(build_dir)

if __name__=="__main__":
  run()#sys.argv[1])
