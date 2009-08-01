import sys, os
from optparse import OptionParser

def run():
  #CM 31jul09
  myparser = OptionParser()
  myparser.add_option("--current_working_directory", dest="cwd", 
                      default=os.getcwd(), help="current working dir")
  myparser.add_option("--build", dest="bld", default="release",
                      help="switch for build mode")
  # store values in myopts.{cwd,bld}:
  (myopts, myargs) = myparser.parse_args()
  # validargs taken from 'cctbx_install_script.csh':
  validargs = ['fftw3tbx', 'rstbx', 'mmtbx', 'smtbx', 'clipper']
  finalargs = list(set(filter((lambda x: x in validargs), myargs)))
  if (not finalargs):
    raise RuntimeError('\n%s\nNo valid command line argument. Valid args are:\
                        \n%s\n%s' % ('*'*78,validargs,'*'*78))
  else:
    newsysargv = [
        '%s'%sys.argv[0],
        '--current_working_directory=%s'%myopts.cwd,
        '--build=%s'%myopts.bld ]
    sys.argv = newsysargv + finalargs
  if (not hasattr(sys, "version_info")
      or sys.version_info[0] < 2
      or (sys.version_info[0] == 2 and sys.version_info[1] < 3)):
    print
    print "*" * 78
    print "FATAL: Python 2.3 or higher is required."
    print "Version currently in use:", sys.version
    print "*" * 78
    print
    return
  if (os.name == "nt"):
    open("shortpath.bat", "w").write("@echo off\necho %~s1\n")
  sys.path.insert(1, os.path.join(sys.path[0], "pythonpath"))
  sys.path[0] = os.path.dirname(sys.path[0])
  import libtbx.env_config
  libtbx.env_config.cold_start(sys.argv)

if (__name__ == "__main__"):
  run()
