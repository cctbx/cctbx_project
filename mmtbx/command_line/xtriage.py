# LIBTBX_SET_DISPATCHER_NAME phenix.xtriage

from mmtbx.scaling import xtriage
import sys

if (__name__ == "__main__"):
  if "--gui" in sys.argv :
    try :
      from phenix.command_line import xtriage_gui
      xtriage_gui.run(args=sys.argv[1:])
    except ImportError :
      print "PHENIX GUI is not installed."
      sys.exit()
    except Exception, e :
      print e
      sys.exit()
  else :
    xtriage.run(args=sys.argv[1:])
