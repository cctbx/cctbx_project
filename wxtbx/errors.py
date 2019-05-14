from __future__ import absolute_import, division, print_function

from libtbx.utils import Abort, Sorry
from traceback import format_exception
import wx

def show_sorry(message):
  wx.MessageBox(message=message,
    caption="Error",
    style=wx.OK|wx.ICON_EXCLAMATION)

def wx_excepthook(type, value, traceback):
  if (type is Abort):
    pass
  elif (type is Sorry):
    show_sorry(str(value))
  else :
    message = process_exception(type, value, traceback)
    if (message is None):
      print(str(type.__name__) + ": " + str(value))
      print("".join(format_exception(type, value, traceback)))
    else :
      assert isinstance(message, str)
      show_sorry(message)

def process_exception(type, value, traceback):
  if (type is IOError):
    if (value.errno == 5) : # FIXME how to handle this?
      return ("Input/output error!  This could be a problem with your "+
        "filesystem and/or network (if an NFS mount is in use).")
    elif (value.errno == 13):
      return (("Can't write file %s - permission denied.  Please " +
        "change the permissions on the enclosing folder, or use a " +
        "different directory.") % value.filename)
    elif (value.errno == 28):
      return ("System error: no space left on device.  Please delete "+
        "or move unused files; you can remove obsolete or failed PHENIX "+
        "jobs from the job history interface, which will delete the "+
        "associated folders (if any).")
    elif (value.errno == 122):
      return ("System error: disk quota exceeded.  Please delete some "+
        "files to avoid this error (or move them to another disk).")
    elif (value.errno == 2001):
      return (("The file '%s' appears to have been moved or deleted - "+
        "can't process data.") % value.filename)
  elif (type is OSError):
    if (value.errno == 13):
      return (("Permission denied: %s.  Please change the permissions "+
        "on the enclosing folder, or use a different directory.") %
        value.filename)
    elif (value.errno == 16):
      return (("Filesystem error: %s.  This usually indicates a problem "+
        "communicating with an NFS server.") % str(value.message))
    elif (value.errno == 17):
      return (("A directory named %s already exists in the destination "+
        "directory; please move or delete it, or select a different "+
        "destination directory.") % value.filename)
  return None
