
from wxtbx import errors
import wx
import sys

class CCTBXApp (wx.App) :
  def __init__ (self, *args, **kwds) :
    wx.App.__init__(self, *args, **kwds)

  def OnInit (self) :
    sys.excepthook = errors.wx_excepthook
    return True
