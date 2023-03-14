from __future__ import absolute_import, division, print_function

from wxtbx import errors
import wx
import sys

class CCTBXApp(wx.App):
  def __init__(self, *args, **kwds):

    # hide GTK messages
    if kwds.pop('gtk_suppress_diagnostics', True):
      try:
        wx.App.GTKSuppressDiagnostics()
      except AttributeError:
        pass

    wx.App.__init__(self, *args, **kwds)

  def OnInit(self):
    sys.excepthook = errors.wx_excepthook
    return True

  def SetTaskbarIcon(self, icon_path, label):
    icon = wx.Icon(icon_path, wx.BITMAP_TYPE_PNG)
    if ((wx.Platform == '__WXMAC__') and (wx.VERSION >= (2,9)) and
        (hasattr(wx, "TBI_DOCK"))):
      self.icon = wx.TaskBarIcon(wx.TBI_DOCK)
    else :
      self.icon = wx.TaskBarIcon()
    self.icon.SetIcon(icon, label)
