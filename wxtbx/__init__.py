from __future__ import division
import wx

global default_font_size

if (wx.Platform == '__WXMSW__') :
  default_font_size = 9
elif (wx.Platform == '__WXMAC__') :
  default_font_size = 12
else :
  default_font_size = 11

class MouseWheelTransparencyMixin (object) :
  """
  This mixin provides an event handler for passing the mouse wheel event to
  the parent, presumably a ScrolledPanel or similar.  For this to happen, the
  actual class must bind wx.EVT_MOUSEWHEEL to self.OnMouseWheel.
  """
  def OnMouseWheel (self, evt) :
    parent = self.GetParent()
    evt.SetId(parent.GetId())
    evt.SetEventObject(parent)
    parent.GetEventHandler().ProcessEvent(evt)

def is_unicode_build () :
  return (wx.PlatformInfo[2] == 'unicode')
