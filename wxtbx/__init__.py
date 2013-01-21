from __future__ import division

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
