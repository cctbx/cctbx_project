
from __future__ import absolute_import, division, print_function
from wxtbx import phil_controls
from libtbx import Auto
import wx

WXTBX_TRIBOOL_TRUE_FALSE = 0
WXTBX_TRIBOOL_YES_NO = 1
WXTBX_TRIBOOL_DISABLE_AUTO = 2
WXTBX_TRIBOOL_AUTO_IS_NONE = 4

class TriBoolCtrl(wx.Choice, phil_controls.PhilCtrl):
  """
  Three-way boolean control: returns True, False, or Auto (optionally None).
  Can be two-way if WXTBX_TRIBOOL_DISABLE_AUTO is specified.
  """
  def __init__(self, *args, **kwds):
    phil_controls.PhilCtrl.__init__(self)
    self.SetOptional(True)
    kwds = dict(kwds)
    self._tribool_style = kwds.get("style", WXTBX_TRIBOOL_TRUE_FALSE)
    kwds['style'] = 0
    wx.Choice.__init__(self, *args, **kwds)
    items = []
    if (self._tribool_style & WXTBX_TRIBOOL_YES_NO):
      items.extend(["Yes","No"])
    else :
      items.extend(["True","False"])
    if (not self._tribool_style & WXTBX_TRIBOOL_DISABLE_AUTO):
      items.append("Auto")
    self.SetItems(items)
    self.Bind(wx.EVT_CHOICE, lambda evt: self.DoSendEvent(), self)

  def SetValue(self, value):
    if (value == True):
      self.SetSelection(0)
    elif (value == False):
      self.SetSelection(1)
    else :
      assert (value in [None, Auto])
      if (self._tribool_style & WXTBX_TRIBOOL_DISABLE_AUTO):
        self.SetSelection(wx.NOT_FOUND)
      else :
        self.SetSelection(2)

  def GetValue(self):
    return self.GetPhilValue()

  def GetPhilValue(self):
    sel = self.GetSelection()
    if (sel < 0):
      return self.ReturnNoneIfOptional()
    vals = [True, False]
    if (self._tribool_style & WXTBX_TRIBOOL_AUTO_IS_NONE):
      vals.append(None)
    else :
      vals.append(Auto)
    return vals[self.GetSelection()]

  def GetStringValue(self):
    return str(self.GetPhilValue())

if (__name__ == "__main__"):
  app = wx.App(0)
  frame = wx.Frame(None, -1, "Tribool test")
  panel = wx.Panel(frame, -1, size=(720,480))
  txt1 = wx.StaticText(panel, -1, "Hydrogens:", pos=(20,180))
  choice1 = TriBoolCtrl(panel, -1, pos=(240,180))
  choice1.SetOptional(True)
  txt2 = wx.StaticText(panel, -1, "Add waters:", pos=(20, 240))
  choice2 = TriBoolCtrl(panel, -1, pos=(240,240),
    style=WXTBX_TRIBOOL_DISABLE_AUTO)
  assert (choice1.GetPhilValue() == True)
  assert (choice2.GetPhilValue() == True)
  choice1.SetValue(None)
  assert (choice1.GetStringValue() == "Auto")
  choice2.SetOptional(False)
  choice2.SetValue(None)
  try :
    print(choice2.GetPhilValue())
  except Exception as e :
    assert (str(e) == "Value required for 'choice'.")
  else :
    assert 0
  frame.Fit()
  frame.Show()
  app.MainLoop()
