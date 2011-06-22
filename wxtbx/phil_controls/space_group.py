
from wxtbx.phil_controls import TextCtrlValidator
from libtbx.utils import Abort
import wx

class SpaceGroupControl (wx.ComboBox) :
  def __init__ (self, *args, **kwds) :
    super(SpaceGroupControl, self).__init__(*args, **kwds)
    self.SetValidator(SpaceGroupValidator())
    self.Bind(wx.EVT_TEXT_ENTER, lambda evt: self.Validate(), self)

  def SetSpaceGroup (self, sg) :
    from cctbx import sgtbx
    assert (isinstance(sg, sgtbx.space_group_info) or isinstance(sg, str) or
            (sg is None))
    items = self.GetItems()
    if (sg is not None) :
      sg = str(sg)
    else :
      sg = ""
    if (not sg in items) :
      self.Append(sg)
    self.SetStringSelection(sg)

  def SetValue (self, sg) :
    self.SetSpaceGroup(sg)

  def GetPhilValue (self) :
    self.Validate()
    sg = str(self.GetValue())
    if (sg == "") :
      return None
    from cctbx import sgtbx
    return sgtbx.space_group_info(symbol=sg)

  def Validate (self) :
    # XXX why doesn't self.Validate() work?
    if self.GetValidator().Validate(self.GetParent()) :
      return True
    else :
      raise Abort()

  if (wx.Platform == "__WXMAC__") and (wx.PlatformInfo[4] != 'wxOSX-cocoa') :
    def SetBackgroundColour (self, color) :
      self.GetChildren()[0].SetBackgroundColour(color)

class SpaceGroupValidator (TextCtrlValidator) :
  def CheckFormat (self, value) :
    if len(value) == 0 :
      return ""
    from cctbx import sgtbx
    sg = str(sgtbx.space_group_info(symbol=value))
    return sg

if (__name__ == "__main__") :
  app = wx.App(0)
  frame = wx.Frame(None, -1, "Space group test")
  panel = wx.Panel(frame, -1, size=(600,400))
  txt1 = wx.StaticText(panel, -1, "Space group:", pos=(100,180))
  sg_ctrl = SpaceGroupControl(panel, -1, pos=(200,180),
    name="Space group")
  btn = wx.Button(panel, -1, "Process input", pos=(400, 360))
  def OnOkay (evt) :
    sg = sg_ctrl.GetPhilValue()
    print type(sg).__name__, str(sg)
  frame.Bind(wx.EVT_BUTTON, OnOkay, btn)
  frame.Fit()
  frame.Show()
  app.MainLoop()
