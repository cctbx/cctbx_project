
from wxtbx import phil_controls
import wx

class ArrayLabelCtrl (wx.Choice, phil_controls.PhilCtrl) :
  def __init__ (self, *args, **kwds) :
    super(ArrayLabelCtrl, self).__init__(*args, **kwds)
    if (wx.Platform == '__WXMAC__') and (wx.VERSION[1] > 8) :
      self._default_value = "---"
    else :
      self._default_value = ""
    self.Bind(wx.EVT_CHOICE, self.OnChoose)

  def SetLabel (self, label) :
    assert isinstance(label, str) or isinstance(label, None)
    if (label is None) :
      self.SetSelection(0)
    else :
      self.SetStringSelection(label)

  def SetLabelChoices (self, labels) :
    old_label = self.GetStringSelection()
    self.SetItems([self._default_value] + labels)
    if (old_label in labels) :
      self.SetLabel(old_label)
    else :
      self.SetLabel(labels[0])
    self.Layout()

  def GetPhilValue (self) :
    value = self.GetStringSelection()
    if (value == self._default_value) or (value == "") :
      return None
    else :
      return value

  def GetStringVAlue (self) :
    return str(self.GetPhilValue())

  def OnChoose (self, event) :
    label = self.GetPhilValue()
    print label

class ArrayLabelsCtrl (ArrayLabelCtrl) :
  def SetLabel (self, label) :
    if (isinstance(label, list)) :
      assert (len(label) == 1)
      ArrayLabelCtrl.SetLabel(self, label[0])
    else :
      ArrayLabelCtrl.SetLabel(self, label)

  def GetPhilValue (self) :
    value = self.GetStringSelection()
    if (value == self._default_value) or (value == "") :
      return None
    else :
      return [value]

  def GetStringVAlue (self) :
    return str(ArrayLabelCtrl.GetPhilValue(self))

# testing
if (__name__ == "__main__") :
  app = wx.App(0)
  frame = wx.Frame(None, -1, "Array label test")
  panel = wx.Panel(frame, -1, size=(640,480))
  txt1 = wx.StaticText(panel, -1, "Data labels:", pos=(20,180))
  choice_ctrl = ArrayLabelCtrl(panel, -1, pos=(240,180), size=(160,-1),
    name="Data labels")
  assert (choice_ctrl.GetPhilValue() is None)
  choices1 = ["F,SIGF", "F(+),SIGF(+),F(-),SIGF(-)", "IMEAN,SIGIMEAN"]
  choice_ctrl.SetLabelChoices(choices1)
  assert (choice_ctrl.GetPhilValue() == "F,SIGF")
  choices2 = ["F(+),SIGF(+),F(-),SIGF(-)", "I(+),SIGI(+),I(-),SIGI(-)",
    "IMEAN,SIGIMEAN", "F,SIGF"]
  choice_ctrl.SetLabelChoices(choices2)
  assert (choice_ctrl.GetPhilValue() == "F,SIGF")
  frame.Fit()
  frame.Show()
  app.MainLoop()
