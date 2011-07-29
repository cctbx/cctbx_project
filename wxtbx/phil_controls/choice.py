
from wxtbx import phil_controls
import wx
import re

class ChoiceCtrl (wx.Choice, phil_controls.PhilCtrl) :
  def __init__ (self, *args, **kwds) :
    super(ChoiceCtrl, self).__init__(*args, **kwds)
    self._options = None
    self.Bind(wx.EVT_CHOICE, lambda evt: self.DoSendEvent(), self)

  def SetChoices (self, choices, captions=None) :
    selection = None
    is_selected = [ ("*" in choice) for choice in choices ]
    if (True in is_selected) :
      selection = is_selected.index(True)
    choices = [ re.sub("\*", "", choice) for choice in choices ]
    if (captions is None) :
      captions = list(choices) # XXX force copy
    if (selection is None) :
      captions.insert(0, "---")
      choices.insert(0, None)
      selection = 0
    assert (len(captions) == len(choices))
    self._options = choices
    self.SetItems(captions)
    self.SetSelection(selection)

  def SetValue (self, value) :
    selection = self._options.index(value)
    self.SetSelection(selection)

  def GetValue (self) :
    raise NotImplementedError("Please use GetPhilValue()")

  def GetPhilValue (self) :
    """Returns a single string."""
    return self._options[self.GetSelection()]

  def GetStringValue (self) :
    """Returns the long format (all choices, '*' denotes selected)."""
    selection = self.GetSelection()
    choices_out = []
    for i, choice in enumerate(self._options) :
      if (choice is None) :
        continue
      elif (i == selection) :
        choices_out.append("*" + choice)
      else :
        choices_out.append(choice)
    return " ".join(choices_out)

if (__name__ == "__main__") :
  app = wx.App(0)
  frame = wx.Frame(None, -1, "Choice test")
  panel = wx.Panel(frame, -1, size=(720,480))
  txt1 = wx.StaticText(panel, -1, "NCS restraint type:", pos=(20,180))
  choice1 = ChoiceCtrl(panel, -1, pos=(240,180))
  choice1.SetChoices(["torsion","dihedral"])
  txt2 = wx.StaticText(panel, -1, "Output map format:", pos=(20, 240))
  choice2 = ChoiceCtrl(panel, -1, pos=(240,240))
  choice2.SetChoices(["*ccp4", "xplor"], ["CCP4", "X-PLOR"])
  choice1.SetValue("torsion")
  assert (choice1.GetPhilValue() == "torsion")
  choice1.SetValue(None)
  assert (choice1.GetPhilValue() == None)
  assert (choice1.GetStringValue() == "torsion dihedral")
  assert (choice2.GetPhilValue() == "ccp4")
  assert (choice2.GetStringValue() == "*ccp4 xplor")
  frame.Fit()
  frame.Show()
  app.MainLoop()
