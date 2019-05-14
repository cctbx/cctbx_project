
from __future__ import absolute_import, division, print_function
from wxtbx import phil_controls
from wxtbx.phil_controls.text_base import ValidatedTextCtrl, TextCtrlValidator
from libtbx.utils import Sorry
from libtbx.phil import strings_from_words, strings_as_words
from libtbx.phil import tokenizer
import wx

class StringsCtrl(ValidatedTextCtrl):
  def __init__(self, *args, **kwds):
    ValidatedTextCtrl.__init__(self, *args, **kwds)

  def GetPhilValue(self):
    self.Validate()
    vals = parse_strings(self.GetValue())
    if (len(vals) == 0):
      return None
    return vals

  def FormatValue(self, value):
    return " ".join([ str(w) for w in strings_as_words(value) ])

  def CreateValidator(self):
    return StringsValidator()

  def SetValue(self, value):
    if (value is None):
      ValidatedTextCtrl.SetValue(self, "")
    else :
      assert isinstance(value, list), value
      ValidatedTextCtrl.SetValue(self, self.FormatValue(value))

class StringsValidator(TextCtrlValidator):
  def CheckFormat(self, value):
    return parse_strings(value)

def parse_strings(value):
  try :
    if (value != "") and (not "\"" in value) and (not "'" in value):
      value = "\"" + "\" \"".join(value.split()) + "\""
    words = list(tokenizer.word_iterator(value))
    string_list = strings_from_words(words)
  except ValueError as e :
    raise Sorry(str(e))
  else :
    return string_list

if (__name__ == "__main__"):
  app = wx.App(0)
  frame = wx.Frame(None, -1, "Strings control test")
  panel = wx.Panel(frame, -1, size=(600,400))
  txt1 = wx.StaticText(panel, -1, "Elements:", pos=(100,180))
  elems_ctrl = StringsCtrl(panel, -1, pos=(300,180), size=(120,-1),
      name="Elements", style=wx.TE_PROCESS_ENTER)
  txt2 = wx.StaticText(panel, -1, "Labels list:", pos=(100,240))
  labels_ctrl = StringsCtrl(panel, -1, pos=(300,240), size=(200,-1),
      name="Labels list")
  txt2 = wx.StaticText(panel, -1, "Map types:", pos=(100,300))
  maps_ctrl = StringsCtrl(panel, -1, pos=(300,300), size=(200,-1),
      name="Map types")
  def OnUpdate(evt):
    elems = elems_ctrl.GetPhilValue()
    print("Current elements:", elems)
  def OnOkay(evt):
    print("elems:", elems_ctrl.GetPhilValue())
    print("elems phil:", elems_ctrl.GetStringValue())
    print("labels:", labels_ctrl.GetPhilValue())
    print("labels phil:", labels_ctrl.GetStringValue())
    print("map types:", maps_ctrl.GetPhilValue())
    print("map types phil:", maps_ctrl.GetStringValue())
  frame.Bind(phil_controls.EVT_PHIL_CONTROL, OnUpdate, elems_ctrl)
  btn = wx.Button(panel, -1, "Process input", pos=(400, 360))
  frame.Bind(wx.EVT_BUTTON, OnOkay, btn)
  frame.Fit()
  frame.Show()
  app.MainLoop()
