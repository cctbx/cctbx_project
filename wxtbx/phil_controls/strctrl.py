
from __future__ import absolute_import, division, print_function
from wxtbx.phil_controls.text_base import ValidatedTextCtrl, TextCtrlValidator
from wxtbx import phil_controls
import wxtbx
from libtbx.phil import tokenizer
from libtbx.utils import to_unicode, to_str
from libtbx import Auto
import libtbx.phil
import wx
import sys

class StrCtrl(ValidatedTextCtrl):
  def __init__(self, *args, **kwds):
    kwds = dict(kwds)
    if (kwds.get("size", wx.DefaultSize) == wx.DefaultSize):
      kwds['size'] = (200,-1)
    super(StrCtrl, self).__init__(*args, **kwds)
    self._min_len = 0
    self._max_len = sys.maxunicode

  def CreateValidator(self):
    return StrValidator()

  def SetValue(self, value):
    if (value in [None, Auto]):
      ValidatedTextCtrl.SetValue(self, "")
    else :
      if isinstance(value, str):
        if wxtbx.is_unicode_build():
          ValidatedTextCtrl.SetValue(self, to_unicode(value))
        else :
          ValidatedTextCtrl.SetValue(self, value)
      else :
        if (sys.version_info.major < 3) and (not isinstance(value, unicode)):
          raise RuntimeError("Improper value (type: %s) for control %s" %
            (type(value).__name__, self.GetName()))
        if (sys.version_info.major >= 3) and (not isinstance(value, str)):
          raise RuntimeError("Improper value (type: %s) for control %s" %
            (type(value).__name__, self.GetName()))
        ValidatedTextCtrl.SetValue(self, value)

  def GetPhilValue(self):
    self.Validate()
    val_str = self.GetValue().strip()
    if (val_str in ["", "none", "None"]):
      return self.ReturnNoneIfOptional()
    return val_str

  def GetStringValue(self):
    value = self.GetPhilValue()
    if (value is None):
      return "None"
    else :
      return parse_str(value)

  def FormatValue(self, value):
    if wxtbx.is_unicode_build():
      return to_str(value)
    else :
      return value

  def SetMinLength(self, n):
    assert (n >= 0)
    self._min_len = n

  def SetMaxLength(self, n):
    assert (n >= 1)
    self._max_len = n

  def GetMinLength(self):
    return self._min_len

  def GetMaxLength(self):
    return self._max_len

class StrValidator(TextCtrlValidator):
  def CheckFormat(self, value):
    window = self.GetWindow()
    if ("$" in value):
      raise ValueError("The dollar symbol ($) may not be used here.")
    elif (len(value) > window.GetMaxLength()):
      raise ValueError("Value must be %d characters or less." %
        window.GetMaxLength())
    elif (len(value) < window.GetMinLength()):
      raise ValueError("Value must be at least %d characters." %
        window.GetMinLength())
    return value # XXX does anything else need to be done here?

def parse_str(value):
  #value = value.decode("utf-8")
  try :
    word = tokenizer.word(value, quote_token='"""')
    phil_string = str(word)
  except ValueError as e :
    raise
  else :
    return phil_string

if (__name__ == "__main__"):
  app = wx.App(0)
  frame = wx.Frame(None, -1, "String parameter test")
  panel = wx.Panel(frame, -1, size=(720,400))
  txt1 = wx.StaticText(panel, -1, "Job title:", pos=(10,100))
  ctrl1 = StrCtrl(panel, -1, value=None, pos=(160, 100), size=(400,-1),
    name="Job title")
  txt2 = wx.StaticText(panel, -1, "Output file prefix:", pos=(10,200))
  ctrl2 = StrCtrl(panel, -1, value="refine", pos=(160,200),
    name="Output file prefix")
  ctrl2.SetOptional(False)
  btn = wx.Button(panel, -1, "Process input", pos=(400, 360))
  btn2 = wx.Button(panel, -1, "Toggle title", pos=(200, 360))
  master_phil = libtbx.phil.parse("""
    title = None
      .type = str
    prefix = None
      .type = str""")
  def OnOkay(evt):
    print("""title = %s""" % ctrl1.GetStringValue())
    title_phil = libtbx.phil.parse("""title = %s""" % ctrl1.GetStringValue())
    prefix_phil = libtbx.phil.parse("""prefix = %s""" % ctrl2.GetStringValue())
    p = master_phil.fetch(sources=[title_phil, prefix_phil]).extract()
    print("title recycled via phil:", p.title)
    print("prefix recycled via phil:", p.prefix)
    value1 = ctrl1.GetPhilValue()
    value2 = ctrl2.GetPhilValue()
    assert (p.title == value1), value1
    assert (p.prefix == value2)
  assert (ctrl1.GetPhilValue() is None)
  assert (ctrl1.GetStringValue() == "None")
  inp_str = """This string has bad; characters f\""""
  ctrl1.SetValue(inp_str)
  assert (ctrl1.GetPhilValue() == inp_str)
  #assert (ctrl1.GetStringValue() == '"This string has bad; characters f\\""')
  title_phil = libtbx.phil.parse("""title = %s""" % ctrl1.GetStringValue())
  p = master_phil.fetch(source=title_phil).extract()
  assert (p.title == inp_str)
  assert (ctrl2.GetPhilValue() == "refine")
  assert (ctrl2.GetStringValue() == '"""refine"""')
  def OnChange(evt):
    pass
  def OnToggle(evt):
    if (ctrl1.IsEnabled()) : ctrl1.Enable(False)
    else : ctrl1.Enable()
  frame.Bind(wx.EVT_BUTTON, OnOkay, btn)
  frame.Bind(phil_controls.EVT_PHIL_CONTROL, OnChange)
  frame.Bind(wx.EVT_BUTTON, OnToggle, btn2)
  frame.Fit()
  frame.Show()
  app.MainLoop()
