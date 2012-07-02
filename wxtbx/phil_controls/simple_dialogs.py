
from wxtbx.phil_controls import intctrl, floatctrl, symop, strctrl
from libtbx.utils import Abort
import wx

class SimpleInputDialog (wx.Dialog) :
  def __init__ (self,
                parent,
                title,
                label,
                value=None,
                caption=None) :
    style = wx.CAPTION|wx.CLOSE_BOX|wx.RAISED_BORDER| \
      wx.WS_EX_VALIDATE_RECURSIVELY
    wx.Dialog.__init__(self,
      parent=parent,
      title=title,
      style=style)
    self.sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.sizer)
    self.inner_sizer = wx.BoxSizer(wx.VERTICAL)
    self.sizer.Add(self.inner_sizer, 1, wx.EXPAND|wx.ALL, 5)
    if (caption is not None) :
      caption_txt = wx.StaticText(self, -1, caption)
      caption_txt.Wrap(480)
      self.inner_sizer.Add(caption_txt, 0, wx.ALL, 5)
    input_szr = wx.BoxSizer(wx.HORIZONTAL)
    self.inner_sizer.Add(input_szr, 0, wx.ALIGN_CENTER)
    label_txt = wx.StaticText(self, -1, label + ":")
    input_szr.Add(label_txt, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.phil_ctrl = self.CreatePhilControl(value)
    input_szr.Add(self.phil_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    ok_btn = wx.Button(self, wx.ID_OK)
    cancel_btn = wx.Button(self, wx.ID_CANCEL)
    btn_szr = wx.StdDialogButtonSizer()
    btn_szr.Add(cancel_btn, 0, wx.ALL, 5)
    btn_szr.Add(ok_btn, 0, wx.ALL, 5)
    ok_btn.SetDefault()
    btn_szr.Realize()
    self.sizer.Add(btn_szr, 0, wx.ALL|wx.ALIGN_RIGHT, 5)
    self.Fit()
    self.Centre(wx.BOTH)

  def CreatePhilControl (self, value) :
    raise NotImplementedError()

  def GetPhilValue (self) :
    return self.phil_ctrl.GetPhilValue()

  def __getattr__ (self, name) :
    return getattr(self.phil_ctrl, name)

class IntegerDialog (SimpleInputDialog) :
  def CreatePhilControl (self, value) :
    return intctrl.IntCtrl(
      parent=self,
      value=value)

class FloatDialog (SimpleInputDialog) :
  def CreatePhilControl (self, value) :
    return floatctrl.FloatCtrl(
      parent=self,
      value=value)

class StringDialog (SimpleInputDialog) :
  def CreatePhilControl (self, value) :
    return strctrl.StrCtrl(
      parent=self,
      value=value)

class SymopDialog (SimpleInputDialog) :
  def CreatePhilControl (self, value) :
    return symop.SymopCtrl(
      parent=self,
      value=value)

def get_phil_value_from_dialog (dlg) :
  abort = False
  if (dlg.ShowModal() == wx.ID_OK) :
    value = dlg.GetPhilValue()
  else :
    abort = True
  wx.CallAfter(dlg.Destroy)
  if (abort) :
    raise Abort()
  return value

def get_float_value (**kwds) :
  dlg = FloatDialog(**kwds)
  return get_phil_value_from_dialog(dlg)

def get_integer_value (**kwds) :
  dlg = IntegerDialog(**kwds)
  return get_phil_value_from_dialog(dlg)

if (__name__ == "__main__") :
  app = wx.App(0)
  value1 = get_float_value(
    parent=None,
    title="Float input",
    label="Bond sigma",
    value=None,
    caption="Please enter a sigma value (in Angstroms) for the selected bond.")
  print value1
  value2 = get_integer_value(
    parent=None,
    title="Integer input",
    label="Number of cycles",
    value=5)
  print value2
  dlg = SymopDialog(
    parent=None,
    title="Symmetry operator input",
    label="Symmetry operator",
    value=None)
  print get_phil_value_from_dialog(dlg)
