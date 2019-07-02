from __future__ import absolute_import, division, print_function

from wxtbx.phil_controls import intctrl, floatctrl, symop, strctrl, ints, choice
from wxtbx.utils import std_sizer_flags, add_ok_cancel_buttons
from libtbx.utils import Abort
import wx
from six.moves import range

class SimpleInputDialog(wx.Dialog):
  def __init__(self,
                parent,
                title,
                label,
                value=None,
                caption=None):
    style = wx.CAPTION|wx.CLOSE_BOX|wx.RAISED_BORDER| \
      wx.WS_EX_VALIDATE_RECURSIVELY | wx.RESIZE_BORDER
    wx.Dialog.__init__(self,
      parent=parent,
      title=title,
      style=style)
    self.sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.sizer)
    self.inner_sizer = wx.BoxSizer(wx.VERTICAL)
    self.sizer.Add(self.inner_sizer, 1, wx.EXPAND|wx.ALL, 5)
    if (caption is not None):
      caption_txt = wx.StaticText(self, -1, caption)
      caption_txt.Wrap(480)
      self.inner_sizer.Add(caption_txt, 0, wx.ALL, 5)
    input_szr = wx.BoxSizer(wx.HORIZONTAL)
    self.inner_sizer.Add(input_szr, 0, wx.ALIGN_CENTER)
    label_txt = wx.StaticText(self, -1, label + ":")
    input_szr.Add(label_txt, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.phil_ctrl = self.CreatePhilControl(value)
    input_szr.Add(self.phil_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    add_ok_cancel_buttons(self, self.sizer)
    self.Fit()
    self.Centre(wx.BOTH)

  def CreatePhilControl(self, value):
    raise NotImplementedError()

  def GetPhilValue(self):
    return self.phil_ctrl.GetPhilValue()

  def __getattr__(self, name):
    return getattr(self.phil_ctrl, name)

class IntegerDialog(SimpleInputDialog):
  def CreatePhilControl(self, value):
    return intctrl.IntCtrl(
      parent=self,
      value=value)

class FloatDialog(SimpleInputDialog):
  def CreatePhilControl(self, value):
    return floatctrl.FloatCtrl(
      parent=self,
      value=value)

class StringDialog(SimpleInputDialog):
  def CreatePhilControl(self, value):
    return strctrl.StrCtrl(
      parent=self,
      value=value)

class SymopDialog(SimpleInputDialog):
  def CreatePhilControl(self, value):
    return symop.SymopCtrl(
      parent=self,
      value=value)

class SymopChoiceDialog(SimpleInputDialog):
  def CreatePhilControl(self, value):
    return symop.SymopChoiceCtrl(
      parent=self)

class IntegersDialog(SimpleInputDialog):
  def CreatePhilControl(self, value):
    return ints.IntsCtrl(
      parent=self,
      value=value)

class ChoiceDialog(SimpleInputDialog):
  def CreatePhilControl(self, value=None):
    return choice.ChoiceCtrl(
      parent=self)

  def SetChoices(self, *args, **kwds):
    self.phil_ctrl.SetChoices(*args, **kwds)
    self.Layout()

def get_phil_value_from_dialog(dlg):
  abort = False
  if (dlg.ShowModal() == wx.ID_OK):
    value = dlg.GetPhilValue()
  else :
    abort = True
  wx.CallAfter(dlg.Destroy)
  if (abort):
    raise Abort()
  return value

def get_float_value(**kwds):
  dlg = FloatDialog(**kwds)
  return get_phil_value_from_dialog(dlg)

def get_integer_value(**kwds):
  dlg = IntegerDialog(**kwds)
  return get_phil_value_from_dialog(dlg)

def get_miller_index(**kwds):
  dlg = IntegersDialog(**kwds)
  dlg.SetSizeMin(3)
  dlg.SetSizeMax(3)
  dlg.SetOptional(False)
  result = get_phil_value_from_dialog(dlg)
  if (isinstance(result, list)):
    return tuple(result)
  return result

RT_DIALOG_ENABLE_FRACTIONAL = 1
class RTDialog(wx.Dialog):
  def __init__(self, *args, **kwds):
    kwds = dict(kwds)
    style = kwds.get('style', 0)
    style |= wx.CAPTION|wx.CLOSE_BOX|wx.RAISED_BORDER| \
      wx.WS_EX_VALIDATE_RECURSIVELY
    kwds['style'] = style
    extra_style = kwds.pop("wxtbxStyle", 0)
    wx.Dialog.__init__(self, *args, **kwds)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    self._r_ctrls = []
    self._t_ctrls = []
    txt = wx.StaticText(self, -1,
      "Please enter rotation and translation operations to apply.  All "+
      "values must be defined.")
    szr.Add(txt, 0, wx.ALL, 5)
    grid = wx.FlexGridSizer(rows=4, cols=4)
    szr.Add(grid, 0, wx.ALL, 5)
    label1 = wx.StaticText(self, -1, "Rotation:")
    grid.Add(label1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    for i in range(3):
      if (i > 0):
        grid.Add((1,1))
      for j in range(3):
        _value = 0.0
        if (j == i):
          _value = 1.0
        ctrl = floatctrl.FloatCtrl(
          parent=self,
          name="Rotation",
          value=_value)
        ctrl.SetOptional(False)
        grid.Add(ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
        self._r_ctrls.append(ctrl)
    label2 = wx.StaticText(self, -1, "Translation:")
    grid.Add(label2, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    for i in range(3):
      ctrl = floatctrl.FloatCtrl(
        parent=self,
        name="Translation",
        value=0)
      ctrl.SetOptional(False)
      grid.Add(ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
      self._t_ctrls.append(ctrl)
    self.frac_box = None
    if (extra_style & RT_DIALOG_ENABLE_FRACTIONAL):
      self.frac_box = wx.CheckBox(self, -1, "Use fractional coordinates")
      szr.Add(self.frac_box, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    add_ok_cancel_buttons(self, szr)
    self.Fit()
    self.Centre(wx.BOTH)

  def GetMatrix(self):
    from scitbx import matrix
    r = [ c.GetPhilValue() for c in self._r_ctrls ]
    t = [ c.GetPhilValue() for c in self._t_ctrls ]
    return matrix.rt((r,t))

  def IsFractional(self):
    if (self.frac_box is not None):
      return self.frac_box.GetValue()
    return False

def get_rt_matrix(parent=None, enable_fractional=False):
  style = 0
  if (enable_fractional):
    style = RT_DIALOG_ENABLE_FRACTIONAL
  dlg = RTDialog(
    parent=parent,
    title="Rotation/translation operator",
    wxtbxStyle=style)
  rt = None
  if (dlg.ShowModal() == wx.ID_OK):
    rt = dlg.GetMatrix()
  wx.CallAfter(dlg.Destroy)
  if (rt is None):
    raise Abort()
  return rt

class HTTPProxyDialog(wx.Dialog):
  def __init__(self, *args, **kwds):
    wx.Dialog.__init__(self, *args, **kwds)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    szr2 = wx.BoxSizer(wx.VERTICAL)
    szr.Add(szr2, 0, wx.ALL|wx.ALIGN_CENTER, 10)
    txt = wx.StaticText(parent=self,
      label="Please enter the information required by your HTTP proxy.  The "+
        "authentication information may not be needed at all sites.",
      size=(400,-1))
    txt.Wrap(400)
    szr2.Add(txt)
    grid = wx.FlexGridSizer(cols=2)
    szr2.Add(grid)
    grid.Add(wx.StaticText(self, label="HTTP proxy server:"), 0,
      std_sizer_flags, 5)
    self.server_ctrl = strctrl.StrCtrl(parent=self, size=(300,-1),
      name="HTTP proxy server")
    self.server_ctrl.SetOptional(False)
    grid.Add(self.server_ctrl, 0, std_sizer_flags, 5)
    grid.Add(wx.StaticText(self, label="Proxy port:"), 0, std_sizer_flags, 5)
    self.port_ctrl = intctrl.IntCtrl(parent=self, size=(100,-1),
      name="Proxy port")
    self.port_ctrl.SetOptional(False)
    grid.Add(self.port_ctrl, 0, std_sizer_flags, 5)
    grid.Add(wx.StaticText(self, label="User name:"), 0, std_sizer_flags, 5)
    self.user_ctrl = strctrl.StrCtrl(parent=self, size=(100,-1),
      name="User name")
    self.user_ctrl.SetOptional(False)
    grid.Add(self.user_ctrl, 0, std_sizer_flags, 5)
    grid.Add(wx.StaticText(self, label="Password:"), 0, std_sizer_flags, 5)
    self.password_ctrl = strctrl.StrCtrl(parent=self, size=(200,-1),
      name="Password", style=wx.TE_PASSWORD)
    self.password_ctrl.SetOptional(False)
    grid.Add(self.password_ctrl, 0, std_sizer_flags, 5)
    add_ok_cancel_buttons(self, szr)
    szr.Layout()
    self.Fit()
    self.Centre(wx.BOTH)

  def SetProxyServer(self, server):
    self.server_ctrl.SetValue(server)

  def SetProxyUser(self, user):
    self.user_ctrl.SetValue(user)

  def SetProxyPort(self, port):
    assert isinstance(port, int)
    self.port_ctrl.SetValue(port)

  def GetProxyServer(self):
    return self.server_ctrl.GetPhilValue()

  def GetProxyPort(self):
    return self.port_ctrl.GetPhilValue()

  def GetProxyPassword(self):
    return self.password_ctrl.GetPhilValue()

  def GetProxyUser(self):
    return self.user_ctrl.GetPhilValue()

  def InstallProxy(self):
    import libtbx.utils
    libtbx.utils.install_urllib_http_proxy(
      server=self.GetProxyServer(),
      port=self.GetProxyPort(),
      user=self.GetProxyUser(),
      password=self.GetProxyPassword())

if (__name__ == "__main__"):
  app = wx.App(0)
  value1 = get_float_value(
    parent=None,
    title="Float input",
    label="Bond sigma",
    value=None,
    caption="Please enter a sigma value (in Angstroms) for the selected bond.")
  print(value1)
  value2 = get_integer_value(
    parent=None,
    title="Integer input",
    label="Number of cycles",
    value=5)
  print(value2)
  dlg = SymopDialog(
    parent=None,
    title="Symmetry operator input",
    label="Symmetry operator",
    value=None)
  print(get_phil_value_from_dialog(dlg))
  dlg = RTDialog(
    parent=None,
    title="Rotation/translation operator",
    wxtbxStyle=RT_DIALOG_ENABLE_FRACTIONAL)
  if (dlg.ShowModal() == wx.ID_OK):
    rt = dlg.GetMatrix()
  dlg = HTTPProxyDialog(None, title="HTTP proxy authentication")
  if (dlg.ShowModal() == wx.ID_OK):
    dlg.InstallProxy()
