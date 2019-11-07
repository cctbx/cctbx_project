
from __future__ import absolute_import, division, print_function
from wxtbx.phil_controls.text_base import ValidatedTextCtrl, TextCtrlValidator
from wxtbx import phil_controls
from libtbx.utils import Sorry
from libtbx import Auto
import wx

class SymopCtrl(ValidatedTextCtrl):
  def CreateValidator(self):
    return SymopValidator()

  def SetSymop(self, value):
    if type(value)==type(u'abc'): value = value.encode("ascii", "ignore")
    if (value is None) or (value is Auto):
      ValidatedTextCtrl.SetValue(self, "")
    elif (isinstance(value, str)):
      try :
        from cctbx import sgtbx
        rt_mx = sgtbx.rt_mx(symbol=value)
      except ValueError :
        raise Sorry("Inappropriate value '%s' for %s." % (value,
          self.GetName()))
      else :
        ValidatedTextCtrl.SetValue(self, str(value))
    else :
      raise TypeError("Type '%s' not allowed!" % type(value).__name__)

  def SetValue(self, value):
    self.SetSymop(value)

  def GetPhilValue(self):
    self.Validate()
    val_str = ValidatedTextCtrl.GetValue(self)
    if (val_str == ""):
      return self.ReturnNoneIfOptional()
    return val_str

  def FormatValue(self, value):
    return str(value)

class SymopValidator(TextCtrlValidator):
  def CheckFormat(self, value):
    if type(value)==type(u'abc'): value = value.encode("ascii", "ignore")
    from cctbx import sgtbx
    rt_mx = sgtbx.rt_mx(symbol=value)
    return value

class SymopChoiceCtrl(wx.Choice, phil_controls.PhilCtrl):
  def __init__(self, *args, **kwds):
    super(SymopChoiceCtrl, self).__init__(*args, **kwds)
    self._space_group = None
    self.Bind(wx.EVT_CHOICE, lambda evt: self.DoSendEvent(), self)

  def SetSpaceGroup(self, space_group):
    if (type(space_group).__name__ == "space_group_info"):
      space_group = space_group.group()
    self._space_group = space_group
    if (space_group is None):
      self.SetItems([])
    else :
      current_symop = self.GetStringSelection()
      symops = [ str(smx) for smx in space_group.smx() ]
      self.SetItems(symops)
      if (current_symop in symops):
        self.SetStringSelection(current_symop)

  def GetValue(self):
    """
    Returns a rotation/translation matrix.
    """
    if (self._space_group is not None):
      selected_symop = self.GetStringSelection()
      for smx in self._space_group.smx():
        if (str(smx) == selected_symop):
          return smx
    return None

  def GetPhilValue(self):
    """
    Returns a string if defined, otherwise None.
    """
    sel = self.GetStringSelection()
    if (sel == ""):
      return None
    return sel

  def SetValue(self, value):
    items = self.GetItems()
    if (str(value) in items):
      self.SetStringSelection(str(value))

if (__name__ == "__main__"):
  app = wx.App(0)
  frame = wx.Frame(None, -1, "Symop control test")
  panel = wx.Panel(frame, -1, size=(600,400))
  txt1 = wx.StaticText(panel, -1, "Symmetry operator:", pos=(100,180))
  sym_ctrl = SymopCtrl(panel, -1, pos=(300,180), size=(80,-1),
      name="Symmetry operator")
  txt2 = wx.StaticText(panel, -1, "Symmetry operator:", pos=(100,240))
  sym_ctrl_2 = SymopChoiceCtrl(panel, -1, pos=(300, 240), size=(160,-1),
      name="Symmetry operator choice")
  from cctbx import sgtbx
  sym_ctrl_2.SetSpaceGroup(sgtbx.space_group_info("P6322"))
  btn = wx.Button(panel, -1, "Process input", pos=(400, 360))
  def OnOkay(evt):
    symop = sym_ctrl.GetPhilValue()
    symop2 = sym_ctrl_2.GetPhilValue()
    print(symop, symop2)
  def OnPhilCtrl(evt):
    symop3 = sym_ctrl_2.GetValue()
    symop3.show_geometrical_elements()
  frame.Bind(wx.EVT_BUTTON, OnOkay, btn)
  frame.Bind(phil_controls.EVT_PHIL_CONTROL, OnPhilCtrl, sym_ctrl_2)
  frame.Fit()
  frame.Show()
  app.MainLoop()
