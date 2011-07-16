
from wxtbx import phil_controls
import wx

class MultiChoiceCtrl (wx.Panel, phil_controls.PhilCtrl) :
  """
  Control for phil type 'choice(multi=True)'

  Composed of a grid of checkboxes, any or none of which may be checked.
  Each choice may be associated with an action upon clicking, such as
  opening a dialog, or with a contextual menu.

  In phil syntax, choices may be specified one of two ways:

    strategy=individual_adp+individual_sites+occupancies

    strategy = *individual_adp *individual_sites group_adp tls \
               *occupancies group_anomalous

  The first is typically used on the command line, the latter is used in
  parameter files.  When extracted as a Python object, the parameter will be
  a list of strings corresponding to the selected choices.
  """
  def __init__ (self, *args, **kwds) :
    super(MultiChoiceCtrl, self).__init__(*args, **kwds)
    self._n_columns = 3
    self._sizer = wx.FlexGridSizer(cols=self._n_columns, hgap=20, vgap=5)
    self.SetSizer(self._sizer)
    self._option_names = []
    self._option_controls = {}

  def SetCols (self, n_columns) :
    self._n_columns = n_columns
    self._sizer.SetCols(n_columns)
    if (len(self._option_controls) > 0) :
      self.Realize()

  def SetToolTip (self, tooltip) :
    super(MultiChoiceCtrl, self).SetToolTip(tooltip)
    for key, ctrl in self._option_controls.iteritems() :
      ctrl.SetToolTip(tooltip)

  def Realize (self) :
    """Finalizes layout and fits the panel."""
    self._sizer.Layout()
    self._sizer.Fit(self)

  def GetChoiceCtrl (self, choice) :
    return self._option_controls.get(choice)

  def AddChoice (self, choice, label=None) :
    """Adds a checkbox to the grid.  If choice starts with '*', the box will
    be checked automatically."""
    auto_enable = False
    if (choice.startswith("*")) :
      auto_enable = True
      choice = choice[1:]
    assert (not choice in self._option_names)
    if (label is None) :
      label = choice
    box = wx.CheckBox(
      parent=self,
      label=label,
      name=choice)
    self.Bind(wx.EVT_CHECKBOX, self.OnCheck, box)
    #self.Bind(wx.EVT_CONTEXT_MENU, self.OnRightClick, box)
    self._sizer.Add(box) #, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self._option_names.append(choice)
    self._option_controls[choice] = box
    if (auto_enable) :
      box.SetValue(True)

  def Enable (self, enable=True) :
    for choice in self._option_names :
      self._option_controls[choice].Enable(enable)

  def EnableChoice (self, choice, enable=True) :
    """Enables or disables the checkbox corresponding to the choice."""
    self._option_controls[choice].Enable(enable)

  def SetChoice (self, choice, selected=True) :
    """Select or deselect the specified choice."""
    self._option_controls[choice].SetValue(selected)

  def GetChoiceIndex (self, choice) :
    if (choice in self._option_names) :
      return self._option_names.index(choice)
    return None

  def GetValue (self) :
    raise NotImplementedError()

  def SetValue (self, value) :
    """
    Check the boxes specified by value.  A list is preferred, but either of
    the two string conventions may be used (selected choices separated by '+',
    or all choices separated by spaces with '*' denoting selected).
    """
    if (isinstance(value, str)) :
      if ("+" in value) :
        value = value.split("+")
      else :
        value_ = value.split()
        value = []
        for choice in value_ :
          if (choice.startswith("*")) :
            value.append(choice[1:])
    assert (isinstance(value, list)) or (isinstance(value, tuple))
    for choice in self._option_names :
      box = self._option_controls[choice]
      if (choice in value) :
        box.SetValue(True)
      else :
        box.SetValue(False)

  def GetPhilValue (self) :
    """Returns a list of strings."""
    values = []
    for choice in self._option_names :
      box = self._option_controls[choice]
      if (box.GetValue() == True) :
        values.append(choice)
    return values

  def GetStringValue (self) :
    """Returns the long format (all choices, '*' denotes selected)."""
    values = []
    for choice in self._option_names :
      box = self._option_controls[choice]
      if (box.GetValue() == True) :
        values.append("*" + choice)
      else :
        values.append(choice)
    return " ".join(values)

  def OnCheck (self, event) :
    box = event.GetEventObject()
    self.DoSendEvent(original_window=box)

# testing
if (__name__ == "__main__") :
  app = wx.App(0)
  frame = wx.Frame(None, -1, "Multi-choice test")
  panel = wx.Panel(frame, -1, size=(720,480))
  txt1 = wx.StaticText(panel, -1, "Refinement strategy:", pos=(20,180))
  choice_ctrl = MultiChoiceCtrl(panel, -1, pos=(240,160),
    name="Refinement strategy")
  choices1 = ["*individual_adp",
                 "*individual_sites",
                 "individual_sites_real_space",
                 "tls",
                 "group_adp",
                 "*occupancies",
                 "group_anomalous"]
  labels = ["Individual B-factors", "XYZ coordinates", "Real-space", "TLS",
    "Group B-factors", "Occupancies", "Anomalous scatterers"]
  for choice, label in zip(choices1, labels) :
    choice_ctrl.AddChoice(choice, label)
  choice_ctrl.Realize()
  txt2 = wx.StaticText(panel, -1, "Output formats:", pos=(20,300))
  choice_ctrl2 = MultiChoiceCtrl(panel, -1, pos=(240,300),
    name="Output formats")
  for choice in ["mtz", "pdb", "ccp4_map", "cif", "pkl"] :
    choice_ctrl2.AddChoice(choice)
  choice_ctrl2.Realize()
  choice_ctrl2.SetCols(5)
  choice_ctrl2.SetValue("mtz+pdb")
  frame.Fit()
  frame.Show()
  assert (choice_ctrl.GetPhilValue() == ["individual_adp","individual_sites",
    "occupancies"])
  assert (choice_ctrl2.GetPhilValue() == ["mtz", "pdb"])
  assert (choice_ctrl2.GetStringValue() == "*mtz *pdb ccp4_map cif pkl")
  choice_ctrl2.EnableChoice("ccp4_map", False)
  choice_ctrl2.SetChoice("cif")
  assert (choice_ctrl2.GetPhilValue() == ["mtz", "pdb", "cif"])
  choice_ctrl2.SetValue(["mtz", "pdb"])
  assert (choice_ctrl2.GetPhilValue() == ["mtz", "pdb"])
  app.MainLoop()
