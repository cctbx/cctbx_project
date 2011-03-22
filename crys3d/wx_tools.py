
import wx

draw_modes = [
  ("trace", "Backbone trace"),
  ("all_atoms", "All atoms"),
  ("bonded_only", "Bonded atoms")
]
draw_flags = [ ("flag_show_hydrogens", "Show hydrogens"),
               ("flag_show_ellipsoids", "Show B-factor ellipsoids"),
               ("flag_show_labels", "Show labels"),
               ("flag_show_noncovalent_bonds", "Show non-covalent bonds"), ]
color_modes = [ ("rainbow", "Sequence (rainbow)"),
                ("b", "Isotropic B-factor"),
                ("element", "Atomic element"),
                ("chain", "Chain ID"), ]

def set_tiny_font (control) :
  font = control.GetFont()
  font.SetPointSize(10)
  control.SetFont(font)

def set_tiny_bold_font (control) :
  font = control.GetFont()
  font.SetPointSize(10)
  font.SetWeight(wx.FONTWEIGHT_BOLD)
  control.SetFont(font)

class ModelControlPanel (wx.MiniFrame) :
  def __init__ (self, *args, **kwds) :
    super(ModelControlPanel, self).__init__(*args, **kwds)
    self.Bind(wx.EVT_CLOSE, self.OnClose)
    self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy)
    self.parent = self.GetParent()
    self.model_id = None
    self.inner_panel = wx.Panel(self, -1)
    self.panel_sizer = wx.BoxSizer(wx.VERTICAL)
    self.inner_panel.SetSizer(self.panel_sizer)
    self.panel1 = wx.Panel(self.inner_panel, -1, style=wx.RAISED_BORDER)
    self.panel_sizer.Add(self.panel1, 0)
    self.sizer = wx.BoxSizer(wx.VERTICAL)
    self.sizer.Add(self.inner_panel, 0, wx.EXPAND)
    szr1 = wx.BoxSizer(wx.HORIZONTAL)
    self.panel1.SetSizer(szr1)
    txt1 = wx.StaticText(self.panel1, -1, "Model:")
    set_tiny_bold_font(txt1)
    chooser = wx.Choice(self.panel1, -1,
      choices=self.parent.model_ids,
      size=(300,-1))
    set_tiny_font(chooser)
    self.Bind(wx.EVT_CHOICE, self.OnChooseModel, chooser)
    self.model_chooser = chooser
    szr1.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    szr1.Add(chooser, 0, wx.RIGHT|wx.TOP|wx.BOTTOM|wx.ALIGN_CENTER_VERTICAL, 5)
    self.panel2 = wx.Panel(self.inner_panel, -1)
    if (len(self.parent.model_ids) > 0) :
      chooser.SetSelection(0)
      self.set_model(self.parent.model_ids[0])
    else :
      self.sizer.Fit(self.inner_panel)
      self.Fit()

  def refresh_model_list (self) :
    models = self.parent.model_ids
    current_model = self.model_chooser.GetStringSelection()
    self.model_chooser.SetItems(models)
    if (current_model in models) :
      self.model_chooser.SetStringSelection(current_model)

  def set_model (self, model_id) :
    self.panel_sizer.Detach(self.panel2)
    self.panel2.Destroy()
    self.model_id = model_id
    self.panel2 = wx.Panel(self.inner_panel, -1)
    self.panel_sizer.Add(self.panel2, 1)
    p = self.panel2
    szr2 = wx.BoxSizer(wx.VERTICAL)
    self.panel2.SetSizer(szr2)
    visibility_ctrl = wx.CheckBox(p, -1,
      label="Model visible")
    if self.parent.show_object.get(model_id, False) :
      visibility_ctrl.SetValue(True)
    self.Bind(wx.EVT_CHECKBOX, self.OnSetVisibility, visibility_ctrl)
    set_tiny_bold_font(visibility_ctrl)
    szr2.Add(visibility_ctrl, 0, wx.ALL, 5)
    szr3 = wx.BoxSizer(wx.HORIZONTAL)
    txt = wx.StaticText(p, -1, "Show atoms:", size=(120,-1))
    set_tiny_bold_font(txt)
    szr3.Add(txt, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    viewer = self.parent
    model = viewer.get_model(model_id)
    modes = [ label for name, label in draw_modes ]
    mode_choice = wx.Choice(p, -1,
      choices=modes)
    set_tiny_font(mode_choice)
    for (name, label) in draw_modes :
      if (model.draw_mode == name) :
        mode_choice.SetStringSelection(label)
        break
    self.Bind(wx.EVT_CHOICE, self.OnChooseDrawMode, mode_choice)
    szr3.Add(mode_choice, 0, wx.RIGHT|wx.TOP|wx.BOTTOM|wx.ALIGN_CENTER_VERTICAL,
      5)
    szr2.Add(szr3)
    szr4 = wx.BoxSizer(wx.HORIZONTAL)
    txt2 = wx.StaticText(p, -1, "Color by:", size=(120,-1))
    set_tiny_bold_font(txt2)
    szr4.Add(txt2, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    color_choice = wx.Choice(p, -1,
      choices=[ label for name, label in color_modes ])
    set_tiny_font(color_choice)
    for name, label in color_modes :
      if (model.color_mode == name) :
        color_choice.SetStringSelection(label)
        break
    self.Bind(wx.EVT_CHOICE, self.OnChooseColorMode, color_choice)
    szr4.Add(color_choice, 0,
      wx.RIGHT|wx.TOP|wx.BOTTOM|wx.ALIGN_CENTER_VERTICAL, 5)
    szr2.Add(szr4)
    txt3 = wx.StaticText(p, -1, "Other display options:")
    set_tiny_bold_font(txt3)
    szr2.Add(txt3, 0, wx.ALL, 5)
    szr5 = wx.BoxSizer(wx.HORIZONTAL)
    szr5.Add((10,-1))
    szr2.Add(szr5)
    szr6 = wx.GridSizer(cols=2)
    szr5.Add(szr6)
    for name, label in draw_flags :
      box = wx.CheckBox(p, -1, label=label)
      self.Bind(wx.EVT_CHECKBOX, self.OnSetFlag, box)
      box.SetValue(getattr(model, name))
      set_tiny_font(box)
      szr6.Add(box)
    szr2.Fit(p)
    delete_btn = wx.Button(p, -1, label="Delete model")
    set_tiny_font(delete_btn)
    self.Bind(wx.EVT_BUTTON, self.OnDeleteModel, delete_btn)
    szr2.Add(delete_btn, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5)
    self.panel_sizer.Layout()
    self.sizer.Fit(self.inner_panel)
    self.Fit()

  def OnChooseModel (self, event) :
    model_id = event.GetEventObject().GetStringSelection()
    self.set_model(model_id)

  def OnSetVisibility (self, event) :
    visible = event.GetEventObject().GetValue()
    self.parent.toggle_visibility(visible, self.model_id)
    self.parent.Refresh()

  def OnChooseDrawMode (self, event) :
    mode_label = event.GetEventObject().GetStringSelection()
    model = self.parent.get_model(self.model_id)
    if (model is not None) :
      for mode, label in draw_modes :
        if (mode_label == label) :
          model.set_draw_mode(mode)
          self.parent.update_scene = True
          self.parent.Refresh()
          break

  def OnChooseColorMode (self, event) :
    mode_label = event.GetEventObject().GetStringSelection()
    model = self.parent.get_model(self.model_id)
    if (model is not None) :
      for mode, label in color_modes :
        if (mode_label == label) :
          model.set_color_mode(mode)
          self.parent.update_scene = True
          self.parent.Refresh()
          break

  def OnSetFlag (self, event) :
    box = event.GetEventObject()
    flag_value = box.GetValue()
    label = box.GetLabel()
    model = self.parent.get_model(self.model_id)
    if (model is not None) :
      for flag_name, flag_label in draw_flags :
        if (label == flag_label) :
          setattr(model, flag_name, flag_value)
          model.refresh()
          self.parent.update_scene = True
          self.parent.Refresh()
          break

  def OnClose (self, event) :
    self.Destroy()

  def OnDestroy (self, event) :
    self.parent.model_panel = None

  def OnDeleteModel (self, event) :
    self.parent.delete_model(self.model_id)
    self.refresh_model_list()
    self.parent.Refresh()
