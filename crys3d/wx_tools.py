from __future__ import absolute_import, division, print_function

import wx
from six.moves import zip

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


class EllipsizePaths(object):
  def __init__(self, parent, pathlst):
    self.csize = (300, -1)
    self.__ellipsispaths__ = [ (p, wx.Control.Ellipsize(p, wx.ClientDC(parent),
                                 wx.ELLIPSIZE_START,
                                  self.csize[0]-60)
                                )
                          for p in pathlst
                          ]
  def length(self):
    return len(self.__ellipsispaths__)
  def getpaths(self):
    return [ep[0] for ep in self.__ellipsispaths__]
  def getpath(self, i):
    return self.__ellipsispaths__[i][0]
  def getellipspaths(self):
    return [ep[1] for ep in self.__ellipsispaths__]
  def getellipspath(self, i):
    return self.__ellipsispaths__[i][1]


class ModelControlPanel (wx.MiniFrame) :
  def __init__ (self, *args, **kwds) :
    super(ModelControlPanel, self).__init__(*args, **kwds)
    self.Bind(wx.EVT_CLOSE, self.OnClose, self)
    self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy, self)
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
    self.ep = EllipsizePaths(self, self.parent.model_ids)
    chooser = wx.Choice(self.panel1, -1,
      choices = self.ep.getellipspaths() ,
      size=self.ep.csize)
    set_tiny_font(chooser)
    self.Bind(wx.EVT_CHOICE, self.OnChooseModel, chooser)
    self.model_chooser = chooser
    szr1.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    szr1.Add(chooser, 0, wx.RIGHT|wx.TOP|wx.BOTTOM|wx.ALIGN_CENTER_VERTICAL, 5)
    self.panel2 = wx.Panel(self.inner_panel, -1)
    if (self.ep.length() > 0) :
      chooser.SetSelection(0)
      self.set_model(self.ep.getpath(0))
    else :
      self.sizer.Fit(self.inner_panel)
      self.Fit()
    self.panel_sizer.Layout()


  def refresh_model_list (self) :
    self.ep = EllipsizePaths(self, self.parent.model_ids)
    models = self.ep.getpaths()
    nsel = min( self.model_chooser.GetSelection(), self.ep.length() -1 )
    current_model = self.ep.getpath(nsel) if nsel >= 0 else None
    self.model_chooser.SetItems( self.ep.getellipspaths() )
    if (current_model in models) :
      self.model_chooser.SetSelection( nsel)
    else :
      self.set_model(None)

  def set_model (self, model_id) :
    self.panel_sizer.Detach(self.panel2)
    self.panel2.Destroy()
    self.model_id = model_id
    self.panel2 = wx.Panel(self.inner_panel, -1)
    self.panel_sizer.Add(self.panel2, 1)
    if (model_id is not None) :
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
    model_id = event.GetEventObject().GetSelection()
    self.set_model(self.ep.getpath(model_id))

  def OnSetVisibility (self, event) :
    visible = event.GetEventObject().GetValue()
    self.parent.toggle_visibility(visible, self.model_id)
    self.parent.Refresh()

  def OnChooseDrawMode (self, event) :
    mode_label = str(event.GetEventObject().GetStringSelection())
    model = self.parent.get_model(self.model_id)
    if (model is not None) :
      for mode, label in draw_modes :
        if (mode_label == label) :
          model.set_draw_mode(mode)
          self.parent.update_scene = True
          self.parent.Refresh()
          break

  def OnChooseColorMode (self, event) :
    mode_label = str(event.GetEventObject().GetStringSelection())
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
    label = str(box.GetLabel())
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
    self.Refresh()

class MapControlPanel (wx.MiniFrame) :
  def __init__ (self, *args, **kwds) :
    super(MapControlPanel, self).__init__(*args, **kwds)
    self.Bind(wx.EVT_CLOSE, self.OnClose, self)
    self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy, self)
    self.parent = self.GetParent()
    self.map_id = None
    self.inner_panel = wx.Panel(self, -1)
    self.panel_sizer = wx.BoxSizer(wx.VERTICAL)
    self.inner_panel.SetSizer(self.panel_sizer)
    self.panel1 = wx.Panel(self.inner_panel, -1, style=wx.RAISED_BORDER)
    self.panel_sizer.Add(self.panel1, 0)
    self.sizer = wx.BoxSizer(wx.VERTICAL)
    self.sizer.Add(self.inner_panel, 0, wx.EXPAND)
    szr1 = wx.BoxSizer(wx.HORIZONTAL)
    self.panel1.SetSizer(szr1)
    txt1 = wx.StaticText(self.panel1, -1, "Map:")
    set_tiny_bold_font(txt1)
    self.ep = EllipsizePaths(self, self.parent.map_ids)
    chooser = wx.Choice(self.panel1, -1,
      choices = self.ep.getellipspaths() ,
      size=self.ep.csize)
    set_tiny_font(chooser)
    self.Bind(wx.EVT_CHOICE, self.OnChooseMap, chooser)
    self.map_chooser = chooser
    szr1.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    szr1.Add(chooser, 0, wx.RIGHT|wx.TOP|wx.BOTTOM|wx.ALIGN_CENTER_VERTICAL, 5)
    self.panel2 = wx.Panel(self.inner_panel, -1)
    if (len(self.parent.map_ids) > 0) :
      chooser.SetSelection(0)
      self.set_map(self.parent.map_ids[0])
    else :
      self.sizer.Fit(self.inner_panel)
      self.Fit()
    self.panel_sizer.Layout()

  def refresh_map_list (self) :
    self.ep = EllipsizePaths(self, self.parent.map_ids)
    maps = self.ep.getpaths()
    nsel = min( self.map_chooser.GetSelection(), self.ep.length() -1 )
    current_map = self.ep.getpath(nsel) if nsel >= 0 else None
    self.map_chooser.SetItems( self.ep.getellipspaths() )
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    if (current_map in maps) :
      self.map_chooser.SetSelection( nsel)
      self.set_map(current_map)
    else :
      self.set_map(None)

  def set_map (self, map_id) :
    self.panel_sizer.Detach(self.panel2)
    self.panel2.Destroy()
    self.map_id = map_id
    if (map_id is not None) :
      map_object = self.parent.get_map(map_id)
      self.panel2 = wx.Panel(self.inner_panel, -1)
      self.panel_sizer.Add(self.panel2, 1)
      p = self.panel2
      szr2 = wx.BoxSizer(wx.VERTICAL)
      self.panel2.SetSizer(szr2)
      visibility_ctrl = wx.CheckBox(p, -1,
        label="Map visible")
      visibility_ctrl.SetValue(self.parent.show_object[map_id])
      set_tiny_bold_font(visibility_ctrl)
      szr2.Add(visibility_ctrl, 0, wx.ALL, 5)
      self.Bind(wx.EVT_CHECKBOX, self.OnSetVisibility, visibility_ctrl)
      panel = self.panel2
      self.ctrls = []
      for iso_level, color in zip(map_object.iso_levels, map_object.colors) :
        level_sizer = wx.BoxSizer(wx.HORIZONTAL)
        slider = wx.Slider(panel,
                           size=(160,-1),
                           minValue=-100,
                           maxValue=100,
                           value=int(iso_level * 10),
                           style=wx.SL_AUTOTICKS)
        level_txt = wx.TextCtrl(panel,
                                size=(64,-1),
                                style=wx.TE_READONLY|wx.TE_RIGHT)
        set_tiny_font(level_txt)
        level_txt.SetValue("%.2f" % iso_level)
        initial_color = [ int(x*255) for x in color ]
        color_ctrl = wx.lib.colourselect.ColourSelect(panel,
                                                      colour=initial_color)
        level_sizer.Add(slider, 0, wx.ALL, 5)
        level_sizer.Add(level_txt, 0, wx.ALL, 5)
        level_sizer.Add(color_ctrl, 0, wx.ALL, 5)
        szr2.Add(level_sizer, 1, wx.ALL, 5)
        self.ctrls.append((slider, level_txt, color_ctrl))
      self.attach_mouse_wheel = wx.CheckBox(panel,
        label="Attach mouse wheel to this map")
      set_tiny_bold_font(self.attach_mouse_wheel)
      self.attach_mouse_wheel.SetValue(map_id == self.parent.selected_map_id)
      szr2.Add(self.attach_mouse_wheel, 0, wx.ALL, 5)
      self.Bind(wx.EVT_CHECKBOX, self.OnAttachMouse, self.attach_mouse_wheel)
      self.Bind(wx.EVT_SLIDER, self.OnUpdate)
      self.Bind(wx.lib.colourselect.EVT_COLOURSELECT, self.OnUpdate)
      delete_btn = wx.Button(p, -1, label="Delete map")
      set_tiny_font(delete_btn)
      self.Bind(wx.EVT_BUTTON, self.OnDeleteMap, delete_btn)
      szr2.Fit(panel)
      szr2.Add(delete_btn, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5)
    self.panel_sizer.Layout()
    self.sizer.Fit(self.inner_panel)
    self.Fit()

  def OnUpdate (self, event) :
    source_ctrl = event.GetEventObject()
    map_object = self.parent.get_map(self.map_id)
    for i, (s, t, c) in enumerate(self.ctrls) :
      if source_ctrl is s :
        new_level = s.GetValue() / 10.0
        map_object.iso_levels[i] = new_level
        t.SetValue("%.2f" % new_level)
        self.parent.update_maps = True
      elif source_ctrl is c :
        new_color = c.GetValue()
        new_color = [ x / 255.0 for x in new_color ]
        map_object.colors[i] = new_color
        self.parent.update_maps = True
    self.parent.OnRedrawGL()

  def refresh_iso_levels (self) :
    map_object = self.parent.get_map(self.map_id)
    levels = map_object.iso_levels
    assert len(levels) == len(self.ctrls)
    for i, (s, t, c) in enumerate(self.ctrls) :
      s.SetValue(int(10 * levels[i]))
      t.SetValue("%.2f" % levels[i])

  def refresh_attach_mouse (self) :
    self.attach_mouse_wheel.SetValue(self.map_id==self.parent.selected_map_id)

  def OnAttachMouse (self, event) :
    if self.attach_mouse_wheel.GetValue() == True :
      self.parent.set_selected_map(self.map_id)
    else :
      self.parent.set_selected_map(None)

  def OnChooseMap (self, event) :
    map_id = event.GetEventObject().GetSelection()
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    self.set_map(self.ep.getpath(map_id))

  def OnSetVisibility (self, event) :
    visible = event.GetEventObject().GetValue()
    self.parent.toggle_visibility(visible, self.map_id)
    self.parent.Refresh()

  def OnClose (self, event) :
    self.Destroy()

  def OnDestroy (self, event) :
    self.parent.map_panel = None

  def OnDeleteMap (self, event) :
    self.parent.delete_map(self.map_id)
    self.refresh_map_list()
    self.parent.Refresh()
