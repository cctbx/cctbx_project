from __future__ import absolute_import, division, print_function

'''
Author      : Lyubimov, A.Y.
Created     : 11/15/2018
Last Changed: 01/30/2019
Description : IOTA GUI base classes (with backwards compatibility for
              wxPython 3)
'''

import os
import wx
from wx.lib.buttons import GenToggleButton
from wx.lib.scrolledpanel import ScrolledPanel

from wxtbx import bitmaps
from iotbx.phil import parse

import iota.components.iota_ui_controls as ct
from iota.components.iota_utils import norm_font_size

wx4 = wx.__version__[0] == '4'

gui_phil = parse('''
gui
  .help = Options for IOTA GUI only
  .alias = GUI Options
{
  image_viewer = *dials.image_viewer cctbx.image_viewer distl.image_viewer cxi.view
    .type = choice
    .help = Select image viewer (GUI only)
    .alias = Image Viewer
  monitor_mode = False
    .type = bool
    .help = Set to true to keep watch for incoming images (GUI only)
    .alias = Process in Monitor Mode
  monitor_mode_timeout = False
    .type = bool
    .help = Set to true to auto-terminate continuous mode (GUI only)
    .alias = Monitor Mode Timeout
  monitor_mode_timeout_length = 0
    .type = int
    .help = Timeout length in seconds (GUI only)
    .alias = Timeout (sec)
}
''')

class IOTAFrameError(Exception):
  def __init__(self, msg):
    Exception.__init__(self, msg)

class IOTABaseFrame(wx.Frame):
  """ New frame that will show processing info """

  def __init__(self, parent, id, title, *args, **kwargs):
    wx.Frame.__init__(self, parent, id, title, *args, **kwargs)
    self.parent = parent

    self.main_sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.main_sizer)

  def initialize_toolbar(self):
    self.toolbar = self.CreateToolBar(style=wx.TB_3DBUTTONS | wx.TB_TEXT)

  def add_tool(self, id=wx.ID_ANY, label=None, bitmap=None,
               kind=wx.ITEM_NORMAL, shortHelp=None, longHelp=None):
    if type(bitmap) in (tuple, list):
      if len(bitmap) == 3:
        bsz = bitmap[2]
        bsc = (bitmap[2], bitmap[2])
      else:
        bsz = bsc = None
      if bitmap[0] == 'custom':
        bmp = bitmaps.fetch_custom_icon_bitmap(bitmap[1], size=bsz)
      else:
        bmp = bitmaps.fetch_icon_bitmap(bitmap[0], bitmap[1], scale=bsc)
    else:
      bmp = bitmap

    short_string = shortHelp if shortHelp else ''
    long_string = longHelp if longHelp else ''

    if wx4:
      return self.toolbar.AddTool(toolId=id,
                                  label=label,
                                  kind=kind,
                                  bitmap=bmp,
                                  shortHelp=short_string)
    else:
      return self.toolbar.AddLabelTool(id=id,
                                       label=label,
                                       kind=kind,
                                       bitmap=bmp,
                                       shortHelp=short_string,
                                       longHelp=long_string)

  def add_toolbar_separator(self, stretch=False):
    if stretch and wx4:
      self.toolbar.AddStretchableSpace()
    else:
      self.toolbar.AddSeparator()

  def realize_toolbar(self):
    self.toolbar.Realize()

  def set_tool_states(self, tools):
    for tool in tools:
      self.set_tool_state(tool=tool[0], enable=tool[1],
                          toggle=tool[2] if len(tool) == 3 else None)

  def set_tool_state(self, tool, enable=None, toggle=None):
    id = tool.GetId()
    if enable is not None:
      self.toolbar.EnableTool(id, enable)
    if toggle is not None:
      self.toolbar.ToggleTool(id, toggle)


  def place_and_size(self, set_size=False, set_by=None, center=False):
    """ Place and size the frame"""

    # Determine effective minimum size
    if set_size:
      self.SetMinSize(self.GetEffectiveMinSize())

    # Find mouse position
    if set_by == 'mouse':
      self.SetPosition(wx.GetMousePosition())
    elif set_by == 'parent':
      self.SetPosition(self.set_relative_position())
    else:
      self.SetPosition((0, 0))

    # Center on display
    if center:
      self.Center()

  def set_relative_position(self):
    """ Determines screen position w/ respect to parent window; will also
    detect if it goes beyond the display edge, and adjust """

    # Position proc window w/ respect to IOTA window
    mx, my = self.parent.GetPosition()
    px = mx + 50
    py = my + 50

    # Calculate if proc window is going out of bounds, and adjust
    disp_idx = wx.Display.GetFromWindow(self.parent)
    disp_geom = wx.Display(disp_idx).GetClientArea()
    dxmin = disp_geom[0]
    dxmax = disp_geom[0] + disp_geom[2]
    dymin = disp_geom[1]
    dymax = disp_geom[1] + disp_geom[3]

    pw, pl = self.GetSize()
    if not (px + pw * 1.1 in range(dxmin, dxmax)):
      px = dxmax - pw * 1.1
    if not (py + pl * 1.1 in range(dymin, dymax)):
      py = dymax - pl * 1.1

    return (px, py)


class IOTABasePanel(wx.Panel):
  def __init__(self, parent, *args, **kwargs):
    wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY, size=(800, 500),
                      *args, **kwargs)

    self.main_sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.main_sizer)

class BaseDialog(wx.Dialog):
  def __init__(self, parent, style=wx.DEFAULT_DIALOG_STYLE,
               label_style='bold',
               content_style='normal',
               *args, **kwargs):
    wx.Dialog.__init__(self, parent, style=style, *args, **kwargs)

    self.envelope = wx.BoxSizer(wx.VERTICAL)
    self.main_sizer = wx.BoxSizer(wx.VERTICAL)
    self.envelope.Add(self.main_sizer, 1, flag=wx.EXPAND | wx.ALL, border=5)
    self.SetSizer(self.envelope)

    if label_style == 'normal':
      self.font = wx.Font(norm_font_size, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL)
    elif label_style == 'bold':
      self.font = wx.Font(norm_font_size, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD)
    elif label_style == 'italic':
      self.font = wx.Font(norm_font_size, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_ITALIC, wx.FONTWEIGHT_NORMAL)
    elif label_style == 'italic_bold':
      self.font = wx.Font(norm_font_size, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_ITALIC, wx.FONTWEIGHT_BOLD)

    if content_style == 'normal':
      self.cfont = wx.Font(norm_font_size, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL)
    elif content_style == 'bold':
      self.cfont = wx.Font(norm_font_size, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD)
    elif content_style == 'italic':
      self.cfont = wx.Font(norm_font_size, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_ITALIC, wx.FONTWEIGHT_NORMAL)
    elif content_style == 'italic_bold':
      self.cfont = wx.Font(norm_font_size, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_ITALIC, wx.FONTWEIGHT_BOLD)


class BaseBackendDialog(BaseDialog):
  def __init__(self, parent, phil,
               backend_name = 'BACKEND',
               target=None,
               content_style='normal',
               label_style='bold',
               opt_size=(500, 500),
               phil_size=(500, 500),
               *args, **kwargs):
    BaseDialog.__init__(self, parent,
                        content_style=content_style,
                        label_style=label_style,
                        *args, **kwargs)

    self.parent = parent
    self.target_phil = target
    self.backend = backend_name
    self.params = phil.extract()
    self.opt_size = opt_size
    self.phil_size = phil_size
    self.sash_position = None

    self.splitter = wx.SplitterWindow(self, style=wx.SP_LIVE_UPDATE |
                                                  wx.SP_3DSASH |
                                                  wx.SP_NOBORDER)

    # Create options panel (all objects should be called as self.options.object)
    self.options = ScrolledPanel(self.splitter, size=self.opt_size)
    self.options_sizer = wx.BoxSizer(wx.VERTICAL)
    self.options.SetSizer(self.options_sizer)

    # Create PHIL panel
    phil_label = "{} Target Settings".format(backend_name)
    self.phil_panel = wx.Panel(self.splitter, size=self.opt_size)
    phil_box = wx.StaticBox(self.phil_panel, label=phil_label)
    self.phil_sizer = wx.StaticBoxSizer(phil_box, wx.VERTICAL)
    self.phil_panel.SetSizer(self.phil_sizer)

    # Dialog control
    self.dlg_ctr = ct.DialogButtonsCtrl(self, preset='PROC_DIALOG')

    # Splitter button
    self.btn_hide_script = GenToggleButton(self, label='Show Script >>>')
    self.show_hide_script()
    self.btn_hide_script.SetValue(False)

    self.main_sizer.Add(self.btn_hide_script, flag=wx.ALIGN_RIGHT | wx.ALL,
                        border=10)
    self.main_sizer.Add(self.splitter, 1, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.dlg_ctr,
                        flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.RIGHT,
                        border=10)

  def show_hide_script(self, initialized=False):
    if self.btn_hide_script.GetValue():
      if initialized:
        h = self.GetSize()[1]
        w = self.GetSize()[0] + self.phil_size[0]
        self.SetSize((w, h))
      self.splitter.SplitVertically(self.options, self.phil_panel)
      self.splitter.SetSashPosition(self.sash_position)
      self.phil_panel.SetSize(self.phil_size)
      self.options.SetSize(self.opt_size)
      self.btn_hide_script.SetLabel('<<< Hide Script')
    else:
      h = self.GetSize()[1]
      w = self.GetSize()[0] - self.phil_size[0]
      self.SetSize((w, h))
      self.splitter.Unsplit()
      self.phil_panel.SetSize(self.phil_size)
      self.options.SetSize(self.opt_size)
      self.btn_hide_script.SetLabel('Show Script >>>')
    self.splitter.UpdateSize()

  def get_target_file(self):
    dlg = wx.FileDialog(
      self, message="Select CCTBX.XFEL target file",
      defaultDir=os.curdir,
      defaultFile="*.phil",
      wildcard="*",
      style=wx.FD_OPEN | wx.FD_CHANGE_DIR
    )
    if dlg.ShowModal() == wx.ID_OK:
      filepath = dlg.GetPaths()[0]

      with open(filepath, 'r') as phil_file:
        phil_content = phil_file.read()
      return phil_content
    else:
      return None

  def write_default_phil(self):
    if str.lower(self.backend) in ('cctbx.xfel', 'dials'):
      method = 'cctbx.xfel'
    elif str.lower(self.backend) in ('cctbx', 'ha14', 'labelit'):
      method = 'ha14'
    else:
      method = 'current'
    from iota.components.iota_input import write_defaults
    default_phil, _ = write_defaults(method=method, write_target_file=False,
                                     write_param_file=False)
    self.target_phil = default_phil.as_str()

class BaseOptionsDialog(BaseDialog):
  ''' Test class to work out AutoPHIL, etc. '''

  def __init__(self, parent, input, *args, **kwargs):
    dlg_style = wx.CAPTION | wx.CLOSE_BOX | wx.RESIZE_BORDER | wx.STAY_ON_TOP
    BaseDialog.__init__(self, parent, style=dlg_style, *args, **kwargs)

    self.parent = parent

    self.phil_panel = PHILPanelFactory.from_scope_objects(self, scope=input)
    self.main_sizer.Add(self.phil_panel, -1, flag=wx.EXPAND)

    # Dialog control
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)
    self.main_sizer.Add(dialog_box,
                        flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                        border=10)

    self.Fit()


class PHILPanelFactory(IOTABasePanel):
  def __init__(self, parent, objects, layers=None, *args, **kwargs):
    IOTABasePanel.__init__(self, parent=parent, *args, **kwargs)

    self.layers = layers

    for obj in objects:
      if type(obj) in (list, tuple):
        pass
      elif obj.is_scope:
        self.add_scope_box(obj=obj)
      elif obj.is_definition:
        self.add_definition_control(self, obj)


  def get_all_path_names(self, phil_object, paths=None):
    if paths is None:
      paths = []
    if phil_object.is_scope:
      for object in phil_object.objects:
        paths = self.get_all_path_names(object, paths)
        paths.extend(paths)
    elif phil_object.is_definition:
      full_path = phil_object.full_path()
      if not full_path in paths:
        paths.append(full_path)
    return paths

  def add_definition_control(self, parent, obj):
    alias = obj.alias_path()
    label = alias if alias else obj.full_path().split('.')[-1]

    wdg = ct.WidgetFactory.make_widget(parent, obj, label)
    sizer = parent.GetSizer()
    sizer.Add(wdg, flag=wx.RIGHT|wx.LEFT|wx.BOTTOM|wx.EXPAND,
                        border=5)
    self.__setattr__(label, wdg)

  def add_scope_box(self, obj):
    obj_name = obj.full_path().split('.')[-1]
    label = obj.alias_path() if obj.alias_path() else obj_name

    # Make scope panel
    panel = wx.Panel(self)
    box = wx.StaticBox(panel, label=label)
    box_sz = wx.StaticBoxSizer(box, wx.VERTICAL)
    panel.SetSizer(box_sz)

    self.main_sizer.Add(panel, flag=wx.ALL|wx.EXPAND, border=10)

    # Add widgets to box (do one layer so far)
    for box_obj in obj.active_objects():
      if box_obj.is_definition:
        self.add_definition_control(panel, box_obj)

  @classmethod
  def from_scope_objects(cls, parent, scope):

    scope_type = type(scope).__name__
    if scope_type in ('list', 'tuple'):
      objects = (r for r in scope)
    elif scope_type == 'scope':
      objects = scope.active_objects()
    else:
      objects = scope

    return cls(parent, objects)

  @classmethod
  def from_filename(cls, filepath):
    pass
