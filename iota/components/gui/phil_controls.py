from __future__ import division, print_function, absolute_import

'''
Author      : Lyubimov, A.Y.
Created     : 04/02/2019
Last Changed: 12/02/2019
Description : IOTA GUI controls for PHIL-formatted settings
'''

import os
import sys
import wx
from wx.lib.scrolledpanel import ScrolledPanel

from wxtbx import bitmaps

from libtbx.utils import Sorry
from libtbx import Auto

from iota.components import gui
from iota.components.gui import base
from iota.components.gui import controls as ct
from iota.components.iota_utils import InputFinder, makenone
from iota.components.gui.dialogs import DirView

ginp = InputFinder()

def get_test_phil():
  from iotbx.phil import parse

  test_phil = """
    string_definition = None
      .help = A string definition in the main scope
      .optional = False
      .type = str
      .alias = String Definition
    path_to_folder = $PWD
      .help = A path to a folder in the main scope
      .optional = False
      .type = path
      .multiple = False
      .alias = Main Folder
      .style = path:folder
    input_file = None
      .help = A path to a file (multiple)
      .optional = True
      .type = path
      .multiple = True
      .alias = Input File
      .style = input_list
    multi_string = None
      .help = A string definition that is multiple
      .type = str
      .multiple = True
      .alias = Multi-string
    child_scope_alpha
      .help = The first child scope
      .alias = Child Alpha
    {
      flag_on = False
        .help = toggle this scope on and off
        .type = bool
        .alias = Use Child Scope Alpha
        .style = scope_switch
      float_1_def = None
        .help = A float definition in the child scope alpha
        .type = float
        .alias = First Float
        .optional = True
      integer_2_def = None
        .help = An integer definition in the child scope alpha
        .type = int
        .alias = First Integer
        .optional = False
      bool_3_def = True
        .help = A boolean definition in the child scope alpha
        .type = bool
        .alias = True or False in Alpha
      choice_4_def = one two three *four five
        .help = A choice definition in the child scope alpha
        .type = choice
        .alias = Number of Items
      }
    child_scope_beta
      .help = The second child scope
      .alias = Child Beta
    {
      string_1_def = None
        .help = A string definition in the child scope beta
        .type = str
        .alias = String Beta
      space_group_def = None
        .help = A space group definition in the child scope beta
        .type = space_group
        .alias = Space Group
      unit_cell_def = None
        .help = A unit cell definition in the child scope beta
        .type = unit_cell
        .alias = Unit Cell
      grandchild_scope_beta
        .help = a multiple grandchild scope in child scope beta
        .alias = Grandchild Beta
        .multiple = True
      {
        flag_on = False
          .help = a flag_on checkbox (should turn whole scope on/off)
          .type = bool
          .alias = Activate Grandchild Beta
          .style = scope_switch
        string_1b_beta = None
          .help = A string in grandchild beta
          .type = str
          .alias = String 1B Beta
        float_2b_beta = None
          .help = A float in grandchild beta
          .type = float
          .optional = False
          .alias = Float 1B Beta
        choice_3b_beta = one *several many
          .help = A choice in grandchild beta
          .type = choice
          .optional = True
          .alias = Choice 3B Beta
      }
    }
    child_scope_gamma
      .help = Third child scope (test multiple panel)
      .multiple = False
      .alias = Child Gamma
    {
      apply_gamma = False
        .help = Boolean definition for child scope gamma
        .type = bool
        .alias = Apply
      string_1_gamma = one some many
        .help = A string definition in child scope gamma
        .type = choice
        .alias = String One Gamma
      string_2_gamma = None
        .help = A string definition in child scope gamma
        .type = str
        .alias = String Two Gamma
      string_3_gamma = None
        .help = A string definition in child scope gamma
        .type = str
        .alias = String Three Gamma
    }


"""
  return parse(test_phil)


# Platform-specific stuff
# TODO: Will need to test this on Windows at some point
if wx.Platform == '__WXGTK__':
  norm_font_size = 10
  button_font_size = 12
  LABEL_SIZE = 14
  CAPTION_SIZE = 12
elif wx.Platform == '__WXMAC__':
  norm_font_size = 12
  button_font_size = 14
  LABEL_SIZE = 14
  CAPTION_SIZE = 12
elif (wx.Platform == '__WXMSW__'):
  norm_font_size = 9
  button_font_size = 11
  LABEL_SIZE = 11
  CAPTION_SIZE = 9

# Metallicbutton globals
GRADIENT_NORMAL = 0
GRADIENT_PRESSED = 1
GRADIENT_HIGHLIGHT = 2

MB_STYLE_DEFAULT = 1
MB_STYLE_BOLD_LABEL = 2
MB_STYLE_DROPARROW = 4

# --------------------------- PHIL-specific Mixins --------------------------- #


class PHILPanelMixin(object):
  """ PHIL-handling mixin for PHIL panels """

  def initialize_phil_panel(self, parent, box=None, direction=wx.VERTICAL):
    # Establish the top parent window
    # self.window = getattr(parent, 'window', None)
    # if not self.window:
    #   self.window = self.GetTopLevelParent()
    self.window = self.GetTopLevelParent()
    self.parent = parent
    self.multi_scope = False
    self.multi_definition = False
    self._multiples = []
    self._input_lists = {}
    self._toggled_scopes = {}
    self.scope_switch = None
    self.control_index = {}

    # Create a box around panel
    if box is not None:
      assert type(box) == str
      panel_box = wx.StaticBox(self, label=box)
      self.main_sizer = PHILBoxSizer(self, panel_box, direction)
    else:
      self.main_sizer = PHILSizer(self, direction)
    self.SetSizer(self.main_sizer)

  def set_phil_index(self, master_phil, working_phil=None, fetch_new=True):
    self.phil_index = gui.PHILIndex(master_phil=master_phil,
                                    working_phil=working_phil,
                                    fetch_new=fetch_new)

  def redraw_by_expert_level(self, expert_level=0):
    self.expert_level = expert_level
    self.show_hide_controls(expert_level=self.expert_level)

  def show_hide_controls(self, expert_level=0):
    if hasattr(self, 'controls') and self.controls:
      for idx, ctrl in self.controls.items():
        if ctrl.expert_level > expert_level:
          ctrl.Hide()
        else:
          ctrl.Show()
          if ctrl.is_scope:
            try:
              ctrl.show_hide_controls(expert_level=expert_level)
            except Exception as e:
              raise e
    else:
      pass
      # raise Sorry('IOTA PHIL error: no controls in {}'.format(self))

  def flatten_scope(self, scope):
    saved_values = {}

    if isinstance(scope, list):
      objects = scope
    else:
      objects = scope.active_objects()

    for obj in objects:
      if obj.is_definition:
        path = obj.full_path()
        value = obj.type.from_words(obj.words, obj)
        saved_values[path] = str(value)
      elif obj.is_scope:
        try:
          from_scope = self.flatten_scope(scope=obj)
        except RuntimeError as e:
          raise e
        saved_values.update(from_scope)
    return saved_values

  def check_full_path(self, full_path=None):
    if not full_path:
      if hasattr(self, 'full_path'):
        return self.full_path
      else:
        return
    else:
      return full_path

  def mark_non_defaults(self, full_path=None):
    """ Loop through controls and for text controls set background to amber
        if the value is not a default """
    full_path = self.check_full_path(full_path)
    if full_path is None:
      return

    # Create a dictionary of default values
    defaults = self.get_default_values(full_path)
    if not defaults:
      return

    # Iterate through controls and set background as appropriate
    for idx, ctrl in self.controls.items():
      if ctrl.is_scope:
        ctrl.mark_non_defaults(full_path=ctrl.full_path)
      elif ctrl.is_definition:
        if ctrl.multi_definition:
          for i, c in ctrl.controls.items():
            path = c.full_path
            if path and path in defaults:
              # c.default_value = defaults[path]
              c.set_background()
        else:
          path = ctrl.full_path
          if path and path in defaults:
            # ctrl.default_value = defaults[path]
            ctrl.set_background()

  def get_default_scope(self, full_path=None):
    """ Extract the master version of the given scope (the full scope
        associated with this panel, or any specified scope) """
    full_path = self.check_full_path(full_path)
    if full_path is not None:
      master_scopes = self.phil_index.get_master_scope(full_path)
      if isinstance(master_scopes, list):
        return master_scopes[0]
      else:
        return master_scopes

  def get_default_values(self, full_path=None):
    """ Obtain a dictionary of path and value from default scope """
    full_path = self.check_full_path(full_path)
    if full_path is not None:
      scope = self.get_default_scope(full_path=full_path)
      if scope:
       dv = self.flatten_scope(scope=[scope])
       return dv

  def get_max_label_size(self, scope=None, scopes=None):
    # Get font info for string-to-pixels conversion
    panel_font = self.GetFont()
    dc = wx.WindowDC(self)
    dc.SetFont(panel_font)

    # Identify the longest label
    max_label_size = 0

    source = scope if scope else scopes
    assert source

    if type(source).__name__ == 'scope':
      active_objects = source.active_objects()
    else:
      try:
        active_objects = list(source)
      except TypeError:
        raise Sorry('IOTA PHIL Error: input is {}, must be scope, list, '
                    'or tuple'.format(type(source).__name__))
      except Exception as e:
        raise e

    for o in active_objects:
      if (
              (o.is_definition and o.type.phil_type != 'bool') or
              (o.is_scope and (o.style and 'grid' in o.style))
      ):
        alias = o.alias_path()
        label = alias if alias else o.full_path().split('.')[-1] + ': '
        label_size = dc.GetTextExtent(label)[0]
        if label_size > max_label_size:
          max_label_size = label_size
    # Pad max label size
    max_label_size += dc.GetTextExtent('   ')[0]
    if max_label_size > 0:
      return wx.Size(max_label_size, -1)
    else:
      return wx.DefaultSize

  # Scope widget tools ------------------------------------------------------- #

  def add_scope_box(self, scope, sizer=None, label=None,
                    index=None, border=10, flag=wx.EXPAND | wx.BOTTOM,
                    ignore_multiple=False):
    btn_scopes = []
    if index is None:
      index = len(self.controls)
    if sizer is None:
      sizer = self.main_sizer

    # Generate name for the scope box
    if not label:
      obj_name = scope.full_path().split('.')[-1]
      box_label = scope.alias_path() if scope.alias_path() else obj_name
    else:
      box_label = label

    # Make scope box, add to sizer, and include in controls dictionary
    if scope.multiple and not ignore_multiple:
      if scope.full_path() not in self._multiples:
        multi_scopes = self.phil_index.get_scopes(scope.full_path())
        panel = PHILMultiScopePanel(self, multi_scopes, box=box_label)
        sizer.Add(panel, 1, flag=flag, border=border)
        self.controls.update({index:panel})
        self.control_index.update({scope.full_path():panel})
        self._multiples.append(scope.full_path())
      return
    else:
      panel = PHILScopePanel(self, scope, box=box_label)
      sizer.Add(panel, flag=flag, border=border)
      self.controls.update({index:panel})
      self.control_index.update({scope.full_path(): panel})

    # Iterate through PHIL objects; if scope is found, it's added to the list
    # of scopes that will become buttons
    for obj in scope.active_objects():
      max_label_size = self.get_max_label_size(scopes=scope)
      if obj.is_scope:
        # only collect a single example of a full path to a multiple object
        if obj.full_path() not in self._multiples:
          # check if this is a grid-style scope
          style = self.phil_index.get_scope_style(obj.full_path())
          if style.grid or len(obj.objects) <= 3:
            self.add_scope_grid(parent=panel, scope=obj, style=style,
                                label_size=max_label_size)
          else:
            btn_scopes.append(obj)
        if obj.multiple:
          self._multiples.append(obj.full_path())
      elif obj.is_definition:
        self.add_definition_control(parent=panel, obj=obj,
                                    label_size=max_label_size)

    # If any scopes were found, generate buttons all at once, and place at
    # bottom of the scope box
    if btn_scopes:
      self.add_scope_buttons(parent=panel, scopes=btn_scopes)

  def add_definition_control(self, parent, obj,
                             sizer=None,
                             proportion=0,
                             border=5,
                             in_grid=False,
                             ignore_multiple=False,
                             label=None,
                             label_size=wx.DefaultSize):
    sizer = sizer if sizer else parent.main_sizer
    value = obj.type.from_words(obj.words, obj)
    style = self.phil_index.style[obj.full_path()]
    idx = len(parent.controls)

    # Handle multiple definitions
    if obj.multiple:
      if ignore_multiple:
        sizer = parent.scope_sizer
      elif style.input_list:
        if obj.full_path() not in self._input_lists:
          sizer = parent.main_sizer
        else:
          self._input_lists[obj.full_path()].add_item(path=value)
          return
      else:
        if obj.full_path() not in self._multiples:
          try:
            multi_defs = self.phil_index.get_scopes(obj.full_path())
            multi_panel = PHILMultiDefPanel(parent, multi_defs, style=style)
            sizer.Add(multi_panel, flag=wx.EXPAND | wx.BOTTOM, border=10)
            parent.controls[idx] = multi_panel
            parent.control_index[obj.full_path()] = multi_panel
          except RuntimeError:
            pass
          self._multiples.append(obj.full_path())
        return

    # Create widget and add to sizer / controls dictionary
    extras = getattr(self, '_{}_extras'.format(obj.type.phil_type), None)
    wdg = WidgetFactory.make_widget(parent, obj, label_size, value=value,
                                    border=border, style=style, label=label,
                                    extras=extras)

    # Set widget_specific formatting
    if obj.type.phil_type in ('str', 'strings', 'unit_cell', 'space_group',
                              'path'):
      expand = True
      if style.input_list:
        proportion = 1
    else:
      expand = False

    # Add widget to sizer and dictionaries
    parent.controls[idx] = wdg
    parent.control_index[obj.full_path()] = wdg
    if style.input_list:
      self._input_lists[obj.full_path()] = wdg
    if style.scope_switch:
      setattr(parent, 'scope_switch', wdg)
      if in_grid:
        parent.setup_switch(switch_ctrl=wdg)
        return
    sizer.add_widget(widget=wdg, expand=expand, proportion=proportion,
                     border=border)

  def add_scope_buttons(self, parent, scopes):
    parent.main_sizer.Add((15, 0))
    btn_sizer = PHILSizer(parent=parent, direction=wx.HORIZONTAL)
    for scope in scopes:
      name = scope.full_path().split('.')[-1]
      btn = PHILDialogButton(parent=parent,
                             scope=scope,
                             phil_index=self.phil_index,
                             name=name,
                             expert_level=scope.expert_level)
      parent.controls.update({len(parent.controls):btn})
      parent.control_index.update({scope.full_path():btn})
      btn_sizer.add_widget(btn)
    parent.main_sizer.add_widget(btn_sizer, expand=True)

  def add_scope_grid(self, parent, scope, style, label_size=wx.DefaultSize):
    one_line = (len(scope.objects) == 2 and style.has_scope_switch) or \
               (len(scope.objects) == 1)

    if one_line:
      grid_label_size = label_size
    else:
      grid_label_size = self.get_max_label_size(scopes=scope)

    panel = PHILGridScopePanel(parent, scope, style, label_size=label_size)
    for obj in scope.active_objects():
      if obj.is_definition:
        obj_style = self.phil_index.get_scope_style(obj.full_path())
        if obj_style.scope_switch:
          label = ''
          dfn_label_size = (0, -1)
        else:
          label = None
          dfn_label_size = grid_label_size
        self.add_definition_control(panel, obj, label=label,
                                    label_size=dfn_label_size,
                                    sizer=panel.main_sizer,
                                    in_grid=True)
      else:
        # todo: handle scopes here!
        pass
    parent.main_sizer.Add(panel, flag=wx.EXPAND)
    idx = len(parent.controls)
    parent.controls[idx] = panel
    parent._toggled_scopes.update(self._toggled_scopes)

  def rebuild_panel(self, scope=None):
    if not scope:
      scope = self.scope
    self.clear_panel()
    self.construct_panel(scope=scope)

  def redraw_panel(self, panel=None, reset=False, exempt=None):
    if not panel:
      panel = self
    if exempt is None:
      exempt = []
    for idx, ctrl in panel.controls.items():
      style = self.phil_index.get_scope_style(scope_name=ctrl.full_path)
      if ctrl.is_definition:
        if reset and ctrl.full_path not in exempt:
          reset_ctrl = True
        else:
          reset_ctrl = False
        if ctrl.multi_definition:
          ctrl.redraw_dialog(reset_to_default=reset_ctrl)
        else:
          if reset_ctrl:
            value = ctrl.default_value
          else:
            value = self.phil_index.get_value(path=ctrl.full_path)
          if style.input_list:
            ctrl.ResetValue(value)
          else:
            ctrl.SetValue(value)
          if style.scope_switch:
            value = bool(value)
            ctrl.ctr.SetValue(value)
            ctrl.parent.check_scope_switches()
      else:
        if ctrl.full_path in exempt:
          continue
        if ctrl.multi_scope:
          ctrl.redraw_dialog(reset_to_default=reset)
        ctrl.redraw_panel(reset=reset)
    self.mark_non_defaults()

  def construct_panel(self, scope=None):
    if not scope:
      scope = self.scope

    if isinstance(scope, list):
      phil_objects = scope
    else:
      try:
        assert scope.is_scope
      except AssertionError:
        raise Sorry('IOTA PHIL ERROR: Scope object required, {} received'
                    ''.format(type(scope).__name__))
      else:
        phil_objects = scope.objects

    for obj in phil_objects:
      max_label_size = self.get_max_label_size(scopes=phil_objects)
      if obj.is_scope:
        style = self.phil_index.get_scope_style(obj.full_path())
        if style.grid:
          self.add_scope_grid(self, scope=obj, style=style,
                              label_size=max_label_size)
        else:
          self.add_scope_box(scope=obj)
      elif obj.is_definition:
        self.add_definition_control(parent=self, obj=obj,
                                    label_size=max_label_size)

    self.check_scope_switches()
    self.Layout()

    # Mark widgets with non-default values
    self.mark_non_defaults(self.full_path)

  def check_scope_switches(self, force_switch=False):
    # Go through all scope boxes and disables ones that are turned off
    if self.scope_switch is not None or force_switch:
      for idx, ctrl in self.controls.items():
        enable = self.scope_switch.ctr.GetValue()

        if ctrl != self.scope_switch:
          if ctrl.is_scope:
            ctrl.check_scope_switches(force_switch=True)
          else:
            ctrl.enable_panel(enable=enable)
    else:
      for idx, ctrl in self.controls.items():
        if ctrl.is_scope:
          ctrl.check_scope_switches()

  def clear_panel(self, panel=None):
    if panel is None:
      panel = self
    self._multiples = []
    self._input_lists = {}
    self._toggled_scopes = {}
    self.controls = {}
    panel.main_sizer.DeleteWindows()

  def enable_panel(self, enable=True, children=None):
    if not children:
      if hasattr(self, 'GetChildren'):
        children = self.GetChildren()
    if not children:
      return
    for child in children:
      if 'PHIL' in child.__class__.__name__:
        self.enable_panel(enable=enable, children=child.GetChildren())
      child.Enable(enable=enable)

  def collect_errors(self, panel=None):
    """ Go through all controls recursively and collect any format errors """

    if panel is None:
      panel = self

    errors = {}
    for idx, ctrl in panel.controls.items():
      if ctrl.is_definition:
        if ctrl.multi_definition:
          multi_def_errors = ctrl.collect_errors()
          if multi_def_errors:
            errors.update(multi_def_errors)
        else:
          if hasattr(ctrl.ctr, 'error_msg') and ctrl.ctr.error_msg:
            errors[ctrl.name] = ctrl.ctr.error_msg
      elif ctrl.is_scope:
        scope_errors = ctrl.collect_errors()
        if scope_errors:
          errors.update(scope_errors)
      else:
        return None

    return errors

  def change_value(self, full_path, value):
    # NOTE: won't work with multiples!
    if full_path in self.control_index:
      if self.control_index[full_path].is_definition:
        self.control_index[full_path].SetValue(value)

  def get_value(self, full_path):
    if full_path in self.control_index and \
            self.control_index[full_path].is_definition:
        return self.control_index[full_path].GetStringValue()
    return None


class MultiObjectPanelMixin(object):
  """ Control-handing mixin for multi-scope and multi-definition panels """

  def redraw_dialog(self, remove_controls=None, reset_to_default=False):
    if reset_to_default:
      self.clear_boxes()
      self.add_default()
      return
    elif remove_controls:
      for ctrl in remove_controls:
        self.main_sizer.Detach(ctrl)
        ctrl.Destroy()
    for idx, ctrl in self.controls.items():
      ctrl.toggle.Hide()
    self.dialog_layout()

  def clear_boxes(self):
    self.scope_sizer.DeleteWindows()

  def reset_boxes(self):
    self._current_objects = {}
    self.controls = {}
    self.redraw_dialog(reset_to_default=True)

  def delete_boxes(self, last_only=False, idx=None):
    if last_only:
      last_index = max([i for i in self._current_objects])
      self._current_objects.pop(last_index)
      selected_controls = [self.controls.pop(last_index)]
    elif idx:
      self._current_objects.pop(idx)
      selected_controls = [self.controls.pop(idx)]
    else:
      selected_controls = []
      for idx, ctrl in self.controls.items():
        if ctrl.selected:
          self._current_objects.pop(idx)
          ctrl = self.controls.pop(idx)
          selected_controls.append(ctrl)
    self.redraw_dialog(remove_controls=selected_controls,
                       reset_to_default=(not self.controls))

# -------------------------- PHIL Sizers and Panels -------------------------- #

class PHILSizer(wx.BoxSizer, gui.WidgetHandlerMixin):
  def __init__(self, parent, direction=wx.VERTICAL):
    super(PHILSizer, self).__init__(direction)
    self.parent = parent


class PHILBoxSizer(wx.StaticBoxSizer, gui.WidgetHandlerMixin):
  def __init__(self, parent, box, direction=wx.VERTICAL):
    super(PHILBoxSizer, self).__init__(box, direction)
    self.parent = parent


class PHILFlexGridSizer(wx.FlexGridSizer, gui.WidgetHandlerMixin):
  def __init__(self, parent, rows, cols, vgap, hgap):
    super(PHILFlexGridSizer, self).__init__(rows, cols, vgap, hgap)
    self.parent = parent
    self.label = None

  def add_growable(self, cols=None, rows=None, proportion=0):
    if cols:
      for col in cols:
        self.AddGrowableCol(idx=col, proportion=proportion)
    if rows:
      for row in rows:
        self.AddGrowableRow(idx=row, proportion=proportion)


class PHILBaseScrolledPanel(ScrolledPanel, PHILPanelMixin):
  def __init__(self, parent, box=None, direction=wx.VERTICAL, *args, **kwargs):
    if 'phil_index' in kwargs:
      phil_index = kwargs.pop('phil_index')
      self.phil_index = phil_index
    else:
      self.phil_index = getattr(parent, 'phil_index', None)
    ScrolledPanel.__init__(self, parent, *args, **kwargs)
    self.initialize_phil_panel(parent, box=box, direction=direction)


class PHILBaseFixedPanel(wx.Panel, PHILPanelMixin):
  def __init__(self, parent, box=None, direction=wx.VERTICAL, *args, **kwargs):
    if 'phil_index' in kwargs:
      phil_index = kwargs.pop('phil_index')
      self.phil_index = phil_index
    else:
      self.phil_index = getattr(parent, 'phil_index', None)
    wx.Panel.__init__(self, parent, *args, **kwargs)
    self.initialize_phil_panel(parent, box=box, direction=direction)


class PHILBaseDialogPanel(PHILBaseScrolledPanel, gui.IOTAScopeCtrl):
  def __init__(self, parent, scope, box=None, direction=wx.VERTICAL,
               *args, **kwargs):
    super(PHILBaseDialogPanel, self).__init__(parent, box=box,
                                              direction=direction,
                                              *args, **kwargs)
    gui.IOTAScopeCtrl.__init__(self, scope)

    # Set expert level
    self.expert_level = self.phil_index.get_min_expert_level(scope)
    if self.expert_level is None:
      self.expert_level = 0


class PHILBaseScopePanel(PHILBaseFixedPanel, gui.IOTAScopeCtrl):
  def __init__(self, parent, scope, box=None, direction=wx.VERTICAL,
               *args, **kwargs):
    super(PHILBaseScopePanel, self).__init__(parent, box=box,
                                             direction=direction,
                                             *args, **kwargs)
    gui.IOTAScopeCtrl.__init__(self, scope)
    if not self.phil_index:
      self.phil_index = self.parent.phil_index

    # Set expert level
    self.expert_level = self.phil_index.get_min_expert_level(scope)
    if self.expert_level is None:
      self.expert_level = 0


class PHILBaseDefPanel(wx.Panel, gui.IOTADefinitionCtrl, PHILPanelMixin):
  def __init__(self, parent, phil_object, box=None, direction=wx.VERTICAL,
               *args, **kwargs):
    wx.Panel.__init__(self, parent, *args, **kwargs)
    gui.IOTADefinitionCtrl.__init__(self, phil_object)
    self.multi_definition = False
    self.color_change = True

    # Initialize panel
    self.initialize_phil_panel(parent, box=box, direction=direction)

  def GetStringValue(self):
    """ Override in subclasses with subclass-specific function """
    raise NotImplementedError()

  def SetValue(self, value):
    """ Overwrite in subclasses with subclass-specific function """
    raise NotImplementedError()


class PHILBaseDialog(base.FormattedDialog):
  """ Base dialog class for PHIL-formatted settings """
  def __init__(self, parent, *args, **kwargs):
    style = wx.CAPTION | wx.CLOSE_BOX | wx.RESIZE_BORDER | wx.STAY_ON_TOP
    super(PHILBaseDialog, self).__init__(parent, style=style, *args, **kwargs)
    self.phil_sizer = PHILSizer(self)
    self.envelope.Add(self.phil_sizer, 1, flag=wx.EXPAND | wx.ALL, border=5)

  def _set_default_size(self, disp_size):
    # Get x, y of Display as well as of Dialog
    sw, sh = disp_size
    dw, dh = self.GetSize()

    # Set dialog height (soft maximum to 2/3 of display height)
    soft_max_height = sh * (2/3)
    dlg_height = dh if dh <= soft_max_height else soft_max_height

    # Set dialog width (soft maximum to 1/3 of display width, pad if height
    # has been truncated, to accommodate the scroll bar)
    soft_max_width = sw / 3
    dw = dw + 20 if dlg_height < dh else dw
    dlg_width = dw if dw <= soft_max_width else soft_max_width

    # Apply size restrictions (wide enough to accommodate the button bar) and
    # set dialog size
    if sw > 500 and sh > 300:
      self.SetMinSize((500, 300))
    self.SetSize(wx.Size(dlg_width, dlg_height))

  def _place_near_parent_window(self, disp_geom):
    # initially set coordinates to be slightly offset from the top window
    wx, wy = self.parent.GetTopLevelParent().GetPosition()
    dx, dy, dw, dh = disp_geom
    x = wx + int(0.025 * dw)
    y = wy + int(0.025 * dh)

    # check for boundary situations and adjust
    w, h = self.GetSize()
    if x + w > dx + dw:
      x = dx + dw - w
    if y + h > dy + dh:
      y = dy + dh - h

    self.SetPosition((x, y))

  def size_and_place(self):
    self.set_size(set_size='dialog')
    x, y = self.set_relative_position()
    self.SetPosition((x, y))

  def OnCancel (self, event):
    self.EndModal(wx.ID_CANCEL)

  def OnExpertLevel(self, event):
    expert_level = self.dlg_ctr.choice.GetSelection()
    self.phil_panel.redraw_by_expert_level(expert_level=expert_level)
    self.Layout()
    self.phil_panel.SetupScrolling()

# ------------------------------- PHIL Widgets ------------------------------- #


class PHILDefPanel(PHILBaseDefPanel):
  """ Base class for the PHIL control, subclassed from wx.Panel and
      IOTADefinitionCtrl. The panel's main sizer is a FlexGridSizer, with a
      customizable grid depending on which kind of control is desired. """

  def __init__(self, parent, phil_object, rows=1, cols=2, vgap=10, hgap=10,
               *args, **kwargs):
    """ Constructor
    :param parent: parent object, typically another panel
    :param rows: Number of rows in the FlexGridSizer
    :param cols: Number of columns in the FlexGridSizer
    :param vgap: gap (in pixels) between rows (set to zero if rows=1)
    :param hgap: gap (in pixels) between columns (set to zero if cols=1)
    :param size: size of the panel
    """

    super(PHILDefPanel, self).__init__(parent, phil_object=phil_object)

    self.error_btn = None
    self.ctr = None
    self.selected = False

    # Set grid
    vgap = vgap if rows > 1 else 0
    hgap = hgap if cols > 1 else 0

    self.ctrl_sizer = PHILFlexGridSizer(self, rows, cols+3, vgap, hgap)
    self.SetSizer(self.ctrl_sizer)

    # Attach and hide a checkbox for use with multi-definition widgets
    self.toggle = wx.CheckBox(self, label='')
    self.ctrl_sizer.Add(self.toggle, flag=wx.ALIGN_CENTER)
    self.Bind(wx.EVT_CHECKBOX, self.onToggle, self.toggle)
    self.toggle.Hide()

    # Attach and hide a small button that would show the format error message
    err_bmp = bitmaps.fetch_icon_bitmap('actions', 'status_unknown', size=16)
    self.error_btn = ct.GradButton(parent=self, bmp=err_bmp, button_margin=1,
                                   size=(22, 22), gradient_percent=0)
    self.ctrl_sizer.add_widget(self.error_btn)
    self.Bind(wx.EVT_BUTTON, self.ShowError, self.error_btn)
    self.error_btn.Hide()

  def onToggle(self, e):
    self.selected = self.toggle.GetValue()

  def GetLabelSize(self):
    size = self.ctrl_sizer.label.GetSize() if self.ctrl_sizer.label else (0, 0)
    return size

  def SetLabelSize(self, size=wx.DefaultSize):
    if self.ctrl_sizer.label:
      self.ctrl_sizer.label.SetSize(size)

  def GetStringValue(self):
    """ Override in subclasses with subclass-specific function """
    raise NotImplementedError()

  def SetValue(self, value):
    """ Overwrite in subclasses with subclass-specific function """
    raise NotImplementedError()

  def SetError(self, err):
    self.error_btn.user_data = err
    self.error_btn.Show()
    self.set_background(is_error=True)
    self.parent.Layout()
    self.Refresh()

  def RemoveError(self):
    self.error_btn.Hide()
    self.set_background(is_error=False)
    self.parent.Layout()
    self.Refresh()

  def ShowError(self, e):
    err = self.error_btn.user_data
    wx.MessageBox(caption='Format Error!',
                  message=str(err),
                  style=wx.ICON_EXCLAMATION|wx.STAY_ON_TOP)

  def set_background(self, is_error=False, force_null=False):
    c_err = (215, 48, 31)
    c_new = (254, 240, 217)

    if force_null:
      color = wx.NullColour
    elif is_error:
      color = c_err
    else:
      if self.is_default():
        color = wx.NullColour
      else:
        color = c_new

    self.SetBackgroundColour(color)
    self.Refresh()

  def is_default(self):
    default_value = str(self.default_value)
    control_value = str(self.GetStringValue())
    return control_value == default_value


class PHILScopePanel(PHILBaseScopePanel):
  """ Based class for a non-scrolled panel to display PHIL scopes; to be
      instantiated for scope boxes """

  def __init__(self, parent, scope, direction=wx.VERTICAL, box=None,
               *args, **kwargs):
    super(PHILScopePanel, self).__init__(parent, scope, box=box,
                                         direction=direction, *args, **kwargs)
    self.multi_scope = False

    self.window = getattr(parent, 'window', None)
    if not self.window:
      self.window = self.GetTopLevelParent()
    self.parent = parent
    self.selected = False

    # Set up panel sizer
    self.panel_sizer = PHILFlexGridSizer(self, 1, 2, 0, 10)
    self.panel_sizer.AddGrowableCol(1)
    self.panel_sizer.AddGrowableRow(0)

    # Set checkbox (hidden for most PHIL Scope Panels)
    self.toggle = wx.CheckBox(self, label='')
    self.panel_sizer.Add(self.toggle, flag=wx.ALIGN_CENTER)
    self.Bind(wx.EVT_CHECKBOX, self.onToggle, self.toggle)
    self.toggle.Hide()

    # Set up main PHIL sizer
    if box:
      assert type(box) == str
      panel_box = wx.StaticBox(self, label=box)
      self.main_sizer = PHILBoxSizer(self, panel_box, direction)
    else:
      self.main_sizer = PHILSizer(self, direction)
    self.panel_sizer.Add(self.main_sizer, 1, flag=wx.EXPAND)
    self.SetSizer(self.panel_sizer)

    self.Layout()

  def onToggle(self, e):
    self.selected = self.toggle.GetValue()


class PHILGridScopePanel(PHILBaseScopePanel):
  def __init__(self, parent, scope, style,
               label_size=wx.DefaultSize,
               *args, **kwargs):
    super(PHILGridScopePanel, self).__init__(parent, scope, *args, **kwargs)
    self.multi_scope = False
    self.window = getattr(parent, 'window', None)
    if not self.window:
      self.window = self.GetTopLevelParent()
    self.parent = parent
    self.selected = False
    self.switch_ctrl = None

    n_items = len(scope.objects)
    has_switch = style.has_scope_switch
    if not has_switch:
      bool_objects = [(o.style and 'scope_switch' in o.style) for o in
                      scope.objects if o.is_definition]
      has_switch = True in bool_objects

    if not style.grid:
      grid_style = 'auto'
    else:
      grid_style = style.grid

    if grid_style in ('none', 'auto'):
      rows = n_items
      cols = 1
    else:
      delimiters = ['x', ':', '-', ',']
      bool_d = [(d in grid_style) for d in delimiters]
      d_idx = bool_d.index(True)
      delimiter = delimiters[d_idx]
      rows, cols = grid_style.split(delimiter)

    # Populate grid
    self.main_sizer = PHILFlexGridSizer(self, int(rows), int(cols), 0, 10)
    self.main_sizer.add_growable(cols=range(int(cols)))

    if (has_switch and n_items == 2) or n_items == 1:     # No label
      self.panel_sizer = wx.GridBagSizer(0, 0)
      self.panel_sizer.Add(self.main_sizer, flag=wx.EXPAND, pos=(0, 1))
      self.panel_sizer.AddGrowableCol(1)
    else:                               # Apply label
      self.panel_sizer = wx.GridBagSizer(10, 10)
      label = scope.alias
      if not label:
        if scope.name:
          label = scope.name.capitalize()
        elif scope.full_path():
          label = scope.full_path().split('.')[-1].replace('_',
                                                           ' ').capitalize()
        else:
          label = 'Options'
      txt_label = wx.StaticText(self, label=label, size=label_size)
      if has_switch:
        self.panel_sizer.Add(txt_label, pos=(0, 1))
      else:
        self.panel_sizer.Add(txt_label, pos=(0, 0), span=(1, 2))
        self.panel_sizer.Add((25, 0), pos=(1, 0))
      self.panel_sizer.Add(self.main_sizer, flag=wx.EXPAND, pos=(1, 1))
      self.panel_sizer.AddGrowableCol(1)

    self.SetSizer(self.panel_sizer)

  def setup_switch(self, switch_ctrl, value=None):
    self.switch_ctrl = switch_ctrl
    self.panel_sizer.Add(self.switch_ctrl, pos=(0, 0))
    if not value:
      value = self.switch_ctrl.ctr.GetValue()
    self.switch_ctrl.ctr.SetValue(state=value)


class PHILMultiScopePanel(PHILBaseDialogPanel, MultiObjectPanelMixin):
  def __init__(self, parent, scope, box=None, direction=wx.VERTICAL,
               *args, **kwargs):
    super(PHILMultiScopePanel, self).__init__(parent, scope, box=box,
                                              direction=direction,
                                              *args, **kwargs)
    self.multi_scope = True
    self._current_objects = {}

    # Add a sizer for multiple scopes, add existing scopes to it
    self.scope_sizer = PHILSizer(self)

    # Create a dictionary of PHIL objects
    for obj in scope:
      idx = scope.index(obj)
      self._current_objects[idx] = obj

    self.add_boxes()
    self.main_sizer.Add(self.scope_sizer, flag=wx.EXPAND)

    # Add button box
    self.btn_box = ct.AddDeleteButtonBox(self, reset_button=True)
    self.main_sizer.Add(self.btn_box, flag=wx.EXPAND | wx.LEFT, border=15)

    # Bindings to buttons
    self.Bind(wx.EVT_BUTTON, self.onAdd, self.btn_box.btn_add)
    self.Bind(wx.EVT_BUTTON, self.onDelete, self.btn_box.btn_del)
    self.Bind(wx.EVT_BUTTON, self.onReset, self.btn_box.btn_rst)
    self.Bind(wx.EVT_BUTTON, self.onDeleteLast, self.btn_box.btn_del_lst)
    self.Bind(wx.EVT_BUTTON, self.onDeleteSelected, self.btn_box.btn_del_sel)
    self.Bind(wx.EVT_BUTTON, self.onCancelDelete, self.btn_box.btn_del_not)

    self.check_scope_switches()

  def add_boxes(self):
    for idx, obj in self._current_objects.items():
      label = '{} {}'.format(obj.full_path().split('.')[-1].replace('_', ' '),
                             idx + 1)
      self.add_scope_box(scope=obj, sizer=self.scope_sizer, index=idx,
                         label=label, ignore_multiple=True)

  def add_default(self, idx=None):
    if not idx:
      indices = [i for i in self._current_objects]
      if indices:
        idx = max(indices) + 1
      else:
        idx = 0
    master_object = self.phil_index.get_master_scope(self.full_path)[0]
    self._current_objects[idx] = master_object
    label = '{} {}'.format(
      master_object.full_path().split('.')[-1].replace('_', ' '),
      idx + 1
    )
    self.add_scope_box(scope=master_object, sizer=self.scope_sizer, index=idx,
                       label=label, ignore_multiple=True)
    self.dialog_layout()

  def dialog_layout(self):
    self.check_scope_switches()
    # self.SetupScrolling()
    self.Layout()
    self.GetTopLevelParent().Layout()

  def onAdd(self, e):
    self.add_default()

  def onDelete(self, e):
    for idx, ctrl in self.controls.items():
      ctrl.toggle.Show()
    self.Layout()

  def onReset(self, e):
    self.reset_boxes()

  def onDeleteLast(self, e):
    self.delete_boxes(last_only=True)

  def onDeleteSelected(self, e):
    self.delete_boxes()

  def onCancelDelete(self, e):
    for idx, ctrl in self.controls.items():
      ctrl.toggle.Hide()
    self.dialog_layout()


class PHILMultiDefPanel(PHILBaseFixedPanel, gui.IOTADefinitionCtrl,
                        MultiObjectPanelMixin):
  def __init__(self, parent, scope, *args, **kwargs):
    style = kwargs.pop('style', None)

    super(PHILMultiDefPanel, self).__init__(parent, box='', *args, **kwargs)
    gui.IOTADefinitionCtrl.__init__(self, scope)
    self.multi_definition = True
    self._current_objects = {}
    self.controls = {}

    # Add a sizer for multiple scopes, add existing scopes to it
    self.scope_sizer = PHILSizer(self)

    # Create a dictionary of PHIL objects
    for obj in scope:
      idx = scope.index(obj)
      self._current_objects[idx] = obj

    self.add_definitions()
    self.main_sizer.Add(self.scope_sizer, flag=wx.EXPAND)

    # Add button box
    self.btn_box = ct.AddDeleteButtonBox(self, reset_button=True)
    self.main_sizer.Add(self.btn_box, flag=wx.EXPAND | wx.LEFT, border=15)

    # Bindings to buttons
    self.Bind(wx.EVT_BUTTON, self.onAdd, self.btn_box.btn_add)
    self.Bind(wx.EVT_BUTTON, self.onDelete, self.btn_box.btn_del)
    self.Bind(wx.EVT_BUTTON, self.onReset, self.btn_box.btn_rst)
    self.Bind(wx.EVT_BUTTON, self.onDeleteLast, self.btn_box.btn_del_lst)
    self.Bind(wx.EVT_BUTTON, self.onDeleteSelected, self.btn_box.btn_del_sel)
    self.Bind(wx.EVT_BUTTON, self.onCancelDelete, self.btn_box.btn_del_not)

    self.Fit()

  def add_definitions(self):
    for idx, obj in self._current_objects.items():
      self.add_definition_control(self, obj=obj, ignore_multiple=True)

  def add_default(self, idx=None):
    if not idx:
      indices = [i for i in self._current_objects]
      if indices:
        idx = max(indices) + 1
      else:
        idx = 0
    master_def = self.phil_index.get_master_scope(self.full_path)[0]
    self._current_objects[idx] = master_def
    self.add_definition_control(self, obj=master_def, ignore_multiple=True)
    self.dialog_layout()

  def dialog_layout(self):
    self.GetTopLevelParent().Layout()

  def onAdd(self, e):
    self.add_default()

  def onDelete(self, e):
    for idx, ctrl in self.controls.items():
      ctrl.toggle.Show()
    self.Layout()

  def onReset(self, e):
    self.reset_boxes()

  def onDeleteLast(self, e):
    self.delete_boxes(last_only=True)

  def onDeleteSelected(self, e):
    self.delete_boxes()

  def onCancelDelete(self, e):
    for idx, ctrl in self.controls.items():
      ctrl.toggle.Hide()
    self.dialog_layout()

  def GetStringValue(self):
    """ Override in subclasses with subclass-specific function """
    raise NotImplementedError()

  def SetValue(self, value):
    """ Overwrite in subclasses with subclass-specific function """
    raise NotImplementedError()


class ValidatedTextCtrl(wx.TextCtrl):
  ''' Base class for a wx.TextCtrl that performs PHIL-specific validation and
      format checking (sub-classes will customize those functions) '''

  def __init__(self, *args, **kwargs):
    self.error_msg = None

    # Intercept a specified value to be set after initialization
    # saved_value = None
    if 'value' in kwargs:
      saved_value = kwargs['value']
      kwargs['value'] = ""
    else:
      saved_value = None

    # Initialize base class
    super(ValidatedTextCtrl, self).__init__(*args, **kwargs)
    self.parent = self.GetParent()

    # Set font style for text control
    font = wx.Font(norm_font_size, wx.FONTFAMILY_MODERN,
                   wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL)
    self.SetFont(font)
    style = self.GetWindowStyle()

    # Enforce "process ENTER" option
    if not style & wx.TE_PROCESS_ENTER:
      style |= wx.TE_PROCESS_ENTER
      self.SetWindowStyle(style)

    # Create appropriate validator (done in subclasses)
    self.SetValidator(self.CreateValidator())

    # Bindings
    self.Bind(wx.EVT_TEXT_ENTER, self.OnEnter, self)
    self.Bind(wx.EVT_KILL_FOCUS, self.OnFocusLost, self)
    self.Bind(wx.EVT_SET_FOCUS, self.OnSetFocus, self)

    # Apply value if one was passed to this control
    saved_value = self.parent.ReformatValue(value=saved_value,
                                            raise_error=False)
    self.SetStringValue(value=saved_value)

  def SetStringValue(self, value=None):
    if type(value) in (list, tuple):
      value = ' '.join([str(v) for v in value])
    self.SetValue(str(value))

  def GetStringValue(self):
    return str(self.GetValue())

  def Validate(self):
    is_good = self.GetValidator().Validate(parent=self.parent)
    if is_good:
      self.parent.RemoveError()
      self.Refresh()
    else:
      if not self.error_msg:
        self.error_msg = 'Unknown format error, please double-check your entry.'
      self.parent.SetError(self.error_msg)
      self.Refresh()

  def OnEnter(self, e=None):
    self.Validate()
    e.Skip()

  def OnFocusLost(self, e=None):
    self.Validate()
    e.Skip()

  def OnSetFocus(self, e=None):
    e.Skip()

  def CreateValidator(self):
    return gui.TextCtrlValidator().Clone()


class ValidatedStringCtrl(ValidatedTextCtrl):
  def __init__(self, *args, **kwargs):
    super(ValidatedStringCtrl, self).__init__(*args, **kwargs)

    self._min_len = 0
    self._max_len = sys.maxint

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

  def CheckFormat(self, value):
    if "$" in value:
      raise ValueError("The dollar symbol ($) may not be used here.")
    elif len(value) > self.GetMaxLength():
      raise ValueError("Value must be {} characters or less."
                       "".format(self.GetMaxLength()))
    elif len(value) < self.GetMinLength():
      raise ValueError("Value must be at least {} characters."
                       "".format(self.GetMinLength()))
    return value


class ValidatedMultiStringCtrl(ValidatedStringCtrl):
  """ A subclass of ValidatedNumberCtrl for multi-number PHIL objects """

  def __init__(self, *args, **kwargs):
    super(ValidatedMultiStringCtrl, self).__init__(*args, **kwargs)

  def CheckFormat(self, value):
    if value in (None, Auto):
      return value

    if isinstance(value, str):
      values = value.split(' ')
    elif type(value) in (list, tuple):
      values = value
    else:
      raise ValueError('Unrecognized format: expecting string or iterable')

    # Iterate through values and check format for each
    errors = []
    for v in values:
      err = None
      idx = values.index(v) + 1
      if v in (None, Auto):
        pass
      if "$" in v:
        err = "Item #{}: The dollar symbol ($) may not be used here." \
              "".format(idx)
      elif len(value) > self.GetMaxLength():
       err = "Item #{}: Value must be {} characters or less." \
             "".format(idx, self.GetMaxLength())
      elif len(value) < self.GetMinLength():
        err = "Item #{}: Value must be {} characters or more." \
              "".format(idx, self.GetMinLength())
      if err:
        errors.append(err)

    # Raise Value Error if any errors are found
    if errors:
      error_msg = 'Error(s) found!\n{}'.format('\n'.join(errors))
      raise ValueError(error_msg)
    return value


class ValidatedPathCtrl(ValidatedTextCtrl):
  def __init__(self, *args, **kwargs):
    super(ValidatedPathCtrl, self).__init__(*args, **kwargs)

    self.read = self.GetParent().read
    self.write = self.GetParent().write
    self.is_file = self.GetParent().is_file
    self.is_folder = self.GetParent().is_folder

  def CheckFormat(self, value=None):
    """ A hacky way to validate path syntax in a platform-independent way
    Args:
      value: path as string
    Returns: value if validated, raises ValueError if not
    """

    import errno

    # Just in case, check that entered path is string or unicode
    if not isinstance(value, str) and not isinstance(value, unicode):
      raise ValueError('Path must be a string!')

    # Check for None or blank space
    if value.lower() == 'none' or value.isspace() or not value:
      if self.parent.IsOptional():
        return value
      else:
        raise ValueError('Path is not optional! Please provide a valid path.')

    # Check each path component for OS errors
    # Strip Windows-specific drive specifier (e.g., `C:\`) if it exists
    _, pathname = os.path.splitdrive(value)

    # Define directory guaranteed to exist
    root_dirname = os.environ.get('HOMEDRIVE', 'C:') \
      if sys.platform == 'win32' else os.path.sep
    assert os.path.isdir(root_dirname)

    # Append a path separator to this directory if needed
    root_dirname = root_dirname.rstrip(os.path.sep) + os.path.sep

    # Validate each path component
    for pathname_part in pathname.split(os.path.sep):
      try:
        os.lstat(root_dirname + pathname_part)
      except OSError as e:
        if hasattr(e, 'winerror'):
          if e.winerror == 123:
            raise ValueError("{} is not a valid Windows path!"
                             "".format(value))
        elif e.errno in {errno.ENAMETOOLONG, errno.ERANGE}:
          raise ValueError("{} is not a valid Posix path!"
                           "".format(value))

    # Check path existence and read/write permissions
    permission_errors = []
    if os.path.isdir(value):
      dirname = value
    else:
      dirname = os.path.dirname(value)
    if self.read:
      if not self.write and not os.path.exists(value):
        raise ValueError('Path {} not found!'.format(value))
      if not os.access(dirname, os.R_OK):
        permission_errors.append('Read permission for {} denied!'
                                 ''.format(dirname))
      elif self.is_file and not os.access(value, os.R_OK):
        permission_errors.append('Read permission for {} denied!'
                                 ''.format(value))
    if self.write:
      if not os.access(dirname, os.W_OK):
        permission_errors.append('Write permission to {} denied!'
                                 ''.format(dirname))
      elif self.is_file and not os.access(value, os.W_OK):
        permission_errors.append('Write permission to {} denied!'
                                 ''.format(value))
    if permission_errors:
      raise ValueError('Permission errors for {}:\n{}'
                       ''.format(dirname, '\n'.join(permission_errors)))
    return value


class ValidatedNumberCtrl(ValidatedTextCtrl):
  def __init__(self, *args, **kwargs):

    # Check for "type" kwarg, for int or float
    self._num_type = kwargs.pop('as_type', 'float')

    # Check for 'min' kwarg to set min value for control
    vmin = kwargs.pop('min', -sys.maxint)
    if vmin is None:
      vmin = -sys.maxint
    self._value_min = int(vmin) if self._num_type == 'int' else float(vmin)

    # Check for 'max' kwarg to set max value for control
    vmax = kwargs.pop('max', sys.maxint)
    if vmax is None:
      vmax = sys.maxint
    self._value_max = int(vmax) if self._num_type == 'int' else float(vmax)

    # Check for 'allow_none' option (won't be used in subclassed multinumber
    # control)
    self._allow_none = kwargs.pop('allow_none', True)

    super(ValidatedNumberCtrl, self).__init__(*args, **kwargs)

  def SetMinValue(self, n):
    assert (n >= 0)
    self._value_min = n

  def SetMaxValue(self, n):
    assert (n >= 1)
    self._value_max = n

  def GetMinValue(self):
    return self._value_min

  def GetMaxValue(self):
    return self._value_max

  def determine_type(self, value, is_none=True, is_auto=True):
    suggested_type = 'a number'
    try:
      if 'int' in self._num_type:
        value = int(value)
        suggested_type = 'an integer'
      else:
        value = float(value)
        suggested_type = 'a float'
    except ValueError:
      value = str(value)
      if value.lower() == 'none':
        value = None
      elif value.lower() == 'auto':
        value = Auto
    if isinstance(value, str):
      if is_auto and is_none:
        suggested_type += ', None, or Auto'
      elif is_none:
        suggested_type += 'or None'
      elif is_auto:
        suggested_type += ', None, or Auto'
    return suggested_type, value

  def CheckFormat(self, value):
    """ Checks that the format of the value is numerical; if string is found,
        only 'none' and 'auto' can be accepted as valid
    :param value: entered value
    :return: checked value or error
    """
    is_none = self._allow_none
    is_auto = self.parent.UseAuto()
    suggested_type, value = self.determine_type(value,
                                                is_none=is_none,
                                                is_auto=is_auto)
    if (value is None and is_none) or (value is Auto and is_auto):
      pass
    elif isinstance(value, str):
      if (value.lower() == 'none' and is_none) or \
              (value.lower() == 'auto' and is_auto):
        pass
      else:
        raise ValueError("String entries are not allowed! Enter {}."
                         "".format(suggested_type))
    else:
      if value > self.GetMaxValue():
        raise ValueError("Value ({}) must be less than the maximum of {}."
                         "".format(value, self.GetMaxValue()))
      elif value < self.GetMinValue():
        raise ValueError("Value ({}) must be more than the minimum of {}."
                         "".format(value, self.GetMinValue()))
    return value


class ValidatedMultiNumberCtrl(ValidatedNumberCtrl):
  """ A subclass of ValidatedNumberCtrl for multi-number PHIL objects """

  def __init__(self, *args, **kwargs):
    self._size_min = kwargs.pop('size_min', 0)
    self._size_max = kwargs.pop('size_max', sys.maxint)
    self._allow_none_elements = kwargs.pop('allow_none_elements', True)
    self._allow_auto_elements = kwargs.pop('allow_auto_elements', True)
    super(ValidatedMultiNumberCtrl, self).__init__(*args, **kwargs)

  def CheckFormat(self, value):
    """ Checks that the value list is the right size; that each of the items
        is numerical; and if string values are found, only 'none' and 'auto' can
        be accepted as valid
    :param value: string or list containing entered value(s)
    :return: checked values or error
    """
    is_none = self._allow_none_elements
    is_auto = self._allow_auto_elements
    if (str(value).lower() == 'none' and is_none) or \
            (str(value).lower() == 'auto' and is_auto):
      return value
    if isinstance(value, str):
      values = value.strip().split(' ')
    else:
      assert type(value) in (list, tuple)
      values = value
    if len(values) < self._size_min:
      raise ValueError("Need a minimum of {} values in this field!"
                       "".format(self._size_min))
    if len(values) > self._size_max:
      raise ValueError("Cannot have more than {} values in this field!"
                       "".format(self._size_max))

    # Iterate through values and check format for each; error message will be
    # a summary of all found errors
    errors = []
    new_values = []
    for item in values:
      idx = values.index(item)
      err = None
      suggested_type, item = self.determine_type(item,
                                                 is_none=is_none,
                                                 is_auto=is_auto)
      new_values.append(item)
      if (item is None and is_none) or (item is Auto and is_auto):
        pass
      if isinstance(item, str):
        if (item.lower() == 'none' and is_none) or\
                (item.lower() == 'auto' and is_auto):
          pass
        else:
          err = "String entries are not allowed! Enter {}."\
                "".format(suggested_type)
      else:
        if item > self.GetMaxValue():
          err = "Value ({}) must be less than the maximum of {}."
          "".format(item, self.GetMaxValue())
        elif item < self.GetMinValue():
          err = "Value ({}) must be more than the minimum of {}."
          "".format(item, self.GetMinValue())
      if err:
        msg = '   Item #{}: {}'.format(idx, err)
        errors.append(msg)

    # Raise Value Error if any errors are found
    if errors:
      error_msg = 'Error(s) found!\n{}'.format('\n'.join(errors))
      raise ValueError(error_msg)

    return new_values


class ValidatedUnitCellCtrl(ValidatedTextCtrl):
  def __init__(self, *args, **kwargs):
    super(ValidatedUnitCellCtrl, self).__init__(*args, **kwargs)

  def CheckFormat(self, value):
    """ Check that the entry is a valid unit cell notation
    Args:
      value: Unit cell parameters
    Returns: value if validated, raises ValueError if not
    """
    value = makenone(value)

    # Break up string by possible delimiters (set up to handle even a mixture
    # of delimiters, because you never know)
    if value:
      uc = value
      for dlm in [',', ';', '|', '-', '/']:
        if uc.count(dlm) > 0:
          uc = uc.replace(dlm, ' ')
      uc_params = uc.rsplit()

      error_msg = 'Invalid unit cell entry {}'.format(value)

      # Check if there are six parameters in unit cell (that's as far as
      # validation will go; the correctness of the unit cell is up to the user)
      if len(uc_params) != 6:
        raise ValueError('{}: should be six parameters!'.format(error_msg))

      # Check that every parameter can be converted to float (i.e. that the
      # parameters are actually all numbers)
      try:
        uc_params = [float(p) for p in uc_params]
      except ValueError as e:
        if 'invalid literal' in e.message.lower():
          raise ValueError('{}: unit cell should only contain '
                           'numbers'.format(error_msg))
        raise ValueError('{}: {}'.format(error_msg, e.message))

      uc = ' '.join([str(p) for p in uc_params])

    else:
      uc = None

    return str(uc)


class ValidatedSpaceGroupCtrl(ValidatedTextCtrl):
  def __init__(self, *args, **kwargs):
    super(ValidatedSpaceGroupCtrl, self).__init__(*args, **kwargs)

  def CheckFormat(self, value):
    """ Check that the entry is a valid space group notation
    Args:
      value: Space group symbol and/or number
    Returns: value if validated, raises ValueError if not
    """

    # Attempt to create a space group object; if symbol or number are
    # invalid, this will fail
    from cctbx import crystal
    value = makenone(value)       # Convert to Nonetype if None or Null string

    if value:
      try:
        sym = crystal.symmetry(space_group_symbol=str(value))
      except RuntimeError as e:
        raise ValueError('Invalid space group entry:\n{}'.format(e))

      # Return space group symbol even if number is entered
      sg = sym.space_group_info()
    else:
      sg = None

    return str(sg)


# ------------------------------- PHIL Buttons ------------------------------- #

class PHILDialogButton(ct.IOTAButton):
  """ Button that launches a wx.Dialog auto-populated with PHIL settings.
      Will also take out a specified scope from the PHILIndex that resides in
      the master window. A selection of scopes to show (while hiding all the
      rest) can also be provided. """

  def __init__(self, parent, scope, phil_index, name=None,
               expert_level=None, title=None, *args, **kwargs):
    """ Constructor
    :param parent: parent control (here, likely a panel)
    :param scope_name: Name of scope that will appear as a dialog
    :param include: A list of paths for scopes or definitions that can be
                      modified in the dialog; the omitted paths will not be
                      displayed or modified. Can be full paths or names.
    :param label: Custom label for the button; if None, a label will be made
                  from the scope name
    :param args: Any other arguments
    :param kwargs: Any other keyword arguments
    """

    self.parent = parent
    self.is_dlg_button = True
    self.is_scope = False
    self.is_definition = False
    self.scope = scope
    self.phil_index = phil_index
    self.name = name
    self.expert_level = expert_level
    self.title = title

    ct.IOTAButton.__init__(self, parent, handler_function=self.open_phil_dialog,
                           *args,  **kwargs)

    # Set attributes from scope(s) and set the button label
    self._set_attributes()
    label = self.title + '...'
    self.SetLabel(label)

  def _set_attributes(self):
    if isinstance(self.scope, list):  # means it's a multiple!
      scp = self.scope[0]
    else:
      scp = self.scope

      # Determine overall expert level
    self.full_path = scp.full_path()
    self.name = scp.name
    if not self.name:
      if self.full_path:
        self.name = self.full_path.split('.')[-1]
      else:
        self.name = 'phil_options'

    # Create button and dialog title from variable name or alias
    if not self.title:
      if scp.short_caption:
        self.title = scp.short_caption
      elif scp.alias:
        self.title = scp.alias.capitalize()
      elif self.name and not self.name.isspace():
        self.title = self.name.replace('_', ' ').capitalize()

  def get_phil_strings(self):
    if isinstance(self.scope, list):
      phil_strings = []
      for scp in self.scope:
        scp_strings = scp.as_str().split('\n')
        phil_strings.extend(scp_strings)
    else:
      phil_strings = self.scope.as_str().split('\n')
    return phil_strings

  def open_phil_dialog(self, e):
    """ Generate a PHIL Dialog from scope; on_OK, the dialog will generate a
        PHIL object that can be used to update the button scope """
    with PHILDialog(parent=self.parent,
                    name=self.name,
                    scope=self.scope,
                    phil_index=self.phil_index,
                    title=self.title) as phil_dlg:
      if phil_dlg.ShowModal() == wx.ID_OK:
        self.scope = phil_dlg.scope

# ------------------------------ PHIL Controls ------------------------------- #


class PHILFileListCtrl(ct.FileListCtrl, gui.IOTADefinitionCtrl):
  def __init__(self, parent, phil_object, value=None, *args, **kwargs):

    extras = kwargs.pop('extras', None)
    if extras:
      file_types = extras.get('file_types', None)
      folder_types = extras.get('folder_types', None)
      data_types = extras.get('data_types', None)
      input_filter = extras.get('input_filter', None)
    else:
      file_types = folder_types = data_types = None

    ct.FileListCtrl.__init__(self, parent=parent,
                             size=(600, 300),
                             file_types=file_types,
                             folder_types=folder_types,
                             data_types=data_types,
                             input_filter=input_filter
                             )
    gui.IOTADefinitionCtrl.__init__(self, phil_object=phil_object)
    self.multi_definition = False
    self.SetValue(value)

  def update_from_phil(self):
    new_values = self.parent.phil_index.get_value(path=self.full_path)
    self.ResetValue(value=new_values)

  def ResetValue(self, value=None):
    self.delete_all()
    self.SetValue(value)

  def SetValue(self, value=None):
    if isinstance(value, list):
      values = [v for v in value if v]
      if values:
        for v in values:
          self.add_item(path=v)
    else:
      if value:
        self.add_item(path=value)

  def set_background(self):
    pass

  def GetPHIL(self, full_path=False, indent_length=2):
    """ Overridden because FileListCtrl is tricky """
    idxs = self.ctr.GetItemCount()
    inputs = [self.ctr.GetItemData(i).path for i in range(idxs)]

    phil_strings = []
    indent = len(self.full_path.split('.')) * indent_length
    for inp in inputs:
      value_string = '{} = {}'.format(self.name, inp)
      phil_string = '{:{ind}}{}'.format(' ', value_string, ind=indent)
      phil_strings.append(phil_string)

    return '\n'.join(phil_strings)


class PHILPathCtrl(PHILDefPanel):
  """ Control for the PHIL path type """

  def __init__(self, parent, phil_object, label='', style=None, value=None,
               label_size=wx.DefaultSize, *args, **kwargs):
    cols = kwargs.pop('cols', 4)
    vgap = kwargs.pop('vgap', 0)
    PHILDefPanel.__init__(self, parent=parent, phil_object=phil_object,
                          label_size=label_size, cols=cols, vgap=vgap,
                          *args, **kwargs)

    # Extract relevant styles
    self.is_file = 'file' in style.path if style.path else False
    self.is_folder = 'folder' in style.path if style.path else False
    self.read = 'read' in style.permissions if style.permissions else True
    self.write = 'write' in style.permissions if style.permissions else True

    # If neither read-only nor write-only is specified, read/write is enabled
    if not self.read and not self.write:
      self.read = self.write = True

    # If file *and* folder are set to True (unlikely), self.is_folder = True
    if not self.is_file and not self.is_folder:
      self.is_folder = True
    elif self.is_file and self.is_folder:
      print ('IOTA PHIL DEBUG: {} specified as both file and folder! Setting '
             'to folder...'.format(phil_object.full_path()))
      self.is_file = False

    # Set defaultfile and wildcard parameters for file dialog
    if self.is_file:
      self.defaultfile = style.defaultfile if style.defaultfile else '*'
      self.wildcard = style.wildcard if style.wildcard else '*'
    else:
      self.defaultfile = '*'
      self.wildcard = '*'

    # Create path control
    self.ctr = ValidatedPathCtrl(self, value=value)
    self.SetValue(value=phil_object)
    self.ctrl_sizer.add_widget_and_label(widget=self.ctr, label=label,
                                         expand=True, label_size=label_size)
    self.ctrl_sizer.add_growable(cols=[3])

    # Create browse and view buttons
    self.btn_browse = wx.Button(self, label='...', style=wx.BU_EXACTFIT)
    viewmag_bmp = bitmaps.fetch_icon_bitmap('actions', 'viewmag', size=16)
    self.btn_mag = wx.BitmapButton(self, bitmap=viewmag_bmp)
    self.ctrl_sizer.add_widget(self.btn_browse)
    self.ctrl_sizer.add_widget(self.btn_mag)

    # Bindings
    self.Bind(wx.EVT_BUTTON, self.OnBrowse, self.btn_browse)

  def SetValue(self, value):
    try:
      self.SetStringValue(phil_object=value)
    except AttributeError:
      self.ctr.SetValue(str(value))

  def SetStringValue(self, phil_object):
    value = self.value_from_words(phil_object=phil_object)[0]
    self.ctr.SetValue(str(value))

  def GetStringValue(self):
    return self.ctr.GetValue()

  def OnBrowse(self, e):
    if self.is_file:
     self._open_file_dialog()
    elif self.is_folder:
      self._open_folder_dialog()
    else:
      command_list = [('Browse files...',
                       lambda evt: self._open_file_dialog()),
                      ('Browse folders...',
                       lambda evt: self._open_folder_dialog())]
      browse_menu = ct.Menu(self)
      browse_menu.add_commands(command_list)
      self.PopupMenu(browse_menu)
      browse_menu.Destroy()

    self.ctr.Validate()
    e.Skip()

  def onMagButton(self, e):
    dirview = DirView(self, title='Current Folder')
    if dirview.ShowModal() == wx.ID_OK:
      dirview.Destroy()
    e.Skip()

  def _open_folder_dialog(self):
    dlg = wx.DirDialog(self, "Choose folder:", style=wx.DD_DEFAULT_STYLE)
    if dlg.ShowModal() == wx.ID_OK:
      self.ctr.SetValue(dlg.GetPath())
    dlg.Destroy()

  def _open_file_dialog(self):
    dlg = wx.FileDialog(
      self, message="Choose file",
      defaultDir=os.curdir,
      defaultFile=self.defaultfile,
      wildcard=self.wildcard,
      style=wx.FD_OPEN | wx.FD_CHANGE_DIR
    )
    if dlg.ShowModal() == wx.ID_OK:
      filepath = dlg.GetPaths()[0]
      self.ctr.SetValue(filepath)


class PHILStringCtrl(PHILDefPanel):
  """ Control for the PHIL string type """

  def __init__(self, parent, phil_object, control=None, label='', value=None,
               label_size=wx.DefaultSize, *args, **kwargs):
    vgap = kwargs.pop('vgap', 0)
    PHILDefPanel.__init__(self,
                          parent=parent,
                          phil_object=phil_object,
                          label_size=label_size,
                          vgap=vgap,
                          *args, **kwargs)

    if not control:
      control = ValidatedStringCtrl

    self._create_control(control=control, value=value)
    self._place_control(label, label_size)

  def _create_control(self, control, value=None):
    self.ctr = control(self, value=value)
    self.ctr.Validate()

  def _place_control(self, label, label_size):
    self.ctrl_sizer.add_widget_and_label(widget=self.ctr, label=label,
                                         expand=True, label_size=label_size)
    self.ctrl_sizer.add_growable(cols=[3])

  def GetStringValue(self):
    return self.ctr.GetValue()

  def SetValue(self, value=None):
    value = self.ReformatValue(value, raise_error=False)
    self.ctr.SetStringValue(value=value)
    self.ctr.Validate()


class PHILMultiStringCtrl(PHILDefPanel):
  """ Control for the PHIL ints and floats (multiple numbers) types """

  def __init__(self, parent, phil_object, label='', value=None,
               label_size=wx.DefaultSize, *args, **kwargs):
    vgap = kwargs.pop('vgap', 0)
    PHILDefPanel.__init__(self, parent=parent, phil_object=phil_object,
                          vgap=vgap, *args, **kwargs)

    # make control
    self.ctr = ValidatedMultiStringCtrl(self, value=value)
    self.ctr.Validate()

    # place control in sizer
    self.ctrl_sizer.add_widget_and_label(self.ctr, label=label,
                                         label_size=label_size, expand=True)
    self.ctrl_sizer.add_growable(cols=[3])

  def SetValue(self, value):
    if isinstance(value, str):
      value = value.split()
    elif type(value) not in (list, tuple):
      raise ValueError('IOTA GUI Error: Multi-string Control: Value should be a '
                       'string or an iterable!')
    value = self.ReformatValue(value, raise_error=False)
    self.ctr.SetStringValue(value=value)
    self.ctr.Validate()

  def GetStringValue(self):
    return self.ctr.GetValue()

  def is_default(self):
    if type(self.default_value) in (list, tuple):
      default_value = ' '.join([str(v).lower() for v in self.default_value])
    else:
      default_value = str(self.default_value).lower()
    control_value = self.GetStringValue().lower()
    return default_value == control_value


class PHILSpaceGroupCtrl(PHILStringCtrl):
  """ Control for the PHIL Space Group, subclassed from PHILStringCtrl """

  def __init__(self, parent, phil_object, label='', value=None,
               label_size=wx.DefaultSize, *args, **kwargs):
    PHILStringCtrl.__init__(self, parent,
                            label=label,
                            control=ValidatedSpaceGroupCtrl,
                            phil_object=phil_object,
                            label_size=label_size,
                            value=value,
                            *args, **kwargs)


class PHILUnitCellCtrl(PHILStringCtrl):
  """ Control for the PHIL Unit Cell, subclassed from PHILStringCtrl """

  def __init__(self, parent, phil_object, label='', value=None,
               label_size=wx.DefaultSize, *args, **kwargs):

    PHILStringCtrl.__init__(self, parent,
                            phil_object=phil_object,
                            control=ValidatedUnitCellCtrl,
                            label=label,
                            label_size=label_size,
                            value=value,
                            *args, **kwargs)


class PHILBaseChoiceCtrl(PHILDefPanel):
  """ Choice control for PHIL choice item, with label """

  def __init__(self, parent,
               phil_object,
               label='',
               label_size=wx.DefaultSize,
               ctrl_size=wx.DefaultSize,
               *args, **kwargs):
    """ Constructor
    :param parent: parent object
    :param label: choice control label
    :param label_size: size of choice control label
    :param label_style: normal, bold, italic, or italic_bold
    :param ctrl_size: size of choice control
    """
    # Initialize the base class
    PHILDefPanel.__init__(self, parent=parent, phil_object=phil_object,
                          label_size=label_size, *args, **kwargs)
    self._options = None

    # Set choice control
    self.ctr = wx.Choice(self, size=ctrl_size)
    self.ctrl_sizer.add_widget_and_label(widget=self.ctr, label=label,
                                         label_size=label_size)
    self.Bind(wx.EVT_CHOICE, self.onChoice, self.ctr)

  def is_default(self):
    default_value = self.default_value
    control_value = self.GetPHILValue()

    # Sometimes None is replaced with '---' in PHIL choice controls
    if control_value == '---':
      control_value = None

    return str(control_value) == str(default_value)

  def onChoice(self, e):
    self.set_background()

  def SetChoices(self, choices, captions=None, value=None, allow_none=True):
    ''' Insert choices into the control
    :param choices: list of choices (must be list or tuple)
    :param value: definition value from PHIL object
    :param captions: list of captions (optional)
    :param allow_none: allow no selection
    '''

    # Determine selection
    selection = None
    if isinstance(value, list):
      value = value[0]
    if value:
      is_selected = [(choice.replace('*', '') == value) for choice in choices]
    else:
      is_selected = [("*" in choice) for choice in choices]
    if True in is_selected:
      selection = is_selected.index(True)

    # Strip asterisk(s) from list of choices
    choices = [choice.replace("*", "") for choice in choices]

    # Apply or create captions and set selection
    if captions is None:
      captions = list(choices)
    if len(captions) != len(choices):
      raise RuntimeError("Wrong number of caption items for {}\n"
                         "Choices: {}\n"
                         "Captions: {}"
                         "".format(self.name,
                                   '\n'.join(choices),
                                   '\n'.join(captions)))

    # Add a dashed line if parameter is optional (or None is an option)
    if allow_none:
      if choices[0] is None or choices[0].lower() == 'none':
        captions[0] = '---'
        choices[0] = None
      elif self.IsOptional():
        captions.insert(0, "---")
        choices.insert(0, None)
        # Increment selection to account for item insertion
        if selection is not None:
          selection += 1

    # Sometimes selection may be None; if so, set it to zero
    if selection is None:
      selection = 0

    # Set options, captions, and selection
    self._options = choices
    self.ctr.SetItems(captions)
    self.ctr.SetSelection(selection)

  def SetValue(self, value):
    ''' Set choice selection to specific value '''
    if value not in self._options:
      raise Sorry('Value {} not found! Available choices are:\n{}'
                  ''.format(value, '\n'.join(self._options)))

    selection = self._options.index(value)
    self.ctr.SetSelection(selection)

  def GetValue(self):
    raise NotImplementedError("Please use GetPhilValue()")

  def GetPHILValue(self):
    """Returns a single string."""
    return self._options[self.ctr.GetSelection()]

  def GetStringValue(self):
    """Returns the long format (all choices, '*' denotes selected)."""
    selection = self.ctr.GetSelection()
    choices_out = []
    for i, choice in enumerate(self._options):
      if choice is None:
        continue
      elif i == selection:
        choices_out.append("*" + choice)
      else:
        choices_out.append(choice)
    return " ".join(choices_out)


class PHILChoiceCtrl(PHILBaseChoiceCtrl):
  def __init__(self, parent, phil_object, label='', captions=None, value=None,
               *args, **kwargs):
    super(PHILChoiceCtrl, self).__init__(parent=parent,
                                          phil_object=phil_object,
                                          label=label,
                                          *args, **kwargs)
    choices = [str(i) for i in phil_object.words]
    self.SetChoices(choices=choices, captions=captions, value=value)

  def GetValue(self):
    raise NotImplementedError("Please use GetPhilValue()")


class PHILTriBoolCtrl(PHILChoiceCtrl):
  """ Three-way boolean control: returns True, False, or None. Currently used as
      the option if a boolean PHIL definition has a default value of None.
      PHIL definitions with default values of Auto or True/False are
      automatically made as wx.CheckBox controls. That can be overridden by
      specifying 'tribool' in a definition's style card.
  """
  def __init__(self, parent, phil_object, label='', value=None,
               *args, **kwargs):
    super(PHILTriBoolCtrl, self).__init__(parent=parent,
                                          phil_object=phil_object,
                                          label=label,
                                          value=value,
                                          *args, **kwargs)
    self._options = None
    self.SetOptional(True)
    self.SetChoices(choices=['None', 'Yes', 'No'], value=value,
                    allow_none=False)

  def SetValue(self, value):
    if value is True:
      self.ctr.SetSelection(1)
    elif value is False:
      self.ctr.SetSelection(2)
    else:
      assert value in [None, Auto]
      self.ctr.SetSelection(0)

  def GetValue(self):
    return self.GetPhilValue()

  def GetPhilValue(self):
    vals = [None, True, False]
    return vals[self.ctr.GetSelection()]

  def GetStringValue(self):
    return str(self.GetPhilValue())


class PHILNumberCtrl(PHILDefPanel):
  """ Control for the PHIL int and float types """

  def __init__(self, parent, phil_object, label='', value=None,
               label_size=wx.DefaultSize, *args, **kwargs):
    vgap = kwargs.pop('vgap', 0)
    PHILDefPanel.__init__(self, parent=parent, phil_object=phil_object,
                          vgap=vgap, *args, **kwargs)

    self.ctr = ValidatedNumberCtrl(self,
                                   as_type=phil_object.type.phil_type,
                                   min=phil_object.type.value_min,
                                   max=phil_object.type.value_max,
                                   value=value)
    self.ctr.Validate()   # Validate to make sure the input value is legit

    if hasattr(phil_object.type, 'allow_none'):
      self.SetOptional(optional=phil_object.type.allow_none)

    self.ctrl_sizer.add_widget_and_label(self.ctr, label=label,
                                         label_size=label_size, expand=True)
    self.ctrl_sizer.add_growable(cols=[2])

  def GetStringValue(self):
    """ Extract value as a string """
    return self.ctr.GetValue()

  def SetValue(self, value):
    value = self.ReformatValue(value, raise_error=False)
    self.ctr.SetStringValue(value=value)
    self.ctr.Validate()

  def is_default(self):
    default_value = str(self.default_value)
    control_value = str(self.GetStringValue())

    # convert to float (unless None, Auto) for accurate comparison
    if default_value.isdigit():
      default_value = float(default_value)
    if control_value.isdigit():
      control_value = float(control_value)
    return control_value == default_value


class PHILMultiNumberCtrl(PHILDefPanel):
  """ Control for the PHIL ints and floats (multiple numbers) types """

  def __init__(self, parent, phil_object, label='', value=None,
               label_size=wx.DefaultSize, *args, **kwargs):
    vgap = kwargs.pop('vgap', 0)
    PHILDefPanel.__init__(self, parent=parent, phil_object=phil_object,
                          vgap=vgap, *args, **kwargs)

    # make control
    self.ctr = ValidatedMultiNumberCtrl(
      self,
      as_type=phil_object.type.phil_type,
      min=phil_object.type.value_min,
      max=phil_object.type.value_max,
      size_min=phil_object.type.size_min,
      size_max=phil_object.type.size_max,
      value=value
    )

    # Validate to make sure the input value is legit
    self.ctr.Validate()

    # place control in sizer
    self.ctrl_sizer.add_widget_and_label(self.ctr, label=label,
                                         label_size=label_size, expand=True)
    self.ctrl_sizer.add_growable(cols=[2])

  def GetStringValue(self):
    """ Extract value as a string """
    return self.ctr.GetValue()

  def SetValue(self, value):
    if isinstance(value, str):
      value = value.split()
    value = self.ReformatValue(value, raise_error=False)
    self.ctr.SetStringValue(value=value)
    self.ctr.Validate()

  def is_default(self):
    # Convert strings into lists for comparison
    if isinstance(self.default_value, str):
      default_values = self.default_value.split()
    else:
      default_values = self.default_value
    control_values = self.GetStringValue().split()

    # Iterate through values and compare (assume None and Auto values will be
    # NoneType and AutoType objects, rather than strings)
    if default_values is not None:
      for dv in default_values:
        idx = default_values.index(dv)
        cv = str(control_values[idx])
        dv = str(dv)
        if dv.isdigit():
          dv = float(dv)
        if cv.isdigit():
          cv = float(cv)
        if dv != cv:
          return False
      return True
    else:
      if len(control_values) == 1:
        control_values = control_values[0]
      elif len(control_values) == 0:
        control_values = None
      return str(default_values) == str(control_values)


class PHILCheckBoxCtrl(PHILDefPanel):
  """ Checkbox control for PHIL bool item (with tribool option as default;
      PHIL style has to be set to noauto in order to make this a regular
      bool) """

  def __init__(self, parent, phil_object, label='', style=None,
               label_size=wx.DefaultSize, value=False, *args, **kwargs):
    """ Constructor """
    vgap = kwargs.pop('vgap', 0)
    PHILDefPanel.__init__(self, parent=parent, phil_object=phil_object,
                          label_size=label_size, vgap=vgap, *args, **kwargs)
    self.style = style
    self.scope_switch = self.style.scope_switch
    if self.default_value is Auto:
      self.ctr = wx.CheckBox(self, label=label,
                             style=wx.CHK_ALLOW_3RD_STATE_FOR_USER |
                                   wx.CHK_3STATE)
    else:
      self.ctr = wx.CheckBox(self, label=label)
    self.SetValue(value)

    if label == '':
      border = 0
    else:
      border = 5
    self.ctrl_sizer.add_widget(widget=self.ctr, border=border)

    self.Bind(wx.EVT_CHECKBOX, self.onChangeValue, self.ctr)

  def GetStringValue(self):
    """ Extract value as a string """
    return str(self.GetValue())

  def onChangeValue(self, e):
    self.set_background()
    if self.scope_switch:
      self.parent.check_scope_switches()

  def SetValue(self, value=None):
    """ Set checkbox state, None is interpreted as Auto if either is allowed """
    if value in (None, Auto):
      assert (self.ctr.Is3State())
      self.ctr.Set3StateValue(wx.CHK_UNDETERMINED)
    else:
      value = bool(value)
      if self.ctr.Is3State():
        if value:
          self.ctr.Set3StateValue(wx.CHK_CHECKED)
        else:
          self.ctr.Set3StateValue(wx.CHK_UNCHECKED)
      else:
        self.ctr.SetValue(value)

  def GetValue(self):
    """ Set checkbox state """
    if self.ctr.Is3State():
      value = self.ctr.Get3StateValue()
      if value == wx.CHK_UNDETERMINED:
        return Auto
      else:
        return value == wx.CHK_CHECKED
    else:
      return self.ctr.GetValue()


class WidgetFactory(object):
  ''' Class that will automatically make widgets for automated dialog making '''
  widget_types = {
    'path'        : PHILPathCtrl        ,
    'str'         : PHILStringCtrl      ,
    'strings'     : PHILMultiStringCtrl ,
    'choice'      : PHILChoiceCtrl      ,
    'number'      : PHILNumberCtrl      ,
    'numbers'     : PHILMultiNumberCtrl ,
    'bool'        : PHILCheckBoxCtrl    ,
    'tribool'     : PHILTriBoolCtrl     ,
    'space_group' : PHILSpaceGroupCtrl  ,
    'unit_cell'   : PHILUnitCellCtrl
  }

  def __init__(self):
    pass

  @staticmethod
  def make_widget(parent,
                  phil_object,
                  label_size=wx.DefaultSize,
                  widget_types=widget_types,
                  *args, **kwargs):

    style = kwargs.pop('style', None)
    value = kwargs.pop('value', None)
    label = kwargs.pop('label', None)
    wtype = phil_object.type.phil_type

    if label is None:
      alias = phil_object.alias_path()
      if alias:
        label = alias
      else:
        label = phil_object.full_path().split('.')[-1]
        label = label.replace('_', ' ').capitalize()

      if wtype == 'bool':
        label_size = wx.DefaultSize
      else:
        label += ": "

    if wtype in ('int', 'float'):
      wtype = 'number'
    elif wtype in ('ints', 'floats'):
      wtype = 'numbers'
    elif style.tribool or (wtype == 'bool' and value is None):
      wtype = 'tribool'
    if wtype == 'path' and (style and style.input_list):
      widget = PHILFileListCtrl(parent=parent,
                                phil_object=phil_object,
                                value=value,
                                *args, **kwargs)
    else:
      if wtype in widget_types:
        widget_ctrl = widget_types[wtype]
      else:
        widget_ctrl = PHILStringCtrl

      widget = widget_ctrl(parent=parent,
                           phil_object=phil_object,
                           label=label,
                           label_size=label_size,
                           value=value,
                           style=style,
                           *args, **kwargs)

    return widget

# ----------------------------- PHIL Controls -------------------------------- #


class PHILDialogPanel(PHILBaseDialogPanel):
  """ Panel automatically created from PHIL settings """

  _str_extras = {}
  _path_extras = {}
  _number_extras = {}
  _unit_cell_extras = {}
  _space_group_extras = {}
  _choice_extras = {}
  _checkbox_extras = {}

  def __init__(self, parent, scope, *args, **kwargs):
    str_extras = kwargs.pop('str_extras', None)
    if str_extras:
      self._str_extras.update(str_extras)

    path_extras = kwargs.pop('path_extras', None)
    if path_extras:
      self._path_extras.update(path_extras)

    number_extras = kwargs.pop('number_extras', None)
    if number_extras:
      self._number_extras.update(number_extras)

    unit_cell_extras = kwargs.pop('unit_cell_extras', None)
    if unit_cell_extras:
      self._unit_cell_extras.update(unit_cell_extras)

    space_group_extras = kwargs.pop('space_group_extras', None)
    if space_group_extras:
      self._space_group_extras.update(space_group_extras)

    choice_extras = kwargs.pop('choice_extras', None)
    if choice_extras:
      self._choice_extras.update(choice_extras)

    checkbox_extras = kwargs.pop('checkbox_extras', None)
    if checkbox_extras:
      self._checkbox_extras.update(checkbox_extras)

    super(PHILDialogPanel, self).__init__(parent, scope, *args, **kwargs)
    self.scope = scope

    if not self.phil_index:
      if hasattr(self.window, 'phil_index'):
        self.phil_index = self.window.phil_index
      else:
        raise Sorry('IOTA PHIL ERROR: PHILIndex not found!')

    # Recurse through scope and create controls and widgets
    self.construct_panel()

    # Redraw the window to show/hide panels based on expert level
    self.redraw_by_expert_level()


class PHILDialog(PHILBaseDialog):
  """ Dialog auto-populated with PHIL settings """

  def __init__(self, parent, scope, phil_index, name=None, title=None,
               *args, **kwargs):
    """ Constructor
    :param parent: parent GUI element
    :param name: name of the PHIL scope rendered by this dialog
    :param scope: PHIL scope from which to build the dialog
    """
    self.phil_index = phil_index
    self.parent = parent
    self.name = name
    if not self.name:
      self.name = parent.name
    if not title:
      title = self.name.replace('_', ' ').capitalize()
    super(PHILDialog, self).__init__(parent, title=title, *args, **kwargs)

    if not isinstance(scope, list) and scope.multiple:
      self.scope = [scope]
    else:
      self.scope = scope

    if isinstance(self.scope, list):
      self._scope_paths = list(set([s.full_path() for s in self.scope]))
      self.phil_panel = PHILMultiScopePanel(self,
                                            scope=self.scope,
                                            phil_index=self.phil_index)
    else:
      self._scope_paths = [self.scope.full_path()]
      self.phil_panel = PHILDialogPanel(self,
                                        scope=self.scope,
                                        phil_index=self.phil_index)
    self.phil_sizer.add_panel(self.phil_panel)

    # Dialog control
    self.dlg_ctr = ct.DialogButtonsCtrl(self, preset='PHIL_DIALOG')
    self.envelope.Add(self.dlg_ctr,
                      flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                      border=10)

    # Set up size and scrolling (adjust size if auto-fit dialog is too big
    # for screen)
    self.Fit()
    self.phil_panel.SetupScrolling()
    self.size_and_place()
    self.Layout()

    # Bindings
    self.Bind(wx.EVT_BUTTON, self.OnOkay, id=wx.ID_OK)
    self.Bind(wx.EVT_BUTTON, self.OnCancel, id=wx.ID_CANCEL)
    self.Bind(wx.EVT_CHOICE, self.OnExpertLevel, self.dlg_ctr.choice)

  def OnOkay (self, event):
    """ Check for saved errors and pop a warning if any are found (user
        cannot continue if any errors are present) """

    all_errors = self.phil_panel.collect_errors()
    if all_errors:
      # Check for errors and pop up a message if any are present
      wx.MessageBox(caption='Errors in Settings!',
                    message='Correct all errors to accept the changes.',
                    style=wx.OK|wx.ICON_EXCLAMATION)

    else:
      phil_string = self.phil_panel.GetPHIL(expand=True)
      if phil_string:
        self.phil_index.update_phil(phil_string=phil_string)
        self.scope = self.phil_index.get_scopes(include=self._scope_paths)
      self.EndModal(wx.ID_OK)


# -- end
