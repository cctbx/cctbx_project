from __future__ import division, print_function, absolute_import

'''
Author      : Lyubimov, A.Y.
Created     : 04/02/2019
Last Changed: 08/01/2019
Description : IOTA GUI initialization module
'''

import wx
from wxtbx import wx4_compatibility as wx4c

from libtbx import Auto, str_utils
from libtbx.utils import Sorry, to_unicode
from libtbx.phil import gui_objects, interface

Validator = wx4c.get_wx_mod(wx, wx.Validator)

additional_style_args = [
  'input_list',
  'scope_switch',
  'grid',
  'has_scope_switch'
]
additional_style_kwds = [
  'path',
  'permissions',
  'defaultfile',
  'wildcard',
  'align_label'
]

style = gui_objects.style
style.style_args.extend(additional_style_args)
style.style_kwds.extend(additional_style_kwds)


def make_phil_index(master_phil, working_phil=None, fetch_new=True):

  # collect scopes of the same name in phil
  merged_scopes = collect_scopes(phil=master_phil)
  updated_master_phil = master_phil.customized_copy(objects=merged_scopes)

  # self-fetch to resolve variables
  updated_master_phil = updated_master_phil.fetch(source=updated_master_phil)

  return PHILIndex(master_phil=updated_master_phil,
                   working_phil=working_phil,
                   fetch_new=fetch_new)


def collect_scopes(phil):
  """ Collect and merge scopes with the same path. Assuming that the only
      groups of same-name objects in a PHIL would be scopes; several
      definitions of the same name that aren't multiple would throw an
      error on fetch. """
  scp_paths = {}
  phil_objects = []

  # Go through all objects and collect ones with same name in a dictionary
  for obj in phil.objects:
    if obj.is_scope:
      full_path = obj.full_path()
      if full_path not in scp_paths:
        scp_paths[full_path] = [obj]
      else:
        scp_paths[full_path].append(obj)
    elif obj.is_definition:
      phil_objects.append(obj)

  # Merge scopes that have the same path
  for path, obj_list in scp_paths.items():
    if len(obj_list) > 1:
      top_object = obj_list.pop(0)
      for obj in obj_list:
        top_object.adopt_scope(obj)
      scp_paths[path] = top_object
    else:
      scp_paths[path] = obj_list[0]

  # Create a new master PHIL with merged scopes
  merged_scopes = [scp_paths[path] for path in scp_paths]
  phil_objects.extend(merged_scopes)
  return phil_objects


class IOTAWindowMixin(object):

  def IsFrame(self):
    is_frame = True in ['Frame' in cls.__name__ for cls in
                        self.__class__.__mro__]
    return is_frame

  def IsDialog(self):
    is_dialog = True in ['Dialog' in cls.__name__ for cls in
                        self.__class__.__mro__]
    return is_dialog

  def place_and_size(self,
                     set_size=None,
                     set_by=None,
                     center=False,
                     position=None):
    """ Place and size the frame"""
    # Find mouse position
    if set_by == 'mouse':
      self.SetPosition(wx.GetMousePosition())
    elif set_by == 'parent':
      self.SetPosition(self.set_relative_position())
    else:
      self.SetPosition((0, 0))

    if set_size:
      self.set_size(set_size)

     # Position on display
    if position:
      disp_geom = self.get_display()
      x = disp_geom[0] + position[0]
      y = disp_geom[1] + position[1]
      self.SetPosition((x, y))
    else:
      if center:
        self.Center()

  def get_display(self, window=None):
    if window is None:
      window = self.GetTopLevelParent()
    disp_idx = wx.Display.GetFromWindow(window)
    try:
      disp = wx.Display(disp_idx)
    except ValueError:
      disp = wx.Display(0)
    disp_geom = disp.GetClientArea()
    return disp_geom

  def set_size(self, set_size):

    # Determine effective minimum size
    if set_size is True:
      self.SetSize(self.GetEffectiveMinSize())
    elif type(set_size).__name__ in ('list', 'tuple'):
      assert len(set_size) == 2
      self.SetSize(wx.Size(set_size[0], set_size[1]))
    elif type(set_size).__name__ == 'Size':
      self.SetSize(set_size)
    elif type(set_size).__name__ == 'str':
      disp_geom = self.get_display()

      if set_size.lower() == 'wide':
        h_percent = 0.75
        w_percent = 0.75
      elif set_size.lower() == 'tall':
        h_percent = 0.5
        w_percent = 0.75
      elif set_size.lower == 'figure':
        h_percent = w_percent = 2 / 3
      else:
        h_percent = w_percent = 0.5

      if disp_geom[2] <= 640:
        win_w = disp_geom[2] * 0.9
      elif disp_geom[2] >= 1600:
        win_w = 800
      else:
        win_w = disp_geom[2] * w_percent

      if disp_geom[3] <= 480:
        win_h = disp_geom[3] * 0.9
      elif disp_geom[3] >= 1200:
        win_h = 600
      else:
        win_h = disp_geom[3] * h_percent

      if set_size.lower() == 'figure':
        # Initial figure should be square, with side = 2/3 of display height
        size = wx.Size(disp_geom[3] * 2/3, disp_geom[3] * 2/3)
      elif set_size.lower() == 'v_default':
        size = wx.Size(win_w, -1)
      elif set_size.lower() == 'h_default':
        size = wx.Size(-1, win_h)
      elif set_size.lower() in ('full', 'fit'):
        win_w, win_h = disp_geom[2:]
        size = wx.Size(win_w, win_h)
      elif set_size.lower() == 'dialog':
        # For dialogs, a) restrict window size to 2/3 screen height and 1/3
        # screen width; b) if height was reduced, pad the width to
        # accommodate the scroll bar.
        win_w = disp_geom[2] * 1/3
        win_h = disp_geom[3] * 2/3
        dw, dh = self.GetSize()
        if dh <= win_h:
          dlg_h = dh
          dw += 20
        else:
          dlg_h = win_h
        dlg_w = dw if dw <= win_w else win_w
        size = wx.Size(dlg_w, dlg_h)
      else:
        size = wx.Size(win_w, win_h)
      self.SetSize(size)

  def set_relative_position(self):
    """ Determines screen position w/ respect to parent window; will also
    detect if it goes beyond the display edge, and adjust """

    # Position window w/ respect to parent window
    mx, my = self.parent.GetTopLevelParent().GetScreenPosition()
    disp_geom = self.get_display(window=self.parent.GetTopLevelParent())
    dw = disp_geom[2]
    dxmin = disp_geom[0]
    dxmax = disp_geom[0] + dw
    dh = disp_geom[3]
    dymin = disp_geom[1]
    dymax = disp_geom[1] + dh

    if self.IsFrame():
      adj_x = -50
      adj_y = -50
    elif self.IsDialog():
      adj_x = int(dw * 0.025)
      adj_y = int(dh * 0.025)
    else:
      adj_x = 50
      adj_y = 50
    px = mx + adj_x
    py = my + adj_y

    ww, wh = self.GetSize()
    pxmax = px + ww
    pymax = py + wh

    # Calculate if window is going out of bounds, and adjust
    if pxmax >= dxmax:
      px = int((dxmin + dw - ww) * 0.9)
    if pymax >= dymax:
      py = int((dymin + dh - wh) * 0.9)
    if px <= dxmin:
      px = 25
    if py <= dymin:
      py = 25

    return px, py


class PHILIndex(interface.index):
  """ An index for PHIL parameters; will allow easy PHIL operations, such as
      updating from GUI elements, finding of specific scopes or definitions,
      extraction of values, etc. """

  def __init__(self, *args, **kwargs):
    super(PHILIndex, self).__init__(*args, **kwargs)

  def get_scope_phil(self, scope_name, as_string=False):
    """ Override with option to get scope from master PHIL or working PHIL
    :param scope_name: scope name
    :param as_string: set to True to return a PHIL string rather than object
    :return: PHIL scope(s)
    """
    scope_name = self.get_full_path(scope_name)
    _phil_string = str_utils.StringIO()
    scope_objects = self.get_scope_by_name(scope_name=scope_name)

    if as_string:
      if isinstance(scope_objects, list):
        for phil_object in scope_objects:
          phil_object.show(out=_phil_string)
      else :
        scope_objects.show(out=_phil_string)
      return _phil_string.getvalue()
    else:
      return scope_objects

  def get_scopes(self, scope_name=None, include=None, exclude=None):
    """ Given a scope name, return an actual scope; given an additional
        selection of paths or names, return the scope with only the specified
        objects included
    :param scope_name: scope name as a string
    :param include: scope names or paths as a list of strings to be included
                    in the returned phil object
    :param exclude: scope names or paths as a list of strings to be excluded
                    from the returned phil object (all other scopes will be
                    included; this is mutually exclusive with 'include')
    :return: scope copy (either the full scope, or included sub-scopes only,
             or all but the excluded sub-scopes)
    """

    assert (len([i for i in [include, exclude] if i]) <= 1)

    # If only a single name is passed on, that can be converted to scope
    # name; this also avoids a "custom_copy problem" for multiple scopes
    if include and len(include) == 1:
      scope_name = include[0]
      include = None

    if not scope_name:
      found_scope = self.get_scope_phil(scope_name=self.working_phil.name)
    else:
      found_scope = self.get_scope_phil(scope_name=scope_name)

    phil_objects = []
    selection = include if include else exclude
    if selection:
      if isinstance(found_scope, list):
        # For now disable scope selection from multiples (tricky to handle)
        raise Sorry('IOTA PHILIndex: Cannot include/exclude scopes, {} is a '
                    'multiple'.format(found_scope[0].full_path()))

      parent_path = found_scope.full_path()
      paths = []
      for sel in selection:
        if '.' in sel:
          paths.append(sel)
        else:
          if parent_path:
            paths.append('{}.{}'.format(parent_path, sel))
          else:
            paths.append(sel)

      if exclude:
        paths = [obj.full_path() for obj in found_scope.objects
                 if obj.full_path() not in paths]
      for path in paths:
        scp = self.get_scope_phil(scope_name=path)
        if scp:
          if isinstance(scp, list):
            phil_objects.extend(scp)
          else:
            phil_objects.append(scp)
      scope = found_scope.customized_copy(objects=phil_objects)
    else:
      scope = found_scope
    return scope

  def get_master_scope(self, path):
    if path:
      return self.master_phil.get_without_substitution(path=path)
    else:
      return self.master_phil

  def get_min_expert_level(self, scope):
    if isinstance(scope, list):
      expert_levels = []
      for scp in scope:
        expert_levels.append(self.get_expert_level(scp.full_path()))
      return min(expert_levels)
    else:
      return self.get_expert_level(scope.full_path())

  def update_phil(self, phil_string, reindex=False):
    self.update(phil_string=phil_string)
    if reindex:
      self.reindex_to_unique_paths()

  def update_from_object(self, python_object=None, reindex=False):
    if python_object:
      self.update_from_python(python_object=python_object)
      if reindex:
        self.reindex_to_unique_paths()

  def get_value(self, path, default=False):
    #fixme: multiple values don't work: defs from multiple scopes only found
    # once!
    if default:
      dfn = self.get_master_scope(path=path)
    else:
      dfn = self.get_scope_by_name(scope_name=path)
    if isinstance(dfn, list):    # is multiple!
      value = [d.type.from_words(d.words, d) for d in dfn]
    else:
      value = dfn.type.from_words(dfn.words, dfn)
    return value

  def reset_phil(self, phil=None, reindex=False):
    if phil is not None:
      self.working_phil = self.master_phil.fetch(phil)
    else:
      self.working_phil = self.master_phil.copy()
    if reindex:
      self.rebuild_index()

  def reindex_to_unique_paths(self, scope=None,
                              in_template=False,
                              in_multiple=False):
    """ If the initial indexing module encounters two scopes of the same
        name, it simply overwrites the first one with the second one. To
        correct that, this function (ran shortly after the parent class is
        initialized) goes through the master PHIL, collects all top objects,
        and concatenates the ones that have the same name. The adopt_scope
        function makes sure that the "multiple scopes with same name"
        situation is also resolved in sub-scopes. Scopes marked 'multiple'
        are tacked on to existing index entries.
    :param scope: specify a scope; otherwise it'll run on self.master_phil
    """

    if scope is None:
      scope = self.working_phil
    same_name_objects = {}

    # Loop through active objects and save in a separate dictionary,
    # grouping them under their shared names
    for obj in scope.objects:
      full_path = obj.full_path()
      if full_path in self._full_path_index:
        if obj.multiple:
          # do not include duplicate objects or objects with negative template
          # values (they would persist in PHIL, but are to be ignored in index)
          if obj.is_template >= 0:
            if obj not in self._full_path_index[full_path]:
              self._full_path_index[full_path].append(obj)
            if obj.is_scope:
              if full_path not in self._multiple_scopes:
                self._multiple_scopes[full_path] = True
              self.reindex_to_unique_paths(scope=obj)
            elif obj.is_definition and full_path not in self._multiple_defs:
                self._multiple_defs[full_path] = True
        else:
          if full_path in same_name_objects:
            same_name_objects[full_path].append(obj)
          else:
            same_name_objects[full_path] = [obj]

    # Loop through the dictionary and concatenate objects under shared names
    for full_path, obj_list in same_name_objects.items():
      top_object = obj_list.pop(0)
      if len(obj_list) > 0:
        for obj in obj_list:
          if obj.is_scope:
            top_object.adopt_scope(obj)
          elif obj.is_definition:
            top_object.adopt(obj)

      if top_object.is_template != 0:
        self._template_index[full_path] = top_object
        if top_object.is_template != -1:
          in_template = True
      elif in_template:
        self._template_index[full_path] = top_object

      # Do this recursively in case of scopes (some sub-scopes can be also
      # split between multiple entries)
      if top_object.is_scope:
        self.reindex_to_unique_paths(scope=top_object,
                                     in_template=in_template,
                                     in_multiple=in_multiple)

      # Update the full path index
      self._full_path_index[full_path] = top_object

  def create_style(self, style_string):
    """ Overrode to use customized style class with additional args """
    return style(style_string)


class IOTAPHILCtrl(object):
  """ Master class for all PHIL-compatible controls used in IOTA UI. Contains
      all the necessary functions for handling PHIL objects, so that GUI
      widgets can return PHIL strings and/or objects in a straightforward
      fashion """

  def __init__(self, phil_object):
    """ Constructor
    :param phil_object: PHIL object
    """

    assert phil_object is not None
    self.full_path = phil_object.full_path()
    self.name = phil_object.name
    self.type = phil_object.type if phil_object.is_definition else 'scope'
    self.style = phil_object.style
    self.alias = phil_object.alias
    self.multiple = phil_object.multiple
    self.helpstring = phil_object.help

    self.is_definition = phil_object.is_definition
    self.is_scope = phil_object.is_scope
    self.is_dlg_button = False

    self.is_primary_parent = phil_object.primary_parent_scope is None
    self.SetOptional(optional=phil_object.optional)

    # 'int' and 'float' types don't seem to have an attribute for whether to
    # use Auto; thus it will be set to True in all such cases. UseNone will
    # not affect the GUI for the time being before I learn more.
    if self.is_definition:
      if self.type.phil_type in ('int', 'float'):
        self.SetUseAuto()
        self.SetUseNone(enable=self.type.allow_none)
      elif self.type.phil_type in ('ints', 'floats'):
        self.SetUseAuto(enable=self.type.allow_auto_elements)
        self.SetUseNone(enable=self.type.allow_none_elements)

    # Set expert level
    self.expert_level = phil_object.expert_level
    if self.expert_level is None:
      self.expert_level = 0

  def __str__(self):
    if hasattr(self, 'name'):
      return type(self).__name__ + (" ({})".format(self.name))
    else:
      return type(self).__name__

  def SetUseAuto(self, enable=True):
    self._blank_is_auto = enable

  def SetUseNone(self, enable=True):
    self._use_none = enable

  def UseAuto(self):
    return getattr(self, "_blank_is_auto", False)

  def UseNone(self):
    return getattr(self, '_use_none', False)

  def ReformatValue(self, value=None, raise_error=True):
    """ Takes value of any format and returns a string. Returns an AutoType
        object if value is Auto or blank; returns a NonType object if value
        is None or blank (blank can be a string of spaces). If value is not
        a blank, None, or Auto, this function will return the original value
        converted to string format. """

    # Convert a "unit_cell" type to a tuple
    if 'unit_cell' in type(value).__name__:
      value = value.parameters()

    # Join a list/tuple into a single string; convert non-strings into strings
    if type(value) in (list, tuple):
      errors = []
      new_values = []
      for v in value:
        fv, e = self.check_value(value=v)
        if e:
          errors.append(e)
        new_values.append(fv)
      if errors and raise_error:
        error_msg = 'The following errors were found: \n{}' \
                    ''.format('\n'.join(errors))
        raise Sorry(error_msg)
      else:
        return ' '.join([str(v) for v in new_values])
    else:
      if type(value) != str:
        value = str(value)
      new_value, error = self.check_value(value)
      if error and raise_error:
        raise Sorry(error)
      return new_value

  def check_value(self, value):
    # Check if value is blank ('' or whitespace of any length), None, or Auto
    error = None
    if str(value).isspace():
      if self.IsOptional():
        if self.UseAuto():
          value = Auto
        else:
          value = None
      else:
        error = "Value required for {}.".format(self.GetPHILName())
        value = None
    elif str(value).lower() == 'auto':
      if self.UseAuto():
        value = Auto
      else:
        value = None
    elif str(value).lower() == 'none':
      if self.IsOptional():
        value = None
      else:
        error = "Value required for {}.".format(self.GetPHILName())
        value = None

    return value, error

  def SetOptional(self, optional=None):
    if optional is None:
      optional = True
    self.optional = optional

  def IsOptional(self):
    return getattr(self, "optional", True)

  def SetPHILName(self, name):
    self.name = name

  def GetPHILName(self):
    return getattr(self, "name", None)

  def GetPHIL(self):
    """ Override in subclasses with subclass-specific function """
    raise NotImplementedError()


class IOTADefinitionCtrl(IOTAPHILCtrl):
  """ Master class for all PHIL definition-compatible controls. """

  def __init__(self, phil_object):
    """ Constructor
    :param phil_object: PHIL definition object
    """
    if isinstance(phil_object, list):
      phil_object = phil_object[0]
    assert phil_object.is_definition
    IOTAPHILCtrl.__init__(self, phil_object=phil_object)
    self.default_value = phil_object.type.from_words(phil_object.words,
                                                     phil_object)

  def GetPHIL(self, full_path=False, indent_length=2):
    """ Get PHIL entry for this control
    :param full_path: full object path, e.g. scope.sub_scope.definition
    :param indent_length: number of spaces per indent
    :return: if full_path=True, will return the full object path; otherwise,
    will return the object name, indented by indent_length * number of levels
    """

    assert self.name
    if full_path:
      return_string = self.full_path
    else:
      # Construct PHIL string from current value
      if self.multi_definition:
        phil_strings = []
        for idx, ctrl in self.controls.items():
          phil_name = ctrl.full_path.split(".")[-1]
          value = ctrl.GetStringValue()
          value_string = '{} = {}'.format(phil_name, value)
          indent = len(ctrl.full_path.split('.')) * indent_length
          phil_string = '{:{ind}}{}'.format(' ', value_string, ind=indent)
          phil_strings.append(phil_string)
        return_string = '\n'.join(phil_strings)
      else:
        phil_name = self.full_path.split(".")[-1]
        value = self.GetStringValue()
        phil_string = '{} = {}'.format(phil_name, value)
        indent = len(self.full_path.split('.')) * indent_length
        return_string = '{:{ind}}{}'.format(' ', phil_string, ind=indent)
    return return_string

  def value_from_words(self, phil_object):
    """ Extract value, type, and whether value is None/Auto from PHIL object """
    value = phil_object.type.from_words(phil_object.words, phil_object)
    phil_type = phil_object.type.phil_type
    is_none_auto = (value is Auto or value is None)
    return value, phil_type, is_none_auto

  def GetStringValue(self):
    """ Override in subclasses with subclass-specific function """
    raise NotImplementedError()

  def SetValue(self, value):
    """ Overwrite in subclasses with subclass-specific function """
    raise NotImplementedError()


class IOTAScopeCtrl(IOTAPHILCtrl):
  """ Master class for all PHIL scope-compatible controls """

  def __init__(self, scope):
    """ Constructor
    :param scope: PHIL scope object
    """
    self.SetPHIL(scope=scope)
    self.controls = {}

  def SetPHIL(self, scope):
    if isinstance(scope, list):
      scope = scope[0]
    assert scope.is_scope
    super(IOTAScopeCtrl, self).__init__(phil_object=scope)

  # def GetFullPHIL(self):
  #   """ Get PHIL entry for this control; unlike GetPHIL, this function will
  #       the parameters in a full_path + value format (e.g.
  #       indexing.known_symmetry.space_group='P 1 21 1') in a combined list,
  #       which can be parsed by libtbx.phil and then fetched to the main
  #       scope. This approach will ensure that parameters returned by sub-scopes
  #       can be matched to their global paths. (The GetPHIL approach loses the
  #       full path information during formatting.)
  #   :return: a string of full_path+value formatted parameters delimited by a
  #            line break
  #   """
  #
  #   paths = self.add_scp_paths(self)
  #   phil_string = '\n'.join(paths)
  #   return phil_string
  #
  #
  # def add_scp_paths(self, phil_ctrl):
  #   """ Recursively go through all controls in this scope control, extract
  #       their PHIL strings (or scopes) and create a list of PHIL strings in a
  #       full_path + value format
  #   :param phil_ctrl: a PHIL control
  #   :return: a list of full paths with values
  #   """
  #   phil_paths = []
  #
  #   for ctrl in phil_ctrl.controls:
  #     if ctrl.is_definition:
  #       dfn_path = self.add_dfn_path(ctrl)
  #       phil_paths.append(dfn_path)
  #     elif ctrl.is_scope:
  #       scp_paths = ctrl.add_scp_paths(ctrl)
  #       phil_paths.extend(scp_paths)
  #     elif ctrl.is_dlg_button:
  #       btn_paths = ctrl.get_phil_strings()
  #       phil_paths.extend(btn_paths)
  #   return phil_paths
  #
  # def add_dfn_path(self, dfn):
  #   """ Add a single PHIL string for a given definition control
  #   :param dfn: definition-type PHIL control
  #   :return: a single PHIL string in a full_path + value format
  #   """
  #   try:
  #     # Construct PHIL string from current value
  #     phil_name = dfn.full_path
  #     value = dfn.GetStringValue()
  #     phil_string = '{}={}'.format(phil_name, value)
  #   except Exception:
  #     # Use the default PHIL string
  #     phil_string = dfn.phil_string.replace('\n', '')
  #   return phil_string

  def GetPHIL(self, path_only=False, expand=False, multiple=False,
              indent_length=2):
    """ Get PHIL entry for this control
    :param path_only: return only the scope path (full path will be used)
    :param indent_length: number of spaces per indent
    :param multiple: mark multiple PHIL scope; requires different handling
    :param expand: if True, expand PHIL string to cover the full path
    :return: if path_only=True, return the scope path; else return the PHIL
    entry for the entire scope as a PHIL-formatted string
    """

    if path_only:
      full_phil_string = self.full_path
    else:
      assert self.controls
      if expand:
        indent = len(self.full_path.split('.')) * indent_length
      else:
        indent = 0
      phil_entries = self.add_scp_string(phil_ctrl=self, indent=indent,
                                         indent_length=indent_length)
      phil_string = '\n'.join(phil_entries)

      # check if this is a sub-scope, and pre-pend elements of the full path
      scope_paths = self.full_path.split('.')
      if expand and len(scope_paths) > 1:
        prepended_strings = []
        appended_strings = []
        level = 0
        close_level = len(scope_paths) - 2
        for path in scope_paths[:-1]:
          indent = level * indent_length
          if indent == 0:
            path_string = '{} {{'.format(path)
          else:
            path_string = '{:{ind}}{} {{'.format(' ', path, ind=indent)
          prepended_strings.append(path_string)

          # create closing bracket strings
          close_indent = close_level * indent_length
          if close_indent == 0:
            close_string = '}'
          else:
            close_string = '{:{ind}}}}'.format(' ', ind=close_indent)
          appended_strings.append(close_string)

          level += 1
          close_level -= 1

        open = '\n'.join(prepended_strings)
        close = '\n'.join(appended_strings)
        full_phil_string = '{}\n{}\n{}'.format(open, phil_string, close)
      else:
        full_phil_string = phil_string

    return full_phil_string

  def add_scp_string(self, phil_ctrl, indent, indent_length):
    """ Recursively go through all controls in this scope control, extract
        their PHIL strings (or scopes) and create a list of PHIL strings
    :param phil_ctrl: a scope-type PHIL control object
    :param indent: total number of spaces by which to offset a PHIL string
    :param indent_length: number of spaces per indent 'tab'
    :return: list of PHIL strings
    """

    if self.multi_scope or self.is_primary_parent:
      phil_strings = []
      end_point = ''
    else:
      if indent == 0:
        phil_strings = ['{} {{'.format(phil_ctrl.name)]
        end_point = '}'
      else:
        phil_strings = ['{:{ind}}{} {{'.format(' ', phil_ctrl.name, ind=indent)]
        end_point = '{:{ind}}}}'.format(' ', ind=indent)
      indent += indent_length

    for idx, ctrl in phil_ctrl.controls.items():
      if ctrl.is_definition:
        dfn_string = self.add_dfn_string(ctrl)
        phil_strings.append(dfn_string)
      elif ctrl.is_scope:
        scp_strings = ctrl.add_scp_string(ctrl, indent, indent_length)
        phil_strings.extend(scp_strings)
      elif ctrl.is_dlg_button:
        btn_strings = ctrl.get_phil_strings()
        indented_strings = ["{:{ind}}{}".format('', s, ind=indent)
                            for s in btn_strings]
        phil_strings.extend(indented_strings)

    if end_point:
      phil_strings.append(end_point)

    return phil_strings

  def add_dfn_string(self, dfn):
    """ Add a single PHIL string for a given definition control """
    # return_string = '{:{ind}}{}'.format(' ', dfn.GetPHIL(), ind=indent)
    return dfn.GetPHIL()

  def _create_controls(self, scope):
    """ Create IOTA PHIL controls for a given scope (mainly for test
        purposes, as actual non-generic controls will be created by IOTA)
    :param scope: scope-type PHIL control
    """
    for obj in scope.active_objects():
      if obj.is_definition:
        self.controls.append(IOTADefinitionCtrl(obj))
      else:
        self.controls.append(IOTAScopeCtrl(obj))

  def paths_from_scope(self, scope):
    """ Generate PHIL strings (full path + value) from any scope, rather than
        controls. This is most useful for generating lists of PHIL strings from
        buttons.
    :param scope: a PHIL scope
    :return: a list of full_path+value formatted PHIL strings for individual
             parameters
    """

    phil_strings = []
    for obj in scope.active_objects():
      if obj.is_definition:
        value = self.ReformatValue(value=obj.type.from_words(obj.words, obj))
        dfn_string = '{}={}'.format(obj.full_path(), value)
        phil_strings.append(dfn_string)
      elif obj.is_scope:
        scope_strings = self.paths_from_scope(scope=obj)
        phil_strings.extend(scope_strings)

    return phil_strings

class WidgetHandlerMixin(object):
  """ A mixin with sizer-specific widget-handling methods """

  def reset_layout(self):
    if isinstance(self, wx.Sizer):
      self.Layout()

  def HideAll(self):
    self.ShowItems(False)

  def ShowAll(self):
    self.ShowItems(True)

  def add_panel(self, panel, border=5):
    self.Add(panel, proportion=1, flag=wx.ALL | wx.EXPAND, border=border)

  def add_widget(self, widget, border=5, proportion=0, expand=False,
                 center=False):
    flags = wx.RIGHT | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL
    if expand:
      flags |= wx.EXPAND
    if center:
      flags |= wx.ALIGN_CENTER_HORIZONTAL

    self.Add(widget, proportion, flags, border)

  def add_row_widget(self, widget, proportion=0, border=5):
    self.Add(widget, proportion, wx.LEFT | wx.BOTTOM | wx.ALIGN_RIGHT,
             border)

  def add_expanding_widget(self, widget, proportion=0, border=5):
    self.Add(widget, proportion, wx.ALL | wx.EXPAND, border)

  def add_widgets(self,
                  widget_list,
                  flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL,
                  border=5):
    for widget in widget_list:
      self.Add(widget, 0, flag, border)

  def add_labeled_widget(self, widget, label, label_size=wx.DefaultSize,
                         expand=False, h_space=0):
    flag = wx.ALIGN_CENTER_HORIZONTAL
    if expand:
      flag |= wx.EXPAND
    w_sizer = wx.FlexGridSizer(1, 2, 0, h_space)
    self.label = wx.StaticText(self.parent, label=label, size=label_size)
    w_sizer.Add(self.label)
    w_sizer.Add(widget, flag=flag)

  def add_widget_and_label(self, widget, label, label_size=wx.DefaultSize,
                           *args, **kwargs):
    self.label = wx.StaticText(self.parent, label=label, size=label_size)
    self.add_widget(widget=self.label)
    self.add_widget(widget=widget, *args, **kwargs)

  def __str__(self):
    return type(self).__module__ + "." + type(self).__name__


class TextCtrlValidator(Validator):
  """ A custom validator for text-based controls (including numbers, paths,
      unit cells, and space groups). Runs each control's CheckFormat function
      and responds accordingly to its results. """

  def __init__(self):
    super(TextCtrlValidator, self).__init__()
    self.Bind(wx.EVT_TEXT_ENTER, self.OnEnter)
    self.ctrl = None

  def Clone(self):
    return self.__class__()

  def TransferToWindow(self):
    return True

  def TransferFromWindow(self):
    return True

  def CheckFormat(self, value):
    if type(value) not in (list, tuple):
      value = str(value)
    return self.ctrl.CheckFormat(value)

  def Validate(self, parent):
    self.ctrl = self.GetWindow()
    self.parent = parent

    try:
      raw_value = self.ctrl.GetValue()
      adj_value = self.parent.ReformatValue(value=raw_value)
      value = self.CheckFormat(value=adj_value)
      if type(value) in (list, tuple):
        reformatted = to_unicode(' '.join([str(v) for v in value]))
      else:
        reformatted = to_unicode(str(value))
    except UnicodeEncodeError:
      self.ctrl.error_msg = "Only the standard UTF-8 character set is allowed."
      return False
    except Exception as e:
      self.ctrl.error_msg = str(e)
      return False
    else:
      self.ctrl.error_msg = None
      self.ctrl.SetValue(reformatted)
      return True

  def OnEnter(self, event):
    event.Skip()

# -- end
