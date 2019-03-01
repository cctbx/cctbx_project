from __future__ import absolute_import, division, print_function

# various containers for gui directives and objects

class menu_hierarchy(object):
  is_menu_item = False
  is_submenu = True

  def __init__(self, menu_name):
    self.menu_name = menu_name
    self.submenus = {}
    self.menu_items = []
    self._used = []
    self.submenu_items = []

  def get_items(self):
    return self.menu_items + self.submenu_items

  def add_menu_item(self, menu_item_name):
    if not menu_item_name in self._used :
      self.menu_items.append( menu_item(menu_item_name) )
      self._used.append(menu_item_name)

  def get_submenu(self, submenu_name):
    return self.submenus[submenu_name]

  # this will simply do nothing if the submenu already exists
  def add_submenu(self, submenu_name):
    if submenu_name in self.submenus :
      pass
    else :
      self.submenus[submenu_name] = menu_hierarchy(submenu_name)
      self.submenu_items.append(self.submenus[submenu_name])

  def __str__(self):
    return self.menu_name

class menu_item(object):
  is_menu_item = True
  is_submenu = False

  def __init__(self, menu_item_name):
    self.menu_item_name = menu_item_name

  def __str__(self):
    return self.menu_item_name

class style(object):
  """
  Container for flags used to alter the appearance of controls in the Phenix
  GUI.  These can either be booleans (style_args) or other values (style_kwds).
  """
  style_args = [
    "auto_align", # left-align controls and labels separately (scope)
    "bold", # bold text label (definition)
    "box",  # display controls in wx.StaticBox (scope)
    "color",
    "date", # control specifies a date (definition)
    "use_list", # use wx.ListCtrl-based widget (atom_selection, multiple=True)
    "menu_item", # add to Settings menu (scope)
    "narrow", # fit controls in 300px (path)
    "noauto", # don't automatically display as part of parent scope (any)
    "noedit", # not editable (definition)
    "none_is_auto", # treat None as Auto
    "scrolled",
    "selection", # tags a str definition as atom_selection
    "spinner", # attach wx.SpinCtrl widget
    "submenu", # add to Settings menu as submenu (scope)
    "tribool", # use wx.Choice with True/False/Auto (bool)
    "hidden",  # never display in GUI (definition)
    "directory", # specifies a directory (path)
    "new_file", # specifies a new file (path)
    "default_cwd", # default to current directory (path)
    "resolution", # treat as resolution limit (float)
    "hide_label", # don't display control label
    "not_none", # don't allow None (definition)
    "fixed",
    "checklist", # display as wx.CheckListBox (definition)
    "set_resolution",
    "force_data",
    "anom",
    "non_anom",
    "no_view",
    "process_hkl", # specifies a reflections file with expected data items,
                   # which will be automatically extracted
    "force_rfree",
    "combo_box", # use a wx.ComboBox control (definition)
    "optional",
    "no_map",
    "file_type_default", # specifies that all files of the stated file_type
                         # should be automatically associated with this path
                         # parameter.  only used in ListCtrl-based file input
                         # fields
    "expand",
    "output_dir", # specifies that the path defines the output directory, which
                  # is determined automatically.  only one of these is allowed
                  # per program, and the GUI will extract this automatically.
    "seq_file",
    "single_input", # disable multiple controls (definition, multiple=True)
  ]
  style_kwds = [
    "auto_launch_dialog", # when clicked, launch a window to edit options in
                          # the named scope.  mostly used in phenix.refine.
                          # (definition, especially bool)
    "columns", # number of columns for controls (scope)
    "dialog_link", # add a button to edit the named scope (definition)
    "extensions",
    "file_type", # specifies the file type filter, as well as the behavior of
                 # the ListCtrl-based file input fields (path)
    "min",
    "max",
    "parent_submenu", # add to the specified submenu (scope)
    "OnUpdate", # specifies event handler in
                # phenix/wxGUI2/Programs/Extensions.py (definition)
    "OnChange", # specifies event handler in
                # phenix/wxGUI2/Programs/Extensions.py (definition)
    "renderer", # specifies custom drawing method in
                # phenix/wxGUI2/Programs/CustomControls.py (definition)
    "caption_img", # image to display next to caption (scope)
    "rlabel",
    "llabel",
    "cols",
    "height",
    "help_page", # file name in Phenix documentation (definition)
    "help_anchor", # anchor on help_page (definition)
    "caption_width",
  ]
  multiple_kwds = ["child", "parent"]
  convert_to_int = ["columns", "min", "max"]
  special_words = []
  #__slots__ = ["_style_words"] + style_args + style_kwds

  def __init__(self, style_string=None):
    self._style_words = []
    self.child = None
    self.parent = None
    for arg in self.style_args :
      setattr(self, arg, False)
    for kwd in self.style_kwds :
      setattr(self, kwd, None)

    if style_string is not None :
      style_words = style_string.split()
      for style_word in style_words :
        if style_word in self.style_args :
          setattr(self, style_word, True)
        elif not style_word in self.special_words :
          try :
            fields = style_word.split(":")
            kwd = fields[0]
            value = ":".join(fields[1:])
            #kwd, value = style_word.split(":")
            if kwd in self.multiple_kwds :
              if getattr(self, kwd, None) is not None :
                getattr(self, kwd).append(value)
              else :
                setattr(self, kwd, [value])
            elif kwd in self.special_words :
              self.process_style_keyword(style_word)
            else :
              value_words = value.split(",")
              if len(value_words) == 1 :
                if kwd in self.convert_to_int :
                  try :
                    int_value = int(value)
                    setattr(self, kwd, int_value)
                  except ValueError :
                    print("Ignoring incorrect .style keyword '%s'"%style_word)
                else :
                  setattr(self, kwd, value)
              else :
                setattr(self, kwd, value_words)
          except Exception as e :
            print(e)
            pass
            #wxGUI2.DEBUG2("Unrecognized style keyword(s): %s" % style_string)
            #wxGUI2.DEBUG2(style_string)
        else :
          self.process_style_keyword(style_word)
      self._style_words = style_words

  def process_style_keyword(self, style_word):
    pass

  def __str__(self):
    return " ".join(self._style_words)

  def get_list(self, kwd):
    val = getattr(self, kwd, None)
    if isinstance(val, list):
      return val
    elif val is not None :
      return [val]
    else :
      return []

  def get_child_params(self):
    params = []
    if self.child is not None :
      assert isinstance(self.child, list)
      for value in self.child :
        child_type, child_name = value.split(":")
        params.append((child_type, child_name))
    return dict(params)

  def get_parent_params(self):
    params = []
    if self.parent is not None :
      assert isinstance(self.parent, list)
      for value in self.parent :
        parent_type, parent_name = value.split(":")
        params.append((parent_type, parent_name))
    return dict(params)

default_style = style(None)

class file_type_map(object):
  def __init__(self, param_info, default_label=None):
    self.param_info = param_info
    self.default_label = default_label
    self._names = {}
    self._labels = {}
    self._counts = {}
    self._max_allowed = None
    for args in param_info :
      name = args[0]
      label = args[1]
      self._names[label] = name
      self._labels[name] = label
      if (len(args) >= 3):
        count = args[2]
      else :
        assert (len(args) == 2)
        count = 1
      self._counts[name] = count
    assert (default_label is None) or (default_label in self._names.keys())

  def get_params_and_labels(self):
    return [ (args[0], args[1]) for args in self.param_info ]

  def get_overall_max_count(self):
    n = 0
    for name, count in self._counts.items():
      if (count is None):
        return None
      n += count
    return n

  def get_max_count(self, param_name):
    return self._counts[param_name]

  def get_default_label(self):
    if (len(self._names) == 1):
      return list(self._names.keys())[0]
    return self.default_label

  def get_default_param(self):
    if (len(self._labels) == 1):
      return list(self._labels.keys())[0]
    if (self.default_label is not None):
      return self._names[self.default_label]
    return None

  def get_param_names(self):
    return [ args[0] for args in self.param_info ]

  def get_multiple_params(self):
    names = []
    for name, count in self._counts.items():
      if count is None :
        names.append(name)
    return names

  def disable_param(self, param_name):
    if (param_name in self._labels):
      label = self._labels[param_name]
      self._names.pop(label)
      self._labels.pop(param_name)
      self._counts.pop(param_name)
