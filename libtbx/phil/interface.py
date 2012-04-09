
# XXX: this module is used exclusively by the Phenix GUI, which needs an
# index of all current phil parameters, and an easy way to change them.

import os, sys, re, string
import libtbx.phil
from libtbx.phil import gui_objects
from libtbx.utils import Sorry
from libtbx import easy_pickle, str_utils, smart_open
from libtbx import adopt_init_args, Auto

tracking_params = libtbx.phil.parse("""
  job_title = None
    .type = str
    .input_size = 400
    .help = Job title in PHENIX GUI, not used on command line
""")

class index (object) :
  def __init__ (self, master_phil, working_phil=None, parse=None,
      fetch_new=False) :
    adopt_init_args(self, locals())
    self._states = []
    self._full_path_index = {}
    self._full_text_index = {}
    self._template_index = {}
    self._multiple_scopes = {}
    self._multiple_defs = {}
    self._expert_levels = {}
    self._input_files = []
    self._hidden = [] # XXX: not implemented here (for phenix GUI)
    self._phil_has_changed = False
    self._log = str_utils.StringIO()
    self._prefix = None
    if parse is None :
      self.parse = libtbx.phil.parse
    self.setup_phil(working_phil, fetch_new)
    self.parse_styles()

  def setup_phil (self, working_phil, fetch_new=False) :
    if working_phil is None :
      self.working_phil = self.master_phil.fetch()
    elif fetch_new :
      self.working_phil = self.master_phil.fetch(source=working_phil)
    else :
      self.working_phil = working_phil
    self.build_index(collect_multiple=True)
    self.params = self.working_phil.extract()

  def set_prefix (self, prefix) :
    assert (prefix is None) or (isinstance(prefix, str))
    self._prefix = prefix

  def get_full_path (self, phil_path) :
    if (self._prefix is None) :
      return phil_path
    elif (phil_path.startswith(".")) :
      return self._prefix + phil_path
    else :
      return phil_path

  def clear_index (self, clear_multiple=True) :
    self._full_path_index = {}
    self._full_text_index = {}
    self._template_index = {}
    if clear_multiple :
      self._multiple_scopes = {}
      self._multiple_defs = {}

  def clear_state_stack (self) :
    while len(self._states) > 0 :
      saved_phil = self._states.pop()
      del saved_phil

  def push_state (self) :
    self._states.append(self.working_phil.fetch())
    return len(self._states) - 1

  def pop_state (self) :
    if len(self._states) > 0 :
      self.working_phil = self._states.pop()
      self.rebuild_index()
      return True
    return False

  def set_state (self, index) :
    assert (index >= 0)
    if (len(self._states) == 0) :
      pass
    else :
      self.working_phil = self._states[index].fetch()
      self.rebuild_index()
      return True
    return False

  def build_index (self, collect_multiple=False) :
    self.clear_index(clear_multiple=collect_multiple)
    index_phil_objects(phil_object=self.working_phil,
                       path_index=self._full_path_index,
                       text_index=self._full_text_index,
                       template_index=self._template_index,
                       multiple_scopes=self._multiple_scopes,
                       multiple_defs=self._multiple_defs,
                       collect_multiple=collect_multiple,
                       expert_levels=self._expert_levels,
                       input_files=self._input_files)

  def rebuild_index (self, only_scope=None) :
    self._full_path_index = {}
    reindex_phil_objects(phil_object=self.working_phil,
                         path_index=self._full_path_index,
                         only_scope=only_scope)

  def save_param_file (self,
                       file_name,
                       sources=None,
                       extra_phil="",
                       diff_only=False,
                       save_state=False,
                       replace_path=None) :
    if sources is None :
      sources = []
    if extra_phil != "" :
      self.merge_phil(phil_string=extra_phil, rebuild_index=False)
    final_phil = self.master_phil.fetch(sources=[self.working_phil] +
      list(sources))
    if diff_only :
      output_phil = self.master_phil.fetch_diff(source=final_phil)
    else :
      output_phil = final_phil
    if (replace_path is not None) :
      substitute_directory_name(
        phil_object=output_phil,
        path_name=replace_path,
        sub_name="LIBTBX_BASE_DIR")
    try :
      f = smart_open.for_writing(file_name, "w")
    except IOError, e :
      raise Sorry(str(e))
    else :
      if (replace_path is not None) :
        f.write("LIBTBX_BASE_DIR = \"%s\"\n" % replace_path)
      output_phil.show(out=f)
      f.close()
    if save_state :
      cache_file = "%s_cache.pkl" % file_name
      easy_pickle.dump(cache_file, self)

  def get_combined_phil (self, sources=None) :
    if sources is None :
      sources = []
    new_phil = self.master_phil.fetch(sources=[self.working_phil]+sources)
    return new_phil

  def save_diff (self, file_name, replace_path=None) :
    self.save_param_file(file_name=file_name,
      diff_only=True,
      replace_path=replace_path)

  def get_diff (self) :
    return self.master_phil.fetch_diff(source=self.working_phil)

  def copy (self, preserve_changes=True) :
    if preserve_changes :
      return easy_pickle.loads(easy_pickle.dumps(self))
    else :
      return self.copy_master()

  def copy_master (self) :
    return index(master_phil=self.master_phil, parse=self.parse)

  #---------------------------------------------------------------------
  # Conversion to/from extracted objects
  def update_from_python (self, python_object=None) :
    if python_object is None :
      if self.params is not None :
        python_object = self.params
      else :
        return False
    self.push_state()
    self.working_phil = self.master_phil.format(python_object=python_object)
    self.rebuild_index()

  def get_python_object (self, make_copy=False) :
    assert self.working_phil is not None
    if make_copy :
      return self.working_phil.extract()
    elif self._phil_has_changed or self.params is None :
      self.params = self.working_phil.extract()
      self._phil_has_changed = False
    return self.params

  def get_python_from_other (self, phil_object=None, file_name=None,
      phil_string=None) :
    assert [phil_object, file_name, phil_string].count(None) < 3
    try :
      if phil_object is None :
        if file_name is not None :
          phil_object = self.parse(file_name=file_name)
        elif phil_string is not None :
          phil_object = self.parse(phil_string)
      new_phil = self.master_phil.fetch(source=phil_object)
    except KeyboardInterrupt :
      raise
    except Exception, e :
      return None
    else :
      return new_phil.extract()

  def get_python_from_params (self, phil_object) :
    return self.get_python_from_params(phil_object=phil_object)

  def get_python_from_file (self, file_name) :
    return self.get_python_from_other(file_name=file_name)

  def get_python_from_string (self, phil_string) :
    return self.get_python_from_other(phil_string=phil_string)

  #---------------------------------------------------------------------
  # Retrieval methods
  # FIXME: not working properly for ncs restraint group phil
  def get_scope_by_name (self, scope_name, phil_parent=None) :
    scope_name = self.get_full_path(scope_name)
    if scope_name in self._full_path_index :
      indexed_phil_objects = self._full_path_index[scope_name]
      if isinstance(indexed_phil_objects, list) and phil_parent is not None :
        child_objects = []
        for object1 in indexed_phil_objects :
          for object2 in phil_parent.objects :
            if object1 is object2 :
              child_objects.append(object1)
        return child_objects
      else :
        return indexed_phil_objects
    else :
      return None

  def get_root_scope_names (self) :
    paths = []
    for phil_object in self.working_phil.objects :
      paths.append(phil_object.full_path())
    return paths

  def get_template_copy (self, phil_name) :
    phil_name = self.get_full_path(phil_name)
    template = self._template_index.get(phil_name)
    new_copy = None
    if template is not None :
      new_copy = template.customized_copy()
    return new_copy

  def get_validated_param (self, def_name) :
    phil_name = self.get_full_path(def_name)
    phil_objects = self.get_scope_by_name(def_name)
    if isinstance(phil_objects, list) :
      vals = []
      for obj in phil_objects :
        python_value = obj.extract()
        phil_value = python_value.format()
        vals.append(phil_value)
      return vals
    elif phil_objects is not None :
      #phil_objects.show()
      python_value = phil_objects.extract()
      return python_value

  # TODO: test this
  def validate_format (self, def_name, str_value) :
    if str_value is None :
      return None
    indexed_def = self.get_scope_by_name(def_name)
    d = indexed_def
    if isinstance(indexed_def, list) :
      d = indexed_def[0]
    elif indexed_def is None :
      return None
    proxy = d.validate_and_format(input_string=str_value, source_info="GUI")
    if proxy.error_message is not None :
      raise Exception(proxy.error_message)
    else :
      return str(proxy.formatted.extract())

  def get_label (self, def_name) :
    phil_text = self._full_text_index.get(def_name, None)
    if (phil_text is not None) :
      label = phil_text[0]
      if (label is not None) :
        return label
    scope = self.get_scope_by_name(def_name)
    if isinstance(scope, list) :
      scope = scope[0]
    return get_standard_phil_label(scope)

  def search_phil_text (self, search_text, match_all=False, labels_only=True) :
    fields = search_text.split()
    for word in fields :
      # this allows matching of phil param paths
      if re.search("[^a-zA-Z\.\_]", word) is not None :
        raise Sorry("Invalid search string '%s'." % word)
    regex_list = [ re.compile(word, re.I) for word in fields ]
    matching_defs = []
    n_words = len(regex_list)
    for phil_name, phil_text in self._full_text_index.iteritems() :
      (label, caption, help, is_def) = phil_text
      if (phil_name in self._hidden) or (not is_def) :
        continue
      n_found = 0
      if labels_only :
        for regex in regex_list :
          if (regex.search(label) is None) :
            if match_all : break
            else :         continue
          n_found += 1
      else :
        for regex in regex_list :
          if ((regex.search(label) is None) and
              (regex.search(phil_name) is None) and
              (regex.search(caption) is None) and
              (regex.search(help) is None)) :
            if match_all : break
            else :         continue
          n_found += 1
      if (match_all and (n_found == n_words)) or (n_found > 0) :
        matching_defs.append(phil_name)
    return matching_defs

  def get_scope_phil (self, scope_name) :
    scope_name = self.get_full_path(scope_name)
    _phil_string = str_utils.StringIO()
    scope_objects = self.get_scope_by_name(scope_name)
    if isinstance(scope_objects, list) :
      for phil_object in scope_objects :
        phil_object.show(out=_phil_string)
    else :
      scope_objects.show(out=_phil_string)
    return _phil_string.getvalue()

  def get_phil_help_string (self, phil_name) :
    phil_name = self.get_full_path(phil_name)
    if phil_name in self._full_text_index :
      (label, caption, help, is_def) = self._full_text_index[phil_name]
      return str(help)
    else :
      return None

  def get_expert_level (self, phil_name) :
    phil_name = self.get_full_path(phil_name)
    return self._expert_levels.get(phil_name, 0)

  def get_input_files (self) :
    files = []
    for def_name in self._input_files :
      phil_object = self.get_scope_by_name(def_name)
      label = self.get_label(def_name)
      #phil_text = self._full_text_index[def_name] #.get(def_name, [None]*4)
      #(label, caption, help, is_def) = phil_text
      if isinstance(phil_object, list) :
        for def_copy in phil_object :
          assert (def_copy.is_definition)
          file_name = def_copy.extract()
          if (file_name is not None) and (file_name is not Auto) :
            if isinstance(file_name, list) :
              for fn in file_name :
                files.append((fn, label, def_name))
            else :
              files.append((file_name, label, def_name))
      else :
        assert phil_object.is_definition
        file_name = phil_object.extract()
        if (file_name is not None) and (file_name is not Auto) :
          if isinstance(file_name, list) :
            for fn in file_name :
              files.append((fn, label, def_name))
          else :
            files.append((file_name, label, def_name))
    return files

  def get_run_title (self) :
    for def_name in self._full_path_index.keys() :
      if ((def_name in ["title", "job_title"]) or
          (def_name.endswith(".title") or def_name.endswith(".job_title"))) :
        phil_def = self._full_path_index[def_name]
        assert phil_def.is_definition
        return phil_def.extract()
    return None

  def is_list_type (self, phil_scope_name) :
    phil_object = self.get_scope_by_name(phil_scope_name)
    if isinstance(phil_object, list) :
      return False
    elif phil_object.type.phil_type in ["strings", "ints", "floats"] :
      return True
    return False

  #---------------------------------------------------------------------
  # EDITING METHODS
  def substitute_directory (self, path_name, sub_name) :
    substitute_directory_name(
      phil_object=self.working_phil,
      path_name=path_name,
      sub_name=sub_name)

  def reset_scope (self, phil_scope_name) :
    old_phil = self.working_phil
    delete_phil_objects(old_phil, [phil_scope_name])
    self.working_phil = self.master_phil.fetch(sources=[old_phil])

  def merge_phil (self,
                  phil_object=None,
                  phil_string=None,
                  phil_file=None,
                  overwrite_params=True,
                  rebuild_index=True,
                  only_scope=None) :
    assert ([phil_object, phil_string, phil_file].count(None) == 2)
    if (phil_string is not None) :
      phil_object = self.parse(phil_string)
    elif (phil_file is not None) :
      phil_object = self.parse(file_name=phil_file)
    if (phil_object is not None) :
      old_phil = self.working_phil
      if overwrite_params :
        new_paths = []
        get_all_path_names(phil_object, new_paths)
        redundant_paths = []
        for path in new_paths :
          if path in self._multiple_scopes or path in self._multiple_defs :
            redundant_paths.append(path)
        if len(redundant_paths) > 0 :
          delete_phil_objects(old_phil, redundant_paths, only_scope=only_scope)
      self.log2("Fetching new working phil")
      new_phil = None
      if False : # XXX this was causing too many problems
      # if (only_scope is not None) :
        new_scope = phil_object.get(only_scope)
        scope_master = self.master_phil.get(only_scope)
        fetched_scope = scope_master.fetch(source=new_scope)
        #fetched_scope.show()
        find_and_replace_scope(
          current_phil=self.working_phil,
          new_scope=fetched_scope,
          scope_name=only_scope)
        new_phil = self.working_phil
      else :
        new_phil = self.master_phil.fetch(sources=[old_phil, phil_object])
      if (new_phil is not None) :
        self.working_phil = new_phil
        if rebuild_index :
          self.log2("rebuilding index")
          self.rebuild_index(only_scope=only_scope)
      else :
        self.log("*** ERROR: new phil object is empty")
      self._phil_has_changed = True
      self.params = None

  def erase_scope (self, phil_scope) :
    delete_phil_objects(self.working_phil, [phil_scope])

  # Safe wrapper of merge_phil for loading parameter files from GUI
  def merge_param_file (self, file_name) :
    if not os.path.isfile(file_name) :
      raise Sorry("The path %s does not exist or is not a file." % file_name)
    try :
      phil_object = self.parse(file_name=file_name)
    except KeyboardInterrupt :
      raise
    except Exception, e :
      self.log(e)
      raise Sorry("This parameter file could not be parsed correctly.")
    try :
      new_phil = self.master_phil.fetch(source=phil_object)
    except KeyboardInterrupt :
      raise
    except Exception, e :
      self.log(e)
      self.log(open(file_name).read())
      raise Sorry("This file contains invalid parameters for this program. "+
                  "Check the manual for a list of allowed parameters "+
                  "for each module.")
    self.merge_phil(phil_object=phil_object)

  # Safe wrapper of merge_phil for phil strings
  def update (self, phil_string, only_scope=None) :
    try :
      phil_object = self.parse(phil_string)
      new_phil = self.master_phil.fetch(source=phil_object)
    except KeyboardInterrupt :
      raise
    except Exception, e :
      self.log(str(e))
      self.log(str(phil_string))
      raise Sorry("An unknown error occurred parsing internal parameters. "+
                  "This is probably a bug; if the program was launched with "+
                  "the argument --debug, further information will be printed "+
                  "to the console.")
    self.merge_phil(phil_object=phil_object, only_scope=only_scope)

  #---------------------------------------------------------------------
  # DEBUG/TEST METHODS
  def check_scopes (self, phil_names) :
    for phil_name in phil_names :
      if self.get_scope_by_name(phil_name) is None :
        raise AttributeError("Scope %s does not exist!" % phil_name)
    return True

  def log (self, message) :
    self._log.write(message + "\n")

  def log2 (self, message) :
    f = sys._getframe(1)
    filename = os.path.basename(f.f_code.co_filename)
    self._log.write("%s (%s:%d): %s\n" %
      (f.f_code.co_name, filename, f.f_lineno, str(message).strip()))

  #---------------------------------------------------------------------
  # GUI style handling
  def parse_styles (self) :
    self.style = {}
    self._event_handlers = {}
    self._update_handlers = {}
    self._renderers = {}
    self._file_type_mappings = {}
    self._menu_tree = gui_objects.menu_hierarchy("settings")
    self.generate_gui_components(self.working_phil)

  def create_style (self, style_string) :
    return gui_objects.style(style_string)

  def generate_gui_components (self, phil_scope, in_submenu=False,
      current_menu=None) :
    use_submenu = in_submenu
    if not current_menu :
      current_menu = self._menu_tree
    if phil_scope.is_template < 0 :
      return
    for object in phil_scope.objects :
      next_menu = None
      full_object_path = object.full_path()
      if object.style is not None and phil_scope.is_template != -1 :
        style = self.create_style(object.style)
        if (style.selection) and (object.type.phil_type == "str") :
          print "WARNING: deprecated 'str' type with 'selection' style"
          print "   name: %s" % full_object_path
        self.style[full_object_path] = style
        if style.hidden :
          self._hidden.append(full_object_path)
        if style.OnUpdate is not None :
          print "OnUpdate is deprecated (%s)" % full_object_path
          self._update_handlers[full_object_path] = style.OnUpdate
        elif style.process_hkl :
          self._event_handlers[full_object_path] = "auto_extract_hkl_params"
        if style.OnChange is not None :
          self._event_handlers[full_object_path] = style.OnChange
        if style.renderer is not None :
          self._renderers[full_object_path] = style.renderer
        if style.menu_item :
          if phil_scope.multiple and phil_scope.is_template == 0 :
            pass
          elif style.parent_submenu :
            current_menu.add_submenu(style.parent_submenu)
            current_menu.get_submenu(style.parent_submenu).add_menu_item(
              full_object_path)
          else :
            current_menu.add_menu_item(full_object_path)
        elif style.submenu :
          if phil_scope.multiple and phil_scope.is_template == 0 :
            pass
          elif style.parent_submenu :
            current_menu.add_submenu(style.parent_submenu)
            parent_submenu = current_menu.get_submenu(style.parent_submenu)
            parent_submenu.add_submenu(full_object_path)
            next_menu = parent_submenu.get_submenu(full_object_path)
          else :
            current_menu.add_submenu(full_object_path)
            next_menu = current_menu.get_submenu(full_object_path)
      else :
        self.style[full_object_path] = gui_objects.style()
      if not object.is_definition :
        self.generate_gui_components(object, use_submenu, next_menu)
      use_submenu = False

  def get_scope_style (self, scope_name=None) :
    if scope_name in self.style :
      return self.style[scope_name]
    else :
      return gui_objects.style()

  def get_menu_db (self) :
    return self._menu_tree

  def get_file_type_map (self, file_type, default_label=None,
      exclude_params=()) :
    if (file_type in self._file_type_mappings) :
      return self._file_type_mappings[file_type]
    param_info = []
    for path_name, def_style in self.style.iteritems() :
      def_types = []
      if (def_style.file_type is not None) :
        def_types = def_style.get_list("file_type")
        if (file_type in def_types) :
          if ((def_style.no_map) or (def_style.new_file) or
              (path_name in exclude_params)) :
            continue
          phil_object = self.get_scope_by_name(path_name)
          if isinstance(phil_object, list) :
            phil_object = phil_object[0]
          label = get_standard_phil_label(phil_object)
          parent_scope = phil_object.primary_parent_scope
          if ((phil_object.multiple) or (parent_scope.multiple) or
              (phil_object.type.phil_type=="strings")) :
            count = None
          else :
            count = 1
          if def_style.file_type_default :
            default_label = label
          param_info.append((phil_object.full_path(), label, count))
    type_map = gui_objects.file_type_map(param_info, default_label)
    self._file_type_mappings[file_type] = type_map
    return type_map

########################################################################
#--- STANDALONE FUNCTIONS
def delete_phil_objects (current_phil, phil_path_list, only_scope=None) :
  assert isinstance(phil_path_list, list)
  i = 0
  while i < len(current_phil.objects) :
    full_path = current_phil.objects[i].full_path()
    if (only_scope is not None) :
      if not ((only_scope == full_path) or
              (only_scope.startswith(full_path + ".")) or
              (full_path.startswith(only_scope + "."))) :
        i += 1
        continue
    if current_phil.objects[i].is_template != 0 :
      i += 1
    elif full_path in phil_path_list :
      del current_phil.objects[i]
    else :
      # XXX: is this always true?
      if hasattr(current_phil.objects[i], "objects") :
        for path_name in phil_path_list :
          if path_name.startswith(full_path) :
            delete_phil_objects(current_phil=current_phil.objects[i],
              phil_path_list=phil_path_list,
              only_scope=only_scope)
      i += 1

def find_and_replace_scope (current_phil, new_scope, scope_name) :
  i = 0
  while (i < len(current_phil.objects)) :
    full_path = current_phil.objects[i].full_path()
    if (full_path == scope_name) :
      #assert (not current_phil.objects[i].multiple)
      new_scope.change_primary_parent_scope(current_phil)
      j = i
      while (j < len(current_phil.objects)) :
        if (current_phil.objects[j].full_path() == scope_name) :
          del current_phil.objects[j]
        else :
          j += 1
      current_phil.objects[i:i] = new_scope.objects
      break
    elif (scope_name.startswith(full_path + ".")) :
      find_and_replace_scope(
        current_phil=current_phil.objects[i],
        new_scope=new_scope,
        scope_name=scope_name)
    i += 1

def collect_redundant_paths (master_phil, new_phil, multiple_only=True) :
  phil_diff = master_phil.fetch_diff(source=new_phil)
  return _collect_unique_paths(phil_diff, multiple_only)

def _collect_unique_paths (phil_object, multiple_only=True) :
  paths = []
  if phil_object.multiple :
    paths.append(phil_object.full_path())
  elif phil_object.is_scope :
    for object in phil_object.objects :
      if not object.full_path() in paths :
        paths.extend(_collect_unique_paths(object, multiple_only))
  elif not multiple_only and phil_object.is_definition :
    paths.append(phil_object.full_path())
  return paths

def get_all_path_names (phil_object, paths=None) :
  if paths is None :
    paths = []
  full_path = phil_object.full_path()
  if not full_path in paths :
    paths.append(full_path)
  if phil_object.is_scope :
    for object in phil_object.objects :
      get_all_path_names(object, paths)

def index_phil_objects (phil_object,
                        path_index,
                        text_index,
                        template_index,
                        multiple_scopes=None,
                        multiple_defs=None,
                        collect_multiple=True,
                        in_template=False,
                        expert_levels=None,
                        input_files=None) :
  full_path = phil_object.full_path()
  if expert_levels is not None :
    if phil_object.expert_level is not None :
      expert_levels[full_path] = phil_object.expert_level
    else :
      parent_scope = ".".join(full_path.split(".")[:-1])
      expert_levels[full_path] = expert_levels.get(parent_scope, 0)
  if phil_object.is_template != 0 :
    template_index[full_path] = phil_object
    if phil_object.is_template == -1 :
      return
    else :
      in_template = True
  elif in_template :
    template_index[full_path] = phil_object
  if (phil_object.multiple == True) :
    if collect_multiple :
      if (phil_object.is_scope) and (multiple_scopes is not None) :
        multiple_scopes[full_path] = True
      elif multiple_defs is not None :
        multiple_defs[full_path] = True
    if full_path in path_index :
      path_index[full_path].append(phil_object)
    else :
      path_index[full_path] = [phil_object]
  else :
    path_index[full_path] = phil_object
  if phil_object.is_definition and phil_object.type is None :
    raise RuntimeError("Type required for parameter '%s'." % full_path)
  label = get_standard_phil_label(phil_object)
  text_fields = (label, str(phil_object.caption), str(phil_object.help),
    phil_object.is_definition)
  text_index[full_path] = text_fields
  if phil_object.is_scope :
    for object in phil_object.objects :
      index_phil_objects(phil_object=object,
                         path_index=path_index,
                         text_index=text_index,
                         template_index=template_index,
                         multiple_scopes=multiple_scopes,
                         multiple_defs=multiple_defs,
                         collect_multiple=collect_multiple,
                         in_template=in_template,
                         expert_levels=expert_levels,
                         input_files=input_files)
  elif (input_files is not None) :
    if (phil_object.type.phil_type in ["path", "strings"]) :
      style = phil_object.style
      if (style is not None) :
        style_words = style.split()
        if ("input_file" in style_words) :
          input_files.append(full_path)

def reindex_phil_objects (phil_object, path_index, only_scope=None) :
  if phil_object.is_template < 0 :
    return
  full_path = phil_object.full_path()
  if phil_object.multiple == True :
    if full_path in path_index :
      path_index[full_path].append(phil_object)
    else :
      path_index[full_path] = [phil_object]
  else :
    path_index[full_path] = phil_object
  if phil_object.is_scope :
    for object in phil_object.objects :
      reindex_phil_objects(object, path_index)

non_alnum = re.compile("[^A-Za-z0-9_]")
def substitute_directory_name (phil_object, path_name, sub_name,
    treat_name_as_var_name=True) :
  assert (not non_alnum.search(sub_name))
  if (treat_name_as_var_name) :
    sub_var = "$(" + sub_name + ")"
  else :
    sub_var = sub_name
  if path_name.endswith("/") :
    path_name = path_name[:-1]
  for object in phil_object.objects :
    if object.is_definition :
      if (object.type is None) :
        raise RuntimeError("Missing type for PHIL parameter %s" %
          object.full_path())
      if (object.type.phil_type == "path") :
        py_object = object.extract()
        if (py_object is None) or (py_object is Auto) : continue
        assert isinstance(py_object, str)
        py_object = py_object.replace(path_name, sub_var)
        new_object = object.format(python_object=py_object)
        object.words = new_object.words
    else :
      substitute_directory_name(object, path_name, sub_name)

def update_phil_file_paths (master_phil, file_name, old_path, new_path,
    use_iotbx_parser=False) :
  if (use_iotbx_parser) :
    import iotbx.phil
    parse = iotbx.phil.parse
  else :
    parse = libtbx.phil.parse
  phil_in = open(file_name).read()
  new_format = False
  out_lines = []
  for line in phil_in.splitlines() :
    if line.startswith("LIBTBX_BASE_DIR") :
      line = line.replace(old_path, new_path)
      new_format = True
      out_lines.append(line)
    else :
      out_lines.append(line)
  if (new_format) :
    open(file_name, "w").write("\n".join(out_lines))
  else :
    file_phil = parse(file_name=file_name)
    working_phil = master_phil.fetch(source=file_phil)
    substitute_directory_name(
      phil_object=working_phil,
      path_name=old_path,
      sub_name="LIBTBX_BASE_DIR")
    f = open(file_name, "w")
    f.write("LIBTBX_BASE_DIR = \"%s\"\n" % new_path)
    working_phil.show(out=f)
    f.close()

def get_standard_phil_label (phil_object=None, phil_name=None, append="") :
  if phil_object is None and phil_name is None :
    raise Exception("No phil object or path name supplied.")
  if phil_object is not None :
    if phil_object.short_caption is None :
      if phil_name is not None :
        phil_object.short_caption = reformat_phil_name(phil_name)
      else :
        phil_object.short_caption = reformat_phil_name(phil_object.name)
    return phil_object.short_caption + append
  else :
    return reformat_phil_name(phil_name)

def reformat_phil_full_name (phil_full_name) :
  phil_name = phil_full_name.split(".")[-1]
  return reformat_phil_name(phil_name)

def reformat_phil_name (phil_name) :
  if phil_name == "" :
    return phil_name
  _name = " ".join(str(phil_name).split("_"))
  name = string.upper(_name[0]) + _name[1:]
  return name

def join_scope_paths (scope1, scope2) :
  if scope1 == "" :
    return scope2
  else :
    return "%s.%s" % (scope1, scope2)

def get_adjoining_phil_path (def_path, def_name) :
  return ".".join(def_path.split(".")[:-1]) + "." + def_name

#---end
