from iotbx.cif import builders, model, errors
import libtbx.load_env
import os
import re
import sys
from urllib2 import urlopen


class ErrorHandler:
  """An error handler for the validator. This class can be subclassed by clients
  that want to use their own error handlers"""

  def __init__(self):
    self.reset()

  def warning(self, w):
    self._add_warning(w)

  def error(self, e):
    self._add_error(e)

  def reset(self):
    self.warning_count = 0
    self.error_count = 0
    self.errors = {}
    self.warnings = {}

  def _add_warning(self, w):
    self.warning_count += 1
    if w.code in self.warnings:
      self.warnings[w.code].append(w)
    else:
      self.warnings.setdefault(w.code, [w])

  def _add_error(self, e):
    self.error_count += 1
    if e.code in self.errors:
      self.errors[e.code].append(e)
    else:
      self.errors.setdefault(e.code, [e])

  def show(self, show_warnings=True, out=None):
    if out is None:
      out = sys.stdout
    codes = self.errors.keys()
    errors = self.errors.values()
    if show_warnings:
      codes.extend(self.warnings.keys())
      errors.extend(self.warnings.values())
    for code, errs in zip(codes, errors):
      printed_messages = set()
      for e in errs:
        if str(e) not in printed_messages: # avoid printing duplicates
          printed_messages.add(str(e))
          print >> out, e


class ValidationError(Exception):
  def __init__(self, code, format_string, **kwds):
    self.code = code
    self.format_string = format_string
    self.kwds = kwds

  def __str__(self):
    return self.format_string %self.kwds



cifdic_register_url = "ftp://ftp.iucr.org/pub/cifdics/cifdic.register"

def smart_load_dictionary(name=None, file_path=None, url=None,
                          registry_location=cifdic_register_url,
                          save_local=False, store_dir=None):
  from iotbx import cif
  assert [name, file_path, url].count(None) < 3
  cif_dic = None
  if store_dir is None:
    store_dir = libtbx.env.under_dist(module_name='iotbx', path='cif')
  if name is not None and [file_path, url].count(None) == 2:
    if file_path is None:
      file_path = os.path.join(store_dir, name)
    if url is None:
      url = locate_dictionary(name, registry_location=registry_location)
  if file_path is not None and os.path.isfile(file_path):
    cif_dic = dictionary(cif.fast_reader(file_path=file_path).model())
  elif url is not None:
    file_object = urlopen(url)
    if save_local:
      if name is None:
        name = os.path.split(url)[-1]
      f = open(os.path.join(store_dir, name), 'wb')
      f.write(file_object.read())
      f.close()
      cif_dic = dictionary(cif.fast_reader(
        file_path=os.path.join(store_dir, name)).model())
    else:
      cif_dic = dictionary(cif.fast_reader(
        file_object=file_object).model())
  assert cif_dic is not None
  return cif_dic

def locate_dictionary(name, version=None, registry_location=cifdic_register_url):
  from iotbx import cif
  cm = cif.fast_reader(file_object=urlopen(registry_location)).model()
  if version is None: version = '.'
  reg = cm["validation_dictionaries"]
  for n, v, url in zip(reg['_cifdic_dictionary.name'],
                       reg['_cifdic_dictionary.version'],
                       reg['_cifdic_dictionary.URL']):
    if n == name and v == str(version):
      return url


class dictionary(model.cif):

  def __init__(self, other):
    model.cif.__init__(self, other.blocks)
    self.item_type_list = {}
    self.child_parent_relations = {}
    self.look_up_table = {} # cached definitions for each data name
    if self.has_key('on_this_dictionary'):
      self.DDL_version = 1
      for key, value in self.blocks.iteritems():
        self[key] = DDL1_definition(value)
      on_this_dict = self['on_this_dictionary']
      self.name = on_this_dict['_dictionary_name']
      self.version = on_this_dict['_dictionary_version']
    else:
      self.DDL_version = 2
      master_block = self.values()[0]
      self.name = master_block['_dictionary.title']
      self.version = master_block['_dictionary.version']
      type_codes = master_block.get('_item_type_list.code')
      type_constructs = master_block.get('_item_type_list.construct')
      for code, construct in zip(type_codes, type_constructs):
        self.item_type_list.setdefault(code, re.compile(construct))
      for key, save in master_block.saves.iteritems():
        master_block[key] = DDL2_definition(save)
        children = save.get('_item_linked.child_name')
        parents = save.get('_item_linked.parent_name')
        if parents is not None and children is not None:
          if not isinstance(parents, basestring):
            for child, parent in zip(children, parents):
              self.child_parent_relations.setdefault(child, parent)
    self.err = ErrorHandler()
    language = "en"
    self.errors = errors.error_dicts[language]

  def set_error_handler(self, handler):
    self.err = handler

  def report_error(self, number, **kwds):
    message = self.errors[number]
    if number < 2000:
      self.err.warning(ValidationError(number, message, **kwds))
    elif number < 3000:
      self.err.error(ValidationError(number, message, **kwds))

  def find_definition(self, key):
    """Returns the name of the data block containing the definition for the
       given key.  Raises a KeyError if item not found."""
    if self.DDL_version == 1:
      if key in self.look_up_table:
        return self.look_up_table[key]
      key_ = key.lstrip("_")
      data = self.get(key_) # try simplest first
      # then try shortening key
      while data is None:
        new_key = key_[:key_.rfind("_", 0, -1)+1]
        if new_key == key_: break
        data = self.get(new_key)
        if data is not None and key not in data['_name']:
          data = None
        key_ = new_key
      if data is not None:
        self.look_up_table.setdefault(key, key_)
        return key_
      # otherwise we have to check every block in turn
      else:
        for k, v in self.iteritems():
          if k == 'on_this_dictionary': continue
          elif key in v['_name']:
            self.look_up_table.setdefault(key, key_)
            return k
        self.report_error(1001, key=key) # item not in dictionary
        raise KeyError, key
    else:
      if key not in self.values()[0]:
        self.report_error(1001, key=key) # item not in dictionary
        raise KeyError, key
      else:
        return key

  def get_definition(self, key):
    if self.DDL_version == 1:
      return self[self.find_definition(key)]
    elif self.DDL_version == 2:
      return self.values()[0][self.find_definition(key)]

  def validate_single_item(self, key, value, block):
    try:
      definition = self.get_definition(key)
    except KeyError: return
    self.validate_type(key, value, definition)
    self.validate_enumeration(key, value, definition)
    self.validate_related(key, block, definition)
    self.validate_dependent(key, block, definition)
    _list = definition.get("_list")
    if _list == 'yes':
      self.report_error(2506, key=key) # must be in looped list

  def validate_type(self, key, value, definition):
    if value in ('?', '.'): return
    item_type = definition.type
    if item_type in self.item_type_list:
      match = re.match(self.item_type_list[item_type], value)
      if match is None:
        self.report_error(2001, key=key, value=value, item_type=item_type)
    elif item_type == 'numb':
      # only for DDL1
      try:
        builders.float_from_string(value)
      except Exception, e:
        # can't interpret as numb
        self.report_error(2001, key=key, value=value, item_type=definition.type)
      else:
        # check any type conditions
        type_condition = definition.type_conditions
        if type_condition not in ('esd', 'su'):
          try:
            float(value)
          except Exception, e:
            # if we have got here, then from the data type checking we can assume
            # that the value is given with an esd, which causes it to be invalid.
            self.report_error(2002, key=key)

  def validate_dependent(self, key, block, definition):
    dependents = definition.dependent
    if dependents is None: return
    elif isinstance(dependents, basestring):
      dependents = [dependents]
    for dependent in dependents:
      if dependent not in block:
        self.report_error(2301, dependent=dependent, key=key)

  def validate_enumeration(self, key, value, definition):
    if isinstance(value, basestring):
      values = [value]
    else:
      values = value
    enum_values = definition.enumeration
    enum_min, enum_max = definition.get_min_max()
    if enum_values is None and enum_max is None and enum_min is None:
      return # nothing to check
    for value in values:
      if value in ('.', '?'): continue
      elif enum_values is not None and value not in enum_values:
        enum_lower = [v.lower() for v in enum_values]
        if value.lower() in enum_lower:
          self.report_error(1002, key=key, value=value) # case sensitive match failure
        else:
          # invalid choice for enumeration
          self.report_error(2102, key=key, value=value, enum=tuple(enum_values))
      if definition.type in ('numb', 'float', 'int'):
        try:
          v = builders.float_from_string(value)
        except: return # this error is handled with elsewhere
        # DDL1 range is inclusive, DDL2 is exclusive
        if self.DDL_version == 1:
          if not ((enum_min is None or v >= float(enum_min)) and
                  (enum_max is None or v <= float(enum_max))):
            self.report_error(
              2101, key=key, value=value, enum="%s:%s" %(enum_min, enum_max))
        else:
          for min, max in zip(enum_min, enum_max):
            if ((min == '.' or v > float(min)) and
                (max == '.' or v < float(max))):
              return # at least one condition was met, so value is inside range
            elif (min == max and v == float(min)):
              return # matched boundary value
          # else value out of range
          self.report_error(2101, key=key, value=value, enum="%s:%s" %(min, max))

  def validate_related(self, key, block, definition):
    related_items = definition.related
    related_functions = definition.related_function
    if related_items is not None and related_functions is not None:
      if (isinstance(related_items, basestring) and
          isinstance(related_functions, basestring)):
        related_items = [related_items]
        related_functions = [related_functions]
      for related_item, related_function in zip(related_items, related_functions):
        if related_function == 'replace':
          if block.get(related_item) is not None:
            self.report_error(2201, key=key, related_item=related_item)
          else: # obsolete definition warning
            self.report_error(1003, key=key, related_item=related_item)
        elif (related_function == 'alternate_exclusive' and
              related_item in block):
          self.report_error(2201, key=key, related_item=related_item)
        elif related_function == 'replacedby': # obsolete definition warning
          self.report_error(1003, key=key, related_item=related_item)
        elif (related_function == 'associated_value' and
              related_item not in block): # missing associated value
          self.report_error(2202, key=key, related_item=related_item)

  def validate_loop(self, loop, block):
    list_category = None
    for key, value in loop.iteritems():
      try:
        definition = self.get_definition(key)
      except KeyError: continue
      self.validate_enumeration(key, value, definition)
      self.validate_dependent(key, block, definition)
      self.validate_related(key, block, definition)
      _list = definition.get("_list")
      if self.DDL_version == 1 and _list in ('no', None):
        self.report_error(2501, key=key) # not allowed in list
      if list_category is None:
        list_category = definition.category
      elif (isinstance(list_category, basestring)
            and definition.category is not None
            and list_category != definition.category):
        self.report_error(2502) # multiple categories in loop
      mandatory = definition.mandatory == 'yes'
      reference = definition.get('_list_reference')
      if reference is not None:
        ref_data = self.get_definition(reference)
        ref_names = ref_data['_name']
        if isinstance(ref_names, basestring):
          ref_names = [ref_names]
        for name in ref_names:
          if name not in loop:
            self.report_error(2505, key=key, reference=name) # missing _list_reference
      elif (self.DDL_version == 2
            and isinstance(definition.category, basestring)):
        category_def = self.get_definition(definition.category)
        if category_def.category_key is not None:
          category_keys = category_def.category_key
          if isinstance(category_keys, basestring):
            category_keys = [category_keys]
          for cat_key in category_keys:
            cat_key_def = self.get_definition(cat_key)
          if (cat_key_def.mandatory == 'yes'
              and isinstance(cat_key_def.mandatory, basestring)
              and cat_key_def.name not in block):
            self.report_error(
              2203, key=cat_key_def.name, category=definition.category)
      #
      link_parent = definition.get(
        '_list_link_parent', self.child_parent_relations.get(key))
      if link_parent is not None:
        parent_values = block.get(link_parent)
        if parent_values is not None:
          for v in loop[key]:
            if v not in parent_values:
              # missing parent value
              self.report_error(2503, value=v, child=key, parent=link_parent)
        else:
          self.report_error(2504, child=key, parent=link_parent) # missing parent

class definition_base:

  def name(self):
    return self.get(self.aliases['name'])

  def type(self):
    return self.get(self.aliases['type'])

  def type_conditions(self):
    return self.get(self.aliases['type_conditions'])

  def category(self):
    return self.get(self.aliases['category'])

  def category_key(self):
    return self.get(self.aliases['category_key'])

  def mandatory(self):
    return self.get(self.aliases['mandatory'])

  def enumeration(self):
    return self.get(self.aliases['enumeration'])

  def dependent(self):
    return self.get(self.aliases['dependent'])

  def related(self):
    return self.get(self.aliases['related'])

  def related_function(self):
    return self.get(self.aliases['related_function'])

  name = property(name)
  type = property(type)
  type_conditions = property(type_conditions)
  category = property(category)
  category_key = property(category_key)
  mandatory = property(mandatory)
  enumeration = property(enumeration)
  dependent = property(dependent)
  related = property(related)
  related_function = property(related_function)

class DDL1_definition(model.block, definition_base):

  aliases = {
  'name': '_name',
  'type': '_type',
  'type_conditions': '_type_conditions',
  'category': '_category',
  'category_key': 'XXX',
  'mandatory': '_list_mandatory',
  'enumeration': '_enumeration',
  'related': '_related_item',
  'related_function': '_related_function',
  }

  def __init__(self, other):
    self._items = other._items
    self.loops = other.loops
    self.saves = other.saves
    self._set = other._set
    self.keys_lower = other.keys_lower

  def dependent(self):
    return None
  dependent = property(dependent)

  def get_min_max(self):
    enum_range = self.get('_enumeration_range')
    if enum_range is not None:
      enum_min, enum_max = enum_range.split(':')
      if enum_min == '': enum_min = None
      if enum_max == '': enum_max = None
      return (enum_min, enum_max)
    else:
      return (None, None)

class DDL2_definition(model.save, definition_base):

  aliases = {
  'name': '_item.name',
  'type': '_item_type.code',
  'type_conditions': '_item_type_conditions.code',
  'category': '_item.category_id',
  'category_key': '_category_key.name',
  'mandatory': '_item.mandatory_code',
  'enumeration': '_item_enumeration.value',
  'dependent': '_item_dependent.dependent_name',
  'related': '_item_related.related_name',
  'related_function': '_item_related.function_code',
  }

  def __init__(self, other):
    self._items = other._items
    self.loops = other.loops
    self._set = other._set
    self.keys_lower = other.keys_lower

  def get_min_max(self):
    return (self.get('_item_range.minimum'), self.get('_item_range.maximum'))
