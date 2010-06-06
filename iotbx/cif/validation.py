from iotbx.cif import builders, model, errors
import sys


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

  def show(self, out=None):
    if out is None:
      out = sys.stdout
    for code, errors in self.errors.iteritems():
      printed_messages = set()
      for e in errors:
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


class dictionary(model.cif):

  def __init__(self, other):
    model.cif.__init__(self, other.blocks)
    if self.has_key('on_this_dictionary'):
      self.DDL_version = 1
      on_this_dict = self['on_this_dictionary']
      self.name = on_this_dict['_dictionary_name']
      self.version = on_this_dict['_dictionary_version']
    else:
      raise NotImplementedError("Only DDL1 is currently supported")
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

  def find_data_block(self, key):
    """Returns the name of the data block containing the definition for the
       given key.  Raises a KeyError if item not found."""

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
      return key_
    # otherwise we have to check every block in turn
    else:
      for k, v in self.iteritems():
        if k == 'on_this_dictionary': continue
        elif key in v['_name']:
          return k
      self.report_error(1001, key=key) # item not in dictionary
      raise KeyError

  def validate_single_item(self, key, value, block):
    try:
      data = self[self.find_data_block(key)]
    except KeyError: return
    self.validate_type(key, value, data)
    self.validate_enumeration(key, value, data)
    self.validate_related(key, block, data)
    _list = data.get("_list")
    if _list == 'yes':
      self.report_error(2506, key=key) # must be in looped list

  def validate_type(self, key, value, data):
    if value in ('?', '.'): return
    data_type = data['_type']
    if data_type == 'numb':
      try:
        builders.float_from_string(value)
      except Exception, e:
        # can't interpret as numb
        self.report_error(2001, key=key, value=value, data_type=data_type)
      else:
        # check any type conditions
        type_condition = data.get('_type_conditions')
        if type_condition not in ('esd', 'su'):
          try:
            float(value)
          except Exception, e:
            # if we have got here, then from the data type checking we can assume
            # that the value is given with an esd, which causes it to be invalid.
            self.report_error(2002, key=key)

  def validate_enumeration(self, key, value, data):
    if isinstance(value, basestring):
      values = [value]
    else:
      values = value
    enum_values = data.get("_enumeration")
    enum_range = data.get("_enumeration_range")
    if enum_values is None and enum_range is None: return # nothing to check
    for value in values:
      if value in ('.', '?'): continue
      if enum_values is not None and value not in enum_values:
        # invalid choice for enumeration
        self.report_error(2102, key=key, value=value, enum=tuple(enum_values))
      if enum_range is not None:
        type_is_numb = data['_type'] == 'numb'
        enum_min, enum_max = enum_range.split(':')
        try:
          if enum_min != '': enum_min = builders.float_from_string(enum_min)
          if enum_max != '': enum_max = builders.float_from_string(enum_max)
        except ValueError:
          pass
        if data['_type'] == 'numb':
          try:
            value = builders.float_from_string(value)
          except: return # this error is handled with elsewhere
        if ((enum_min != '' and value < enum_min) or
            (enum_max != '' and value > enum_max)):
          # value outside range
          self.report_error(2101, key=key, value=value, enum=enum_range)

  def validate_related(self, key, block, data=None):
    related_item = data.get("_related_item")
    related_function = data.get("_related_function")
    if related_item is not None and related_function is not None:
      if related_function == 'replace' and block.get(related_item) is not None:
        self.report_error(2201, key=key, related_item=related_item)

  def validate_loop(self, loop, block):
    list_category = None
    for key, value in loop.iteritems():
      try:
        data = self[self.find_data_block(key)]
      except KeyError: continue
      self.validate_enumeration(key, value, data)
      _list = data.get("_list")
      if _list in ('no', None):
        self.report_error(2501, key=key) # not allowed in list
      if list_category is None:
        list_category = data.get('_category')
      elif list_category != data.get('_category'):
        self.report_error(2502) # multiple categories in loop
      mandatory = data.get('_list_mandatory') == 'yes'
      if mandatory: list_reference = key
      reference = data.get('_list_reference')
      if reference is not None:
        ref_data = self[self.find_data_block(reference)]
        ref_names = ref_data['_name']
        if isinstance(ref_names, basestring):
          ref_names = [ref_names]
        for name in ref_names:
          if name not in loop:
            self.report_error(2505, key=key, reference=name) # missing _list_reference
      #
      link_parent = data.get('_list_link_parent')
      if link_parent is not None:
        parent_values = block.get(link_parent)
        if parent_values is not None:
          for v in loop[key]:
            if v not in parent_values:
              # missing parent value
              self.report_error(2503, value=v, child=key, parent=link_parent)
        else:
          self.report_error(2504, child=key, parent=link_parent) # missing parent
