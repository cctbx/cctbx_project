from __future__ import division
from __future__ import generators
from iotbx import simple_tokenizer
from scitbx.python_utils.str_utils import line_breaker
from libtbx.itertbx import count
from libtbx import introspection
from cStringIO import StringIO
import copy
import math
import sys, os

default_print_width = 79

standard_identifier_start_characters = {}
for c in "_ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz":
  standard_identifier_start_characters[c] = None
standard_identifier_continuation_characters = dict(
  standard_identifier_start_characters)
for c in ".0123456789":
  standard_identifier_continuation_characters[c] = None

def is_standard_identifier(string):
  if (len(string) == 0): return False
  if (string[0] not in standard_identifier_start_characters): return False
  for c in string[1:]:
    if (c not in standard_identifier_continuation_characters): return False
  sub_strings = string.split(".")
  if (len(sub_strings) > 1):
    for sub in sub_strings:
      if (not is_standard_identifier(sub)): return False
  return True

def str_from_assigned_words(assigned_words):
  if (len(assigned_words) == 1 and assigned_words[0].value.lower() == "none"):
    return None
  return " ".join([word.value for word in assigned_words])

def bool_from_assigned_words(assigned_words):
  value_string = str_from_assigned_words(assigned_words)
  if (value_string is None): return None
  word_lower = assigned_words[0].value.lower()
  if (word_lower == "none"): return None
  if (word_lower in ["false", "no", "off", "0"]): return False
  if (word_lower in ["true", "yes", "on", "1"]): return True
  assert len(assigned_words) > 0
  raise RuntimeError(
    'One True of False value expected, "%s" found%s' % (
      value_string, assigned_words[0].where_str()))

def number_from_assigned_words(assigned_words):
  value_string = str_from_assigned_words(assigned_words)
  if (value_string is None): return None
  try: return eval(value_string, math.__dict__, {})
  except Exception, e:
    raise RuntimeError(
      'Error interpreting "%s" as a numeric expression: %s%s' % (
        value_string, str(e), assigned_words[0].where_str()))

def int_from_assigned_words(assigned_words):
  result = number_from_assigned_words(assigned_words)
  if (result is not None):
    if (isinstance(result, float)
        and round(result) == result):
      result = int(result)
    elif (not isinstance(result, int)):
      raise RuntimeError(
        'Integer expression expected, "%s" found%s' % (
          str_from_assigned_words(assigned_words),
          assigned_words[0].where_str()))
  return result

def float_from_assigned_words(assigned_words):
  result = number_from_assigned_words(assigned_words)
  if (result is not None):
    if (isinstance(result, int)):
      result = float(result)
    elif (not isinstance(result, float)):
      raise RuntimeError(
        'Floating-point expression expected, "%s" found%s' % (
          str_from_assigned_words(assigned_words),
          assigned_words[0].where_str()))
  return result

def choice_from_assigned_words(optional, assigned_words):
  result = None
  for word in assigned_words:
    if (word.value.startswith("*")):
      if (result is not None):
        raise RuntimeError("Multiple choices where only one is possible%s" %
          assigned_words[0].where_str())
      result = word.value[1:]
  if (result is None and not optional):
    raise RuntimeError("Unspecified choice%s" % assigned_words[0].where_str())
  return result

def multi_choice_from_assigned_words(assigned_words):
  result = []
  for word in assigned_words:
    if (word.value.startswith("*")):
      result.append(word.value[1:])
  return result

def unit_cell_from_assigned_words(assigned_words):
  from cctbx import uctbx
  return uctbx.unit_cell(str_from_assigned_words(assigned_words))

def space_group_info_from_assigned_words(assigned_words):
  from cctbx import sgtbx
  return sgtbx.space_group_info(symbol=str_from_assigned_words(assigned_words))

default_definition_type_names = [
  "str", "bool", "int", "float",
  "choice", "multi_choice",
  "path", "key",
  "unit_cell", "space_group"]

def definition_type_from_assigned_words(assigned_words, type_names):
  if (len(assigned_words) == 1):
    word_lower = assigned_words[0].value.lower()
    if (word_lower == "none"): return None
    if (word_lower in type_names): return word_lower
  assert len(assigned_words) > 0
  raise RuntimeError(
    'Unexpected definition type: "%s"%s' % (
      assigned_words[0].value, assigned_words[0].where_str()))

def show_attributes(self, out, prefix, attributes_level, print_width):
  if (attributes_level <= 0): return
  for name in self.attribute_names:
    value = getattr(self, name)
    if ((name == "help" and value is not None)
        or (value is not None and attributes_level > 1)
        or attributes_level > 2):
      if (not isinstance(value, str)):
        print >> out, prefix+"  ."+name, "=", value
      else:
        indent = " " * (len(prefix) + 3 + len(name) + 3)
        fits_on_one_line = len(indent+value) < print_width
        if (not is_standard_identifier(value) or not fits_on_one_line):
          value = str(simple_tokenizer.word(value=value, quote_token='"'))
          fits_on_one_line = len(indent+value) < print_width
        if (fits_on_one_line):
          print >> out, prefix+"  ."+name, "=", value
        else:
          is_first = True
          for block in line_breaker(value[1:-1], print_width-2-len(indent)):
            if (is_first):
              print >> out, prefix+"  ."+name, "=", '"'+block+'"'
              is_first = False
            else:
              print >> out, indent+'"'+block+'"'

class object_locator:

  def __init__(self, parent, path, object):
    introspection.adopt_init_args()

class definition: # FUTURE definition(object)

  attribute_names = [
    "help", "caption", "short_caption", "optional",
    "type", "multiple", "input_size", "expert_level"]

  __slots__ = ["name", "words", "primary_id", "primary_parent_scope",
               "is_disabled", "where_str", "merge_names"] + attribute_names

  def __init__(self,
        name,
        words,
        primary_id=None,
        primary_parent_scope=None,
        is_disabled=False,
        where_str="",
        merge_names=False,
        help=None,
        caption=None,
        short_caption=None,
        optional=None,
        type=None,
        multiple=None,
        input_size=None,
        expert_level=None):
    if (name != "include" and "include" in name.split(".")):
      raise RuntimeError('Reserved identifier: "include"%s' % where_str)
    self.name = name
    self.words = words
    self.primary_id = primary_id
    self.primary_parent_scope = primary_parent_scope
    self.is_disabled = is_disabled
    self.where_str = where_str
    self.merge_names = merge_names
    self.help = help
    self.caption = caption
    self.short_caption = short_caption
    self.optional = optional
    self.type = type
    self.multiple = multiple
    self.input_size = input_size
    self.expert_level = expert_level

  def copy(self, name=None, words=None):
    keyword_args = {}
    for keyword in self.__slots__:
      keyword_args[keyword] = getattr(self, keyword)
    if (name is not None):
      keyword_args["name"] = name
    if (words is not None):
      keyword_args["words"] = words
    return definition(**keyword_args)

  def fetch(self, source):
    if (not isinstance(source, definition)):
      raise RuntimeError('Incompatible parameter objects "%s"%s and "%s"%s' %
        (self.name, self.where_str, source.name, source.where_str))
    source = source.resolve_variables()
    if (self.type not in ["choice", "multi_choice"]):
      return self.copy(words=source.words)
    flags = {}
    for word in self.words:
      if (word.value.startswith("*")): value = word.value[1:]
      else: value = word.value
      flags[value] = False
    for word in source.words:
      if (word.value.startswith("*")):
        value = word.value[1:]
        flag = True
      else:
        value = word.value
        if (len(source.words) == 1):
          flag = True
        else:
          flag = False
      if (flag and value not in flags):
        raise RuntimeError("Not a possible choice: %s%s" % (
          str(word), word.where_str()))
      flags[value] = flag
    words = []
    for word in self.words:
      if (word.value.startswith("*")): value = word.value[1:]
      else: value = word.value
      if (flags[value]): value = "*" + value
      words.append(simple_tokenizer.word(
        value=value,
        line_number=word.line_number,
        source_info=word.source_info))
    return self.copy(words=words)

  def has_attribute_with_name(self, name):
    return name in self.attribute_names

  def assign_attribute(self, name, assigned_words, type_names):
    assert self.has_attribute_with_name(name)
    if (name in ["optional", "multiple"]):
      value = bool_from_assigned_words(assigned_words)
    elif (name == "type"):
      value = definition_type_from_assigned_words(assigned_words, type_names)
    elif (name in ["input_size", "expert_level"]):
      value = int_from_assigned_words(assigned_words)
    else:
      value = str_from_assigned_words(assigned_words)
    setattr(self, name, value)

  def show(self,
        out=None,
        merged_names=[],
        prefix="",
        attributes_level=0,
        print_width=None):
    if (out is None): out = sys.stdout
    if (print_width is None): print_width = default_print_width
    if (self.is_disabled): hash = "#"
    else:                  hash = ""
    line = prefix + hash + ".".join(merged_names + [self.name])
    if (self.name != "include"): line += " ="
    indent = " " * len(line)
    for word in self.words:
      line_plus = line + " " + str(word)
      if (len(line_plus) > print_width-2 and len(line) > len(indent)):
        print >> out, line + " \\"
        line = indent + " " + str(word)
      else:
        line = line_plus
    print >> out, line
    show_attributes(
      self=self,
      out=out,
      prefix=prefix,
      attributes_level=attributes_level,
      print_width=print_width)

  def as_str(self,
        prefix="",
        attributes_level=0,
        print_width=None):
    out = StringIO()
    self.show(
      out=out,
      prefix=prefix,
      attributes_level=attributes_level,
      print_width=print_width)
    return out.getvalue()

  def has_same_definitions(self, other):
    return self.as_str() == other.as_str()

  def _all_definitions(self, parent, parent_path, result):
    result.append(object_locator(
      parent=parent, path=parent_path+self.name, object=self))

  def get_without_substitution(self, path):
    if (self.is_disabled or self.name != path): return []
    return [self]

  def automatic_type(self):
    types = {}
    for word in self.words:
      if (word.quote_token is not None):
        types["str"] = None
        continue
      word_lower = word.value.lower()
      if (word_lower in ["false", "no", "off"
                         "true", "yes", "on"]):
        types["bool"] = None
        continue
      try: py_value = eval(word.value, {}, {})
      except:
        if (word.value[0] in standard_identifier_start_characters):
          types["str"] = None
        else:
          types["unknown"] = None
        continue
      if (isinstance(py_value, float)):
        types["float"] = None
        continue
      if (isinstance(py_value, int)):
        types["int"] = None
        continue
    types = types.keys()
    types.sort()
    if (types == ["int", "float"]): return "float"
    if (len(types) == 1): return types[0]
    return None

  def automatic_type_assignment(self, assignment_if_unknown=None):
    if (self.type is None):
      self.type = self.automatic_type()
      if (self.type is None):
        self.type = assignment_if_unknown

  def extract(self, custom_converters=None):
    if (self.type in ["str", "path", "key"]):
      return str_from_assigned_words(self.words)
    if (self.type == "bool"):
      return bool_from_assigned_words(self.words)
    if (self.type == "int"):
      return int_from_assigned_words(self.words)
    if (self.type == "float"):
      return float_from_assigned_words(self.words)
    if (self.type == "choice"):
      return choice_from_assigned_words(self.optional, self.words)
    if (self.type == "multi_choice"):
      return multi_choice_from_assigned_words(self.words)
    if (self.type == "unit_cell"):
      return unit_cell_from_assigned_words(self.words)
    if (self.type == "space_group"):
      return space_group_info_from_assigned_words(self.words)
    if (custom_converters is not None):
      converter = custom_converters.get(self.type, None)
      if (converter is not None):
        return converter.process_assigned_words(self.words)
    if (self.type is None):
      return [word.value for word in self.words]
    raise RuntimeError(
       ('No converter for parameter definition type "%s"'
      + ' required for converting words assigned to "%s"%s') % (
        self.type, self.name, self.where_str))

  def format(self, python_object, custom_converters=None):
    words = None
    if (python_object is None):
      words = [simple_tokenizer.word(value="None")]
    elif (self.type in ["str", "path", "key"]):
      words = [simple_tokenizer.word(value=python_object, quote_token='"')]
    elif (self.type == "bool"):
      if (python_object):
        words = [simple_tokenizer.word(value="True")]
      else:
        words = [simple_tokenizer.word(value="False")]
    elif (self.type == "int"):
      words = [simple_tokenizer.word(value=str(python_object))]
    elif (self.type == "float"):
      words = [simple_tokenizer.word(value="%.10g" % python_object)]
    elif (self.type in ["choice", "multi_choice"]):
      words = []
      for word in self.words:
        if (word.value.startswith("*")): value = word.value[1:]
        else: value = word.value
        if (self.type == "choice"):
          if (value == python_object):
            value = "*" + value
        else:
          if (value in python_object):
            value = "*" + value
        words.append(simple_tokenizer.word(
          value=value, quote_token=word.quote_token))
    elif (self.type == "unit_cell"):
      words = [simple_tokenizer.word(value="%.10g" % v)
        for v in python_object.parameters()]
    elif (self.type == "space_group"):
      words = [simple_tokenizer.word(value=str(python_object),quote_token='"')]
    elif (custom_converters is not None):
      converter = custom_converters.get(self.type, None)
      if (converter is not None):
        words = converter.format_as_assigned_words(self, python_object)
    elif (self.type is None):
      words = [simple_tokenizer.word(value=value, quote_token='"')
        for value in python_object]
    if (words is None):
      raise RuntimeError(
         ('No converter for parameter definition type "%s"'
        + ' required for converting values for "%s"%s') % (
          self.type, self.name, self.where_str))
    return self.copy(words=words)

  def unique(self):
    return self

  def resolve_variables(self):
    new_words = []
    for word in self.words:
      if (word.quote_token == "'"):
        new_words.append(word)
        continue
      substitution_proxy = variable_substitution_proxy(word)
      for fragment in substitution_proxy.fragments:
        if (not fragment.is_variable):
          fragment.result = simple_tokenizer.word(
            value=fragment.value, quote_token='"')
          continue
        variable_words = None
        if (self.primary_parent_scope is not None):
          substitution_source = self.primary_parent_scope.lexical_get(
            path=fragment.value, stop_id=self.primary_id)
          if (substitution_source is not None):
            assert isinstance(substitution_source, definition)
            variable_words = substitution_source.resolve_variables().words
        if (variable_words is None):
          env_var = os.environ.get(fragment.value, None)
          if (env_var is not None):
            variable_words = [simple_tokenizer.word(
              value=env_var,
              quote_token='"',
              source_info='environment: "%s"'%fragment.value)]
        if (variable_words is None):
          raise RuntimeError("Undefined variable: $%s%s" % (
            fragment.value, word.where_str()))
        if (not substitution_proxy.force_string):
          fragment.result = variable_words
        else:
          fragment.result = simple_tokenizer.word(
            value=" ".join([str(v) for v in variable_words]),
            quote_token='"')
      new_words.extend(substitution_proxy.get_new_words())
    return self.copy(words=new_words)

class scope_extract_attribute_error: pass
class scope_extract_is_disabled: pass

class scope_extract:

  def __set__(self, name, multiple, value):
    assert not "." in name
    node = getattr(self, name, scope_extract_attribute_error)
    if (not multiple):
      if (value is scope_extract_is_disabled):
        value = None
      if (node is scope_extract_attribute_error
          or not isinstance(value, scope_extract)
          or not isinstance(node, scope_extract)):
        setattr(self, name, value)
      else:
        node.__dict__.update(value.__dict__)
    else:
      if (node is scope_extract_attribute_error):
        if (value is scope_extract_is_disabled):
          setattr(self, name, [])
        else:
          setattr(self, name, [value])
      elif (not value is scope_extract_is_disabled):
        node.append(value)

  def __get__(self, name):
    assert not "." in name
    return getattr(self, name, scope_extract_attribute_error)

class scope:

  def __init__(self,
        name,
        objects=None,
        primary_id=None,
        primary_parent_scope=None,
        is_disabled=False,
        where_str="",
        merge_names=False,
        style=None,
        help=None,
        caption=None,
        short_caption=None,
        optional=None,
        multiple=None,
        sequential_format=None,
        disable_add=None,
        disable_delete=None,
        expert_level=None):
    introspection.adopt_init_args()
    self.attribute_names = self.__init__varnames__[8:]
    if (objects is None):
      self.objects = []
    assert style in [None, "row", "column", "block", "page"]
    if ("include" in name.split(".")):
      raise RuntimeError('Reserved identifier: "include"%s' % where_str)
    if (sequential_format is not None):
      assert isinstance(sequential_format % 0, str)

  def copy(self, name=None, objects=None):
    keyword_args = {}
    for keyword in self.__init__varnames__[1:]:
      keyword_args[keyword] = getattr(self, keyword)
    if (name is not None):
      keyword_args["name"] = name
    if (objects is not None):
      keyword_args["objects"] = objects
    return scope(**keyword_args)

  def adopt(self, object):
    assert len(object.name) > 0
    primary_parent_scope = self
    name_components = object.name.split(".")
    merge_names = False
    for name in name_components[:-1]:
      child_scope = scope(name=name)
      child_scope.merge_names = merge_names
      child_scope.optional = object.optional
      child_scope.multiple = object.multiple
      primary_parent_scope.adopt(child_scope)
      primary_parent_scope = child_scope
      merge_names = True
    if (len(name_components) > 1):
      object.name = name_components[-1]
      object.merge_names = True
    object.primary_parent_scope = primary_parent_scope
    primary_parent_scope.objects.append(object)

  def has_attribute_with_name(self, name):
    return name in self.attribute_names

  def assign_attribute(self, name, assigned_words):
    assert self.has_attribute_with_name(name)
    if (name in ["optional", "multiple", "disable_add", "disable_delete"]):
      value = bool_from_assigned_words(assigned_words)
    elif (name in ["expert_level"]):
      value = int_from_assigned_words(assigned_words)
    else:
      value = str_from_assigned_words(assigned_words)
      if (name == "style"):
        style = value
        assert style in [None, "row", "column", "block", "page"]
      elif (name == "sequential_format"):
        sequential_format = value
        if (sequential_format is not None):
          assert isinstance(sequential_format % 0, str)
    setattr(self, name, value)

  def active_objects(self):
    for object in self.objects:
      if (object.is_disabled): continue
      yield object

  def show(self,
        out=None,
        merged_names=[],
        prefix="",
        attributes_level=0,
        print_width=None):
    if (out is None): out = sys.stdout
    if (print_width is None): print_width = default_print_width
    is_proper_scope = False
    if (len(self.name) == 0):
      assert len(merged_names) == 0
    elif (len(self.objects) == 1 and self.objects[0].merge_names):
      merged_names = merged_names + [self.name]
    else:
      if (len(self.objects) > 1): assert not self.objects[0].merge_names
      is_proper_scope = True
      if (self.is_disabled): hash = "#"
      else:                  hash = ""
      out_attributes = StringIO()
      show_attributes(
        self=self,
        out=out_attributes,
        prefix=prefix,
        attributes_level=attributes_level,
        print_width=print_width)
      out_attributes = out_attributes.getvalue()
      merged_name = ".".join(merged_names + [self.name])
      merged_names = []
      if (len(out_attributes) == 0):
        print >> out, prefix + hash + merged_name, "{"
      else:
        print >> out, prefix + hash + merged_name
        out.write(out_attributes)
        print >> out, prefix+"{"
      prefix += "  "
    for object in self.objects:
      object.show(
        out=out,
        merged_names=merged_names,
        prefix=prefix,
        attributes_level=attributes_level,
        print_width=print_width)
    if (is_proper_scope):
      print >> out, prefix[:-2] + "}"

  def as_str(self,
        prefix="",
        attributes_level=0,
        print_width=None):
    out = StringIO()
    self.show(
      out=out,
      prefix=prefix,
      attributes_level=attributes_level,
      print_width=print_width)
    return out.getvalue()

  def has_same_definitions(self, other):
    return self.as_str() == other.as_str()

  def _all_definitions(self, parent, parent_path, result):
    parent_path += self.name+"."
    for object in self.active_objects():
      object._all_definitions(
        parent=self, parent_path=parent_path, result=result)

  def all_definitions(self):
    result = []
    for object in self.active_objects():
      object._all_definitions(parent=self, parent_path="", result=result)
    return result

  def automatic_type_assignment(self, assignment_if_unknown=None):
    for item in self.all_definitions():
      item.object.automatic_type_assignment(
        assignment_if_unknown=assignment_if_unknown)

  def get_without_substitution(self, path):
    if (self.is_disabled): return []
    if (len(self.name) == 0):
      if (len(path) == 0): return self.objects
    elif (self.name == path):
      return [self]
    elif (path.startswith(self.name+".")):
      path = path[len(self.name)+1:]
    else:
      return []
    result = []
    for object in self.active_objects():
      result.extend(object.get_without_substitution(path=path))
    return result

  def get(self, path, with_substitution=True):
    result = scope(name="", objects=self.get_without_substitution(path=path))
    if (not with_substitution): return result
    return result.resolve_variables()

  def resolve_variables(self):
    result = []
    for object in self.active_objects():
      result.append(object.resolve_variables())
    return self.copy(objects=result)

  def lexical_get(self, path, stop_id, search_up=True):
    if (path.startswith(".")):
      while (self.primary_parent_scope is not None):
        self = self.primary_parent_scope
      path = path[1:]
    candidates = []
    for object in self.objects:
      if (object.primary_id >= stop_id): break
      if (isinstance(object, definition)):
        if (object.name == path):
          candidates.append(object)
      elif (object.name == path
            or path.startswith(object.name+".")):
        candidates.append(object)
    while (len(candidates) > 0):
      object = candidates.pop()
      if (object.name == path): return object
      object = object.lexical_get(
        path=path[len(object.name)+1:], stop_id=stop_id, search_up=False)
      if (object is not None): return object
    if (not search_up): return None
    if (self.primary_parent_scope is None): return None
    return self.primary_parent_scope.lexical_get(path=path, stop_id=stop_id)

  def extract(self, custom_converters=None):
    result = scope_extract()
    for object in self.objects:
      if (object.is_disabled):
        value = scope_extract_is_disabled
      else:
        value = object.extract(custom_converters=custom_converters)
      result.__set__(name=object.name, multiple=object.multiple, value=value)
    return result

  def format(self, python_object, custom_converters=None):
    result = []
    for object in self.active_objects():
      if (isinstance(python_object, scope_extract)):
        python_object = [python_object]
      for python_object_i in python_object:
        sub_python_object = python_object_i.__get__(object.name)
        if (sub_python_object is not scope_extract_attribute_error):
          if (not object.multiple):
            result.append(object.format(
              sub_python_object, custom_converters))
          else:
            for sub_python_object_i in sub_python_object:
              result.append(object.format(
                sub_python_object_i, custom_converters))
    return self.copy(objects=result)

  def _fetch(self, source):
    assert source.name == self.name
    if (not isinstance(source, scope)):
      raise RuntimeError('Incompatible parameter objects "%s"%s and "%s"%s' %
        (self.name, self.where_str, source.name, source.where_str))
    result_objects = []
    for master_object in self.active_objects():
      if (len(self.name) == 0):
        path = master_object.name
      else:
        path = self.name + "." + master_object.name
      matching_sources = source.get(path=path, with_substitution=False)
      fetch_count = 0
      if (master_object.multiple):
        for matching_source in matching_sources.active_objects():
          fetch_count += 1
          result_objects.append(master_object.fetch(source=matching_source))
        if (fetch_count == 0):
          result_objects.append(copy.deepcopy(master_object))
          if (master_object.optional):
            result_objects[-1].is_disabled = True
        elif (fetch_count == 1
              and master_object.multiple and master_object.optional
              and master_object.has_same_definitions(result_objects[-1])):
          result_objects[-1].is_disabled = True
      else:
        result_object = master_object
        for matching_source in matching_sources.active_objects():
          fetch_count += 1
          result_object = result_object.fetch(source=matching_source)
        if (fetch_count == 0):
          result_objects.append(copy.deepcopy(master_object))
        else:
          result_objects.append(result_object)
    return result_objects

  def fetch(self, source=None, sources=None):
    assert [source, sources].count(None) == 1
    if (source is not None):
      result = objects=self._fetch(source=source)
    elif (len(sources) == 0):
      return self
    else:
      result = []
      for source in sources:
        result.extend(self._fetch(source=source))
    return self.copy(objects=clean_fetched_scope(fetched_scope=result))

  def process_includes(self,
        definition_type_names,
        reference_directory,
        include_stack=None):
    if (definition_type_names is None):
      definition_type_names = default_definition_type_names
    if (include_stack is None): include_stack = []
    result = []
    for object in self.objects:
      if (object.is_disabled):
        result.append(object)
      elif (isinstance(object, definition)):
        if (object.name != "include"):
          result.append(object)
        else:
          object_sub = object.resolve_variables()
          for file_name in [word.value for word in object_sub.words]:
            if (reference_directory is not None
                and not os.path.isabs(file_name)):
              file_name = os.path.join(reference_directory, file_name)
            result.extend(parse(
              file_name=file_name,
              definition_type_names=definition_type_names,
              process_includes=True,
              include_stack=include_stack).objects)
      else:
        result.append(object.process_includes(
          definition_type_names=definition_type_names,
          reference_directory=reference_directory,
          include_stack=include_stack))
    return self.copy(objects=result)

  def unique(self):
    selection = {}
    result = []
    for i_object,object in enumerate(self.active_objects()):
      selection[object.name] = i_object
    for i_object,object in enumerate(self.active_objects()):
      if (selection[object.name] == i_object):
        result.append(object.unique())
    return self.copy(objects=result)

def clean_fetched_scope(fetched_scope):
  result = []
  for object in fetched_scope:
    if (not isinstance(object, scope) or len(object.objects) < 2):
      result.append(object)
    else:
      child_group = []
      for child in object.objects:
        if (not child.merge_names):
          child_group.append(child)
        else:
          if (len(child_group) > 0):
            result.append(object.copy(objects=child_group))
            child_groups = []
          result.append(object.copy(objects=[child]))
      if (len(child_group) > 0):
        result.append(object.copy(objects=child_group))
  return result

class variable_substitution_fragment(object):

  __slots__ = ["is_variable", "value", "result"]

  def __init__(self, is_variable, value):
    self.is_variable = is_variable
    self.value = value

class variable_substitution_proxy(object):

  __slots__ = ["word", "force_string", "have_variables", "fragments"]

  def __init__(self, word):
    self.word = word
    self.force_string = word.quote_token is not None
    self.have_variables = False
    self.fragments = []
    fragment_value = ""
    char_iter = simple_tokenizer.character_iterator(word.value)
    c = char_iter.next()
    while (c is not None):
      if (c != "$"):
        fragment_value += c
        if (c == "\\" and char_iter.look_ahead_1() == "$"):
          fragment_value += char_iter.next()
        c = char_iter.next()
      else:
        self.have_variables = True
        if (len(fragment_value) > 0):
          self.fragments.append(variable_substitution_fragment(
            is_variable=False,
            value=fragment_value))
          fragment_value = ""
        c = char_iter.next()
        if (c is None):
          word.raise_syntax_error("$ must be followed by an identifier: ")
        if (c == "{"):
          while True:
            c = char_iter.next()
            if (c is None):
              word.raise_syntax_error('missing "}": ')
            if (c == "}"):
              c = char_iter.next()
              break
            fragment_value += c
          offs = int(fragment_value.startswith("."))
          if (not is_standard_identifier(fragment_value[offs:])):
            word.raise_syntax_error("improper variable name ")
          self.fragments.append(variable_substitution_fragment(
            is_variable=True,
            value=fragment_value))
        else:
          if (c not in standard_identifier_start_characters):
            word.raise_syntax_error("improper variable name ")
          fragment_value = c
          while True:
            c = char_iter.next()
            if (c is None): break
            if (c == "."): break
            if (c not in standard_identifier_continuation_characters): break
            fragment_value += c
          self.fragments.append(variable_substitution_fragment(
            is_variable=True,
            value=fragment_value))
        fragment_value = ""
    if (len(fragment_value) > 0):
      self.fragments.append(variable_substitution_fragment(
        is_variable=False,
        value=fragment_value))
    if (len(self.fragments) > 1):
      self.force_string = True

  def get_new_words(self):
    if (not self.have_variables):
      return [self.word]
    if (not self.force_string):
      return self.fragments[0].result
    return [simple_tokenizer.word(
      value="".join([fragment.result.value for fragment in self.fragments]),
      quote_token='"')]

def parse(
      input_string=None,
      source_info=None,
      file_name=None,
      definition_type_names=None,
      process_includes=False,
      include_stack=None):
  from iotbx.parameters import parser
  assert source_info is None or file_name is None
  if (input_string is None):
    assert file_name is not None
    input_string = open(file_name).read()
  if (definition_type_names is None):
    definition_type_names = default_definition_type_names
  result = scope(name="", primary_id=0)
  parser.collect_objects(
    word_iterator=simple_tokenizer.word_iterator(
      input_string=input_string,
      source_info=source_info,
      file_name=file_name,
      list_of_settings=[
        simple_tokenizer.settings(
          unquoted_single_character_words="{}=",
          contiguous_word_characters="",
          comment_characters="#"),
        simple_tokenizer.settings(
          unquoted_single_character_words="",
          contiguous_word_characters="")]),
    definition_type_names=definition_type_names,
    primary_id_generator=count(1),
    primary_parent_scope=result)
  if (process_includes):
    if (file_name is None):
      file_name_normalized = None
      reference_directory = None
    else:
      file_name_normalized = os.path.normpath(os.path.abspath(file_name))
      reference_directory = os.path.dirname(file_name_normalized)
      if (include_stack is None):
        include_stack = []
      elif (file_name_normalized in include_stack):
        raise RuntimeError("Include dependency cycle: %s"
          % ", ".join(include_stack+[file_name_normalized]))
      include_stack.append(file_name_normalized)
    result = result.process_includes(
      definition_type_names=definition_type_names,
      reference_directory=reference_directory,
      include_stack=include_stack)
    if (include_stack is not None):
      include_stack.pop()
  return result

def read_default(
      caller_file_name,
      params_extension=".params",
      definition_type_names=None,
      process_includes=True):
  params_file_name = os.path.splitext(caller_file_name)[0] + params_extension
  if (not os.path.isfile(params_file_name)):
    raise RuntimeError("Missing parameter file: %s" % params_file_name)
  return parse(
    file_name=params_file_name,
    definition_type_names=definition_type_names,
    process_includes=process_includes)
