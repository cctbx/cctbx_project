from __future__ import division
from iotbx import simple_tokenizer
from scitbx.python_utils.str_utils import line_breaker
from libtbx.itertbx import count
from libtbx import introspection
import copy
import math
import sys, os

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

def str_from_value_words(value_words):
  if (len(value_words) == 1 and value_words[0].value.lower() == "none"):
    return None
  return " ".join([word.value for word in value_words])

def bool_from_value_words(value_words):
  value_string = str_from_value_words(value_words)
  if (value_string is None): return None
  word_lower = value_words[0].value.lower()
  if (word_lower == "none"): return None
  if (word_lower in ["false", "no", "off", "0"]): return False
  if (word_lower in ["true", "yes", "on", "1"]): return True
  assert len(value_words) > 0
  raise RuntimeError(
    'One True of False value expected, "%s" found%s' % (
      value_string, value_words[0].where_str()))

def number_from_value_words(value_words):
  value_string = str_from_value_words(value_words)
  if (value_string is None): return None
  try: return eval(value_string, math.__dict__, {})
  except Exception, e:
    raise RuntimeError(
      'Error interpreting "%s" as a numeric expression: %s%s' % (
        value_string, str(e), value_words[0].where_str()))

def int_from_value_words(value_words):
  result = number_from_value_words(value_words)
  if (result is not None):
    if (isinstance(result, float)
        and round(result) == result):
      result = int(result)
    elif (not isinstance(result, int)):
      raise RuntimeError(
        'Integer expression expected, "%s" found%s' % (
          str_from_value_words(value_words), value_words[0].where_str()))
  return result

def float_from_value_words(value_words):
  result = number_from_value_words(value_words)
  if (result is not None):
    if (isinstance(result, int)):
      result = float(result)
    elif (not isinstance(result, float)):
      raise RuntimeError(
        'Floating-point expression expected, "%s" found%s' % (
          str_from_value_words(value_words), value_words[0].where_str()))
  return result

def choice_from_value_words(is_required, value_words):
  result = None
  for word in value_words:
    if (word.value.startswith("*")):
      if (result is not None):
        raise RuntimeError("Multiple choices where only one is possible%s" %
          value_words[0].where_str())
      result = word.value[1:]
  if (result is None and is_required):
    raise RuntimeError("Unspecified choice%s" % value_words[0].where_str())
  return result

def multi_choice_from_value_words(value_words):
  result = []
  for word in value_words:
    if (word.value.startswith("*")):
      result.append(word.value[1:])
  return result

def unit_cell_from_value_words(value_words):
  from cctbx import uctbx
  return uctbx.unit_cell(str_from_value_words(value_words))

def space_group_info_from_value_words(value_words):
  from cctbx import sgtbx
  return sgtbx.space_group_info(symbol=str_from_value_words(value_words))

default_definition_type_names = [
  "str", "bool", "int", "float",
  "choice", "multi_choice",
  "path", "key",
  "unit_cell", "space_group"]

def definition_type_from_value_words(value_words, type_names):
  if (len(value_words) == 1):
    word_lower = value_words[0].value.lower()
    if (word_lower == "none"): return None
    if (word_lower in type_names): return word_lower
  assert len(value_words) > 0
  raise RuntimeError(
    'Unexpected definition type: "%s"%s' % (
      value_words[0].value, value_words[0].where_str()))

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
        value = str(simple_tokenizer.word(value=value, quote_token='"'))
        indent = " " * (len(prefix) + 3 + len(name) + 3)
        if (len(indent+value) < print_width):
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
    "help", "caption", "short_caption", "required",
    "type", "input_size", "expert_level"]

  __slots__ = ["name", "values", "is_disabled"] + attribute_names

  def __init__(self,
        name,
        values,
        is_disabled=False,
        help=None,
        caption=None,
        short_caption=None,
        required=None,
        type=None,
        input_size=None,
        expert_level=None):
    self.name = name
    self.values = values
    self.is_disabled = is_disabled
    self.help = help
    self.caption = caption
    self.short_caption = short_caption
    self.required = required
    self.type = type
    self.input_size = input_size
    self.expert_level = expert_level

  def copy(self, values):
    keyword_args = {}
    for keyword in self.__slots__:
      keyword_args[keyword] = getattr(self, keyword)
    keyword_args["values"] = values
    return definition(**keyword_args)

  def merge(self, source):
    if (self.type not in ["choice", "multi_choice"]):
      return self.copy(source.values)
    flags = {}
    for word in self.values:
      if (word.value.startswith("*")): value = word.value[1:]
      else: value = word.value
      flags[value] = False
    for word in source.values:
      if (word.value.startswith("*")):
        value = word.value[1:]
        flag = True
      else:
        value = word.value
        flag = False
      if (flag and value not in flags):
        raise RuntimeError("Not a possible choice: %s%s" % (
          str(word), word.where_str()))
      flags[value] = flag
    values = []
    for word in self.values:
      if (word.value.startswith("*")): value = word.value[1:]
      else: value = word.value
      if (flags[value]): value = "*" + value
      values.append(simple_tokenizer.word(
        value=value,
        line_number=word.line_number,
        file_name=word.file_name))
    return self.copy(values)

  def has_attribute_with_name(self, name):
    return name in self.attribute_names

  def assign_attribute(self, name, value_words, type_names):
    assert self.has_attribute_with_name(name)
    if (name in ["required"]):
      value = bool_from_value_words(value_words)
    elif (name == "type"):
      value = definition_type_from_value_words(value_words, type_names)
    elif (name in ["input_size", "expert_level"]):
      value = int_from_value_words(value_words)
    else:
      value = str_from_value_words(value_words)
    setattr(self, name, value)

  def show(self, out, prefix="", attributes_level=0, print_width=79,
                 previous_object=None):
    if (previous_object is not None
        and not isinstance(previous_object, definition)):
      print >> out, prefix.rstrip()
    if (self.is_disabled): hash = "#"
    else:                  hash = ""
    line = prefix + hash + self.name
    if (self.name != "include"): line += " ="
    indent = " " * len(line)
    for word in self.values:
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

  def all_definitions(self, parent, parent_path, result):
    result.append(object_locator(
      parent=parent, path=parent_path+self.name, object=self))

  def get_without_substitution(self, path):
    if (self.name == path): return [self]
    return []

  def automatic_type(self):
    types = {}
    for word in self.values:
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
      return str_from_value_words(self.values)
    if (self.type == "bool"):
      return bool_from_value_words(self.values)
    if (self.type == "int"):
      return int_from_value_words(self.values)
    if (self.type == "float"):
      return float_from_value_words(self.values)
    if (self.type == "choice"):
      return choice_from_value_words(self.required, self.values)
    if (self.type == "multi_choice"):
      return multi_choice_from_value_words(self.values)
    if (self.type == "unit_cell"):
      return unit_cell_from_value_words(self.values)
    if (self.type == "space_group"):
      return space_group_info_from_value_words(self.values)
    if (custom_converters is not None):
      converter = custom_converters.get(self.type, None)
      if (converter is not None):
        return converter.process_value_words(self.values)
    if (self.type is None):
      return [word.value for word in self.values]
    raise RuntimeError(
       ('No converter for parameter definition type "%s"'
      + ' required for converting values of "%s"%s') % (
        self.type, self.name, self.values[0].where_str()))

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
      for word in self.values:
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
        words = converter.format_as_value_words(self, python_object)
    elif (self.type is None):
      words = [simple_tokenizer.word(value=value, quote_token='"')
        for value in python_object]
    if (words is None):
      raise RuntimeError(
         ('No converter for parameter definition type "%s"'
        + ' required for converting values of "%s"%s') % (
          self.type, self.name, self.values[0].where_str()))
    return self.copy(values=words)

def scope_merge(master, source):
  result_objects = []
  master_lookup_dict = {}
  for master_object in master.objects:
    master_lookup_dict[master_object.name] = master_object
  master_use_flags = {}
  for master_object in master.objects:
    master_use_flags[master_object.name] = False
  for source_object in source.objects:
    master_object = master_lookup_dict.get(source_object.name, None)
    if (master_object is None):
      result_objects.append(source_object)
    else:
      result_objects.append(master_object.merge(source=source_object))
      master_use_flags[source_object.name] = True
  for master_object in master.objects:
    if (not master_use_flags[master_object.name]
        and master_object.required is not False):
      result_objects.append(copy.copy(master_object))
  return result_objects

class scope_class: pass

class scope:

  def __init__(self,
        name,
        objects,
        is_disabled=False,
        style=None,
        help=None,
        caption=None,
        short_caption=None,
        required=None,
        sequential_format=None,
        disable_add=None,
        disable_delete=None,
        expert_level=None):
    introspection.adopt_init_args()
    self.attribute_names = self.__init__varnames__[4:]
    assert style in [None, "row", "column", "block", "page"]
    if (sequential_format is not None):
      assert isinstance(sequential_format % 0, str)

  def copy(self, objects):
    keyword_args = {}
    for keyword in self.__init__varnames__[1:]:
      keyword_args[keyword] = getattr(self, keyword)
    keyword_args["objects"] = objects
    return scope(**keyword_args)

  def merge(self, source):
    return self.copy(objects=scope_merge(self, source=source))

  def has_attribute_with_name(self, name):
    return name in self.attribute_names

  def assign_attribute(self, name, value_words):
    assert self.has_attribute_with_name(name)
    if (name in ["required", "disable_add", "disable_delete"]):
      value = bool_from_value_words(value_words)
    elif (name in ["expert_level"]):
      value = int_from_value_words(value_words)
    else:
      value = str_from_value_words(value_words)
      if (name == "style"):
        style = value
        assert style in [None, "row", "column", "block", "page"]
      elif (name == "sequential_format"):
        sequential_format = value
        if (sequential_format is not None):
          assert isinstance(sequential_format % 0, str)
    setattr(self, name, value)

  def show(self, out, prefix="", attributes_level=0, print_width=79,
                 previous_object=None):
    if (previous_object is not None):
      print >> out, prefix.rstrip()
    if (self.is_disabled): hash = "#"
    else:                  hash = ""
    print >> out, prefix + hash + self.name
    show_attributes(
      self=self,
      out=out,
      prefix=prefix,
      attributes_level=attributes_level,
      print_width=print_width)
    print >> out, prefix+"{"
    previous_object = None
    for object in self.objects:
      object.show(
        out=out,
        prefix=prefix+"  ",
        attributes_level=attributes_level,
        print_width=print_width,
        previous_object=previous_object)
      previous_object = object
    print >> out, prefix + "}"

  def all_definitions(self, parent, parent_path, result):
    parent_path += self.name+"."
    for object in self.objects:
      object.all_definitions(
        parent=self, parent_path=parent_path, result=result)

  def get_without_substitution(self, path):
    if (self.name == path): return [self]
    if (not path.startswith(self.name+".")): return []
    path = path[len(self.name)+1:]
    result = []
    for object in self.objects:
      result.extend(object.get_without_substitution(path=path))
    return result

  def extract(self, custom_converters=None):
    result = scope_class()
    for object in self.objects:
      setattr(result, object.name, object.extract(
        custom_converters=custom_converters))
    return result

class object_list:

  def __init__(self, objects):
    self.objects = objects

  def merge(self, source=None, sources=None):
    assert [source, sources].count(None) == 1
    if (source is not None):
      return object_list(objects=scope_merge(self, source=source))
    elif (len(sources) == 0):
      return self
    else:
      objects = []
      for source in sources:
        objects.extend(scope_merge(self, source=source))
      return object_list(objects=objects)

  def show(self, out=None, prefix="", attributes_level=0, print_width=None):
    if (out is None): out = sys.stdout
    if (print_width is None):
      print_width = 79
    previous_object = None
    for object in self.objects:
      object.show(
        out=out,
        prefix=prefix,
        attributes_level=attributes_level,
        print_width=print_width,
        previous_object=previous_object)
      previous_object = object
    return self

  def all_definitions(self):
    result = []
    for object in self.objects:
      object.all_definitions(parent=self, parent_path="", result=result)
    return result

  def get_without_substitution(self, path):
    if (len(path) == 0):
      result = [self]
    else:
      result = []
      for object in self.objects:
        result.extend(object.get_without_substitution(path))
    return object_list(objects=result)

  def get(self, path, with_substitution=True, path_memory=None):
    result_raw = self.get_without_substitution(path=path)
    if (not with_substitution):
      return result_raw
    if (path_memory is None):
      path_memory = {path: None}
    elif (path not in path_memory):
      path_memory[path] = None
    else:
      raise RuntimeError("Dependency cycle in variable substitution: $%s" % (
        path))
    result_sub = []
    for object in result_raw.objects:
      if (not isinstance(object, definition)):
        result_sub.append(object)
      else:
        result_sub.append(self.variable_substitution(
          object=object, path_memory=path_memory))
    del path_memory[path]
    return object_list(objects=result_sub)

  def get_last(self, path, with_substitution=True):
    all = self.get(path=path, with_substitution=with_substitution)
    if (len(all.objects) == 0): return None
    return all.objects[-1]

  def extract_last(self, path, with_substitution=True, custom_converters=None):
    raw = self.get_last(path=path, with_substitution=with_substitution)
    if (raw is None): return None
    return raw.extract(custom_converters=custom_converters)

  def extract(self, master, custom_converters=None):
    return object_list_extract(
      object_list=self, master=master, custom_converters=custom_converters)

  def merge_and_extract(self,
        source=None,
        sources=None,
        custom_converters=None):
    return self.merge(source=source, sources=sources).extract(
      master=self, custom_converters=custom_converters)

  def variable_substitution(self, object, path_memory):
    new_values = []
    for word in object.values:
      if (word.quote_token == "'"):
        new_values.append(word)
        continue
      substitution_proxy = variable_substitution_proxy(word)
      for fragment in substitution_proxy.fragments:
        if (not fragment.is_variable):
          fragment.result = simple_tokenizer.word(
            value=fragment.value, quote_token='"')
          continue
        variable_values = None
        for variable_object in self.get(
                                 path=fragment.value,
                                 path_memory=path_memory).objects:
          if (isinstance(variable_object, definition)):
            variable_values = variable_object.values
        if (variable_values is None):
          env_var = os.environ.get(fragment.value, None)
          if (env_var is not None):
            variable_values = [simple_tokenizer.word(
              value=env_var,
              quote_token='"')]
        if (variable_values is None):
          raise RuntimeError("Undefined variable: $%s%s" % (
            fragment.value, word.where_str()))
        if (not substitution_proxy.force_string):
          fragment.result = variable_values
        else:
          fragment.result = simple_tokenizer.word(
            value=" ".join([str(v) for v in variable_values]),
            quote_token='"')
      new_values.extend(substitution_proxy.get_new_values())
    return object.copy(values=new_values)

  def process_includes(self,
        definition_type_names,
        reference_directory,
        include_memory=None):
    if (definition_type_names is None):
      definition_type_names = default_definition_type_names
    if (include_memory is None): include_memory = {}
    result = []
    for object in self.objects:
      if (not isinstance(object, definition)
          or object.name != "include"
          or object.is_disabled):
        result.append(object)
      else:
        object_sub = self.variable_substitution(object=object, path_memory={})
        for file_name in [word.value for word in object_sub.values]:
          if (reference_directory is not None
              and not os.path.isabs(file_name)):
            file_name = os.path.join(reference_directory, file_name)
          file_name_normalized = os.path.normpath(os.path.abspath(file_name))
          if (file_name_normalized in include_memory): continue
          include_memory[file_name_normalized] = None
          result.extend(parse(
            file_name=file_name,
            definition_type_names=definition_type_names).process_includes(
              definition_type_names=definition_type_names,
              reference_directory=os.path.dirname(file_name_normalized),
              include_memory=include_memory).objects)
    return object_list(objects=result)

  def automatic_type_assignment(self, assignment_if_unknown=None):
    for item in self.all_definitions():
      item.object.automatic_type_assignment(
        assignment_if_unknown=assignment_if_unknown)

class object_list_extract:

  def __init__(self, object_list, master, custom_converters):
    self.__object_list__ = object_list
    self.__master__ = master
    for object in master.objects:
      attr_name = object.name.split(".")[-1]
      assert attr_name not in ["__object_list__", "__master__"]
      setattr(self, attr_name, object_list.extract_last(
        path=object.name, custom_converters=custom_converters))

  def format(self):
    scopes = []
    for scope_object in self.__master__.objects:
      assert isinstance(scope_object, scope) # not implemented
      attr_name = scope_object.name.split(".")[-1]
      py_scope = getattr(self, attr_name)
      definitions = []
      for definition_object in scope_object.objects:
        assert isinstance(definition_object, definition) # not implemented
        value = getattr(py_scope, definition_object.name)
        definitions.append(definition_object.format(value))
      scopes.append(scope_object.copy(objects=definitions))
    return object_list(objects=scopes)

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
          if (not is_standard_identifier(fragment_value)):
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

  def get_new_values(self):
    if (not self.have_variables):
      return [self.word]
    if (not self.force_string):
      return self.fragments[0].result
    return [simple_tokenizer.word(
      value="".join([fragment.result.value for fragment in self.fragments]),
      quote_token='"')]

def parse(
      input_string=None,
      file_name=None,
      definition_type_names=None,
      process_includes=False):
  from iotbx.parameters import parser
  if (input_string is None):
    assert file_name is not None
    input_string = open(file_name).read()
  if (definition_type_names is None):
    definition_type_names = default_definition_type_names
  result = object_list(objects=parser.collect_objects(
    word_iterator=simple_tokenizer.word_iterator(
      input_string=input_string,
      file_name=file_name,
      list_of_settings=[
        simple_tokenizer.settings(
          unquoted_single_character_words="{}=",
          contiguous_word_characters="",
          comment_characters="#"),
        simple_tokenizer.settings(
          unquoted_single_character_words="",
          contiguous_word_characters="")]),
    definition_type_names=definition_type_names))
  if (process_includes):
    if (file_name is None):
      reference_directory = None
    else:
      reference_directory = os.path.dirname(os.path.abspath(file_name))
    result = result.process_includes(
      definition_type_names=definition_type_names,
      reference_directory=reference_directory)
  return result

def read_default(
      params_extension=".params",
      prefix=None,
      definition_type_names=None,
      process_includes=True):
  caller_file_name = introspection.caller_location(frames_back=1).file_name
  assert os.path.isfile(caller_file_name)
  assert caller_file_name.endswith(".py")
  params_file_name = caller_file_name[:-3] + params_extension
  if (not os.path.isfile(params_file_name)):
    raise RuntimeError("Missing parameter file: %s" % params_file_name)
  result = parse(
    file_name=params_file_name,
    definition_type_names=definition_type_names,
    process_includes=process_includes)
  if (prefix is not None):
    for object in result.objects:
      object.name = prefix + "." + object.name
  return result
