from __future__ import generators
from iotbx import simple_tokenizer
from scitbx.python_utils.str_utils import line_breaker
from libtbx import introspection

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
    'One True of False value expected, "%s" found (input line %d)' % (
      value_string, value_words[0].line_number))

def number_from_value_words(value_words, number_type, number_type_name):
  value_string = str_from_value_words(value_words)
  if (value_string is None): return None
  try: value = number_type(value_string)
  except:
    raise RuntimeError(
      '%s value expected, "%s" found (input line %d)' % (
        number_type_name, value_string, value_words[0].line_number))

def int_from_value_words(value_words):
  return number_from_value_words(value_words, int, "Integer")

def float_from_value_words(value_words):
  return number_from_value_words(value_words, float, "Floating-point")

default_definition_type_names = [
  "bool", "int", "float", "str",
  "choice", "multi_choice",
  "path", "key",
  "unit_cell", "space_group"]

def definition_type_from_value_words(value_words, type_names):
  if (len(value_words) == 1):
    word_lower = value_words[0].value.lower()
    if (word_lower == "none"): return None
    if (word_lower in ["bool", "int", "float", "str",
                       "choice", "multi_choice",
                       "path", "key"]):
      return word_lower
  assert len(value_words) > 0
  raise RuntimeError(
    'Unexpected definition type: "%s" (input line %d)' % (
      value_words[0].value, value_words[0].line_number))

class definition:

  def __init__(self,
        name=None,
        values=None,
        help=None,
        caption=None,
        short_caption=None,
        required=None,
        type=None,
        input_size=None,
        expert_level=None):
    introspection.adopt_init_args()
    self.attribute_names = self.__init__varnames__[3:]

  def has_attribute_with_name(self, name):
    return name in self.attribute_names

  def assign_attribute(self, name, value_words, type_names):
    assert self.has_attribute_with_name(name)
    if (name == "required"):
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
    line = prefix+self.name
    indent = " " * len(line)
    for word in self.values:
      line_plus = line + " " + str(word)
      if (len(line_plus) > print_width-2 and len(line) > len(indent)):
        print >> out, line + " \\"
        line = indent + " " + str(word)
      else:
        line = line_plus
    print >> out, line
    if (attributes_level > 0):
      for name in self.attribute_names:
        value = getattr(self, name)
        if ((name == "help" and value is not None)
            or (value is not None and attributes_level > 1)
            or attributes_level > 2):
          if (not isinstance(value, str)):
            print >> out, prefix+"  ."+name, value
          else:
            value = '"' + value.replace('"', '\\"') + '"'
            indent = " " * (len(prefix) + 3 + len(name) + 1)
            if (len(indent+value) < print_width):
              print >> out, prefix+"  ."+name, value
            else:
              is_first = True
              for block in line_breaker(value[1:-1],print_width-2-len(indent)):
                if (is_first):
                  print >> out, prefix+"  ."+name, '"'+block+'"'
                  is_first = False
                else:
                  print >> out, indent+'"'+block+'"'

class scope:

  def __init__(self,
        name=None,
        objects=[]):
    introspection.adopt_init_args()

  def show(self, out, prefix="", attributes_level=0, print_width=79,
                 previous_object=None):
    if (previous_object is not None):
      print >> out, prefix.rstrip()
    print >> out, prefix + self.name + " {"
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

class table:

  def __init__(self,
        name=None,
        row_names=[],
        row_objects=[],
        style=None,
        help=None,
        caption=None,
        short_caption=None,
        sequential_format=None,
        disable_add=None,
        disable_delete=None):
    introspection.adopt_init_args()
    self.attribute_names = self.__init__varnames__[4:]
    assert style in [None, "row", "column", "block", "page"]
    if (sequential_format is not None):
      assert isinstance(sequential_format % 0, str)

  def has_attribute_with_name(self, name):
    return name in self.attribute_names

  def assign_attribute(self, name, value_words):
    assert self.has_attribute_with_name(name)
    if (name in ["disable_add", "disable_delete"]):
      value = bool_from_value_words(value_words)
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

  def add_row(self, name, objects):
    self.row_names.append(name)
    self.row_objects.append(objects)

  def show(self, out, prefix="", attributes_level=0, print_width=79,
                 previous_object=None):
    if (previous_object is not None):
      print >> out, prefix.rstrip()
    print >> out, "%stable %s" % (prefix, self.name)
    if (attributes_level > 0):
      for name in self.attribute_names:
        value = getattr(self, name)
        if ((name == "help" and value is not None)
            or (value is not None and attributes_level > 1)
            or attributes_level > 2):
          print >> out, prefix+"  ."+name, value
    print >> out, prefix+"{"
    assert len(self.row_names) == len(self.row_objects)
    for name,objects in zip(self.row_names, self.row_objects):
      s = prefix + "  "
      if (name is not None):
        s += name + " "
      print >> out, s+"{"
      previous_object = None
      for object in objects:
        object.show(
          out=out,
          prefix=prefix+"    ",
          attributes_level=attributes_level,
          print_width=print_width,
          previous_object=previous_object)
        previous_object = object
      print >> out, prefix+"  }"
    print >> out, prefix+"}"

class parse:

  def __init__(self, input_string, definition_type_names=None):
    from iotbx.parameters import parser
    if (definition_type_names is None):
      definition_type_names = default_definition_type_names
    self.objects = parser.collect_objects(
      word_stack=simple_tokenizer.as_word_stack(
        input_string=input_string,
        contiguous_word_characters=""),
      definition_type_names=definition_type_names)

  def show(self, out, prefix="", attributes_level=0, print_width=None):
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
