"Documentation: http://cctbx.sourceforge.net/libtbx_phil.html"

from __future__ import division
from libtbx.phil import tokenizer
from libtbx.str_utils import line_breaker
from libtbx.utils import Sorry, format_exception, import_python_object
from itertools import count
from libtbx import Auto, slots_getstate_setstate
from cStringIO import StringIO
import tokenize as python_tokenize
import warnings
import math
import weakref
import sys, os

default_print_width = 79

class PhilDeprecationWarning (DeprecationWarning) :
  pass
warnings.filterwarnings("always", category=PhilDeprecationWarning)

def is_reserved_identifier(string):
  if (len(string) < 5): return False
  return (string.startswith("__") and string.endswith("__"))

standard_identifier_start_characters = set()
for c in "_ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz":
  standard_identifier_start_characters.add(c)
standard_identifier_continuation_characters = set(
  standard_identifier_start_characters)
for c in ".0123456789":
  standard_identifier_continuation_characters.add(c)

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

def is_plain_none(words):
  return (len(words) == 1
          and words[0].quote_token is None
          and words[0].value.lower() == "none")

def is_plain_auto(words):
  return (len(words) == 1
          and words[0].quote_token is None
          and words[0].value.lower() == "auto")

def tokenize_value_literal(input_string, source_info):
  return list(tokenizer.word_iterator(
    input_string=input_string,
    source_info=source_info,
    list_of_settings=[
      tokenizer.settings(contiguous_word_characters="")]))

class words_converters(object):

  phil_type = "words"

  def __str__(self): return self.phil_type

  def from_words(self, words, master):
    if (is_plain_none(words=words)): return None
    if (is_plain_auto(words=words)): return Auto
    return words

  def as_words(self, python_object, master):
    if (python_object is None):
      return [tokenizer.word(value="None")]
    if (python_object is Auto):
      return [tokenizer.word(value="Auto")]
    for word in python_object:
      assert isinstance(word, tokenizer.word)
    return python_object

def strings_from_words(words):
  if (is_plain_none(words=words)): return None
  if (is_plain_auto(words=words)): return Auto
  return [word.value for word in words]

def strings_as_words(python_object):
  if (python_object is None):
    return [tokenizer.word(value="None")]
  if (python_object is Auto):
    return [tokenizer.word(value="Auto")]
  words = []
  for value in python_object:
    if (is_standard_identifier(value)):
      words.append(tokenizer.word(value=value))
    else:
      words.append(tokenizer.word(value=value, quote_token='"'))
  return words

class strings_converters(object):

  phil_type = "strings"

  def __str__(self): return self.phil_type

  def from_words(self, words, master):
    return strings_from_words(words)

  def as_words(self, python_object, master):
    return strings_as_words(python_object)

def str_from_words(words):
  if (is_plain_none(words=words)): return None
  if (is_plain_auto(words=words)): return Auto
  return " ".join([word.value for word in words])

class str_converters(object):

  phil_type = "str"

  def __str__(self): return self.phil_type

  def from_words(self, words, master):
    return str_from_words(words=words)

  def as_words(self, python_object, master):
    if (python_object is None):
      return [tokenizer.word(value="None")]
    if (python_object is Auto):
      return [tokenizer.word(value="Auto")]
    return [tokenizer.word(value=python_object, quote_token='"')]

class qstr_converters(object):

  phil_type = "qstr"

  def __str__(self): return self.phil_type

  def from_words(self, words, master):
    if (is_plain_none(words=words)): return None
    if (is_plain_auto(words=words)): return Auto
    return " ".join([str(word) for word in words])

  def as_words(self, python_object, master):
    if (python_object is None):
      return [tokenizer.word(value="None")]
    if (python_object is Auto):
      return [tokenizer.word(value="Auto")]
    return tokenize_value_literal(
      input_string=python_object,
      source_info="python_object")

class path_converters(str_converters):

  phil_type = "path"

  def __str__(self): return self.phil_type

class key_converters(str_converters):

  phil_type = "key"

  def __str__(self): return self.phil_type

def bool_from_words(words, path):
  value_string = str_from_words(words)
  if (value_string is None): return None
  if (value_string is Auto): return Auto
  value_lower = value_string.lower()
  if (value_lower in ["false", "no", "off", "0"]): return False
  if (value_lower in ["true", "yes", "on", "1"]): return True
  assert len(words) > 0
  raise RuntimeError(
    'One True or False value expected, %s="%s" found%s' % (
      path, value_string, words[0].where_str()))

class bool_converters(object):

  phil_type = "bool"

  def __str__(self): return self.phil_type

  def from_words(self, words, master):
    return bool_from_words(words=words, path=master.full_path())

  def as_words(self, python_object, master):
    if (python_object is None):
      return [tokenizer.word(value="None")]
    if (python_object is Auto):
      return [tokenizer.word(value="Auto")]
    if (python_object):
      return [tokenizer.word(value="True")]
    else:
      return [tokenizer.word(value="False")]

def number_from_value_string(value_string, words, path):
  if (value_string is None): return None
  if (value_string is Auto): return Auto
  value_string_lower = value_string.lower()
  # similar to libtbx.utils.number_from_string
  # (please review if making changes here)
  value_string_lower_strip = value_string.lower().strip()
  if (value_string_lower_strip in ["true", "false"]):
    raise RuntimeError(
      'Error interpreting %s="%s" as a numeric expression%s' % (
        path, value_string, words[0].where_str()))
  if (value_string_lower_strip == "none"): return None
  if (value_string_lower_strip == "auto"): return Auto
  try: return int(value_string)
  except KeyboardInterrupt: raise
  except Exception: pass
  try: return eval(value_string, math.__dict__, {})
  except KeyboardInterrupt: raise
  except Exception:
    raise RuntimeError(
      'Error interpreting %s="%s" as a numeric expression: %s%s' % (
        path, value_string, format_exception(), words[0].where_str()))

def number_from_words(words, path):
  return number_from_value_string(
    value_string=str_from_words(words), words=words, path=path)

def numbers_from_words(words, path):
  all_values_string = str_from_words(words)
  if (all_values_string is None or all_values_string is Auto):
    return all_values_string
  while True:
    have_changes = False
    for o,c in ["()", "[]"]:
      while (    all_values_string.startswith(o)
             and all_values_string.endswith(c)):
        all_values_string = all_values_string[1:-1].strip()
        have_changes = True
    if (not have_changes):
      break
  result = []
  for value_string in all_values_string \
                        .replace(",", " ") \
                        .replace(";", " ") \
                        .split():
    result.append(number_from_value_string(
      value_string=value_string, words=words, path=path))
  return result

def int_from_number(number, words, path):
  if (isinstance(number, int)): return number
  if (isinstance(number, float)
      and round(number) == number):
    return int(number)
  raise RuntimeError(
    'Error interpreting %s="%s" as an integer expression%s' % (
      path, str_from_words(words), words[0].where_str()))

def float_from_number(number, words, path):
  if (isinstance(number, float)): return number
  if (isinstance(number, int)): return float(number)
  raise RuntimeError(
    'Error interpreting %s="%s" as a floating-point expression%s' % (
      path, str_from_words(words), words[0].where_str()))

def int_from_words(words, path):
  result = number_from_words(words=words, path=path)
  if (result is None or result is Auto):
    return result
  return int_from_number(number=result, words=words, path=path)

def float_from_words(words, path):
  result = number_from_words(words=words, path=path)
  if (result is None or result is Auto):
    return result
  return float_from_number(number=result, words=words, path=path)

class _check_value_base(object):

  def _check_value(self, value, path_producer, words=None):
    def where_str():
      if (words is None): return ""
      return words[0].where_str()
    if (self.value_min is not None and value < self.value_min):
      raise RuntimeError(
        "%s element is less than the minimum allowed value:"
        " %s < %s%s"
          % (path_producer(), self._value_as_str(value=value),
             self._value_as_str(value=self.value_min), where_str()))
    if (self.value_max is not None and value > self.value_max):
      raise RuntimeError(
        "%s element is greater than the maximum allowed value:"
        " %s > %s%s"
          % (path_producer(), self._value_as_str(value=value),
             self._value_as_str(value=self.value_max), where_str()))

class number_converters_base(_check_value_base):

  def __init__(self,
      value_min=None,
      value_max=None):
    if (value_min is not None and value_max is not None):
      assert value_min <= value_max
    self.value_min = value_min
    self.value_max = value_max

  def __str__(self):
    kwds = []
    if (self.value_min is not None):
      kwds.append("value_min=" + self._value_as_str(value=self.value_min))
    if (self.value_max is not None):
      kwds.append("value_max=" + self._value_as_str(value=self.value_max))
    if (len(kwds) != 0):
      return self.phil_type + "(" + ", ".join(kwds) + ")"
    return self.phil_type

  def from_words(self, words, master):
    path = master.full_path()
    value = self._value_from_words(words=words, path=master.full_path())
    if (value is None or value is Auto): return value
    self._check_value(
      value=value, path_producer=master.full_path, words=words)
    return value

  def as_words(self, python_object, master):
    if (python_object is None):
      return [tokenizer.word(value="None")]
    if (python_object is Auto):
      return [tokenizer.word(value="Auto")]
    return [tokenizer.word(value=self._value_as_str(value=python_object))]

class int_converters(number_converters_base):

  phil_type = "int"

  def _value_from_words(self, words, path):
    return int_from_words(words=words, path=path)

  def _value_as_str(self, value):
    return "%d" % value

class float_converters(number_converters_base):

  phil_type = "float"

  def _value_from_words(self, words, path):
    return float_from_words(words=words, path=path)

  def _value_as_str(self, value):
    return "%.10g" % value

class numbers_converters_base(_check_value_base):

  def __init__(self,
      size=None,
      size_min=None,
      size_max=None,
      value_min=None,
      value_max=None,
      allow_none_elements=False,
      allow_auto_elements=False):
    assert size is None or (size_min is None and size_max is None)
    if (size is not None):
      assert size > 0
      size_min = size
      size_max = size
    else:
      if (size_min is not None):
        assert size_min > 0
      if (size_max is not None):
        assert size_max > 0
        if (size_min is not None):
          assert size_min <= size_max
    if (value_min is not None and value_max is not None):
      assert value_min <= value_max
    self.size_min = size_min
    self.size_max = size_max
    self.value_min = value_min
    self.value_max = value_max
    self.allow_none_elements = allow_none_elements
    self.allow_auto_elements = allow_auto_elements

  def __str__(self):
    kwds = []
    if (self.size_min == self.size_max):
      if (self.size_min is not None):
        kwds.append("size=%d" % self.size_min)
    else:
      if (self.size_min is not None):
        kwds.append("size_min=%d" % self.size_min)
      if (self.size_max is not None):
        kwds.append("size_max=%d" % self.size_max)
    if (self.value_min is not None):
      kwds.append("value_min=" + self._value_as_str(value=self.value_min))
    if (self.value_max is not None):
      kwds.append("value_max=" + self._value_as_str(value=self.value_max))
    if (self.allow_none_elements):
      kwds.append("allow_none_elements=True")
    if (self.allow_auto_elements):
      kwds.append("allow_auto_elements=True")
    if (len(kwds) != 0):
      return self.phil_type + "(" + ", ".join(kwds) + ")"
    return self.phil_type

  def _check_size(self, size, path_producer, words=None):
    def where_str():
      if (words is None): return ""
      return words[0].where_str()
    if (self.size_max is not None and size > self.size_max):
      if (self.size_max == self.size_min):
        precise = "exactly %d required"
      else:
        precise = "%d allowed at most"
      raise RuntimeError(
        "Too many values for %s: %d given, %s%s"
          % (path_producer(), size, (precise%self.size_max), where_str()))
    if (self.size_min is not None and size < self.size_min):
      if (self.size_max == self.size_min):
        precise = "exactly"
      else:
        precise = "at least"
      raise RuntimeError(
        "Not enough values for %s: %d given, %s %d required%s"
          % (path_producer(), size, precise, self.size_min, where_str()))

  def from_words(self, words, master):
    path = master.full_path()
    numbers = numbers_from_words(words=words, path=path)
    if (numbers is None or numbers is Auto): return numbers
    self._check_size(
      size=len(numbers), path_producer=master.full_path, words=words)
    def where_str():
      if (words is None): return ""
      return words[0].where_str()
    result = []
    for number in numbers:
      if   (number is None):
        if (self.allow_none_elements):
          value = number
        else:
          raise RuntimeError(
            "%s element cannot be None%s" % (path, where_str()))
      elif (number is Auto):
        if (self.allow_auto_elements):
          value = number
        else:
          raise RuntimeError(
            "%s element cannot be Auto%s" % (path, where_str()))
      else:
        value = self._value_from_number(number=number, words=words, path=path)
        self._check_value(
          value=value, path_producer=master.full_path, words=words)
      result.append(value)
    return result

  def as_words(self, python_object, master):
    if (python_object is None):
      return [tokenizer.word(value="None")]
    if (python_object is Auto):
      return [tokenizer.word(value="Auto")]
    self._check_size(
      size=len(python_object), path_producer=master.full_path)
    result = []
    for value in python_object:
      self._check_value(value=value, path_producer=master.full_path)
      if (value is None):
        if (self.allow_none_elements):
          result.append(tokenizer.word(value="None"))
        else:
          raise RuntimeError(
            "%s element cannot be None" % master.full_path())
      elif (value is Auto):
        if (self.allow_auto_elements):
          result.append(tokenizer.word(value="Auto"))
        else:
          raise RuntimeError(
            "%s element cannot be Auto" % master.full_path())
      else:
        result.append(tokenizer.word(value=self._value_as_str(value=value)))
    return result

class ints_converters(numbers_converters_base):

  phil_type = "ints"

  def _value_from_number(self, number, words, path):
    return int_from_number(number=number, words=words, path=path)

  def _value_as_str(self, value):
    return "%d" % value

class floats_converters(numbers_converters_base):

  phil_type = "floats"

  def _value_from_number(self, number, words, path):
    return float_from_number(number=number, words=words, path=path)

  def _value_as_str(self, value):
    return "%.10g" % value

class choice_converters(object):

  phil_type = "choice"

  def __init__(self, multi=False):
    self.multi = multi

  def __str__(self):
    if (self.multi): return self.phil_type+"(multi=True)"
    return self.phil_type

  def from_words(self, words, master):
    if (is_plain_auto(words=words)):
      result = Auto
    elif (self.multi):
      result = []
      for word in words:
        if (word.value.startswith("*")):
          result.append(word.value[1:])
      if (len(result) == 0
          and master.optional is not None and not master.optional):
        raise RuntimeError(
          "Unspecified choice for %s:"
          " at least one choice must be selected%s" % (
            master.full_path(), words[0].where_str()))
    else:
      result = None
      for word in words:
        if (word.value.startswith("*")):
          if (result is not None):
            raise RuntimeError(
              "Multiple choices for %s;"
              " only one choice can be selected%s" % (
                master.full_path(), words[0].where_str()))
          result = word.value[1:]
      if (result is None
          and (master.optional is not None and not master.optional)):
        raise RuntimeError(
          "Unspecified choice for %s:"
          " exactly one choice must be selected%s" % (
            master.full_path(), words[0].where_str()))
    return result

  def as_words(self, python_object, master):
    if (python_object is Auto):
      return [tokenizer.word(value="Auto")]
    assert not self.multi or python_object is not None
    if (self.multi):
      use_flags = dict([(value, False) for value in python_object])
    n_choices = 0
    def raise_improper_master():
      raise RuntimeError("Improper master choice definition: %s%s" % (
        master.as_str().rstrip(), master.words[0].where_str()))
    words = []
    for word in master.words:
      if (word.value.startswith("*")): value = word.value[1:]
      else: value = word.value
      if (python_object is not None):
        if (not self.multi):
          if (value == python_object):
            value = "*" + value
            n_choices += 1
            if (n_choices > 1): raise_improper_master()
        else:
          if (value in use_flags):
            if (use_flags[value]): raise_improper_master()
            use_flags[value] = True
            value = "*" + value
            n_choices += 1
      words.append(tokenizer.word(
        value=value, quote_token=word.quote_token))
    if (not self.multi):
      if (n_choices == 0
          and ((master.optional is not None and not master.optional)
               or python_object is not None)):
        raise RuntimeError("Invalid choice: %s=%s" % (
          master.full_path(), str(python_object)))
    else:
      unused = []
      for value,use_flag in use_flags.items():
        if (not use_flag): unused.append(value)
      n = len(unused)
      if (n != 0):
        raise RuntimeError("Invalid %s: %s=%s" % (
          str(self), master.full_path(), str(unused)))
      if (n_choices == 0
          and (master.optional is not None and not master.optional)):
        raise RuntimeError(
          "Empty list for mandatory %s: %s" % (
            str(self), master.full_path()))
    return words

  def fetch(self, source_words, master):
    assert not is_plain_none(words=master.words)
    assert not is_plain_auto(words=master.words)
    if (is_plain_auto(words=source_words)):
      return master.customized_copy(words=[tokenizer.word(value="Auto")])
    flags = {}
    for word in master.words:
      if (word.value.startswith("*")): value = word.value[1:]
      else: value = word.value
      flags[value.lower()] = False
    if (   (master.optional is not None and not master.optional)
        or not is_plain_none(words=source_words)):
      have_quote_or_star = False
      have_plus = False
      for word in source_words:
        if (word.quote_token is not None or word.value.startswith("*")):
          have_quote_or_star = True
          break
        if (word.value.find("+") >= 0):
          have_plus = True
      process_plus = False
      if (not have_quote_or_star and have_plus):
        values = "".join([word.value for word in source_words]).split("+")
        for value in values[1:]:
          if (len(value.strip()) == 0):
            break
        else:
          process_plus = True
      def raise_not_a_possible_choice(value):
        raise Sorry(
          "Not a possible choice for %s: %s%s\n" % (
            master.full_path(), value, word.where_str())
          + "  Possible choices are:\n"
          + "    " + "\n    ".join([w.value for w in master.words]))
      if (process_plus):
        for word in source_words:
          for value in word.value.split("+"):
            if (len(value) == 0): continue
            if (value not in flags):
              raise_not_a_possible_choice(value)
            flags[value.lower()] = True
      else:
        for word in source_words:
          if (word.value.startswith("*")):
            value = word.value[1:]
            flag = True
          else:
            value = word.value
            if (len(source_words) == 1):
              flag = True
            else:
              flag = False
          if (flag and value.lower() not in flags):
            raise_not_a_possible_choice(value)
          flags[value.lower()] = flag
    words = []
    for word in master.words:
      if (word.value.startswith("*")): value = word.value[1:]
      else: value = word.value
      if (flags[value.lower()]): value = "*" + value
      words.append(tokenizer.word(
        value=value,
        quote_token=word.quote_token,
        line_number=word.line_number,
        source_info=word.source_info))
    return master.customized_copy(words=words)

def get_converters_phil_type(converters):
  result = getattr(converters, "phil_type", None)
  if (result is None):
    result = str(converters()) # backward compatibility
  return result

def extended_converter_registry(additional_converters, base_registry=None):
  if (base_registry is None): base_registry = default_converter_registry
  result = dict(base_registry)
  for converters in additional_converters:
    result[get_converters_phil_type(converters)] = converters
  return result

default_converter_registry = extended_converter_registry(
  additional_converters = [
    words_converters,
    strings_converters,
    str_converters,
    qstr_converters,
    path_converters,
    key_converters,
    bool_converters,
    int_converters,
    float_converters,
    ints_converters,
    floats_converters,
    choice_converters],
  base_registry={})

def extract_args(*args, **keyword_args):
  return args, keyword_args

def normalize_call_expression(expression):
  result = []
  p = ""
  for info in python_tokenize.generate_tokens(StringIO(expression).readline):
    t = info[1]
    if (len(t) == 0): continue
    if (    t != "."
        and t[0] in standard_identifier_start_characters
        and len(p) > 0
        and p != "."
        and p[-1] in standard_identifier_continuation_characters):
      result.append(" ")
    result.append(t)
    if (t[0] == ","):
      result.append(" ")
    p = t
  return "".join(result)

def definition_converters_from_words(
      words,
      converter_registry,
      converter_cache):
  if (is_plain_none(words=words)): return None
  if (is_plain_auto(words=words)): return Auto
  call_expression_raw = str_from_words(words).strip()
  try:
    call_expression = normalize_call_expression(expression=call_expression_raw)
  except python_tokenize.TokenError, e:
    raise RuntimeError(
      'Error evaluating definition type "%s": %s%s' % (
        call_expression_raw, str(e), words[0].where_str()))
  converters_weakref = converter_cache.get(call_expression, None)
  if (converters_weakref is not None):
    converters_instance = converters_weakref()
    if (converters_instance is not None):
      return converters_instance
  flds = call_expression.split("(", 1)
  converters = converter_registry.get(flds[0], None)
  if (converters is not None):
    if (len(flds) == 1): parens = "()"
    else:                parens = ""
    try:
      converters_instance = eval(
        call_expression+parens, math.__dict__, {flds[0]: converters})
    except KeyboardInterrupt: raise
    except Exception:
      raise RuntimeError(
        'Error constructing definition type "%s": %s%s' % (
        call_expression, format_exception(), words[0].where_str()))
  else:
    import_path = flds[0] + "_phil_converters"
    if (len(flds) == 1):
      keyword_args = {}
    else:
      extractor = "__extract_args__(" + flds[1]
      try:
        args, keyword_args = eval(
          extractor, math.__dict__, {"__extract_args__": extract_args})
      except KeyboardInterrupt: raise
      except Exception:
        raise RuntimeError(
          'Error evaluating definition type "%s": %s%s' % (
          call_expression, format_exception(), words[0].where_str()))
    try:
      imported = import_python_object(
        import_path=import_path,
        error_prefix='.type=%s: ' % call_expression,
        target_must_be="; target must be a callable Python object",
        where_str=words[0].where_str())
    except (ValueError, ImportError):
      raise RuntimeError(
        'Unexpected definition type: "%s"%s' % (
          call_expression, words[0].where_str()))
    if (not callable(imported.object)):
      raise TypeError(
        '"%s" is not a callable Python object%s' % (
          import_path, words[0].where_str()))
    try:
      converters_instance = imported.object(**keyword_args)
    except KeyboardInterrupt: raise
    except Exception:
      raise RuntimeError(
        'Error constructing definition type "%s": %s%s' % (
        call_expression, format_exception(), words[0].where_str()))
  converter_cache[call_expression] = weakref.ref(converters_instance)
  return converters_instance

def full_path(self):
  result = [self.name]
  pps = self.primary_parent_scope
  while (pps is not None):
    if (pps.name == ""): break
    result.append(pps.name)
    pps = pps.primary_parent_scope
  result.reverse()
  return ".".join(result)

def show_attributes(self, out, prefix, attributes_level, print_width):
  if (attributes_level <= 0): return
  for name in self.attribute_names:
    value = getattr(self, name)
    if (name == "deprecated") and (not value) :
      continue # only show .deprecated if True
    if ((name == "help" and value is not None)
        or (value is not None and attributes_level > 1)
        or attributes_level > 2):
      if (not isinstance(value, str)):
        # Python 2.2 workaround
        if (name in ["optional", "multiple", "disable_add", "disable_delete"]):
          if   (value is False): value = "False"
          elif (value is True):  value = "True"
        print >> out, prefix+"  ."+name, "=", value
      else:
        indent = " " * (len(prefix) + 3 + len(name) + 3)
        fits_on_one_line = len(indent+value) < print_width
        if (not is_standard_identifier(value) or not fits_on_one_line):
          value = str(tokenizer.word(value=value, quote_token='"'))
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

class object_locator(object):

  def __init__(self, parent, path, object):
    self.parent = parent
    self.path = path
    self.object = object

  def __str__(self):
    return "%s%s" % (self.path, self.object.where_str)

# is_template (set by .fetch() and .format() methods of definition or scope):
#   0: not a template
#  -1: template but there are other copies
#   1: template and there are no copies

class try_tokenize_proxy(object):

  def __init__(self, error_message, tokenized):
    self.error_message = error_message
    self.tokenized = tokenized

class try_extract_proxy(object):

  def __init__(self, error_message, extracted):
    self.error_message = error_message
    self.extracted = extracted

class try_format_proxy(object):

  def __init__(self, error_message, formatted):
    self.error_message = error_message
    self.formatted = formatted

class definition(slots_getstate_setstate):

  is_definition = True
  is_scope = False

  attribute_names = [
    "help", "caption", "short_caption", "optional",
    "type", "multiple", "input_size", "style", "expert_level", "deprecated"]

  __slots__ = ["name", "words", "primary_id", "primary_parent_scope",
               "is_disabled", "is_template", "where_str", "merge_names",
               "tmp"] + attribute_names

  def __init__(self,
        name,
        words,
        primary_id=None,
        primary_parent_scope=None,
        is_disabled=False,
        is_template=0,
        where_str="",
        merge_names=False,
        tmp=None,
        help=None,
        caption=None,
        short_caption=None,
        optional=None,
        type=None,
        multiple=None,
        input_size=None,
        style=None,
        expert_level=None,
        deprecated=None) :
    if (is_reserved_identifier(name)):
      raise RuntimeError('Reserved identifier: "%s"%s' % (name, where_str))
    if (name != "include" and "include" in name.split(".")):
      raise RuntimeError('Reserved identifier: "include"%s' % where_str)
    self.name = name
    self.words = words
    self.primary_id = primary_id
    self.primary_parent_scope = primary_parent_scope
    self.is_disabled = is_disabled
    self.is_template = is_template
    self.where_str = where_str
    self.merge_names = merge_names
    self.tmp = tmp
    self.help = help
    self.caption = caption
    self.short_caption = short_caption
    self.optional = optional
    self.type = type
    self.multiple = multiple
    self.input_size = input_size
    self.style = style
    self.expert_level = expert_level
    self.deprecated = deprecated

  def __setstate__ (self, *args, **kwds) :
    slots_getstate_setstate.__setstate__(self, *args, **kwds)
    # XXX backwards compatibility 2012-03-27
    if (not hasattr(self, "deprecated")) : setattr(self, "deprecated", None)

  def copy(self):
    keyword_args = {}
    for keyword in self.__slots__:
      keyword_args[keyword] = getattr(self, keyword)
    return definition(**keyword_args)

  def customized_copy(self, name=None, words=None):
    result = self.copy()
    if (name is not None): result.name = name
    if (words is not None): result.words = words
    result.is_template = 0
    return result

  def full_path(self):
    return full_path(self)

  def assign_tmp(self, value, active_only=False):
    if (not active_only or not self.is_disabled):
      self.tmp = value

  def fetch_value(self, source, diff_mode=False,
      skip_incompatible_objects=False):
    if (source.is_scope):
      if (skip_incompatible_objects) :
        return self.copy()
      raise RuntimeError(
        'Incompatible parameter objects: definition "%s"%s vs. scope "%s"%s' %
          (self.name, self.where_str, source.name, source.where_str))
    source.tmp = True
    source = source.resolve_variables(diff_mode=diff_mode)
    type_fetch = getattr(self.type, "fetch", None)
    if (self.deprecated) :
      # issue warning if value is not the default, otherwise return None so
      # this parameter stays invisible to users
      result_as_str = strings_from_words(source.words)
      self_as_str = strings_from_words(self.words)
      if (result_as_str != self_as_str) :
        warnings.warn("%s is deprecated - not recommended for use." % \
          self.full_path(), PhilDeprecationWarning)
      else :
        return None
    if (type_fetch is None):
      return self.customized_copy(words=source.words)
    return type_fetch(source_words=source.words, master=self)

  def fetch_diff(self, source, skip_incompatible_objects=False):
    result = self.fetch_value(source=source, diff_mode=True,
      skip_incompatible_objects=skip_incompatible_objects)
    result_as_str = self.extract_format(source=result).as_str()
    self_as_str = self.extract_format().as_str()
    if (result_as_str == self_as_str): result = None
    return result

  def fetch(self, source, diff=False, skip_incompatible_objects=False):
    if (diff): return self.fetch_diff(source=source,
      skip_incompatible_objects=skip_incompatible_objects)
    return self.fetch_value(source=source,
      skip_incompatible_objects=skip_incompatible_objects)

  def has_attribute_with_name(self, name):
    return name in self.attribute_names

  def assign_attribute(self, name, words, converter_registry, converter_cache):
    assert self.has_attribute_with_name(name)
    if (name in ["optional", "multiple"]):
      value = bool_from_words(words=words, path="."+name)
    elif (name == "type"):
      value = definition_converters_from_words(
        words=words,
        converter_registry=converter_registry,
        converter_cache=converter_cache)
    elif (name in ["input_size", "expert_level"]):
      value = int_from_words(words=words, path="."+name)
    else:
      value = str_from_words(words)
    setattr(self, name, value)

  def show(self,
        out=None,
        merged_names=[],
        prefix="",
        expert_level=None,
        attributes_level=0,
        print_width=None):
    if (self.is_template < 0 and attributes_level < 2): return
    elif (self.deprecated and attributes_level < 3) : return
    if (self.expert_level is not None
        and expert_level is not None
        and expert_level >= 0
        and self.expert_level > expert_level): return
    if (out is None): out = sys.stdout
    if (print_width is None): print_width = default_print_width
    if (self.is_disabled): hash = "!"
    else:                  hash = ""
    line = prefix + hash + ".".join(merged_names + [self.name])
    if (self.name != "include"): line += " ="
    indent = " " * len(line)
    if (self.deprecated) :
      print >> out, prefix + "# WARNING: deprecated parameter"
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
        expert_level=None,
        attributes_level=0,
        print_width=None):
    out = StringIO()
    self.show(
      out=out,
      prefix=prefix,
      expert_level=expert_level,
      attributes_level=attributes_level,
      print_width=print_width)
    return out.getvalue()

  def _all_definitions(self,
        suppress_multiple,
        select_tmp,
        parent,
        parent_path,
        result):
    if (suppress_multiple and self.multiple): return
    if (select_tmp is not None and not (self.tmp == select_tmp)): return
    if (self.name == "include"): return
    result.append(object_locator(
      parent=parent, path=parent_path+self.name, object=self))

  def get_without_substitution(self, path):
    if (self.is_disabled or self.name != path): return []
    return [self]

  def _type_from_words(self):
    try: return self.type.from_words
    except AttributeError:
      raise RuntimeError('.type=%s does not have a from_words method%s: %s' %
        (str(self.type), self.where_str, format_exception()))

  def try_extract(self):
    if (self.type is None):
      return try_extract_proxy(
        error_message=None,
        extracted=strings_from_words(words=self.words))
    type_from_words = self._type_from_words()
    try:
      return try_extract_proxy(
        error_message=None,
        extracted=type_from_words(self.words, master=self))
    except RuntimeError, e:
      return try_extract_proxy(error_message=str(e), extracted=None)

  def extract(self, parent=None):
    if (self.type is None):
      return strings_from_words(words=self.words)
    return self._type_from_words()(self.words, master=self)

  def format(self, python_object):
    if (self.type is None):
      words = strings_as_words(python_object=python_object)
    else:
      try: type_as_words = self.type.as_words
      except AttributeError:
        raise RuntimeError('.type=%s does not have an as_words method%s: %s' %
          (str(self.type), self.where_str, format_exception()))
      words = type_as_words(python_object=python_object, master=self)
    return self.customized_copy(words=words)

  def extract_format(self, source=None):
    if (source is None): source = self
    return self.format(python_object=source.extract())

  def try_extract_format(self):
    proxy = self.try_extract()
    if (proxy.error_message is not None):
      return try_format_proxy(
        error_message=proxy.error_message, formatted=None)
    return try_format_proxy(
      error_message=None,
      formatted=self.format(python_object=proxy.extracted))

  def try_tokenize(self, input_string, source_info=None):
    try:
      words = tokenize_value_literal(
        input_string=input_string,
        source_info=source_info)
    except RuntimeError, e:
      return try_tokenize_proxy(
        error_message=str(e),
        tokenized=None)
    if (len(words) == 0):
      words = [tokenizer.word(value="None")]
    return try_tokenize_proxy(
      error_message=None,
      tokenized=self.customized_copy(words=words))

  def _validate(self, input_string, source_info, call):
    proxy = self.try_tokenize(
      input_string=input_string, source_info=source_info)
    if (proxy.error_message is not None):
      return proxy
    return getattr(proxy.tokenized, call)()

  def validate(self, input_string, source_info=None):
    return self._validate(input_string=input_string, source_info=source_info,
      call="try_extract")

  def validate_and_format(self, input_string, source_info=None):
    return self._validate(input_string=input_string, source_info=source_info,
      call="try_extract_format")

  def unique(self):
    return self

  def resolve_variables(self, diff_mode=False):
    new_words = []
    for word in self.words:
      if (word.quote_token == "'"):
        new_words.append(word)
        continue
      substitution_proxy = variable_substitution_proxy(word)
      for fragment in substitution_proxy.fragments:
        if (not fragment.is_variable):
          fragment.result = tokenizer.word(
            value=fragment.value, quote_token='"')
          continue
        variable_words = None
        if (self.primary_parent_scope is not None):
          substitution_source = self.primary_parent_scope.lexical_get(
            path=fragment.value, stop_id=self.primary_id)
          if (substitution_source is not None):
            if (not substitution_source.is_definition):
              raise RuntimeError("Not a definition: $%s%s" % (
                fragment.value, word.where_str()))
            substitution_source.tmp = True
            variable_words = substitution_source.resolve_variables().words
        if (variable_words is None):
          if (diff_mode):
            env_var = "$"+fragment.value
          else:
            env_var = os.environ.get(fragment.value, None)
          if (env_var is not None):
            variable_words = [tokenizer.word(
              value=env_var,
              quote_token='"',
              source_info='environment: "%s"'%fragment.value)]
        if (variable_words is None):
          raise RuntimeError("Undefined variable: $%s%s" % (
            fragment.value, word.where_str()))
        if (not substitution_proxy.force_string):
          fragment.result = variable_words
        else:
          fragment.result = tokenizer.word(
            value=" ".join([word.value for word in variable_words]),
            quote_token='"')
      new_words.extend(substitution_proxy.get_new_words())
    return self.customized_copy(words=new_words)

class scope_extract_call_proxy_object(object):

  def __init__(self, where_str, expression, callable, keyword_args):
    self.where_str = where_str
    self.expression = expression
    self.callable = callable
    self.keyword_args = keyword_args

  def __str__(self):
    return self.expression

def scope_extract_call_proxy(full_path, words, cache):
  if (is_plain_none(words=words)): return None
  if (is_plain_auto(words=words)): return Auto
  call_expression_raw = str_from_words(words).strip()
  try:
    call_expression = normalize_call_expression(expression=call_expression_raw)
  except python_tokenize.TokenError, e:
    raise RuntimeError('scope "%s" .call=%s: %s%s' % (
      full_path, call_expression_raw, str(e), words[0].where_str()))
  call_proxy = cache.get(call_expression, None)
  if (call_proxy is None):
    where_str = words[0].where_str()
    flds = call_expression.split("(", 1)
    import_path = flds[0]
    if (len(flds) == 1):
      keyword_args = {}
    else:
      extractor = "__extract_args__(" + flds[1]
      try:
        args, keyword_args = eval(
          extractor, math.__dict__, {"__extract_args__": extract_args})
      except KeyboardInterrupt: raise
      except Exception:
        raise RuntimeError('scope "%s" .call=%s: %s%s' % (
          full_path, call_expression, format_exception(), where_str))
    imported = import_python_object(
      import_path=import_path,
      error_prefix='scope "%s" .call: ' % full_path,
      target_must_be="; target must be a callable Python object",
      where_str=where_str)
    if (not callable(imported.object)):
      raise TypeError(
        'scope "%s" .call: "%s" is not a callable Python object%s' % (
          full_path, import_path, where_str))
    call_proxy = scope_extract_call_proxy_object(
      where_str=where_str,
      expression=call_expression,
      callable=imported.object,
      keyword_args=keyword_args)
    cache[call_expression] = call_proxy
  return call_proxy

class scope_extract_attribute_error(object): pass
class scope_extract_is_disabled(object): pass

class scope_extract_list(list):

  def __init__(self, optional):
    self.__phil_optional__ = optional
    list.__init__(self)

class scope_extract(object):

  def __init__(self, name, parent, call):
    object.__setattr__(self, "__phil_name__", name)
    object.__setattr__(self, "__phil_parent__", parent)
    object.__setattr__(self, "__phil_call__", call)

  def __phil_path__(self, object_name=None):
    if (   self.__phil_parent__ is None
        or self.__phil_parent__.__phil_name__ is None
        or self.__phil_parent__.__phil_name__ == ""):
      if (object_name is None):
        return self.__phil_name__
      elif (   self.__phil_name__ is None
            or self.__phil_name__ == ""):
        return object_name
      return self.__phil_name__ + "." + object_name
    result = [
      self.__phil_parent__.__phil_path__(),
      self.__phil_name__]
    if (object_name is not None):
      result.append(object_name)
    return ".".join(result)

  def __phil_path_and_value__(self, object_name):
    return (
      self.__phil_path__(object_name=object_name),
      getattr(self, object_name))

  def __setattr__(self, name, value):
    if (getattr(self, name, scope_extract_attribute_error)
          is scope_extract_attribute_error):
      pp = self.__phil_path__()
      if (pp == ""): pp = name
      else:          pp += "." + name
      raise AttributeError(
        'Assignment to non-existing attribute "%s"\n' % pp
          + '  Please correct the attribute name, or to create\n'
          + '  a new attribute use: obj.__inject__(name, value)')
    object.__setattr__(self, name, value)

  def __inject__(self, name, value):
    if (getattr(self, name, scope_extract_attribute_error)
          is not scope_extract_attribute_error):
      pp = self.__phil_path__()
      if (pp == ""): pp = name
      else:          pp += "." + name
      raise AttributeError(
        'Attribute "%s" exists already.' % pp)
    object.__setattr__(self, name, value)

  def __phil_join__(self, other):
    for key,other_value in other.__dict__.items():
      if (is_reserved_identifier(key)): continue
      self_value = self.__dict__.get(key, None)
      if (self_value is None):
        self.__dict__[key] = other_value
      elif (isinstance(self_value, scope_extract_list)):
        assert isinstance(other_value, scope_extract_list)
        for item in other_value:
          if (item is not None):
            self_value.append(item)
        if (len(self_value) > 1 and self_value[0] is None):
          del self_value[0]
      else:
        self_value_phil_join = getattr(self_value, "__phil_join__", None)
        if (self_value_phil_join is None):
          self.__dict__[key] = other_value
        else:
          self_value_phil_join(other_value)

  def __phil_set__(self, name, optional, multiple, value):
    assert not "." in name
    node = getattr(self, name, scope_extract_attribute_error)
    if (not multiple):
      if (value is scope_extract_is_disabled):
        value = None
      if (node is scope_extract_attribute_error
          or not isinstance(value, scope_extract)
          or not isinstance(node, scope_extract)):
        object.__setattr__(self, name, value)
      else:
        node.__phil_join__(value)
    else:
      if (node is scope_extract_attribute_error):
        node = scope_extract_list(optional=optional)
        object.__setattr__(self, name, node)
      if (not value is scope_extract_is_disabled
          and (value is not None or optional is not True)):
        node.append(value)

  def __phil_get__(self, name):
    assert not "." in name
    return getattr(self, name, scope_extract_attribute_error)

  def __call__(self, **keyword_args):
    call_proxy = self.__phil_call__
    if (call_proxy is None):
      raise RuntimeError('scope "%s" is not callable.' % self.__phil_path__())
    if (len(keyword_args) == 0):
      return call_proxy.callable(self, **call_proxy.keyword_args)
    effective_keyword_args = dict(call_proxy.keyword_args)
    effective_keyword_args.update(keyword_args)
    try:
      return call_proxy.callable(self, **effective_keyword_args)
    except KeyboardInterrupt: raise
    except Exception:
      raise RuntimeError('scope "%s" .call=%s execution: %s%s' % (
        self.__phil_path__(), call_proxy.expression, format_exception(),
        call_proxy.where_str))

class scope(slots_getstate_setstate):

  is_definition = False
  is_scope = True
  deprecated = False

  attribute_names = [
    "style",
    "help",
    "caption",
    "short_caption",
    "optional",
    "call",
    "multiple",
    "sequential_format",
    "disable_add",
    "disable_delete",
    "expert_level"]

  __slots__ = [
    "name",
    "objects",
    "primary_id",
    "primary_parent_scope",
    "is_disabled",
    "is_template",
    "where_str",
    "merge_names"] + attribute_names

  def __init__(self,
        name,
        objects=None,
        primary_id=None,
        primary_parent_scope=None,
        is_disabled=False,
        is_template=0,
        where_str="",
        merge_names=False,
        style=None,
        help=None,
        caption=None,
        short_caption=None,
        optional=None,
        call=None,
        multiple=None,
        sequential_format=None,
        disable_add=None,
        disable_delete=None,
        expert_level=None):
    self.name = name
    self.objects = objects
    self.primary_id = primary_id
    self.primary_parent_scope = primary_parent_scope
    self.is_disabled = is_disabled
    self.is_template = is_template
    self.where_str = where_str
    self.merge_names = merge_names
    self.style = style
    self.help = help
    self.caption = caption
    self.short_caption = short_caption
    self.optional = optional
    self.call = call
    self.multiple = multiple
    self.sequential_format = sequential_format
    self.disable_add = disable_add
    self.disable_delete = disable_delete
    self.expert_level = expert_level
    if (objects is None):
      self.objects = []
    if (is_reserved_identifier(name)):
      raise RuntimeError('Reserved identifier: "%s"%s' % (name, where_str))
    if ("include" in name.split(".")):
      raise RuntimeError('Reserved identifier: "include"%s' % where_str)
    if (sequential_format is not None):
      assert isinstance(sequential_format % 0, str)

  def copy(self):
    keyword_args = {}
    for keyword in self.__slots__:
      keyword_args[keyword] = getattr(self, keyword)
    return scope(**keyword_args)

  def customized_copy(self, name=None, objects=None):
    result = self.copy()
    if (name is not None): result.name = name
    if (objects is not None): result.objects = objects
    result.is_template = 0
    return result

  def is_empty(self):
    return len(self.objects) == 0

  def full_path(self):
    return full_path(self)

  def assign_tmp(self, value, active_only=False):
    if (not active_only):
      for object in self.objects:
        object.assign_tmp(value=value)
    else:
      for object in self.objects:
        if (self.is_disabled): continue
        object.assign_tmp(value=value, active_only=True)

  def adopt(self, object):
    assert len(object.name) > 0
    primary_parent_scope = self
    name_components = object.name.split(".")
    merge_names = False
    for name in name_components[:-1]:
      child_scope = scope(name=name)
      child_scope.merge_names = merge_names
      primary_parent_scope.adopt(child_scope)
      primary_parent_scope = child_scope
      merge_names = True
    if (len(name_components) > 1):
      object.name = name_components[-1]
      object.merge_names = True
    object.primary_parent_scope = primary_parent_scope
    primary_parent_scope.objects.append(object)

  def change_primary_parent_scope(self, new_value):
    objects = []
    for object in self.objects:
      obj = object.copy()
      obj.primary_parent_scope = new_value
      if obj.is_scope:
        obj = obj.change_primary_parent_scope(obj)
      objects.append(obj)
    return self.customized_copy(objects=objects)

  def has_attribute_with_name(self, name):
    return name in self.attribute_names

  def assign_attribute(self, name, words, scope_extract_call_proxy_cache):
    assert self.has_attribute_with_name(name)
    if (name in ["optional", "multiple", "disable_add", "disable_delete"]):
      value = bool_from_words(words, path="."+name)
    elif (name == "expert_level"):
      value = int_from_words(words=words, path="."+name)
    elif (name == "call"):
      value = scope_extract_call_proxy(
        full_path=self.full_path(),
        words=words,
        cache=scope_extract_call_proxy_cache)
    else:
      value = str_from_words(words)
      if (name == "style"):
        style = value
      elif (name == "sequential_format"):
        sequential_format = value
        if (sequential_format is not None):
          assert isinstance(sequential_format % 0, str)
    setattr(self, name, value)

  def active_objects(self):
    for object in self.objects:
      if (object.is_disabled): continue
      yield object

  def master_active_objects(self):
    names_object = {}
    for object in self.objects:
      if (object.is_disabled): continue
      master = names_object.setdefault(object.name, object)
      if (master is not object):
        if (master.multiple): continue
        if (object.is_definition):
          raise RuntimeError(
            "Duplicate definitions in master"
            " (first not marked with .multiple=True):\n"
            "  %s%s\n"
            "  %s%s" % (
              master.full_path(), master.where_str,
              object.full_path(), object.where_str))
      yield object

  def show(self,
        out=None,
        merged_names=[],
        prefix="",
        expert_level=None,
        attributes_level=0,
        print_width=None):
    if (self.is_template < 0 and attributes_level < 2): return
    if (self.expert_level is not None
        and expert_level is not None
        and expert_level >= 0
        and self.expert_level > expert_level): return
    if (out is None): out = sys.stdout
    if (print_width is None): print_width = default_print_width
    is_proper_scope = False
    if (len(self.name) == 0):
      assert len(merged_names) == 0
    elif (len(self.objects) > 0 and self.objects[0].merge_names):
      merged_names = merged_names + [self.name]
    else:
      is_proper_scope = True
      if (self.is_disabled): hash = "!"
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
        expert_level=expert_level,
        attributes_level=attributes_level,
        print_width=print_width)
    if (is_proper_scope):
      print >> out, prefix[:-2] + "}"

  def as_str(self,
        prefix="",
        expert_level=None,
        attributes_level=0,
        print_width=None):
    out = StringIO()
    self.show(
      out=out,
      prefix=prefix,
      expert_level=expert_level,
      attributes_level=attributes_level,
      print_width=print_width)
    return out.getvalue()

  def _all_definitions(self,
        suppress_multiple,
        select_tmp,
        parent,
        parent_path,
        result):
    parent_path += self.name+"."
    for object in self.active_objects():
      if (suppress_multiple and object.multiple): continue
      object._all_definitions(
        suppress_multiple=suppress_multiple,
        select_tmp=select_tmp,
        parent=self,
        parent_path=parent_path,
        result=result)

  def all_definitions(self, suppress_multiple=False, select_tmp=None):
    result = []
    for object in self.active_objects():
      if (suppress_multiple and object.multiple): continue
      object._all_definitions(
        suppress_multiple=suppress_multiple,
        select_tmp=select_tmp,
        parent=self,
        parent_path="",
        result=result)
    return result

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
    return self.customized_copy(objects=result)

  def lexical_get(self, path, stop_id, search_up=True):
    if (path.startswith(".")):
      while (self.primary_parent_scope is not None):
        self = self.primary_parent_scope
      path = path[1:]
    candidates = []
    for object in self.objects:
      if (object.primary_id >= stop_id): break
      if (object.is_definition):
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

  def extract(self, parent=None):
    result = scope_extract(name=self.name, parent=parent, call=self.call)
    for object in self.objects:
      if (object.is_template < 0): continue
      if (object.is_disabled or object.is_template > 0):
        value = scope_extract_is_disabled
      else:
        value = object.extract(parent=result)
      result.__phil_set__(
        name=object.name,
        optional=object.optional,
        multiple=object.multiple,
        value=value)
    return result

  def format(self, python_object):
    multiple_scopes_done = {}
    result = []
    for object in self.master_active_objects():
      if (object.multiple and object.is_scope):
        if (object.name in multiple_scopes_done): continue
        multiple_scopes_done[object.name] = False
      if (python_object is None):
        result.append(object.format(None))
      elif (python_object is Auto):
        result.append(object.format(Auto))
      else:
        if (isinstance(python_object, scope_extract)):
          python_object = [python_object]
        for python_object_i in python_object:
          sub_python_object = python_object_i.__phil_get__(object.name)
          if (sub_python_object is not scope_extract_attribute_error):
            if (not object.multiple):
              result.append(object.format(sub_python_object))
            else:
              if (len(sub_python_object) == 0):
                obj = object.copy()
                obj.is_template = 1
                result.append(obj)
              else:
                if (not multiple_scopes_done.get(object.name, True)):
                  multiple_scopes_done[object.name] = True
                  obj = object.copy()
                  obj.is_template = -1
                  result.append(obj)
                for sub_python_object_i in sub_python_object:
                  result.append(object.format(sub_python_object_i))
    return self.customized_copy(objects=result)

  def extract_format(self, source=None):
    if (source is None): source = self
    return self.format(source.extract())

  def clone(self, python_object, converter_registry=None):
    return parse(
      input_string=self.format(python_object=python_object)
        .as_str(attributes_level=3),
      converter_registry=converter_registry).extract()

  def fetch(self,
        source=None,
        sources=None,
        track_unused_definitions=False,
        diff=False,
        skip_incompatible_objects=False,
        warn_deprecated=True):
    combined_objects = []
    if (source is not None or sources is not None):
      assert source is None or sources is None
      combined_objects = []
      if (sources is None): sources = [source]
      for source in sources:
        assert source.name == self.name
        if (source.is_definition):
          if (skip_incompatible_objects) :
            continue
          raise RuntimeError(
            'Incompatible parameter objects:'
            ' scope "%s"%s vs. definition "%s"%s' %
              (self.name, self.where_str, source.name, source.where_str))
        combined_objects.extend(source.objects)
    source = self.customized_copy(objects=combined_objects)
    del sources
    if (track_unused_definitions):
      source.assign_tmp(value=False, active_only=True)
    result_objects = []
    for master_object in self.master_active_objects():
      if (len(self.name) == 0):
        path = master_object.name
      else:
        path = self.name + "." + master_object.name
      matching_sources = source.get(path=path, with_substitution=False)
      if (not master_object.multiple):
        if (master_object.is_definition):
          # loop over all matching_sources to support track_unused_definitions
          result_object = None
          for matching_source in matching_sources.active_objects():
            result_object = master_object.fetch(
              source=matching_source,
              diff=diff,
              skip_incompatible_objects=skip_incompatible_objects)
        else:
          result_object = master_object.fetch(
            sources=matching_sources.active_objects(),
            diff=diff,
            skip_incompatible_objects=skip_incompatible_objects)
          if (diff and len(result_object.objects) == 0):
            result_object = None
        if (result_object is not None):
          result_objects.append(result_object)
        elif (not diff) and (not master_object.deprecated) :
          result_objects.append(master_object.copy())
      else:
        processed_as_str = {}
        result_objs = []
        master_as_str = master_object.extract_format().as_str()
        for from_master,matching in [
              (True, self.get(path=path, with_substitution=False)),
              (False, matching_sources)]:
          for matching_source in matching.active_objects():
            if (matching_source is master_object): continue
            candidate = master_object.fetch(source=matching_source, diff=diff,
              skip_incompatible_objects=skip_incompatible_objects)
            if (diff):
              if (master_object.is_scope):
                if (len(candidate.objects) == 0): continue
              elif (candidate is None):
                continue
            candidate_as_str = master_object.extract_format(
              source=candidate).as_str()
            if (candidate_as_str == master_as_str): continue
            prev_index = processed_as_str.get(candidate_as_str, None)
            if (prev_index is not None):
              if (prev_index == -1): continue
              result_objs[prev_index] = None
            if (diff and from_master):
              processed_as_str[candidate_as_str] = -1
            else:
              processed_as_str[candidate_as_str] = len(result_objs)
              result_objs.append(candidate)
        if (not diff):
          obj = master_object.copy()
          if (    not master_object.optional is None
              and not master_object.optional):
            obj.is_template = 0
          elif (len(processed_as_str) == 0):
            obj.is_template = 1
          else:
            obj.is_template = -1
          result_objects.append(obj)
        del processed_as_str
        for obj in result_objs:
          if (obj is not None):
            result_objects.append(obj)
        del result_objs
    result = self.customized_copy(objects=result_objects)
    if (track_unused_definitions):
      return result, source.all_definitions(select_tmp=False)
    return result

  def fetch_diff(self,
        source=None,
        sources=None,
        track_unused_definitions=False,
        skip_incompatible_objects=False):
    return self.fetch(
      source=source,
      sources=sources,
      track_unused_definitions=track_unused_definitions,
      diff=True,
      skip_incompatible_objects=skip_incompatible_objects)

  def process_includes(self,
        converter_registry,
        reference_directory,
        include_stack=None):
    if (converter_registry is None):
      converter_registry = default_converter_registry
    if (include_stack is None): include_stack = []
    result = []
    for object in self.objects:
      if (object.is_disabled):
        result.append(object)
      elif (object.is_definition):
        if (object.name != "include"):
          result.append(object)
        else:
          object_sub = object.resolve_variables()
          if (len(object_sub.words) < 2):
            raise RuntimeError(
              '"include" must be followed by at least two arguments%s' % (
                object.where_str))
          include_type = object_sub.words[0].value.lower()
          if (include_type == "file"):
            if (len(object_sub.words) != 2):
              raise RuntimeError(
                '"include file" must be followed exactly one argument%s' % (
                  object.where_str))
            file_name = object_sub.words[1].value
            if (reference_directory is not None
                and not os.path.isabs(file_name)):
              file_name = os.path.join(reference_directory, file_name)
            result.extend(parse(
              file_name=file_name,
              converter_registry=converter_registry,
              process_includes=True,
              include_stack=include_stack).objects)
          elif (include_type == "scope"):
            if (len(object_sub.words) > 3):
              raise RuntimeError(
                '"include scope" must be followed one or two arguments,'
                ' i.e. an import path and optionally a phil path%s' % (
                  object.where_str))
            import_path = object_sub.words[1].value
            if (len(object_sub.words) > 2):
              phil_path = object_sub.words[2].value
            else:
              phil_path = None
            result.extend(process_include_scope(
              converter_registry=converter_registry,
              include_stack=include_stack,
              object=object,
              import_path=import_path,
              phil_path=phil_path).objects)
          else:
            raise RuntimeError("Unknown include type: %s%s" % (
              include_type, object.where_str))
      else:
        result.append(object.process_includes(
          converter_registry=converter_registry,
          reference_directory=reference_directory,
          include_stack=include_stack))
    return self.customized_copy(objects=result)

  def unique(self):
    selection = {}
    result = []
    for i_object,object in enumerate(self.active_objects()):
      selection[object.name] = i_object
    for i_object,object in enumerate(self.active_objects()):
      if (selection[object.name] == i_object):
        result.append(object.unique())
    return self.customized_copy(objects=result)

  def command_line_argument_interpreter(self,
        home_scope=None,
        argument_description=None):
    from libtbx.phil.command_line import argument_interpreter as _
    return _(
      master_phil=self,
      home_scope=home_scope,
      argument_description=argument_description)

def process_include_scope(
      converter_registry,
      include_stack,
      object,
      import_path,
      phil_path):
  imported = import_python_object(
    import_path=import_path,
    error_prefix="include scope: ",
    target_must_be="; target must be a phil scope object or phil string",
    where_str=object.where_str)
  source_scope = imported.object
  if (isinstance(source_scope, str)):
    source_scope = parse(
      input_string=source_scope,
      converter_registry=converter_registry)
  elif (source_scope is None or not isinstance(source_scope, scope)):
    raise RuntimeError(
      'include scope: python object "%s" in module "%s" is not a'
      ' libtbx.phil.scope instance%s' % (
        imported.path_elements[-1], imported.module_path, object.where_str))
  source_scope = source_scope.process_includes(
    converter_registry=converter_registry,
    reference_directory=None,
    include_stack=include_stack)
  if (phil_path is None):
    result = source_scope
  else:
    result = source_scope.get(path=phil_path)
    if (len(result.objects) == 0):
      raise RuntimeError(
        'include scope: path "%s" not found in phil scope object "%s"' \
        ' in module "%s"%s' % (
          phil_path, imported.path_elements[-1], imported.module_path,
          object.where_str))
  return result.change_primary_parent_scope(object.primary_parent_scope)

class variable_substitution_fragment(slots_getstate_setstate):

  __slots__ = ["is_variable", "value", "result"]

  def __init__(self, is_variable, value):
    self.is_variable = is_variable
    self.value = value

class variable_substitution_proxy(slots_getstate_setstate):

  __slots__ = ["word", "force_string", "have_variables", "fragments"]

  def __init__(self, word):
    self.word = word
    self.force_string = word.quote_token is not None
    self.have_variables = False
    self.fragments = []
    fragment_value = ""
    char_iter = tokenizer.character_iterator(word.value)
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
        if (c == "("):
          while True:
            c = char_iter.next()
            if (c is None):
              word.raise_syntax_error('missing ")": ')
            if (c == ")"):
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
    return [tokenizer.word(
      value="".join([fragment.result.value for fragment in self.fragments]),
      quote_token='"')]

def parse(
      input_string=None,
      source_info=None,
      file_name=None,
      converter_registry=None,
      process_includes=False,
      include_stack=None):
  from libtbx.phil import parser
  assert source_info is None or file_name is None
  if (input_string is None):
    assert file_name is not None
    input_string = open(file_name).read()
  if (converter_registry is None):
    converter_registry = default_converter_registry
  result = scope(name="", primary_id=0)
  parser.collect_objects(
    word_iterator=tokenizer.word_iterator(
      input_string=input_string,
      source_info=source_info,
      file_name=file_name,
      list_of_settings=[
        tokenizer.settings(
          unquoted_single_character_words="{}=",
          contiguous_word_characters="",
          comment_characters="#",
          meta_comment="phil"),
        tokenizer.settings(
          unquoted_single_character_words="{};",
          contiguous_word_characters="")]),
    converter_registry=converter_registry,
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
      converter_registry=converter_registry,
      reference_directory=reference_directory,
      include_stack=include_stack)
    if (include_stack is not None):
      include_stack.pop()
  return result

def read_default(
      caller_file_name,
      params_extension=".params",
      converter_registry=None,
      process_includes=True):
  params_file_name = os.path.splitext(caller_file_name)[0] + params_extension
  if (not os.path.isfile(params_file_name)):
    raise RuntimeError("Missing parameter file: %s" % params_file_name)
  return parse(
    file_name=params_file_name,
    converter_registry=converter_registry,
    process_includes=process_includes)

def process_command_line(args, master_string, parse=None):
  from libtbx.phil import command_line
  return command_line.process(
    args=args, master_string=master_string, parse=parse)
