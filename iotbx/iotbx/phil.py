from cctbx import sgtbx
from cctbx import uctbx
import libtbx.phil
from libtbx.phil import tokenizer

class unit_cell_converters:

  def __str__(self): return "unit_cell"

  def from_words(self, words, master):
    s = libtbx.phil.str_from_words(words=words)
    if (s is None): return None
    return uctbx.unit_cell(s)

  def as_words(self, python_object, master):
    if (python_object is None):
      return [tokenizer.word(value="None")]
    return [tokenizer.word(value="%.10g" % v)
      for v in python_object.parameters()]

class space_group_converters:

  def __str__(self): return "space_group"

  def from_words(self, words, master):
    symbol = libtbx.phil.str_from_words(words)
    if (symbol is None): return None
    return sgtbx.space_group_info(symbol=symbol)

  def as_words(self, python_object, master):
    if (python_object is None):
      return [tokenizer.word(value="None")]
    return [tokenizer.word(value=str(python_object), quote_token='"')]

default_converter_registry = dict(libtbx.phil.default_converter_registry)
for converters in [unit_cell_converters, space_group_converters]:
  default_converter_registry[str(converters())] = converters

def parse(
      input_string=None,
      source_info=None,
      file_name=None,
      converter_registry=None,
      process_includes=False):
  if (converter_registry is None):
    converter_registry = default_converter_registry
  return libtbx.phil.parse(
    input_string=input_string,
    source_info=source_info,
    file_name=file_name,
    converter_registry=converter_registry,
    process_includes=process_includes)

def read_default(
      caller_file_name,
      params_extension=".params",
      converter_registry=None,
      process_includes=True):
  if (converter_registry is None):
    converter_registry = default_converter_registry
  return libtbx.phil.read_default(
    caller_file_name=caller_file_name,
    params_extension=params_extension,
    converter_registry=converter_registry,
    process_includes=process_includes)
