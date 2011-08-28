from cctbx import sgtbx
from cctbx import uctbx
import libtbx.phil.command_line
from libtbx.utils import Sorry, import_python_object
from libtbx.phil import tokenizer
from libtbx import Auto
import sys, os

class unit_cell_converters(object):

  phil_type = "unit_cell"

  def __str__(self): return self.phil_type

  def from_words(self, words, master):
    s = libtbx.phil.str_from_words(words=words)
    if (s is None): return None
    if (s is Auto): return Auto
    return uctbx.unit_cell(s)

  def as_words(self, python_object, master):
    if (python_object is None):
      return [tokenizer.word(value="None")]
    if (python_object is Auto):
      return [tokenizer.word(value="Auto")]
    return [tokenizer.word(value="%.10g" % v)
      for v in python_object.parameters()]

class space_group_converters(object):

  phil_type = "space_group"

  def __str__(self): return self.phil_type

  def from_words(self, words, master):
    symbol = libtbx.phil.str_from_words(words)
    if (symbol is None): return None
    if (symbol is Auto): return Auto
    return sgtbx.space_group_info(symbol=symbol)

  def as_words(self, python_object, master):
    if (python_object is None):
      return [tokenizer.word(value="None")]
    if (python_object is Auto):
      return [tokenizer.word(value="Auto")]
    return [tokenizer.word(value=str(python_object), quote_token='"')]

class atom_selection_converters(libtbx.phil.qstr_converters):

  phil_type = "atom_selection"

  def __str__(self): return self.phil_type

  def from_words(self, words, master):
    if (len(words) == 1):
      word = words[0]
      if (word.quote_token is not None):
        return word.value # mainly for backward compatibility
    return libtbx.phil.qstr_converters.from_words(self, words, master)

default_converter_registry = libtbx.phil.extended_converter_registry(
  additional_converters=[
    unit_cell_converters,
    space_group_converters,
    atom_selection_converters])

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

def process_command_line(args, master_string, parse=None):
  if (parse is None):
    import iotbx.phil
    parse = iotbx.phil.parse
  return libtbx.phil.process_command_line(
    args=args, master_string=master_string, parse=parse)

class process_command_line_with_files (object) :
  def __init__ (self,
                args,
                master_phil=None,
                master_phil_string=None,
                pdb_file_def=None,
                reflection_file_def=None,
                cif_file_def=None,
                seq_file_def=None) :
    assert (master_phil is not None) or (master_phil_string is not None)
    if (master_phil_string is not None) :
      assert (master_phil is None)
      import iotbx.phil
      master_phil = iotbx.phil.parse(input_string=master_phil_string,
        process_includes=True)
    self.master = master_phil
    self.pdb_file_def = pdb_file_def
    self.reflection_file_def = reflection_file_def
    self.cif_file_def = cif_file_def
    self.seq_file_def = seq_file_def
    cai=libtbx.phil.command_line.argument_interpreter(master_phil=self.master)
    self.work = cai.process_and_fetch(
       args=args,
       custom_processor=self)

  def __call__ (self, arg) :
    if (os.path.isfile(arg)) :
      from iotbx import file_reader
      f = file_reader.any_file(os.path.abspath(arg))
      file_def_name = None
      if (f.file_type == "pdb") and (self.pdb_file_def is not None) :
        file_def_name = self.pdb_file_def
      elif (f.file_type == "hkl") and (self.reflection_file_def is not None) :
        file_def_name = self.reflection_file_def
      elif (f.file_type == "cif") and (self.cif_file_def is not None) :
        file_def_name = self.cif_file_def
      elif (f.file_type == "seq") and (self.seq_file_def is not None) :
        file_def_name = self.seq_file_def
      if (file_def_name is not None) :
        return libtbx.phil.parse("%s=%s" % (file_def_name, f.file_name))
      else :
        return False
    else :
      return self.process_non_file(arg)

  def process_non_file (self, arg) :
    return False

# Utilities for Phenix GUI
class setup_app_generic (object) :
  def __init__ (self, master_phil_path) :
    container = import_python_object(
      import_path=master_phil_path,
      error_prefix="",
      target_must_be="",
      where_str="")
    master_phil = container.object
    if master_phil is None :
      raise Sorry("Couldn't start program using specified phil object (%s)!" %
        master_phil_path)
    if hasattr(master_phil, "__call__") :
      master_phil = master_phil()
    elif isinstance(master_phil, str) :
      master_phil = parse(master_phil, process_includes=True)
    else :
      assert type(master_phil).__name__ == "scope"
    self.master_phil = master_phil

  def __call__ (self, args) :
    (working_phil, options, unused_args) = parse_command_line_phil_args(
      args=args,
      master_phil=self.master_phil,
      command_name="phenix",
      usage_opts=["[model.pdb]", "[data.mtz]"],
      app_options=None,
      home_scope="")
    return (self.master_phil,working_phil,options, unused_args)

# TODO probably redundant, replace with process_command_line or similar?
def parse_command_line_phil_args (args, master_phil, command_name, usage_opts,
    app_options, home_scope, log=sys.stdout) :
  sources = []
  unused_args = []
  interpreter = master_phil.command_line_argument_interpreter(
    home_scope=home_scope)
  for arg in args :
    if os.path.isfile(arg) :
      try :
        user_phil = parse(file_name=arg)
        sources.append(user_phil)
      except Exception, e :
        unused_args.append(os.path.abspath(arg))
    elif arg != "" and not arg.startswith("-") :
      try :
        params = interpreter.process(arg=arg)
      except RuntimeError :
        print >> log, "%s does not appear to be a parameter definition" % arg
        unused_args.append(arg)
      else :
        sources.append(params)
    elif arg != "" :
      unused_args.append(arg)
  try :
    working_phil = master_phil.fetch(sources=sources,
       skip_incompatible_objects=True)
  except Exception, e :
    print >> log, "Error incorporating parameters from user-specified file(s):"
    print >> log, str(e)
    print >> log, "Will revert to default parameters for this run."
    working_phil = master_phil.fetch()

  assert working_phil is not None
  if "--debug" in args :
    working_phil.show()
  cmdline_opts = None
  return (working_phil, cmdline_opts, unused_args)

def load_from_cache_if_possible (phil_path) :
  import libtbx.load_env
  full_path = os.path.join(libtbx.env.build_path, "phil_cache", "%s.phil" %
    phil_path)
  if (os.path.exists(full_path)) :
    return parse(file_name=full_path)
  else :
    return import_python_object(
      import_path=phil_path,
      error_prefix="",
      target_must_be="",
      where_str="")
