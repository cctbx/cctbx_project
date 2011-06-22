class try_read_file(object):

  def __init__(O, file_name, input_types):
    def O_assign(file_type, file_content):
      O.file_name = file_name
      O.file_type = file_type
      O.file_content = file_content
    if ("directory" in input_types):
      import os.path as op
      if (op.isdir(file_name)):
        O_assign(file_type="directory", file_content=None)
        return
    if ("mtz" in input_types):
      lead = open(file_name, "rb").read(3)
      if (lead == "MTZ"):
        import iotbx.mtz
        mtz_obj = iotbx.mtz.object(file_name=file_name)
        O_assign(file_type="mtz", file_content=mtz_obj)
        return
    if ("pdb" in input_types):
      import iotbx.pdb
      try:
        pdb_inp = iotbx.pdb.input(file_name=file_name)
      except KeyboardInterrupt: raise
      except Exception:
        if (iotbx.pdb.is_pdb_file(file_name=file_name)):
          raise
        pdb_inp = None
      else:
        if (pdb_inp.atoms().size() != 0):
          O_assign(file_type="pdb", file_content=pdb_inp)
          return
    if ("phil" in input_types):
      import iotbx.phil
      try:
        phil_obj = iotbx.phil.parse(file_name=file_name)
      except KeyboardInterrupt: raise
      except Exception: pass
      else:
        O_assign(file_type="phil", file_content=phil_obj)
        return
    if ("cif" in input_types):
      import mmtbx.monomer_library.server
      try:
        cif_obj = mmtbx.monomer_library.server.read_cif(file_name=file_name)
      except KeyboardInterrupt: raise
      except Exception: pass
      else:
        if (len(cif_obj) != 0):
          O_assign(file_type="cif", file_content=cif_obj)
          return
    if ("pdb" in input_types and pdb_inp is not None):
      if (pdb_inp.unknown_section().size() != 0):
        O_assign(file_type=None, file_content=None)
        return
      if (pdb_inp.header_section().size() != 0):
        O_assign(file_type="pdb", file_content=pdb_inp)
    O_assign(file_type=None, file_content=None)

def process_command_line_inputs(args, master_phil, input_types):
  assert set(("directory", "mtz", "pdb", "cif")).issuperset(set(input_types))
  input_objects = {}
  for key in input_types:
    input_objects[key] = []
  if (master_phil is None):
    argument_interpreter = None
  else:
    argument_interpreter = master_phil.command_line_argument_interpreter()
    input_objects["phil"] = []
  from libtbx.utils import Sorry
  from libtbx.str_utils import show_string
  import os.path as op
  for arg in args:
    if (len(arg) == 0): continue
    def try_as_file():
      if (not op.isfile(arg)): return False
      obj = try_read_file(file_name=arg, input_types=input_objects)
      if (obj.file_type is None): return False
      if (obj.file_type == "phil"):
        input_objects["phil"].append(obj.file_content)
      else:
        input_objects[obj.file_type].append(obj)
      return True
    def try_as_command_line_params():
      try: command_line_params = argument_interpreter.process(arg=arg)
      except KeyboardInterrupt: raise
      except Exception:
        if (op.isfile(arg)):
          raise Sorry(
            "Error processing file: %s" % show_string(arg))
        raise Sorry(
          "Command-line argument not recognized: %s" % show_string(arg))
      input_objects["phil"].append(command_line_params)
    if (not try_as_file() and argument_interpreter is not None):
      try_as_command_line_params()
  return input_objects
