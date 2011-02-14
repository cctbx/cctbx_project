from __future__ import division
import traceback
import sys, os
op = os.path

cci_pdbmtz_path = os.environ.get("CCI_PDBMTZ")
pdb_mirror_pdb = os.environ.get("PDB_MIRROR_PDB")

assert cci_pdbmtz_path is None or op.isdir(cci_pdbmtz_path)
assert pdb_mirror_pdb is None or op.isdir(pdb_mirror_pdb)

__pdbmtz_codes = None
def pdbmtz_codes():
  global __pdbmtz_codes
  if (__pdbmtz_codes is None):
    __pdbmtz_codes = set()
    for node in os.listdir(cci_pdbmtz_path):
      if (not node.endswith(".mtz")): continue
      assert len(node) == 8
      assert node.lower() == node
      __pdbmtz_codes.add(node[:4])
  return __pdbmtz_codes

def extract_best_resolution(remark_section):
  result = None
  for line in remark_section:
    if (line.startswith("REMARK   2 ")):
      i_resolution = line.find("RESOLUTION")
      if (i_resolution > 0 and line[i_resolution+10:].find("ANGSTROM") > 0):
        flds = line.split()
        if (len(flds) >= 5):
          try: resolution = float(flds[3])
          except ValueError: pass
          else:
            if (result is None or result > resolution):
              result = resolution
  return result

def report_exception(file_name):
  sys.stderr.flush()
  print ">Begin exception"
  print "%s: Exception" % file_name
  sys.stdout.flush()
  traceback.print_exc()
  sys.stderr.flush()
  print ">End exception"
  print
  sys.stdout.flush()

def pdb_code_from_file_name(file_name):
  bn = op.basename(file_name)
  if (bn.endswith(".ent.gz")):
    return bn[-11:-7]
  if (bn.endswith(".ent")):
    return bn[-8:-4]
  return bn.split(".")[0]

class file_info(object):

  __slots__ = ["name", "atom_selection_string"]

  def __init__(O, name, atom_selection_string=None):
    if (not op.isfile(name)):
      from libtbx.str_utils import show_string
      raise RuntimeError("Missing pdb file: %s" % show_string(name))
    O.name = name
    O.atom_selection_string = atom_selection_string

class pdb_info(object):

  __slots__ = [
    "chunk_n",
    "chunk_i",
    "file_name",
    "atom_selection_string",
    "pdb_code",
    "pdb_inp"]

  def __init__(O, chunk_n, chunk_i, file_info, pdb_code, pdb_inp):
    O.chunk_n = chunk_n
    O.chunk_i = chunk_i
    O.file_name = file_info.name
    O.atom_selection_string = file_info.atom_selection_string
    O.pdb_code = pdb_code
    O.pdb_inp = pdb_inp

def pdb_inp_generator(file_infos, chunk_n, chunk_i):
  import iotbx.pdb
  from libtbx.utils import format_cpu_times
  import time
  print "len(file_infos):", len(file_infos)
  sys.stdout.flush()
  t0_total = time.time()
  try:
    for i_file,file_info in enumerate(file_infos):
      if (i_file % chunk_n != chunk_i): continue
      print "i_file:", i_file, file_info.name
      sys.stdout.flush()
      pdb_code = pdb_code_from_file_name(file_name=file_info.name)
      try:
        yield pdb_info(
          chunk_n=chunk_n,
          chunk_i=chunk_i,
          file_info=file_info,
          pdb_code=pdb_code,
          pdb_inp=iotbx.pdb.input(file_name=file_info.name))
      except KeyboardInterrupt: raise
      except Exception:
        report_exception(file_name=file_name)
  finally:
    sys.stderr.flush()
    print "total time: %.2f" % (time.time() - t0_total)
    print format_cpu_times()
    sys.stdout.flush()

def null_generator():
  for never in []: yield never

def run(args, command_call, command_line_add_options=None):
  from iotbx.option_parser import option_parser as iotbx_option_parser
  import libtbx.utils
  show_times = libtbx.utils.show_times(time_start="now")
  command_line = (iotbx_option_parser(
    usage=" ".join(command_call) + " [options] file|directory...")
    .enable_chunk(easy_all=True)
    .enable_multiprocessing()
    .call_with_self_as_first_argument(callable=command_line_add_options)
  ).process(args=args)
  if (command_line.run_multiprocessing_chunks_if_applicable(
        command_call=command_call)):
    show_times()
    return command_line, null_generator()
  #
  from libtbx.str_utils import show_string
  ca = command_line.args
  if (len(ca) == 0 and pdb_mirror_pdb is not None):
    ca = [pdb_mirror_pdb]
  file_infos = []
  for arg in ca:
    if (op.isfile(arg)):
      bn = op.basename(arg)
      if (bn.startswith("pdb_codes_")):
        assert pdb_mirror_pdb is not None
        for i_line,line in enumerate(open(arg).read().splitlines()):
          flds = line.split(None, 1)
          if (len(flds) == 0):
            raise RuntimeError(
              "Error interpreting pdb_codes file:\n"
              "  %s"
              "  line number: %d"
              "  line: %s" % (
                show_string(arg), i_line+1, show_string(line)))
          pdb_code = flds[0]
          atom_selection_string = None
          if (len(flds) > 1):
            atom_selection_string = flds[1]
          file_name = op.join(
            pdb_mirror_pdb, pdb_code[1:3], "pdb%s.ent.gz" % pdb_code)
          file_infos.append(
            file_info(
              name=file_name,
              atom_selection_string=atom_selection_string))
      elif (bn.startswith("file_names_")):
        for file_name in open(arg).read().splitlines():
          file_infos.append(file_info(name=file_name))
      else:
        file_infos.append(file_info(name=arg))
    elif (op.isdir(arg)):
      file_name_index = op.join(arg, "INDEX")
      if (op.isfile(file_name_index)):
        for relative_path in open(file_name_index).read().splitlines():
          file_name = op.join(arg, relative_path)
          file_infos.append(file_info(name=file_name))
      else:
        prev_len_file_infos = len(file_infos)
        for relative_path in os.listdir(arg):
          if (relative_path.endswith((".ent", ".ent.gz", "ent.Z",
                                      ".pdb", ".pdb.gz", "pdb.Z"))):
            file_name = op.join(arg, relative_path)
            file_infos.append(file_info(name=file_name))
        if (len(file_infos) == prev_len_file_infos):
          raise RuntimeError(
            "No INDEX file and no pdb files found in directory: %s" %
              show_string(arg))
    else:
      raise RuntimeError(
        "Not a file or directory: %s" % show_string(arg))
  #
  return command_line, pdb_inp_generator(
    file_infos=file_infos,
    chunk_n=command_line.chunk.n,
    chunk_i=command_line.chunk.i)
