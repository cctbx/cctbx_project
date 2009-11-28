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

def pdb_inp_generator(file_names, chunk_n, chunk_i):
  import iotbx.pdb
  from libtbx.utils import format_cpu_times
  import time
  print "len(file_names):", len(file_names)
  sys.stdout.flush()
  t0_total = time.time()
  try:
    for i_file,file_name in enumerate(file_names):
      if (i_file % chunk_n != chunk_i): continue
      print "i_file:", i_file, file_name
      sys.stdout.flush()
      pdb_code = pdb_code_from_file_name(file_name=file_name)
      try:
        yield pdb_code, iotbx.pdb.input(file_name=file_name)
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

def run_multi(cmd):
  from libtbx import easy_run
  print cmd
  easy_run.call(command=cmd)

def run(args, command_call, command_line_add_options=None):
  from iotbx.option_parser import option_parser as iotbx_option_parser
  import libtbx.utils
  from libtbx.utils import escape_sh_double_quoted
  from libtbx.str_utils import show_string
  show_times = libtbx.utils.show_times(time_start="now")
  command_line = (iotbx_option_parser(
    usage=" ".join(command_call) + " [options] file|directory...")
    .enable_chunk(easy_all=True)
    .option(None, "--multiprocessing",
      action="store",
      type="int",
      default=None,
      help="Use multiprocessing module.",
      metavar="INT")
    .call_with_self_as_first_argument(callable=command_line_add_options)
  ).process(args=args)
  n = command_line.options.multiprocessing
  if (n is not None and n > 1):
    if (command_line.chunk.n == 1):
      cmds = []
      for i in xrange(n):
        cmd = command_call + args + ["--chunk=%d,%d" % (n,i)]
        cmd = " ".join(['"'+escape_sh_double_quoted(s=arg)+'"' for arg in cmd])
        cmds.append(cmd)
      import multiprocessing
      mp_pool = multiprocessing.Pool(processes=n)
      mp_pool.map(run_multi, cmds)
      show_times()
      return command_line, null_generator()
    command_line.chunk.redirect_chunk_stdout_and_stderr(have_array=True)
  #
  ca = command_line.args
  if (len(ca) == 0 and pdb_mirror_pdb is not None):
    ca = [pdb_mirror_pdb]
  file_names = []
  for arg in ca:
    if (op.isfile(arg)):
      bn = op.basename(arg)
      if (bn.startswith("pdb_codes_")):
        assert pdb_mirror_pdb is not None
        for pdb_code in open(arg).read().splitlines():
          file_name = op.join(
            pdb_mirror_pdb, pdb_code[1:3], "pdb%s.ent.gz" % pdb_code)
          assert op.isfile(file_name)
          file_names.append(file_name)
      elif (bn.startswith("file_names_")):
        for file_name in open(arg).read().splitlines():
          assert op.isfile(file_name)
          file_names.append(file_name)
      else:
        file_names.append(arg)
    elif (op.isdir(arg)):
      file_name_index = op.join(arg, "INDEX")
      if (op.isfile(file_name_index)):
        for relative_path in open(file_name_index).read().splitlines():
          file_name = op.join(arg, relative_path)
          assert op.isfile(file_name)
          file_names.append(file_name)
      else:
        prev_len_file_names = len(file_names)
        for relative_path in os.listdir(arg):
          if (relative_path.endswith((".ent", ".ent.gz", "ent.Z",
                                      ".pdb", ".pdb.gz", "pdb.Z"))):
            file_name = op.join(arg, relative_path)
            assert op.isfile(file_name)
            file_names.append(file_name)
        if (len(file_names) == prev_len_file_names):
          raise RuntimeError(
            "No INDEX file and no pdb files found in directory: %s" %
              show_string(arg))
    else:
      raise RuntimeError(
        "Not a file or directory: %s" % show_string(arg))
  #
  return command_line, pdb_inp_generator(
    file_names=file_names,
    chunk_n=command_line.chunk.n,
    chunk_i=command_line.chunk.i)
