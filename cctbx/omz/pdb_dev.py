from __future__ import absolute_import, division, print_function
from libtbx.utils import date_and_time, user_plus_sys_time
from libtbx.str_utils import show_string
import traceback
import os
op = os.path

def process(params, mtz_pdb_pair):
  fn_mtz, fn_pdb = mtz_pdb_pair
  #
  import iotbx.pdb
  xs = iotbx.pdb.input(file_name=fn_pdb).xray_structure_simple()
  xs.show_summary()
  xs.scattering_type_registry(table="it1992").show()
  print()
  #
  f_obs = None
  i_obs = None
  import iotbx.mtz
  mtz_obj = iotbx.mtz.object(file_name=fn_mtz)
  for ma in mtz_obj.as_miller_arrays(crystal_symmetry=xs):
    print(ma.info())
    if (f_obs is None and ma.is_xray_amplitude_array()):
      f_obs = ma
    elif (i_obs is None and ma.is_xray_intensity_array()):
      i_obs = ma
  print()
  if (i_obs is None):
    assert f_obs is not None
    c_obs = f_obs
    f_obs.show_comprehensive_summary()
    i_obs = f_obs.f_as_f_sq(algorithm="shelxl")
  else:
    c_obs = i_obs
    i_obs.show_comprehensive_summary()
    f_obs = i_obs.f_sq_as_f(algorithm="xtal_3_7")
  print()
  #
  if (len(params.optimizers) == 0):
    return
  #
  from cctbx.omz.cod_refine import process_continue
  process_continue(
    params=params,
    cod_id=op.basename(fn_pdb),
    c_obs=c_obs, i_obs=i_obs, f_obs=f_obs,
    structure_prep=xs)

def run(args):
  from iotbx.option_parser import option_parser as iotbx_option_parser
  import libtbx.utils
  show_times = libtbx.utils.show_times(time_start="now")
  command_call = ["iotbx.python", __file__]
  command_line = (iotbx_option_parser(
    usage=" ".join(command_call) + " [options] directory|file...")
    .enable_chunk(easy_all=True)
    .enable_multiprocessing()
  ).process(args=args, min_nargs=1)
  if (command_line.run_multiprocessing_chunks_if_applicable(
        command_call=command_call)):
    show_times()
    return
  co = command_line.options
  #
  print("TIME BEGIN pdb_dev:", date_and_time())
  print()
  libtbx.utils.host_and_user().show()
  print()
  sys.stdout.flush()
  #
  from cctbx.omz import cod_refine
  master_phil = cod_refine.get_master_phil(
    max_atoms=None,
    f_calc_options_algorithm="direct *fft",
    bulk_solvent_correction=True)
  argument_interpreter = master_phil.command_line_argument_interpreter()
  phil_objects = []
  remaining_args = []
  for arg in command_line.args:
    if (arg.find("=") >= 0):
      phil_objects.append(argument_interpreter.process(arg=arg))
    else:
      remaining_args.append(arg)
  work_phil = master_phil.fetch(sources=phil_objects)
  work_phil.show()
  print()
  params = work_phil.extract()
  #
  mtz_pdb_pairs = []
  arg_iter = iter(remaining_args)
  pdb_v3_mirror_dir = os.environ.get("PDB_MIRROR_PDB")
  assert pdb_v3_mirror_dir is None or op.isdir(pdb_v3_mirror_dir)
  cci_pdbmtz_path = os.environ.get("CCI_PDBMTZ")
  assert cci_pdbmtz_path is None or op.isdir(cci_pdbmtz_path)
  for arg in arg_iter:
    def get_next(expected_exts):
      def raise_bad_file(what, fn=None):
        msg = "%s file name (%s expected)" % (what, " or ".join(expected_exts))
        if (fn is None):
          msg += "."
        else:
          msg += ": " + show_string(fn)
        raise RuntimeError(msg)
      try:
        arg = next(arg_iter)
      except StopIteration:
        raise_bad_file("Missing")
      if (not arg.endswith(tuple(expected_exts))):
        raise_bad_file("Unexpected", arg)
      return arg
    if (op.isfile(arg) and arg.endswith((".mtz", ".pdb", ".ent"))):
      if (arg.endswith(".mtz")):
        fn_mtz = arg
        fn_pdb = get_next([".pdb", ".ent"])
      else:
        fn_pdb = arg
        fn_mtz = get_next([".mtz"])
    else:
      fn_mtz = arg+".mtz"
      def raise_mtz_but_no_pdb():
        raise RuntimeError(
          "MTZ file found but no PDB file: %s" % show_string(fn_mtz))
      if (op.isfile(fn_mtz)):
        for ext in [".pdb", ".ent"]:
          fn_pdb = arg+ext
          if (op.isfile(fn_pdb)):
            break
        else:
          raise_mtz_but_no_pdb()
      else:
        fn_mtz = op.join(cci_pdbmtz_path, arg+".mtz")
        if (not op.isfile(fn_mtz)):
          raise RuntimeError(
            "MTZ file not found: %s" % show_string(fn_mtz))
        fn_pdb = op.join(pdb_v3_mirror_dir, arg[1:3], "pdb"+arg+".ent.gz")
        if (not op.isfile(fn_pdb)):
          raise_mtz_but_no_pdb()
    mtz_pdb_pairs.append((fn_mtz, fn_pdb))
  #
  n_caught = 0
  for i_pair,mtz_pdb_pair in enumerate(mtz_pdb_pairs):
    if (i_pair % command_line.chunk.n != command_line.chunk.i): continue
    tm = user_plus_sys_time()
    try:
      process(params, mtz_pdb_pair)
    except KeyboardInterrupt:
      print("CAUGHT EXCEPTION: KeyboardInterrupt", file=sys.stderr)
      traceback.print_exc()
      print(file=sys.stderr)
      sys.stderr.flush()
      return
    except Exception:
      sys.stdout.flush()
      print("CAUGHT EXCEPTION: %s" % ", ".join(mtz_pdb_pair), file=sys.stderr)
      traceback.print_exc()
      print(file=sys.stderr)
      sys.stderr.flush()
      n_caught += 1
    else:
      print("done_with: %s, %s (%.2f seconds)" % (
        mtz_pdb_pair + (tm.elapsed(),)))
      print()
      sys.stdout.flush()
  print()
  print("Number of exceptions caught:", n_caught)
  #
  show_times()
  print()
  print("TIME END pdb_dev:", date_and_time())
  sys.stdout.flush()

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
