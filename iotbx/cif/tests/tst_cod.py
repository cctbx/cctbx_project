from iotbx.cif import cod_tools
import iotbx.cif
from iotbx.cif.builders import CifBuilderError
import os, traceback

def run(args, command_name):
  from iotbx.option_parser import option_parser as iotbx_option_parser
  import libtbx.utils
  show_times = libtbx.utils.show_times(time_start="now")
  command_line = (iotbx_option_parser(
    usage=command_name+" [options] [cod_id...]")
    .enable_chunk(easy_all=True)
    .enable_multiprocessing()
    .option(None, "--parse_only",
            action="store_true")
    .option(None, "--cif_only",
            action="store_true")
    .option(None, "--hkl_only",
            action="store_true")
    .option("-v", "--verbose",
            action="store_true")
  ).process(args=args)
  if (command_line.run_multiprocessing_chunks_if_applicable(
        command_call=[command_name, __file__])):
    show_times()
    return
  co = command_line.options
  assert [co.cif_only, co.hkl_only].count(True) <= 1
  if co.cif_only: ext = "cif"
  elif co.hkl_only: ext = "hkl"
  else: ext = None
  verbose = co.verbose
  parse_only = co.parse_only
  #
  cod_hkl_cif = cod_tools.build_hkl_cif(cod_ids=None, ext=ext)
  cod_hkl_cif.show_summary()
  hkl_files = cod_hkl_cif.hkl
  cif_files = cod_hkl_cif.cif
  #
  n_caught = 0
  n_ignored = 0
  n_success = 0
  n_skipped = 0
  files_to_parse = []
  files_to_parse.extend(hkl_files.values())
  files_to_parse.extend(cif_files.values())
  for i, path in enumerate(files_to_parse):
    if (i % command_line.chunk.n != command_line.chunk.i): continue
    try:
      cif_obj = iotbx.cif.reader(file_path=path)
      if parse_only: continue
      skip_file = False
      for cif_block in cif_obj.model().values():
        value = cif_block.get("_cod_error_flag")
        keys = set(cif_block.keys())
        if (value in ["errors", "retracted"]):
          skip_file = True
          cod_id = os.path.basename(path)
          n_skipped += 1
          if verbose:
            print "SKIPPING: _cod_error_flag %s: %s" % (value, cod_id)
        elif (len(set([
          "_space_group_symop_ssg_operation_algebraic",
          "_space_group_ssg_name"]).intersection(keys)) != 0):
          n_skipped += 1
          if verbose:
            print "SKIPPING: COD entry with super-space group:", cod_id
        elif (len(set([
              "_refln_index_m",
              "_refln_index_m_1"]).intersection(keys)) != 0):
          if verbose:
            print "SKIPPING: COD entry with _refln_index_m:", cod_id
          n_skipped += 1
      if skip_file: continue
      if path.endswith('.cif'):
        cif_obj.build_crystal_structures()
      elif path.endswith('.hkl'):
        cif_obj.build_miller_arrays()
      else:
        iotbx.cif.cctbx_data_structures_from_cif(cif_model=cif_obj.model())
    except KeyboardInterrupt:
      print "CAUGHT EXCEPTION: KeyboardInterrupt"
      return
    except CifBuilderError, e:
      e_str = str(e)
      if not verbose and e_str.startswith(
        "No atomic coordinates could be found"):
        n_ignored += 1
        continue
      elif not verbose and e_str.startswith(
        "No symmetry instructions could be extracted from the cif block"):
        n_ignored += 1
        continue
      sys.stdout.flush()
      cod_id = os.path.basename(path)
      print >> sys.stderr, \
        "CAUGHT EXCEPTION: %s: %s: %s" % (command_name, cod_id, str(e))
      if verbose:
        traceback.print_exc()
        print >> sys.stderr
      sys.stderr.flush()
      n_caught += 1
    except Exception, e:
      sys.stdout.flush()
      cod_id = os.path.basename(path)
      print >> sys.stderr, \
        "CAUGHT EXCEPTION: %s: %s: %s" % (command_name, cod_id, str(e))
      if verbose:
        traceback.print_exc()
        print >> sys.stderr
      sys.stderr.flush()
      n_caught += 1
    else:
      n_success += 1
  print
  print "Number successfully parsed: ", n_success
  if not parse_only:
    print "Number skipped:", n_skipped
    print "Number of exceptions caught:", n_caught
    print "Number of exceptions ignored:", n_ignored
  print
  #
  show_times()
  print

if __name__ == '__main__':
  import libtbx.load_env
  import sys
  run(args=sys.argv[1:], command_name=libtbx.env.dispatcher_name)
