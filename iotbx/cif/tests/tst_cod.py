from iotbx.cif import cod_tools
import iotbx.cif
from iotbx.cif.builders import CifBuilderError
from iotbx.cif import CifParserError
from libtbx import easy_pickle, group_args
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
  cod_ids = command_line.args
  assert [co.cif_only, co.hkl_only].count(True) <= 1
  if co.cif_only: ext = "cif"
  elif co.hkl_only: ext = "hkl"
  else: ext = None
  verbose = co.verbose
  parse_only = co.parse_only
  #
  cod_hkl_cif = cod_tools.build_hkl_cif(cod_ids=cod_ids, ext=ext)
  cod_hkl_cif.show_summary()
  hkl_files = cod_hkl_cif.hkl
  cif_files = cod_hkl_cif.cif
  #
  n_total = 0
  #
  parsing_errors = {}
  build_errors = {}
  ignored_errors = {}
  skipped = set()
  #
  files_to_parse = []
  files_to_parse.extend(hkl_files.values())
  files_to_parse.extend(cif_files.values())
  for i, path in enumerate(files_to_parse):
    n_total += 1
    if (i % command_line.chunk.n != command_line.chunk.i): continue
    try:
      cod_id = os.path.basename(path)
      cif_obj = iotbx.cif.reader(file_path=path)
      if parse_only: continue
      skip_file = False
      for cif_block in cif_obj.model().values():
        value = cif_block.get("_cod_error_flag")
        keys = set(cif_block.keys())
        if (value in ["errors", "retracted"]):
          skip_file = True
          skipped.add(cod_id)
          if verbose:
            print "SKIPPING: _cod_error_flag %s: %s" % (value, cod_id)
        elif (len(set([
          "_space_group_symop_ssg_operation_algebraic",
          "_space_group_ssg_name"]).intersection(keys)) != 0):
          skipped.add(cod_id)
          if verbose:
            print "SKIPPING: COD entry with super-space group:", cod_id
        elif (len(set([
              "_refln_index_m",
              "_refln_index_m_1"]).intersection(keys)) != 0):
          if verbose:
            print "SKIPPING: COD entry with _refln_index_m:", cod_id
          skipped.add(cod_id)
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
      if not verbose and (
        e_str.startswith("No atomic coordinates could be found") or
        e_str.startswith(
          "No symmetry instructions could be extracted from the cif block")):
        ignored_errors.setdefault(cod_id, e_str)
        continue
      sys.stdout.flush()
      print >> sys.stderr, \
        "CAUGHT EXCEPTION: %s: %s: %s" % (command_name, cod_id, str(e))
      if verbose:
        traceback.print_exc()
        print >> sys.stderr
      build_errors.setdefault(cod_id, e_str)
      sys.stderr.flush()
    except CifParserError, e:
      sys.stdout.flush()
      e_str = str(e)
      parsing_errors.setdefault(cod_id, e_str)
      print >> sys.stderr, \
        "PARSING ERROR: %s: %s: %s" % (command_name, cod_id, e_str)
      if verbose:
        traceback.print_exc()
        print >> sys.stderr
      sys.stderr.flush()
    except Exception, e:
      sys.stdout.flush()
      e_str = str(e)
      build_errors.setdefault(cod_id, e_str)
      print >> sys.stderr, \
        "CAUGHT EXCEPTION: %s: %s: %s" % (command_name, cod_id, e_str)
      if verbose:
        traceback.print_exc()
        print >> sys.stderr
      sys.stderr.flush()
  print

  print "Number successfully parsed: %i/%i" \
        % (n_total-len(parsing_errors),n_total)
  if not parse_only:
    print "Number skipped:", len(skipped)
    print "Number of exceptions caught:", len(build_errors)
    print "Number of exceptions ignored:", len(ignored_errors)
  print
  #
  show_times()
  result = group_args(
    n_hkl=len(hkl_files),
    n_cif=len(cif_files),
    n_hkl_cif_pairs=len(cod_hkl_cif.hkl_cif_pairs),
    parsing_errors=parsing_errors,
    build_errors=build_errors,
    ignored_errors=ignored_errors,
    skipped=skipped)
  easy_pickle.dump("result_%03i.pickle" %command_line.chunk.i, result)
  print

if __name__ == '__main__':
  import libtbx.load_env
  import sys
  run(args=sys.argv[1:], command_name=libtbx.env.dispatcher_name)
