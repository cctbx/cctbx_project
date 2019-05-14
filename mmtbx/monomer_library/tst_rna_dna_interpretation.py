from __future__ import absolute_import, division, print_function
from mmtbx.monomer_library import pdb_interpretation
from mmtbx.monomer_library import server
from libtbx.utils import null_out
import libtbx.load_env
import sys, os

def run(args, residue_type, expected_results):
  verbose = "--verbose" in args
  debug = "--debug" in args
  pdb_files = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/%s_pdb_files" % residue_type,
    test=os.path.isdir)
  if (pdb_files is None):
    print("Skipping tst_%s_interpretation: input files not available" \
      % residue_type)
    return
  mon_lib_srv = server.server()
  ener_lib = server.ener_lib()
  for file_name in os.listdir(pdb_files):
    if (file_name[-4:] not in [".ent", ".pdb"]): continue
    if (verbose):
      log = sys.stdout
    else:
      log = null_out()
    processed_pdb_file = pdb_interpretation.process(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      file_name=os.path.join(pdb_files, file_name),
      log=log)
    for mm in processed_pdb_file.all_chain_proxies.all_monomer_mappings:
      residue_name = mm.residue_name
      if (mm.classification == "RNA"):
        for c in ["2", "3"]:
          expected_mod = "%rna"+c+"p_pur"
          if (mm.residue_name.find(expected_mod) > 0):
            residue_name = mm.residue_name.replace(expected_mod, "")
            break
          expected_mod = "%rna"+c+"p_pyr"
          if (mm.residue_name.find(expected_mod) > 0):
            residue_name = mm.residue_name.replace(expected_mod, "")
            break
      assert len(mm.duplicate_atoms) == 0
      assert len(mm.ignored_atoms) == 0
      assert not mm.is_unusual
      unexpected_names = [atom.name for atom in mm.unexpected_atoms]
      result = [
        residue_name,
        len(mm.expected_atoms),
        unexpected_names,
        mm.classification,
        mm.is_terminus,
        mm.incomplete_info]
      key = file_name[:-4]
      if (debug):
        print('"%s":' % key)
        print(" ", str(result)+",")
        print(" ", str(expected_results[key])+",")
      if (expected_results is not None):
        assert result == expected_results[key], "%s not %s" % (
          result,
          expected_results[key],
          )
  print("OK")

expected_results = {
"a_v2": ['  A', 33, [], 'RNA', None, None],
"a_v3": ['  A', 33, [], 'RNA', None, None],
"c_3ter_v2": ['  C%3*END', 32, [], 'RNA', True, None],
"c_3ter_v3": ['  C%3*END', 32, [], 'RNA', True, None],
"c_5ter_v2": ['  C%5*END', 28, [], 'DNA', True, None],
"c_5ter_v3": [' DC%5*END', 28, [], 'DNA', True, None],
"c_p_only_v3": ['  C', 1, [], 'DNA', None, 'p_only'],
"c_v2": ['  C', 31, [], 'RNA', None, None],
"c_v3": ['  C', 31, [], 'RNA', None, None],
"da_v2": ['  A', 32, [], 'DNA', None, None],
"da_v3": [' DA', 32, [], 'DNA', None, None],
"dc_v2": ['  C', 30, [], 'DNA', None, None],
"dc_v3": [' DC', 30, [], 'DNA', None, None],
"dg_v2": ['  G', 33, [], 'DNA', None, None],
"dg_v3": [' DG', 33, [], 'DNA', None, None],
"dt_v2": ['  T', 32, [], 'DNA', None, None],
"dt_v3": [' DT', 32, [], 'DNA', None, None],
"g_5pho_v2": ['  G%p5*END', 36, [], 'RNA', True, None],
"g_5pho_v3": ['  G%p5*END', 36, [], 'RNA', True, None],
"g_5ter_3ter_v2": ['  G%5*END%3*END', 33, [], 'RNA', True, None],
"g_5ter_3ter_v3": ['  G%5*END%3*END', 33, [], 'RNA', True, None],
"g_5ter_v2": ['  G%5*END', 32, [], 'RNA', True, None],
"g_5ter_v3": ['  G%5*END', 32, [], 'RNA', True, None],
"g_v2": ['  G', 34, [], 'RNA', None, None],
"g_v3": ['  G', 34, [], 'RNA', None, None],
"u_v2": ['  U', 30, [], 'RNA', None, None],
"u_v3": ['  U', 30, [], 'RNA', None, None],
"thy_cns_dna_rna_top": ['THY', 20, [], 'DNA', None, None],
}

if (__name__ == "__main__"):
  run(sys.argv[1:], residue_type="rna_dna", expected_results=expected_results)
