from mmtbx.monomer_library import cif_triage
import iotbx.pdb.amino_acid_codes
import libtbx.load_env
import sys, os
op = os.path

def run(args):
  assert len(args) == 0
  amino_acid_resnames = sorted(
    iotbx.pdb.amino_acid_codes.one_letter_given_three_letter.keys())
  geostd_path = libtbx.env.find_in_repositories(
    relative_path="chem_data/geostd", optional=False)
  for resname in amino_acid_resnames:
    file_name = op.join(
      geostd_path, resname[0].lower(), "data_"+resname+".cif")
    assert cif_triage.check_comp(file_name=file_name) == 1
  print "OK"

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
