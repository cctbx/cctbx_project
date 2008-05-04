import iotbx.pdb
import iotbx.pdb.amino_acid_codes
import sys

def run(args):
  if (len(args) == 0):
    raise RuntimeError("Please specify one or more pdb file names.")
  aa_resnames = iotbx.pdb.amino_acid_codes.one_letter_given_three_letter
  ala_atom_names = set([" N  ", " CA ", " C  ", " O  ", " CB "])
  for file_name in args:
    pdb_obj = iotbx.pdb.hierarchy.input(file_name=file_name)
    pdb_obj.hierarchy.overall_counts().show()
    for model in pdb_obj.hierarchy.models():
      for chain in model.chains():
        for rg in chain.residue_groups():
          for ag in rg.atom_groups():
            if (ag.resname in aa_resnames):
              for atom in ag.atoms():
                if (atom.name not in ala_atom_names):
                  ag.remove_atom(atom=atom)
    output_pdb = "v2_truncated_to_ala_"+file_name
    pdb_obj.hierarchy.write_pdb_file(file_name=output_pdb)

if (__name__ == "__main__"):
  run(sys.argv[1:])
