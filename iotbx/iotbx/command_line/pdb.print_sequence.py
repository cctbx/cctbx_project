import os, sys
try:
  from iotbx import pdb
  from iotbx.pdb.hybrid_36 import hy36decode
except ImportError, e:
  print "iotbx not loaded"
  sys.exit()

protein_sequence_to_three = {
  "B" : "ACE",
  "A" : "ALA",
  "C" : "CYS",
  "D" : "ASP",
  "E" : "GLU",
  "F" : "PHE",
  "G" : "GLY",
  "H" : "HIS",
  "I" : "ILE",
  "K" : "LYS",
  "L" : "LEU",
  "M" : "MET",
  "N" : "ASN",
  "P" : "PRO",
  "Q" : "GLN",
  "R" : "ARG",
  "S" : "SER",
  "T" : "THR",
  "V" : "VAL",
  "W" : "TRP",
  "Y" : "TYR",
  "Z" : "NME",
  }

dna_rna_sequence_to_three = {
  "A" : "ADE",
  "C" : "CYT",
  "G" : "GUA",
  "N" : None, # any base
  "T" : "THY",
  "U" : "URI",
  }

ps_lookup = {}
for key in protein_sequence_to_three:
  ps_lookup[protein_sequence_to_three[key]] = key
ds_lookup = {}
for key in dna_rna_sequence_to_three:
  if dna_rna_sequence_to_three[key]:
    ds_lookup[dna_rna_sequence_to_three[key]] = key

def run(filename,
        ignore_residues=[],
        print_unknown=False,
        ):
  protein_only = True
  if os.path.exists(filename):
    pdb_io = pdb.input(filename)
  elif filename.find("\n")>-1:
    pdb_io = pdb.input(source_info=None,
                       lines=flex.split_lines(filename))
  else:
    assert 0
  hierarchy      = pdb_io.construct_hierarchy()
  overall_counts = hierarchy.overall_counts()

  outl = ""
  unk = ""

  for model in hierarchy.models():
    for chain in model.chains():
      for i_conformer, conformer in enumerate(chain.conformers()):
        for residue in conformer.residues():
          atoms = residue.atoms()
          def have_alternative():
            for atom in atoms:
              if(atom.is_alternative()):
                return True
            return False
          if(i_conformer!=0 and not have_alternative()): continue
          if residue.name.strip() in ps_lookup:
            outl += "%s" % ps_lookup[residue.name.strip()]
          elif residue.name.strip() in ds_lookup:
            outl += "%s" % ds_lookup[residue.name.strip()]
            protein_only = False
          elif residue.name.strip() in ignore_residues:
            if unk.find(residue.name.strip())==-1:
              unk += "%s " % residue.name.strip()
          elif residue.name.strip() in dna_rna_sequence_to_three.keys():
            outl += "%s" % residue.name.strip()
            protein_only = False
          else:
            if not protein_only:
              if not residue.atoms()[0].hetero:
                print "This residue is not converted to letter code %s" % (
                  residue.name.strip(),
                  )
                sys.exit()
      outl+="\n"
  if unk and print_unknown:
    print "Unconverted residues",unk
  return outl

def exercise():
  for pdb_file in ["1zap",
                   "1a00",
                   ]:
    print run(os.path.join(os.environ["PDB_MIRROR_UNCOMPRESSED"],
                           "pdb%s.ent",
                           )
              )

if __name__=="__main__":
  print run(sys.argv[1])
