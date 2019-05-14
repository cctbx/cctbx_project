# see iotbx/pdb/common_residue_names.h; additionally here only: U I
from __future__ import absolute_import, division, print_function

ad_hoc_single_metal_residue_element_types = """\
ZN CA MG NA MN K FE CU CD HG NI CO SR CS PT BA TL PB SM AU RB YB LI
MO LU CR OS GD TB LA AG HO GA CE W RU RE PR IR EU AL V PD U
SB SE
TE
""".split()
ad_hoc_non_linking_elements = """\
H D F Cl Br I At
He Ne Ar Kr Xe
""".upper().split()
ad_hoc_first_row = "Li Be B C N O F".upper().split()
#
# most be alpha order
#
amino_acid_bond_cutoff = 1.9
rna_dna_bond_cutoff = 3.4
intra_residue_bond_cutoff = 1.99
saccharide_bond_cutoff = 3.
metal_coordination_cutoff = 3.5
sulfur_bond_cutoff = 2.5
other_bond_cutoff = 2.
#
# defaults to 1
# alpha ordered key
maximum_inter_residue_links = {
  ("common_element",    "other") : 8,
  ("metal",             "other") : 6,
  ("common_amino_acid", "metal") : 2,
  ("common_saccharide", "metal") : 3,
  #("common_element",    "common_saccharide") : 3,
  ("common_rna_dna",    "metal") : 2,
  ("common_rna_dna",    "common_rna_dna") : 5  # for basepairing
  }
maximum_per_atom_links = {
  "common_saccharide" : 1,
  "common_rna_dna"    : 1,
  "common_amino_acid" : 1,
  "other"             : 1,
  }
#
skip_if_both = [
  ["common_water",      "common_water"],
  ["common_amino_acid", "common_water"],
  ["common_saccharide", "common_water"],
  ["common_water",      "other"],
  ]
#
# maximum valence allowed
simple_valence = {
  "O" : 2,
}
# must be sorted pairs
ad_hoc_non_linking_pairs = [
  ["O", "O"], # occurs in nonH models
]
#
def adjust_class(atom, atom_class):
  if atom_class in ["common_element"]:
    if(atom.element.strip().upper() in ad_hoc_single_metal_residue_element_types):
      return "metal"
    elif (atom.element.strip().upper() in ad_hoc_non_linking_elements):
      return "ion" # lone ion
  elif atom_class in ["other"]: # HEM
    if(atom.element.strip().upper() in ad_hoc_single_metal_residue_element_types):
      return "metal"
  return atom_class
#
def sulfur_class(atom, atom_class):
  if atom_class in ["common_amino_acid", "other"]:
    if(atom.element.strip().upper() in ["S"]):
      return "sulfur"
  return atom_class
#
def update_skip_if_longer(amino_acid_bond_cutoff,
                          rna_dna_bond_cutoff,
                          intra_residue_bond_cutoff,
                          saccharide_bond_cutoff,
                          metal_coordination_cutoff,
                          sulfur_bond_cutoff,
                          other_bond_cutoff,
                          ):

  skip_if_longer = {
    ("common_amino_acid", "common_amino_acid") : amino_acid_bond_cutoff**2,
    ("common_amino_acid", "other")             : intra_residue_bond_cutoff**2,
    #
    ("common_rna_dna",       "common_rna_dna") : rna_dna_bond_cutoff**2,
    ("common_rna_dna",       "metal")          : metal_coordination_cutoff**2,
    ("common_rna_dna",       "other")          : other_bond_cutoff**2,
    ("ccp4_mon_lib_rna_dna", "other")          : other_bond_cutoff**2,
    #
    ("common_amino_acid", "common_saccharide") : saccharide_bond_cutoff**2,
    ("common_saccharide", "common_saccharide") : saccharide_bond_cutoff**2,
    #
    ("common_element", "common_water")         : metal_coordination_cutoff**2,
    #
    ("other", "other") : other_bond_cutoff**2,
    #
    ("sulfur", "sulfur") : sulfur_bond_cutoff**2, # is there another place?
    #
    #
    ("common_amino_acid", "common_rna_dna") : amino_acid_bond_cutoff**2,
    ("common_rna_dna",    "common_small_molecule") : other_bond_cutoff**2,
    ("common_amino_acid", "common_small_molecule") : amino_acid_bond_cutoff*other_bond_cutoff,
    ("common_small_molecule", 'other') : other_bond_cutoff*other_bond_cutoff,
    #
    ("metal", "metal") : metal_coordination_cutoff**2,
    }
  for c in ["other",
            "common_water",
            "common_amino_acid",
            ]:
    key = ["metal", c]
    key.sort()
    skip_if_longer[tuple(key)] = metal_coordination_cutoff**2

  ## for key in sorted(skip_if_longer):
  ##   print key, skip_if_longer[key]
  return skip_if_longer

skip_if_longer = update_skip_if_longer(amino_acid_bond_cutoff,
                                       rna_dna_bond_cutoff,
                                       intra_residue_bond_cutoff,
                                       saccharide_bond_cutoff,
                                       metal_coordination_cutoff,
                                       sulfur_bond_cutoff,
                                       other_bond_cutoff,
                                       )

# explicit residue atom names to exclude
ad_hoc_non_linking_atom_names = {'metal': [('HIS', 'CE1'),
                                           ('HIS', 'CD2'),
                                           ('HIS', 'CB'),
                                          ]
                                }
def skip_if_non_linking(lookup, atom1, atom2):
  if 'metal' in lookup:
    excludes = ad_hoc_non_linking_atom_names.get('metal', [])
    test = (atom1.parent().resname.strip(), atom1.name.strip())
    if test in excludes: return False
    test = (atom2.parent().resname.strip(), atom2.name.strip())
    if test in excludes: return False
    return True
  else:
    assert 0

if __name__=="__main__":
  print(skip_if_both)
  print(skip_if_longer)
