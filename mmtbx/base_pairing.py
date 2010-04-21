
dna_rna_params_str = """
base_pair
  .multiple = True
  .optional = True
  .style = noauto
{
  base1 = None
    .type = str
  base2 = None
    .type = str
  pair_type = *wwt
    .type = choice
    .caption = Watson-Crick
# TODO
#  planar = False
#    .type = bool
}
"""
def get_h_bond_dict():
  h_bond_dict = {}
  h_bond_dict['CGWWT'] = [('H42','O6'),('N3','H1'),('O2','H22')]
  h_bond_dict['AUWWT'] = [('H61','O4'),('N1','H3')]
  return h_bond_dict

def get_h_bond_heavy_dict():
  h_bond_dict = {}
  h_bond_dict['CGWWT'] = [('N4','O6'),('N3','N1'),('O2','N2')]
  h_bond_dict['AUWWT'] = [('N6','O4'),('N1','N3')]
  return h_bond_dict

def invert_pairs(h_bond_atoms):
  inverted_h_bond_atoms = []
  for pair in h_bond_atoms:
    inverted_h_bond_atoms.append((pair[1],pair[0]))
  return inverted_h_bond_atoms

def get_h_bond_atoms(residues, pair_type, use_hydrogens=False):
  if use_hydrogens:
    h_bond_dict = get_h_bond_dict()
  else:
    h_bond_dict = get_h_bond_heavy_dict()
  key = residues[0].strip()+residues[1].strip()+pair_type
  try:
    h_bond_atoms = h_bond_dict[key]
  except:
    try:
      key = residues[1].strip()+residues[0].strip()+pair_type
      h_bond_atoms = invert_pairs(h_bond_dict[key])
    except:
      print "Unknown key!"
      return None
  return h_bond_atoms
