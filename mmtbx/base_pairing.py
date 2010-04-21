
import os

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

class pair_database (object) :
  def __init__ (self) :
    self._h_bond_pairs = {}
    self._pseudo_bond_pairs = {}
    base, ext = os.path.splitext(__file__)
    dblines = open("%s.data" % base).readlines()
    for i, line in enumerate(dblines) :
      if line.startswith("#") :
        continue
      fields = line.strip().split()
      assert len(fields) >= 4
      pair_type = fields[0]
      paired_bases = fields[1]
      hydrogen_flag = fields[2]
      atom_pairs = [ (p.split(",")[0], p.split(",")[1]) for p in fields[3:] ]
      if hydrogen_flag == '+' :
        db = self._h_bond_pairs
      else :
        db = self._pseudo_bond_pairs
      if pair_type in db :
        if paired_bases in db[pair_type] :
          raise RuntimeError("Duplicate entry in base pair dictionary, line %d"
            % i)
      else :
        db[pair_type] = {}
      db[pair_type][paired_bases] = atom_pairs

  def get_atoms (self, base_pair, pair_type, use_hydrogens=False) :
    if use_hydrogens :
      db = self._h_bond_pairs
    else :
      db = self._pseudo_bond_pairs
    pair_rules = db[pair_type]
    if base_pair in pair_rules :
      return pair_rules[base_pair]
    elif base_pair[::-1] in pair_rules :
      return invert_pairs(pair_rules[base_pair[::-1]])
    else :
      raise RuntimeError("No entry for base pair %s with type %s (H=%s)." %
        (base_pair, pair_type, str(use_hydrogens)))

  def get_pair_type (self, base_pair, atom_pairs, use_hydrogens=False) :
    return_pair_type = None
    if use_hydrogens :
      db = self._h_bond_pairs
    else:
      db = self._pseudo_bond_pairs
    for pair_type in db:
      if base_pair in db[pair_type]:
        if db[pair_type][base_pair].sort(key=sort_tuple) == atom_pairs.sort(key=sort_tuple):
          if return_pair_type is None:
            return_pair_type = pair_type
          else:
            raise RuntimeError("Redundant entries found for base pair %s." % base_pair)
      elif base_pair[::-1] in db[pair_type]:
        if db[pair_type][base_pair[::-1]].sort(key=sort_tuple) == \
           invert_pairs(atom_pairs).sort(key=sort_tuple):
          if return_pair_type is None:
            return_pair_type = pair_type
          else:
            raise RuntimeError("Redundant entries found for base pair %s." % base_pair)
    return return_pair_type

def invert_pairs(h_bond_atoms):
  inverted_h_bond_atoms = []
  for pair in h_bond_atoms:
    inverted_h_bond_atoms.append((pair[1],pair[0]))
  return inverted_h_bond_atoms

def sort_tuple(tuple):
  return tuple[0]

db = pair_database()

def get_h_bond_atoms(residues, pair_type, use_hydrogens=False):
  base_pair = residues[0].strip()[0] + residues[1].strip()[0]
  return db.get_atoms(base_pair, pair_type, use_hydrogens)

def exercise () :
  assert (db.get_atoms("AU", "WWT", True) == [('H61', 'O4'), ('N1', 'H3')])
  assert (db.get_atoms("GC", "WWT", False) == [('O6', 'N4'), ('N1', 'N3'),
    ('N2', 'O2')])
  assert db.get_pair_type("AU", [('H61', 'O4'), ('N1', 'H3')], True) == "WWT"
  assert db.get_pair_type("AU", [('N1', 'H3'), ('H61', 'O4')], True) == "WWT"
  assert db.get_pair_type("CG", [('N4', 'O6'), ('N3', 'N1'), ('O2', 'N2')], False) == "WWT"
  print "OK"

if __name__ == "__main__" :
  exercise()
