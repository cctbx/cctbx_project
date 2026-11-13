"""
root.atoms_with_labels() returns the atoms of root.atoms(), with their
labels, as an af_shared_atom_with_labels array; every label and flag equals
that of the Python generator it replaces (root._atoms_with_labels_python),
and each entry shares its atom with the hierarchy.
"""
import iotbx.pdb
from iotbx.pdb import hierarchy

# two models; a BREAK, altlocs, an insertion code, a segid, a HETATM, and
# chain A twice in model 1 (two chains)
pdb_str = """\
MODEL        1
ATOM      1  N   MET A   1       1.000   2.000   3.000  1.00 10.00      SEGA N
ATOM      2  CA AMET A   1       2.000   2.000   3.000  0.50 10.00      SEGA C
ATOM      3  CA BMET A   1       2.100   2.000   3.000  0.50 10.00      SEGA C
ATOM      4  N   GLY A   1A      3.000   2.000   3.000  1.00 10.00      SEGA N
BREAK
ATOM      5  N   ALA A   5       4.000   2.000   3.000  1.00 10.00           N
HETATM    6  O   HOH B   1       5.000   2.000   3.000  1.00 10.00           O
ATOM      7  N   LYS A   6       6.000   2.000   3.000  1.00 10.00           N
ENDMDL
MODEL        2
ATOM      8  N   MET A   1       1.000   2.000   3.000  1.00 10.00           N
HETATM    9  O   HOH A   2       5.000   2.000   3.000  1.00 10.00           O
ENDMDL
"""


def labels(awl):
  return (awl.model_id, awl.chain_id, awl.resseq, awl.icode, awl.altloc,
    awl.resname, awl.is_first_in_chain, awl.is_first_after_break,
    awl.id_str(), awl.memory_id())


def compare_with_python(root):
  result = root.atoms_with_labels()
  assert type(result) is hierarchy.af_shared_atom_with_labels, type(result)
  expected = [labels(awl) for awl in root._atoms_with_labels_python()]
  assert len(result) == len(expected), (len(result), len(expected))
  for i, awl in enumerate(result):
    assert labels(awl) == expected[i], (i, labels(awl), expected[i])
  atoms = root.atoms()
  assert len(result) == atoms.size(), (len(result), atoms.size())
  for awl, atom in zip(result, atoms):
    assert awl.memory_id() == atom.memory_id(), awl.id_str()
  return result


def exercise_pdb():
  """Labels and flags of every atom of pdb_str equal the generator's."""
  root = iotbx.pdb.input(source_info=None, lines=pdb_str).construct_hierarchy()
  result = compare_with_python(root)
  n_breaks = [awl.is_first_after_break for awl in result].count(True)
  assert n_breaks == 1, n_breaks
  n_chains = [awl.is_first_in_chain for awl in result].count(True)
  assert n_chains == 4, n_chains


def exercise_empty_levels():
  """
  A residue group without atoms leaves is_first_in_chain to the next atom,
  and empty models, chains and atom groups yield nothing.
  """
  root = hierarchy.root()
  root.append_model(hierarchy.model(id="1"))
  model = hierarchy.model(id="2")
  root.append_model(model)
  model.append_chain(hierarchy.chain(id="A"))
  chain = hierarchy.chain(id="B")
  model.append_chain(chain)
  rg = hierarchy.residue_group(resseq="   1", icode=" ")
  rg.append_atom_group(hierarchy.atom_group(altloc="", resname="GLY"))
  chain.append_residue_group(rg)
  rg = hierarchy.residue_group(resseq="   3", icode=" ",
    link_to_previous=False)
  ag = hierarchy.atom_group(altloc="", resname="ALA")
  ag.append_atom(hierarchy.atom().set_name(" N  "))
  ag.append_atom(hierarchy.atom().set_name(" CA "))
  rg.append_atom_group(ag)
  chain.append_residue_group(rg)
  rg = hierarchy.residue_group(resseq="   4", icode=" ",
    link_to_previous=False)
  ag = hierarchy.atom_group(altloc="", resname="SER")
  ag.append_atom(hierarchy.atom().set_name(" N  "))
  rg.append_atom_group(ag)
  chain.append_residue_group(rg)
  result = compare_with_python(root)
  flags = [(awl.is_first_in_chain, awl.is_first_after_break)
    for awl in result]
  assert flags == [(True, False), (False, False), (False, True)], flags
  assert len(compare_with_python(hierarchy.root())) == 0


def run():
  exercise_pdb()
  exercise_empty_levels()


if __name__ == "__main__":
  run()
  print("OK")
