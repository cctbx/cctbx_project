"""hierarchy.exchangeable_hd_selections() gives the same list, in the same
order, as the implementation it replaced (reference() below)."""
from __future__ import absolute_import, division, print_function
import iotbx.pdb
from iotbx.cli_parser import run_program
from libtbx.utils import null_out

def reference(hierarchy):
  """the implementation before the speed-up, unchanged"""
  result = []
  for model in hierarchy.models():
    for chain in model.chains():
      for residue_group in chain.residue_groups():
        for i_gr1, atom_group_1 in enumerate(residue_group.atom_groups()):
          elements_group1 = [
            e.strip() for e in atom_group_1.atoms().extract_element()]
          non_H_atoms_group1 = list(set(elements_group1) - set(['H','D']))
          for i_gr2, atom_group_2 in enumerate(residue_group.atom_groups()):
            elements_group2 = [
              e.strip() for e in atom_group_2.atoms().extract_element()]
            non_H_atoms_group2 = list(set(elements_group2) - set(['H','D']))
            if non_H_atoms_group1 and non_H_atoms_group2: continue
            if(atom_group_1.altloc != atom_group_2.altloc and i_gr2 > i_gr1):
              for atom1 in atom_group_1.atoms():
                e1 = atom1.element.strip()
                n1 = atom1.name.strip()[1:]
                for atom2 in atom_group_2.atoms():
                  e2 = atom2.element.strip()
                  n2 = atom2.name.strip()[1:]
                  if(e1 in ["H","D"] and e2 in ["H","D"] and e1 != e2 and
                     n1 == n2):
                    result.append([[int(atom1.i_seq)], [int(atom2.i_seq)]])
  return result

def atoms(records):
  """records: (name, altloc, resname, resseq, element); 4-character names"""
  lines = []
  for i, (name, altloc, resname, resseq, element) in enumerate(records):
    lines.append("ATOM  %5d %4s%1s%3s A%4d    %8.3f%8.3f%8.3f  0.50 20.00"
      "          %2s" % (i + 1, name, altloc, resname, resseq, i * 0.7, 0, 0,
      element))
  return "\n".join(lines) + "\n"

def check(text, n_expected=None):
  h = iotbx.pdb.input(source_info=None, lines=text).construct_hierarchy()
  result = h.exchangeable_hd_selections()
  expected = reference(h)
  assert result == expected, (result, expected)
  if n_expected is not None:
    assert len(result) == n_expected, (len(result), n_expected)
  return result

# Ser: heavy-atom altlocs A/B, plus H/D-only groups C and D
ser = atoms([
  (" N  ", " ", "SER", 1, "N"), (" CA ", " ", "SER", 1, "C"),
  (" H  ", "C", "SER", 1, "H"), (" CB ", "A", "SER", 1, "C"),
  (" OG ", "A", "SER", 1, "O"), (" HG ", "A", "SER", 1, "H"),
  (" CB ", "B", "SER", 1, "C"), (" OG ", "B", "SER", 1, "O"),
  (" DG ", "B", "SER", 1, "D"), (" HG ", "C", "SER", 1, "H"),
  (" D  ", "D", "SER", 1, "D"), (" DG ", "D", "SER", 1, "D")])

# waters: HOH A / DOD B (both with O), and a blank O with H in A, D in B
waters = atoms([
  (" O  ", "A", "HOH", 2, "O"), (" H1 ", "A", "HOH", 2, "H"),
  (" H2 ", "A", "HOH", 2, "H"), (" O  ", "B", "DOD", 2, "O"),
  (" D1 ", "B", "DOD", 2, "D"), (" D2 ", "B", "DOD", 2, "D"),
  (" O  ", " ", "HOH", 3, "O"), (" H1 ", "A", "HOH", 3, "H"),
  (" H2 ", "A", "HOH", 3, "H"), (" D1 ", "B", "HOH", 3, "D"),
  (" D2 ", "B", "HOH", 3, "D")])

# digit-first names give no pair; one H with two D (groups B and C, and two
# D of the same name in C); a D in the H group too
names = atoms([
  (" CB ", " ", "LYS", 4, "C"), ("1HB ", "A", "LYS", 4, "H"),
  ("1DB ", "B", "LYS", 4, "D"), (" HZ1", "A", "LYS", 4, "H"),
  (" DZ1", "A", "LYS", 4, "D"), (" DZ1", "B", "LYS", 4, "D"),
  (" DZ1", "C", "LYS", 4, "D"), (" DZ1", "C", "LYS", 4, "D")])

# blank element in an otherwise H/D-only group: the group counts as heavy,
# so H in A and D in C (both groups heavy) give no pair
blank = atoms([
  (" N  ", " ", "GLY", 5, "N"), (" H  ", "A", "GLY", 5, "H"),
  (" X  ", "A", "GLY", 5, "  "), (" D  ", "B", "GLY", 5, "D"),
  (" OG ", "C", "GLY", 5, "O"), (" D  ", "C", "GLY", 5, "D"),
  (" H  ", "D", "GLY", 5, "H"), (" XG ", "E", "GLY", 5, "  ")])

def test_constructed():
  check(ser, 4)
  check(waters, 2)
  check(names, 3)
  check(blank, 3)
  check(ser + waters + names + blank, 12)
  two_models = "MODEL        1\n" + ser + "ENDMDL\nMODEL        2\n" + \
    ser.replace(" HG ", " HX ") + "ENDMDL\n"
  check(two_models, 5)

def test_empty():
  # H only (with altlocs), D only, and H and D only in single-group residues
  h_only = ser.replace(" D", " H")
  d_only = ser.replace(" H", " D")
  assert "D" not in [l[76:78].strip() for l in h_only.splitlines()]
  assert "H" not in [l[76:78].strip() for l in d_only.splitlines()]
  check(h_only, 0)
  check(d_only, 0)
  single = atoms([
    (" N  ", " ", "SER", 1, "N"), (" H  ", " ", "SER", 1, "H"),
    (" N  ", " ", "SER", 2, "N"), (" D  ", " ", "SER", 2, "D")])
  check(single, 0)

def test_reduce2_outputs():
  # reduce2 neutron_option outputs of a small model
  from mmtbx.programs import reduce2
  from mmtbx.regression.tst_reduce2_neutron import pdb_1yjp
  with open("tst_exchangeable_hd_selections.pdb", "w") as f:
    f.write(pdb_1yjp)
  for option in ("hd_and_h", "hd_and_d", "all_hd"):
    prefix = "tst_exchangeable_hd_selections_" + option
    run_program(program_class=reduce2.Program, logger=null_out(), args=[
      "tst_exchangeable_hd_selections.pdb", "neutron_option=" + option,
      "exclude_water=False", "overwrite=True",
      "output.filename=" + prefix + ".pdb",
      "output.description_file=" + prefix + ".txt"])
    with open(prefix + ".pdb") as f:
      result = check(f.read())
    assert len(result) > 0, option

if __name__ == "__main__":
  test_constructed()
  test_empty()
  test_reduce2_outputs()
  print("OK")
