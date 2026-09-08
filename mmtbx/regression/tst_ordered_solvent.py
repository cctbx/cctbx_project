from __future__ import absolute_import, division, print_function
import time
from cctbx.array_family import flex
import iotbx.pdb
import mmtbx.model
import mmtbx.solvent.ordered_solvent as ordered_solvent

pdb_str = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
ATOM      1  N   GLY A   1       5.000   5.000   5.000  1.00 10.00           N
ATOM      2  CA  GLY A   1       6.000   5.000   5.000  1.00 10.00           C
ATOM      3  C   GLY A   1       7.000   5.000   5.000  1.00 10.00           C
ATOM      4  O   GLY A   1       7.500   6.000   5.000  1.00 10.00           O
TER
END
"""

def exercise_added_solvent_atom_name():
  """
  Water added by the ordered-solvent machinery must have the oxygen
  right-justified in the PDB atom-name field, i.e. the element in column 14
  (" O  "), as the PDB spec requires for single-character elements. This calls
  add_solvent_to_model_inplace directly -- the exact call the ordered_solvent
  manager (phenix.refine / ligand pipeline) makes in _add_new_solvent, and
  which model.add_solvent also routes through. Regression test for the oxygen
  ending up in column 13 ("O   "), which coot and other parsers reject.
  """
  model = mmtbx.model.manager(
    model_input = iotbx.pdb.input(source_info=None, lines=pdb_str.splitlines()))
  model.setup_scattering_dictionaries(scattering_table="wk1995")
  sites_frac = flex.vec3_double([(0.40, 0.40, 0.40), (0.60, 0.60, 0.60)])
  params = ordered_solvent.master_params().extract()
  ordered_solvent.add_solvent_to_model_inplace(
    sites=sites_frac, model=model, params=params)
  # print(model.model_as_pdb())
  water_lines = [l for l in model.model_as_pdb().splitlines()
                 if l.startswith(("ATOM", "HETATM")) and l[17:20] == "HOH"]
  assert len(water_lines) == 2, water_lines
  for line in water_lines:
    # PDB columns 13-16 (0-based 12:16) are the atom-name field; a
    # single-character element must sit in column 14 with a leading blank.
    assert line[12:16] == " O  ", \
      "water oxygen must be in column 14, got name field %r in:\n%s" % (
        line[12:16], line)

pdb_str_with_water = pdb_str.replace("END\n", """\
HETATM    5  O   HOH S   1      10.000  10.000  10.000  1.00 30.00           O
HETATM    6  O   HOH S   2      12.000  10.000  10.000  1.00 30.00           O
HETATM    7  O   HOH S   3      14.000  10.000  10.000  1.00 30.00           O
HETATM    8  O   HOH S   4      16.000  10.000  10.000  1.00 30.00           O
END
""")

def exercise_added_solvent_unique_serials():
  """
  Atom serial numbers must stay unique after the filter-then-add sequence the
  ordered_solvent manager runs every macro-cycle. Removing waters with
  model.select() leaves the survivors with their old serials, so new waters
  must not be numbered from model.size(). Regression test for duplicate
  PDB serials / _atom_site.id on waters in phenix.refine output (UR-16031).
  """
  model = mmtbx.model.manager(
    model_input = iotbx.pdb.input(
      source_info=None, lines=pdb_str_with_water.splitlines()))
  model.setup_scattering_dictionaries(scattering_table="wk1995")
  # Emulate the "Filter" steps: drop two of the four waters.
  model = model.select(model.selection("not (chain S and resseq 2:3)"))
  # Emulate "Add new water".
  sites_frac = flex.vec3_double([(0.40, 0.40, 0.40), (0.60, 0.60, 0.60)])
  params = ordered_solvent.master_params().extract()
  ordered_solvent.add_solvent_to_model_inplace(
    sites=sites_frac, model=model, params=params)
  serials = [a.serial.strip() for a in model.get_hierarchy().atoms()]
  assert len(set(serials)) == len(serials), \
    "duplicate atom serials after add_solvent_to_model_inplace: %s" % serials
  pdb_serials = [l[6:11] for l in model.model_as_pdb().splitlines()
                 if l.startswith(("ATOM", "HETATM"))]
  assert len(set(pdb_serials)) == len(pdb_serials), \
    "duplicate serials in PDB output: %s" % pdb_serials

def run():
  exercise_added_solvent_atom_name()
  exercise_added_solvent_unique_serials()

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %6.2f" % (time.time() - t0))
