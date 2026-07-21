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

def run():
  exercise_added_solvent_atom_name()

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %6.2f" % (time.time() - t0))
