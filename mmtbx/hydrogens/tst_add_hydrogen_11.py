from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
import iotbx.pdb
from mmtbx.hydrogens import reduce_hydrogen
from libtbx.utils import null_out

# ------------------------------------------------------------------------------

def run():
  test_000()

# ------------------------------------------------------------------------------

def test_000():
  '''
    An unknown ligand (FCO, no CCD/GeoStd/user restraints) gets a throwaway
    restraint dictionary auto-generated during H placement so that pdb
    interpretation can build it. That dictionary must not leak out as a
    permanent restraint object on the returned model: it is purpose-built for
    placement (bonds shortened to 0.9x, idealized esd=1/period=1 torsions from
    the CCD conformer) and would produce bogus geometry outliers if re-used by
    downstream validation (e.g. mmtbx.development.validate_ligands).
  '''
  pdb_inp = iotbx.pdb.input(lines=pdb_str_000.split("\n"), source_info=None)
  model_initial = mmtbx.model.manager(model_input = pdb_inp, log = null_out())
  assert(model_initial.get_hd_selection().count(True) == 0)

  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(
    model = model_initial,
    stop_for_unknowns = False)
  reduce_add_h_obj.run()
  model_h_added = reduce_add_h_obj.get_model()

  # Sanity: H placement ran fully (H atoms placed on the LEU residue).
  assert(model_h_added.get_hd_selection().count(True) > 0)

  # The auto-generated FCO placement restraints must not persist on the model.
  ro = model_h_added.get_restraint_objects()
  names = [] if ro is None else [name for name, _ in ro]
  assert('auto_FCO' not in names), \
    'auto-generated placement restraints leaked: %s' % names

  # Behavioral check: a downstream consumer that re-derives restraints from the
  # returned model the way validate_ligands does (harvest restraint objects,
  # set them, re-process) must not synthesize bogus FCO geometry. Re-processing
  # rebuilds mon_lib_srv from the restraint objects only, so with the leak fixed
  # the unknown ligand stays unknown and gets no restraints.
  ro = model_h_added.get_restraint_objects()
  model_h_added.set_restraint_objects([] if ro is None else list(ro))
  model_h_added.set_stop_for_unknowns(False)
  model_h_added.process(make_restraints=True)
  isel = model_h_added.iselection('resname FCO and not (element H or element D)')
  fco = model_h_added.select(isel)
  grm = fco.get_restraints_manager().geometry
  assert(grm.dihedral_proxies.size() == 0), \
    'bogus FCO dihedral restraints present: %d' % grm.dihedral_proxies.size()

# ------------------------------------------------------------------------------

pdb_str_000 = """
CRYST1  100.667  101.210  170.826  90.00  90.00  90.00 P 21 21 21
ATOM      1  N   LEU A 482     114.924  99.962 -27.431  1.00 10.00           N
ATOM      2  CA  LEU A 482     114.188 101.177 -27.843  1.00 10.00           C
ATOM      3  C   LEU A 482     115.223 102.183 -28.350  1.00 10.00           C
ATOM      4  O   LEU A 482     116.313 102.353 -27.708  1.00 10.00           O
ATOM      5  CB  LEU A 482     113.411 101.726 -26.625  1.00 10.00           C
ATOM      6  CG  LEU A 482     112.764 103.127 -26.762  1.00 10.00           C
ATOM      7  CD1 LEU A 482     111.736 103.154 -27.867  1.00 10.00           C
ATOM      8  CD2 LEU A 482     112.154 103.534 -25.429  1.00 10.00           C
TER
HETATM    9  C1  FCO A 601     107.099  98.547 -26.445  1.00 10.00           C
HETATM   10  C2  FCO A 601     107.980  98.742 -23.921  1.00 10.00           C
HETATM   11  C3  FCO A 601     108.613 100.597 -25.504  1.00 10.00           C
HETATM   12  N1  FCO A 601     107.159  97.752 -27.314  1.00 10.00           N
HETATM   13  N2  FCO A 601     108.687  98.106 -23.190  1.00 10.00           N
HETATM   14  O3  FCO A 601     109.639 101.050 -25.723  1.00 10.00           O
HETATM   15 FE   FCO A 601     107.039  99.845 -25.099  1.00 10.00          FE
END
"""

# ------------------------------------------------------------------------------

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
