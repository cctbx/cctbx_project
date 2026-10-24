"""Water O-H in model.idealize_h_riding(use_neutron_distances=...).

None (default) keeps the scattering-table rule (0.85, or 0.98 for neutron and
electron tables) and only touches an atom_group with O and two H/D. True uses
the HOH restraint target, False 0.85 A; both find water H/D through their
bonded O, so incomplete waters and blank-altloc O with altloc H/D are included.
"""
from __future__ import absolute_import, division, print_function
import iotbx.pdb
import iotbx.cif
import mmtbx.model
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal
from scitbx import matrix
from scitbx.array_family import flex
from mmtbx.hydrogens import reduce_hydrogen

# one water with distorted H (0.90 A, 95 deg)
pdb_water = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
HETATM    1  O   HOH A   1      10.000  10.000  10.000  1.00 20.00           O
HETATM    2  H1  HOH A   1      10.900  10.000  10.000  1.00 20.00           H
HETATM    3  H2  HOH A   1       9.922  10.897  10.000  1.00 20.00           H
END
"""

# water layouts, O-H 0.90 A
layouts = {
# (a) O + one H
"a": """\
HETATM    1  O   HOH A   1      10.000  10.000  10.000  1.00 20.00           O
HETATM    2  H1  HOH A   1      10.900  10.000  10.000  1.00 20.00           H
""",
# (b) blank-altloc O, H1/H2 in A
"b": """\
HETATM    1  O   HOH A   1      10.000  10.000  10.000  1.00 20.00           O
HETATM    2  H1 AHOH A   1      10.900  10.000  10.000  0.50 20.00           H
HETATM    3  H2 AHOH A   1       9.922  10.897  10.000  0.50 20.00           H
""",
# (c) blank-altloc O, different H1/H2 in A and B
"c": """\
HETATM    1  O   HOH A   1      10.000  10.000  10.000  1.00 20.00           O
HETATM    2  H1 AHOH A   1      10.900  10.000  10.000  0.50 20.00           H
HETATM    3  H2 AHOH A   1       9.922  10.897  10.000  0.50 20.00           H
HETATM    4  H1 BHOH A   1       9.100  10.000  10.000  0.50 20.00           H
HETATM    5  H2 BHOH A   1      10.000   9.100  10.000  0.50 20.00           H
""",
# (d) O, H1, H2 all in A and in B
"d": """\
HETATM    1  O  AHOH A   1      10.000  10.000  10.000  0.50 20.00           O
HETATM    2  H1 AHOH A   1      10.900  10.000  10.000  0.50 20.00           H
HETATM    3  H2 AHOH A   1       9.922  10.897  10.000  0.50 20.00           H
HETATM    4  O  BHOH A   1      10.300  10.000  10.000  0.50 20.00           O
HETATM    5  H1 BHOH A   1       9.400  10.000  10.000  0.50 20.00           H
HETATM    6  H2 BHOH A   1      10.300   9.100  10.000  0.50 20.00           H
""",
# (e) DOD with D
"e": """\
HETATM    1  O   DOD A   1      10.000  10.000  10.000  1.00 20.00           O
HETATM    2  D1  DOD A   1      10.900  10.000  10.000  1.00 20.00           D
HETATM    3  D2  DOD A   1       9.922  10.897  10.000  1.00 20.00           D
""",
}

# user HOH restraints with a neutron O-H target that differs from geostd's
hoh_cif_096 = """\
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 HOH  HOH  'water'  ligand  3  1  .

data_comp_HOH
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
 HOH  O   O  OH2   0  -23.1198  18.3926  -21.5996
 HOH  H1  H  HOH2  0  -22.1601  18.3725  -21.6036
 HOH  H2  H  HOH2  0  -23.3469  18.3903  -20.6669

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
_chem_comp_bond.value_dist_neutron
 HOH  O  H1  single  0.850  0.020  0.960
 HOH  O  H2  single  0.850  0.020  0.960

loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
 HOH  H2  O  H1  103.91  3.000
"""

def make_model(pdb_str, use_neutron_distances=False, scattering_table=None,
               cif=None, process=True):
  inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  restraint_objects = None
  if cif is not None:
    restraint_objects = [("hoh.cif", iotbx.cif.reader(input_string=cif).model())]
  model = mmtbx.model.manager(model_input=inp, log=null_out(),
    restraint_objects=restraint_objects)
  if scattering_table is not None:
    model.setup_scattering_dictionaries(scattering_table=scattering_table)
  if process:
    p = reduce_hydrogen.get_reduce_pdb_interpretation_params(
      use_neutron_distances)
    model.process(make_restraints=True, pdb_interpretation_params=p)
  assert model.riding_h_manager is None
  return model

def layout_model(key, **kwargs):
  cryst1 = "CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1\n"
  return make_model(cryst1 + layouts[key] + "END\n", **kwargs)

def water_oh(model):
  sites = model.get_sites_cart()
  o, h1, h2 = [matrix.col(sites[i]) for i in range(3)]
  return (h1 - o).length(), (h2 - o).length(), (h1 - o).angle(h2 - o, deg=True)

def check_oh(model, d_ideal):
  d1, d2, angle = water_oh(model)
  assert approx_equal(d1, d_ideal, eps=1e-3, out=None), (d1, d_ideal)
  assert approx_equal(d2, d_ideal, eps=1e-3, out=None), (d2, d_ideal)
  assert approx_equal(angle, 103.91, eps=0.01, out=None), angle

def water_groups(model):
  """per water O and conformer: (O i_seq, [H/D i_seq])"""
  h = model.get_hierarchy()
  result = set()
  for chain in h.only_model().chains():
    for conformer in chain.conformers():
      for residue in conformer.residues():
        ra = residue.atoms()
        o = [a.i_seq for a in ra if a.element.strip().upper() == "O"]
        hs = [a.i_seq for a in ra if a.element.strip().upper() in ("H", "D")]
        if len(o) == 1 and hs:
          result.add((o[0], tuple(sorted(hs))))
  return sorted(result)

def check_layout(model, d_ideal):
  sites = model.get_sites_cart()
  groups = water_groups(model)
  assert groups
  for o, hs in groups:
    O = matrix.col(sites[o])
    vs = [matrix.col(sites[h]) - O for h in hs]
    for v in vs:
      assert approx_equal(v.length(), d_ideal, eps=1e-3, out=None), v.length()
    if len(vs) == 2:
      angle = vs[0].angle(vs[1], deg=True)
      assert approx_equal(angle, 103.91, eps=0.01, out=None), angle

def test_default_unchanged():
  # None: today's scattering-table rule, whatever the restraints
  for neutron_restraints in (False, True):
    for table, d in (("n_gaussian", 0.85), ("electron", 0.98),
                     ("neutron", 0.98)):
      model = make_model(pdb_water, neutron_restraints, table)
      model.idealize_h_riding()
      check_oh(model, d)

def test_neutron():
  model = make_model(pdb_water, True, "n_gaussian")
  for m in (model, model.deep_copy(),
            model.select(flex.bool(model.size(), True))):
    m.idealize_h_riding(use_neutron_distances=True)
    check_oh(m, 0.98)
  # the restraint target, not the 0.98 fallback
  model = make_model(pdb_water, True, "n_gaussian", cif=hoh_cif_096)
  model.idealize_h_riding(use_neutron_distances=True)
  check_oh(model, 0.96)
  # False: X-ray length whatever the table
  model = make_model(pdb_water, False, "electron")
  model.idealize_h_riding(use_neutron_distances=False)
  check_oh(model, 0.85)

def test_setup_riding_h_manager_passes_flag():
  model = make_model(pdb_water, True, "n_gaussian", cif=hoh_cif_096)
  model.setup_riding_h_manager(use_neutron_distances=True)
  check_oh(model, 0.96)
  model = make_model(pdb_water, True, "n_gaussian", cif=hoh_cif_096)
  model.setup_riding_h_manager()
  check_oh(model, 0.85)

def test_layouts_default_unchanged():
  # None: only an atom_group with O + two H/D is corrected (d, e); the others
  # keep their coordinates exactly
  for key in sorted(layouts):
    model = layout_model(key, scattering_table="n_gaussian")
    before = model.get_sites_cart().deep_copy()
    model.idealize_h_riding()
    if key in ("d", "e"):
      check_layout(model, 0.85)
    else:
      assert list(model.get_sites_cart()) == list(before), key

def test_layouts_opt_in():
  for key in sorted(layouts):
    # restraint targets
    model = layout_model(key, use_neutron_distances=True)
    model.idealize_h_riding(use_neutron_distances=True)
    check_layout(model, 0.98)
    model = layout_model(key, use_neutron_distances=False)
    model.idealize_h_riding(use_neutron_distances=False)
    check_layout(model, 0.85)
    # no restraints manager: through the hierarchy conformers
    model = layout_model(key, process=False)
    model.idealize_h_riding(use_neutron_distances=True)
    check_layout(model, 0.98)
  # (c): A and B keep distinct H positions
  model = layout_model("c", use_neutron_distances=True)
  model.idealize_h_riding(use_neutron_distances=True)
  sites = model.get_sites_cart()
  assert (matrix.col(sites[1]) - matrix.col(sites[3])).length() > 0.5
  assert (matrix.col(sites[2]) - matrix.col(sites[4])).length() > 0.5

if __name__ == "__main__":
  test_default_unchanged()
  test_neutron()
  test_setup_riding_h_manager_passes_flag()
  test_layouts_default_unchanged()
  test_layouts_opt_in()
  print("OK")
