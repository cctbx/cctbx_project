from __future__ import absolute_import, division, print_function

import iotbx.pdb
import mmtbx.model
import iotbx.cif
from cctbx import geometry_restraints
from cctbx.array_family import flex
from libtbx.utils import null_out

# A tiny custom monomer with a single H whose only third-neighbor source is
# a CONST torsion. The CONST torsion fixes HX-X-Y-Z = 60.0 deg.
custom_cif = """\
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 ZZZ ZZZ test_const non-polymer 4 3 .

data_comp_ZZZ
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
 ZZZ X  C  CH1  0.00
 ZZZ Y  C  C    0.00
 ZZZ Z  C  CH3  0.00
 ZZZ HX H  H    0.00

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 ZZZ X  Y   single 1.520 0.020
 ZZZ Y  Z   single 1.520 0.020
 ZZZ X  HX  single 0.970 0.020

loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
 ZZZ Y  X  HX  109.500 1.500
 ZZZ X  Y  Z   111.000 1.500

loop_
_chem_comp_tor.comp_id
_chem_comp_tor.id
_chem_comp_tor.atom_id_1
_chem_comp_tor.atom_id_2
_chem_comp_tor.atom_id_3
_chem_comp_tor.atom_id_4
_chem_comp_tor.value_angle
_chem_comp_tor.value_angle_esd
_chem_comp_tor.period
 ZZZ CONST_01 HX X Y Z 60.000 0.000 0
"""

pdb_str = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
HETATM    1  X   ZZZ A   1       0.000   0.000   0.000  1.00 20.00           C
HETATM    2  Y   ZZZ A   1       1.520   0.000   0.000  1.00 20.00           C
HETATM    3  Z   ZZZ A   1       2.064   1.434   0.000  1.00 20.00           C
HETATM    4  HX  ZZZ A   1      -0.314   0.769   0.583  1.00 20.00           H
END
"""


def _build_model_with_custom_cif(pdb_str, cif_str):
  """Build an mmtbx.model with a custom monomer cif loaded via restraint_objects.
  The cif is written to a tempfile and parsed via iotbx.cif.reader (matching
  the pattern used in mmtbx.geometry_restraints.torsion_restraints.utils)."""
  import tempfile, os
  fd, cif_path = tempfile.mkstemp(suffix=".cif")
  os.close(fd)
  with open(cif_path, "w") as f:
    f.write(cif_str)
  cif_object = iotbx.cif.reader(file_path=cif_path).model()
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(
    model_input=pdb_inp,
    restraint_objects=[(os.path.basename(cif_path), cif_object)],
    log=null_out())
  model.process(make_restraints=True)
  os.unlink(cif_path)
  return model


def test_const_proxy_populates_b1():
  """An H whose only third-neighbor source is CONST gets b1.dihedral_ideal
  set to the CONST value_angle."""
  model = _build_model_with_custom_cif(pdb_str, custom_cif)
  g = model.get_restraints_manager().geometry
  assert len(g.const_dihedral_proxies) >= 1, \
    "CONST torsion should land in const_dihedral_proxies (got %d)" \
    % len(g.const_dihedral_proxies)

  # Verify canonicalization: the registry stores HX at i4 (outer swap from cif's
  # 'HX X Y Z' to '(Z, X, Y, HX)'). This is what makes the test exercise the
  # negation branch in _set_b1_from_const_proxy. If this assert fails, the test
  # below is not testing what it claims.
  atoms = model.get_hierarchy().atoms()
  name_by_iseq = {a.i_seq: a.name.strip() for a in atoms}
  dp = g.const_dihedral_proxies[0]
  assert name_by_iseq[dp.i_seqs[3]] == "HX", \
    "expected HX at i4 after registry canonicalization; got %s" \
    % name_by_iseq[dp.i_seqs[3]]

  # Build connectivity manager and verify b1
  from mmtbx.hydrogens import connectivity
  cm = connectivity.determine_connectivity(
    pdb_hierarchy=model.get_hierarchy(),
    geometry_restraints=g)
  hx_iseq = next(a.i_seq for a in atoms if a.name.strip() == "HX")
  nb = cm.h_connectivity[hx_iseq]
  assert nb is not None, "HX should have a connectivity entry"
  assert "iseq" in nb.b1, "b1 should have iseq populated"
  assert "dihedral_ideal" in nb.b1, \
    "CONST should populate b1.dihedral_ideal"
  assert abs(nb.b1["dihedral_ideal"] - 60.0) < 0.01, \
    "b1.dihedral_ideal should equal the CONST value_angle (60.0); got %s" \
    % nb.b1["dihedral_ideal"]


# Same monomer but with both a Var and a CONST torsion for the same H.
# Var should win (last-write semantics).
custom_cif_var_and_const = """\
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 ZZW ZZW test_var_const non-polymer 4 3 .

data_comp_ZZW
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
 ZZW X  C  CH1  0.00
 ZZW Y  C  C    0.00
 ZZW Z  C  CH3  0.00
 ZZW HX H  H    0.00

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 ZZW X  Y   single 1.520 0.020
 ZZW Y  Z   single 1.520 0.020
 ZZW X  HX  single 0.970 0.020

loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
 ZZW Y  X  HX  109.500 1.500
 ZZW X  Y  Z   111.000 1.500

loop_
_chem_comp_tor.comp_id
_chem_comp_tor.id
_chem_comp_tor.atom_id_1
_chem_comp_tor.atom_id_2
_chem_comp_tor.atom_id_3
_chem_comp_tor.atom_id_4
_chem_comp_tor.value_angle
_chem_comp_tor.value_angle_esd
_chem_comp_tor.period
 ZZW CONST_01 HX X Y Z  60.000 0.000  0
 ZZW Var_01   HX X Y Z 180.000 10.0   1
"""

pdb_str_zzw = pdb_str.replace("ZZZ", "ZZW")


def test_var_overrides_const():
  """When both a Var and a CONST torsion reference the same H, the Var
  value wins (last-write semantics in find_third_neighbors)."""
  model = _build_model_with_custom_cif(pdb_str_zzw, custom_cif_var_and_const)
  g = model.get_restraints_manager().geometry
  assert len(g.const_dihedral_proxies) >= 1, \
    "CONST torsion should still register (got %d)" \
    % len(g.const_dihedral_proxies)
  assert len(g.dihedral_proxies) >= 1, \
    "Var torsion should register (got %d)" % len(g.dihedral_proxies)

  from mmtbx.hydrogens import connectivity
  cm = connectivity.determine_connectivity(
    pdb_hierarchy=model.get_hierarchy(),
    geometry_restraints=g)
  atoms = model.get_hierarchy().atoms()
  hx_iseq = next(a.i_seq for a in atoms if a.name.strip() == "HX")
  nb = cm.h_connectivity[hx_iseq]
  assert nb is not None, "HX should have a connectivity entry"
  assert "dihedral_ideal" in nb.b1, "b1 must have dihedral_ideal populated"

  # CONST value is 60.0 (period=0). Var has angle_ideal=180.0, period=1, so
  # its path produces a value near +/-180 (the unique periodicity-1 minimum).
  # If the Var override fired, dihedral_ideal lands near +/-180 — clearly
  # distinct from CONST's 60.0. (Period 1 chosen deliberately: with period 3,
  # 60 is itself a periodicity-equivalent of 180 and the values would alias.)
  assert abs(abs(nb.b1["dihedral_ideal"]) - 180.0) < 5.0, \
    "Var torsion should override CONST and produce ~+/-180 deg; got %s" \
    % nb.b1["dihedral_ideal"]


# Snapshot of chem_data/geostd/e/data_EKB.cif, pinned here so this test
# always exercises the CONST pathway. Without pinning, a future geostd
# update that flips CONST_31..34 to Var torsions would silently make the
# test pass via the Var path while no longer testing what it claims.
# CONST_31..34 fix HN1E/HN1A/HN1F/HN1B (the four primary-amine H atoms)
# with period=0, ESD=0 torsions — the only third-neighbor source these
# atoms have, so pre-fix they get dropped by reduce_hydrogen.
EKB_CIF = """\
# electronic Ligand Builder and Optimisation Workbench (eLBOW)
#   - a module of PHENIX version dev-svn-
#   - file written: Mon Apr 27 16:32:00 2020
#   Inital geometry file: a 117 line input string
#   Ligand name: 6-ethyl-5-[3-(3,4,5-trimethoxyphenyl)prop-1-yn-1-yl]pyrimidine-2,4-diamine
#   Quantum optimisation: True
#   Method: PBEh-3c
#   Random seed: 3628800
#   SMILES string: CCc1nc(N)nc(N)c1C#CCc2cc(OC)c(OC)c(OC)c2
#
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
_chem_comp.initial_date
_chem_comp.modified_date
_chem_comp.source
 EKB  EKB  6-ethyl-5-[3-(3,4,5-trimethoxyphenyl)prop-1-yn-1-yl]pyrimidine-2,4-diamine  ligand  47  25  .  2021-08-08  2023-10-17
;
Directly from eLBOW using geometry from QM method PBEh-3c with CPCM solvent model
Validated by Mogul as OK
Added dihedrals for adding hydrogens : 2023-10-17
;

data_comp_EKB
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
 EKB  N1    N  N     0  -0.507  -3.9495   1.2492  -1.5622
 EKB  C2    C  CR6   0   0.542  -5.1272   0.8154  -1.1135
 EKB  N3    N  N     0  -0.501  -5.3324   0.0478  -0.0358
 EKB  C4    C  CR6   0   0.271  -4.2568  -0.3174   0.6502
 EKB  C5    C  CR6   0  -0.320  -2.9708   0.0723   0.2831
 EKB  C6    C  CR6   0   0.468  -2.8717   0.8863  -0.8776
 EKB  C1A   C  CH3   0  -0.666  -4.3162  -2.6880   1.4375
 EKB  C1B   C  CH3   0  -0.438   3.5604   3.1529   1.6379
 EKB  C1C   C  CH3   0  -0.438   3.0103  -3.1025  -1.7846
 EKB  C1D   C  CH3   0  -0.408   6.0311   0.1887  -0.8508
 EKB  N1E   N  NH2   0  -0.822  -6.2142   1.1871  -1.8062
 EKB  N1F   N  NH2   0  -0.824  -1.6807   1.3037  -1.3215
 EKB  C1G   C  CSP   0   0.028  -0.7675  -0.6011   1.5369
 EKB  C1H   C  CSP   0   0.147  -1.8032  -0.3099   0.9880
 EKB  C1I   C  CR16  0  -0.459   2.2638   0.6743   1.4256
 EKB  C1J   C  CR16  0  -0.456   2.0817  -1.4380   0.2744
 EKB  C1K   C  CH2   0  -0.521  -4.4780  -1.2184   1.8274
 EKB  C1L   C  CH2   0  -0.684   0.5030  -0.9525   2.1609
 EKB  O1O   O  O2    0  -0.464   3.9593   2.2208   0.6524
 EKB  O1P   O  O2    0  -0.463   3.6075  -1.8403  -1.5625
 EKB  O1Q   O  O2    0  -0.482   4.7267   0.5586  -1.2770
 EKB  C1S   C  CR6   0   0.133   1.6606  -0.5650   1.2676
 EKB  C1V   C  CR6   0   0.342   3.3069   1.0480   0.5795
 EKB  C1W   C  CR6   0   0.342   3.1225  -1.0683  -0.5749
 EKB  C1Y   C  CR6   0   0.163   3.7323   0.1811  -0.4298
 EKB  H1A   H  HCH3  0   0.238  -4.5157  -3.3335   2.2927
 EKB  H1AA  H  HCH3  0   0.229  -5.0097  -2.9618   0.6422
 EKB  H1AB  H  HCH3  0   0.231  -3.3041  -2.8958   1.0900
 EKB  H1B   H  HCH3  0   0.229   2.5206   3.4652   1.5099
 EKB  H1BA  H  HCH3  0   0.229   3.6883   2.7594   2.6497
 EKB  H1BB  H  HCH3  0   0.261   4.2031   4.0213   1.5169
 EKB  H1C   H  HCH3  0   0.229   1.9565  -3.0122  -2.0600
 EKB  H1CA  H  HCH3  0   0.261   3.5508  -3.5569  -2.6110
 EKB  H1CB  H  HCH3  0   0.229   3.0908  -3.7543  -0.9111
 EKB  H1D   H  HCH3  0   0.220   6.2945   0.6564   0.1011
 EKB  H1DA  H  HCH3  0   0.221   6.1341  -0.8944  -0.7505
 EKB  H1DB  H  HCH3  0   0.236   6.7284   0.5356  -1.6118
 EKB  HN1E  H  HNH2  0   0.397  -6.1256   1.7517  -2.6314
 EKB  HN1A  H  HNH2  0   0.397  -7.1273   0.8825  -1.5213
 EKB  HN1F  H  HNH2  0   0.400  -0.8378   1.1133  -0.8085
 EKB  HN1B  H  HNH2  0   0.398  -1.6271   1.9183  -2.1147
 EKB  H1I   H  HCR6  0   0.276   1.9131   1.3348   2.2075
 EKB  H1J   H  HCR6  0   0.276   1.5959  -2.3996   0.1755
 EKB  H1K   H  HCH2  0   0.244  -3.7729  -0.9738   2.6238
 EKB  H1KA  H  HCH2  0   0.249  -5.4843  -1.0485   2.2116
 EKB  H1L   H  HCH2  0   0.283   0.5941  -0.4546   3.1297
 EKB  H1LA  H  HCH2  0   0.286   0.5268  -2.0269   2.3615

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
_chem_comp_bond.value_dist_neutron
 EKB  C2   N1    aromatic  1.333  0.020  1.333
 EKB  C6   N1    aromatic  1.327  0.020  1.327
 EKB  N3   C2    aromatic  1.339  0.020  1.339
 EKB  C2   N1E   single    1.341  0.020  1.341
 EKB  C4   N3    aromatic  1.327  0.020  1.327
 EKB  C1K  C4    single    1.499  0.020  1.499
 EKB  C4   C5    aromatic  1.393  0.020  1.393
 EKB  C1H  C5    single    1.416  0.020  1.416
 EKB  C5   C6    aromatic  1.421  0.020  1.421
 EKB  C6   N1F   single    1.338  0.020  1.338
 EKB  C1A  C1K   single    1.529  0.020  1.529
 EKB  C1A  H1A   single    0.970  0.020  1.090
 EKB  C1A  H1AA  single    0.970  0.020  1.090
 EKB  C1A  H1AB  single    0.970  0.020  1.090
 EKB  O1O  C1B   single    1.414  0.020  1.414
 EKB  C1B  H1B   single    0.970  0.020  1.090
 EKB  C1B  H1BA  single    0.970  0.020  1.090
 EKB  C1B  H1BB  single    0.970  0.020  1.090
 EKB  O1P  C1C   single    1.414  0.020  1.414
 EKB  C1C  H1C   single    0.970  0.020  1.090
 EKB  C1C  H1CA  single    0.970  0.020  1.090
 EKB  C1C  H1CB  single    0.970  0.020  1.090
 EKB  C1D  O1Q   single    1.421  0.020  1.421
 EKB  C1D  H1D   single    0.970  0.020  1.090
 EKB  C1D  H1DA  single    0.970  0.020  1.090
 EKB  C1D  H1DB  single    0.970  0.020  1.090
 EKB  N1E  HN1E  single    0.860  0.020  1.020
 EKB  N1E  HN1A  single    0.860  0.020  1.020
 EKB  N1F  HN1F  single    0.860  0.020  1.020
 EKB  N1F  HN1B  single    0.860  0.020  1.020
 EKB  C1L  C1G   single    1.458  0.020  1.458
 EKB  C1G  C1H   triple    1.208  0.020  1.208
 EKB  C1V  C1I   aromatic  1.394  0.020  1.394
 EKB  C1I  C1S   aromatic  1.387  0.020  1.387
 EKB  C1I  H1I   single    0.930  0.020  1.080
 EKB  C1W  C1J   aromatic  1.393  0.020  1.393
 EKB  C1J  C1S   aromatic  1.388  0.020  1.388
 EKB  C1J  H1J   single    0.930  0.020  1.080
 EKB  C1K  H1K   single    0.970  0.020  1.090
 EKB  C1K  H1KA  single    0.970  0.020  1.090
 EKB  C1S  C1L   single    1.513  0.020  1.513
 EKB  C1L  H1L   single    0.970  0.020  1.090
 EKB  C1L  H1LA  single    0.970  0.020  1.090
 EKB  O1O  C1V   single    1.344  0.020  1.344
 EKB  O1P  C1W   single    1.344  0.020  1.344
 EKB  O1Q  C1Y   single    1.360  0.020  1.360
 EKB  C1Y  C1V   aromatic  1.397  0.020  1.397
 EKB  C1Y  C1W   aromatic  1.398  0.020  1.398

loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
 EKB  C6    N1   C2    117.05  3.000
 EKB  N1E   C2   N3    116.76  3.000
 EKB  N1E   C2   N1    116.87  3.000
 EKB  N3    C2   N1    126.37  3.000
 EKB  C4    N3   C2    116.72  3.000
 EKB  C1K   C4   C5    120.76  3.000
 EKB  C5    C4   N3    122.35  3.000
 EKB  C1K   C4   N3    116.86  3.000
 EKB  C1H   C5   C6    120.24  3.000
 EKB  C6    C5   C4    116.09  3.000
 EKB  C1H   C5   C4    123.66  3.000
 EKB  N1F   C6   C5    120.78  3.000
 EKB  N1F   C6   N1    117.80  3.000
 EKB  C5    C6   N1    121.41  3.000
 EKB  H1AB  C1A  H1AA  108.07  3.000
 EKB  H1AB  C1A  H1A   107.89  3.000
 EKB  H1AA  C1A  H1A   107.90  3.000
 EKB  H1AB  C1A  C1K   111.27  3.000
 EKB  H1AA  C1A  C1K   111.10  3.000
 EKB  H1A   C1A  C1K   110.48  3.000
 EKB  H1BB  C1B  H1BA  108.74  3.000
 EKB  H1BB  C1B  H1B   108.73  3.000
 EKB  H1BA  C1B  H1B   108.81  3.000
 EKB  H1BB  C1B  O1O   106.40  3.000
 EKB  H1BA  C1B  O1O   112.01  3.000
 EKB  H1B   C1B  O1O   112.03  3.000
 EKB  H1CB  C1C  H1CA  108.77  3.000
 EKB  H1CB  C1C  H1C   108.77  3.000
 EKB  H1CA  C1C  H1C   108.80  3.000
 EKB  H1CB  C1C  O1P   112.07  3.000
 EKB  H1CA  C1C  O1P   106.41  3.000
 EKB  H1C   C1C  O1P   111.90  3.000
 EKB  H1DB  C1D  H1DA  108.63  3.000
 EKB  H1DB  C1D  H1D   108.55  3.000
 EKB  H1DA  C1D  H1D   108.76  3.000
 EKB  H1DB  C1D  O1Q   107.18  3.000
 EKB  H1DA  C1D  O1Q   111.84  3.000
 EKB  H1D   C1D  O1Q   111.78  3.000
 EKB  HN1A  N1E  HN1E  118.97  3.000
 EKB  HN1A  N1E  C2    120.42  3.000
 EKB  HN1E  N1E  C2    120.58  3.000
 EKB  HN1B  N1F  HN1F  118.30  3.000
 EKB  HN1B  N1F  C6    120.01  3.000
 EKB  HN1F  N1F  C6    121.21  3.000
 EKB  C1L   C1G  C1H   180.00  3.000
 EKB  C1G   C1H  C5    180.00  3.000
 EKB  H1I   C1I  C1V   121.16  3.000
 EKB  H1I   C1I  C1S   119.12  3.000
 EKB  C1V   C1I  C1S   119.71  3.000
 EKB  H1J   C1J  C1W   121.03  3.000
 EKB  H1J   C1J  C1S   119.23  3.000
 EKB  C1W   C1J  C1S   119.74  3.000
 EKB  H1KA  C1K  H1K   107.70  3.000
 EKB  H1KA  C1K  C1A   109.71  3.000
 EKB  H1K   C1K  C1A   109.45  3.000
 EKB  H1KA  C1K  C4    108.62  3.000
 EKB  H1K   C1K  C4    110.07  3.000
 EKB  C1A   C1K  C4    111.22  3.000
 EKB  H1LA  C1L  H1L   106.45  3.000
 EKB  H1LA  C1L  C1S   110.09  3.000
 EKB  H1L   C1L  C1S   110.05  3.000
 EKB  H1LA  C1L  C1G   109.53  3.000
 EKB  H1L   C1L  C1G   110.00  3.000
 EKB  C1S   C1L  C1G   110.63  3.000
 EKB  C1V   O1O  C1B   118.43  3.000
 EKB  C1W   O1P  C1C   118.41  3.000
 EKB  C1Y   O1Q  C1D   114.34  3.000
 EKB  C1L   C1S  C1J   119.58  3.000
 EKB  C1L   C1S  C1I   119.63  3.000
 EKB  C1J   C1S  C1I   120.76  3.000
 EKB  C1Y   C1V  O1O   115.65  3.000
 EKB  C1Y   C1V  C1I   120.00  3.000
 EKB  O1O   C1V  C1I   124.35  3.000
 EKB  C1Y   C1W  O1P   115.61  3.000
 EKB  C1Y   C1W  C1J   119.98  3.000
 EKB  O1P   C1W  C1J   124.41  3.000
 EKB  C1W   C1Y  C1V   119.79  3.000
 EKB  C1W   C1Y  O1Q   120.16  3.000
 EKB  C1V   C1Y  O1Q   120.04  3.000

loop_
_chem_comp_tor.comp_id
_chem_comp_tor.id
_chem_comp_tor.atom_id_1
_chem_comp_tor.atom_id_2
_chem_comp_tor.atom_id_3
_chem_comp_tor.atom_id_4
_chem_comp_tor.value_angle
_chem_comp_tor.value_angle_esd
_chem_comp_tor.period
 EKB  CONST_01  C4    N3   C2   N1      0.00   0.0  0
 EKB  CONST_02  C4    C5   C6   N1      0.00   0.0  0
 EKB  CONST_03  C5    C6   N1   C2      0.00   0.0  0
 EKB  CONST_04  C5    C4   N3   C2      0.00   0.0  0
 EKB  CONST_05  C6    N1   C2   N3      0.00   0.0  0
 EKB  CONST_06  C6    C5   C4   N3      0.00   0.0  0
 EKB  CONST_07  C1W   C1Y  C1V  C1I     0.00   0.0  0
 EKB  CONST_08  C1W   C1J  C1S  C1I     0.00   0.0  0
 EKB  CONST_09  C1V   C1Y  C1W  C1J     0.00   0.0  0
 EKB  CONST_10  C1V   C1I  C1S  C1J     0.00   0.0  0
 EKB  CONST_11  C1Y   C1V  C1I  C1S     0.00   0.0  0
 EKB  CONST_12  C1Y   C1W  C1J  C1S     0.00   0.0  0
 EKB  CONST_13  C1H   C5   C6   N1    180.00   0.0  0
 EKB  CONST_14  N1F   C6   N1   C2    180.00   0.0  0
 EKB  CONST_15  C1K   C4   N3   C2    180.00   0.0  0
 EKB  CONST_16  C1H   C5   C4   N3    180.00   0.0  0
 EKB  CONST_17  N1E   C2   N3   C4    180.00   0.0  0
 EKB  CONST_18  N1F   C6   C5   C4    180.00   0.0  0
 EKB  CONST_19  N1E   C2   N1   C6    180.00   0.0  0
 EKB  CONST_20  C1K   C4   C5   C6    180.00   0.0  0
 EKB  CONST_21  O1Q   C1Y  C1V  C1I   180.00   0.0  0
 EKB  CONST_22  O1Q   C1Y  C1W  C1J   180.00   0.0  0
 EKB  CONST_23  C1V   C1I  C1S  C1L   180.00   0.0  0
 EKB  CONST_24  C1W   C1J  C1S  C1L   180.00   0.0  0
 EKB  CONST_25  C1S   C1I  C1V  O1O   180.00   0.0  0
 EKB  CONST_26  C1W   C1Y  C1V  O1O   180.00   0.0  0
 EKB  CONST_27  C1S   C1J  C1W  O1P   180.00   0.0  0
 EKB  CONST_28  C1V   C1Y  C1W  O1P   180.00   0.0  0
 EKB  CONST_29  H1J   C1J  C1S  C1I   180.00   0.0  0
 EKB  CONST_30  H1I   C1I  C1S  C1J   180.00   0.0  0
 EKB  CONST_31  HN1E  N1E  C2   N1      0.00   0.0  0
 EKB  CONST_32  HN1A  N1E  C2   N1    180.00   0.0  0
 EKB  CONST_33  HN1F  N1F  C6   N1    180.00   0.0  0
 EKB  CONST_34  HN1B  N1F  C6   N1      0.00   0.0  0
 EKB  Var_01    C1A   C1K  C4   N3     93.79  30.0  1
 EKB  Var_02    C1I   C1V  O1O  C1B    -0.80  30.0  2
 EKB  Var_03    C1J   C1W  O1P  C1C     1.50  30.0  2
 EKB  Var_04    C1V   C1Y  O1Q  C1D    91.31  30.0  2
 EKB  Var_05    C1I   C1S  C1L  C1G    92.69  30.0  2
 EKB  Var_06    C1S   C1L  C1G  C1H    21.63  30.0  2
 EKB  Var_07    H1A   C1A  C1K  C4   -177.90  30.0  3
 EKB  Var_08    H1B   C1B  O1O  C1V   -60.83  30.0  3
 EKB  Var_09    H1C   C1C  O1P  C1W    61.63  30.0  3
 EKB  Var_10    H1D   C1D  O1Q  C1Y   -61.75  30.0  3
 EKB  Var_11    C1L   C1G  C1H  C5    -29.15  30.0  1

loop_
_chem_comp_plane_atom.comp_id
_chem_comp_plane_atom.plane_id
_chem_comp_plane_atom.atom_id
_chem_comp_plane_atom.dist_esd
 EKB  plan-1  N1    0.020
 EKB  plan-1  C2    0.020
 EKB  plan-1  N3    0.020
 EKB  plan-1  C4    0.020
 EKB  plan-1  C5    0.020
 EKB  plan-1  C6    0.020
 EKB  plan-1  N1E   0.020
 EKB  plan-1  N1F   0.020
 EKB  plan-1  C1G   0.020
 EKB  plan-1  C1H   0.020
 EKB  plan-1  C1K   0.020
 EKB  plan-2  C1I   0.020
 EKB  plan-2  C1J   0.020
 EKB  plan-2  C1L   0.020
 EKB  plan-2  O1O   0.020
 EKB  plan-2  O1P   0.020
 EKB  plan-2  O1Q   0.020
 EKB  plan-2  C1S   0.020
 EKB  plan-2  C1V   0.020
 EKB  plan-2  C1W   0.020
 EKB  plan-2  C1Y   0.020
 EKB  plan-2  H1I   0.020
 EKB  plan-2  H1J   0.020
 EKB  plan-3  C2    0.020
 EKB  plan-3  N1E   0.020
 EKB  plan-3  HN1E  0.020
 EKB  plan-3  HN1A  0.020
 EKB  plan-4  C6    0.020
 EKB  plan-4  N1F   0.020
 EKB  plan-4  HN1F  0.020
 EKB  plan-4  HN1B  0.020
"""

# EKB heavy atoms (coordinates extracted from EKB_CIF above). Hydrogens are
# omitted so reduce_hydrogen.place_hydrogens has to add them.
EKB_HEAVY_PDB = """\
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1
HETATM    1  N1  EKB A   1      -3.949   1.249  -1.562  1.00 20.00           N
HETATM    2  C2  EKB A   1      -5.127   0.815  -1.114  1.00 20.00           C
HETATM    3  N3  EKB A   1      -5.332   0.048  -0.036  1.00 20.00           N
HETATM    4  C4  EKB A   1      -4.257  -0.317   0.650  1.00 20.00           C
HETATM    5  C5  EKB A   1      -2.971   0.072   0.283  1.00 20.00           C
HETATM    6  C6  EKB A   1      -2.872   0.886  -0.878  1.00 20.00           C
HETATM    7  C1A EKB A   1      -4.316  -2.688   1.438  1.00 20.00           C
HETATM    8  C1B EKB A   1       3.560   3.153   1.638  1.00 20.00           C
HETATM    9  C1C EKB A   1       3.010  -3.103  -1.785  1.00 20.00           C
HETATM   10  C1D EKB A   1       6.031   0.189  -0.851  1.00 20.00           C
HETATM   11  N1E EKB A   1      -6.214   1.187  -1.806  1.00 20.00           N
HETATM   12  N1F EKB A   1      -1.681   1.304  -1.322  1.00 20.00           N
HETATM   13  C1G EKB A   1      -0.768  -0.601   1.537  1.00 20.00           C
HETATM   14  C1H EKB A   1      -1.803  -0.310   0.988  1.00 20.00           C
HETATM   15  C1I EKB A   1       2.264   0.674   1.426  1.00 20.00           C
HETATM   16  C1J EKB A   1       2.082  -1.438   0.274  1.00 20.00           C
HETATM   17  C1K EKB A   1      -4.478  -1.218   1.827  1.00 20.00           C
HETATM   18  C1L EKB A   1       0.503  -0.953   2.161  1.00 20.00           C
HETATM   19  O1O EKB A   1       3.959   2.221   0.652  1.00 20.00           O
HETATM   20  O1P EKB A   1       3.608  -1.840  -1.563  1.00 20.00           O
HETATM   21  O1Q EKB A   1       4.727   0.559  -1.277  1.00 20.00           O
HETATM   22  C1S EKB A   1       1.661  -0.565   1.268  1.00 20.00           C
HETATM   23  C1V EKB A   1       3.307   1.048   0.580  1.00 20.00           C
HETATM   24  C1W EKB A   1       3.123  -1.068  -0.575  1.00 20.00           C
HETATM   25  C1Y EKB A   1       3.732   0.181  -0.430  1.00 20.00           C
END
"""


def test_ekb_const_h_placed():
  """End-to-end: HN1A, HN1B, HN1E, HN1F on EKB are placed by
  reduce_hydrogen.place_hydrogens.

  The geostd data_EKB.cif declares these four H atoms via CONST_31..34
  (period=0, ESD=0) with no Var torsion alternative. Pre-fix, they are
  dropped because pdb_interpretation filters out ESD=0 torsions and the
  connectivity manager's no-dihedral fallback wipes b1.iseq. With the
  CONST plumbing in place, find_third_neighbors._set_b1_from_const_proxy
  populates b1 from the CONST proxy and the riding-H parameterizer
  succeeds.
  """
  # Use the pinned EKB_CIF rather than relying on geostd auto-discovery.
  # _build_model_with_custom_cif attaches the cif via restraint_objects;
  # place_hydrogens.run() propagates restraint_objects through its internal
  # model rebuild, so the same cif governs both restraint passes.
  model = _build_model_with_custom_cif(EKB_HEAVY_PDB, EKB_CIF)

  from mmtbx.hydrogens.reduce_hydrogen import place_hydrogens
  placer = place_hydrogens(model=model)
  placer.run()
  out_model = placer.get_model()

  out_atoms = out_model.get_hierarchy().atoms()
  out_names = set(a.name.strip() for a in out_atoms)
  expected_h = ("HN1A", "HN1B", "HN1E", "HN1F")
  for required in expected_h:
    assert required in out_names, \
      "%s missing from reduce_hydrogen output (CONST torsion regression)" \
      % required

  # Verify that this test is *actually* exercising the CONST pathway — i.e.
  # the geostd cif still routes these four H atoms through CONST (period=0,
  # ESD=0) torsions. If a future geostd update flips them to Var, placement
  # would still succeed (via the Var pathway) but the test would silently
  # stop testing what it claims to test. This assertion fires loudly so the
  # test owner knows to refresh the input or pin a local cif copy.
  out_g = out_model.get_restraints_manager().geometry
  name_by_iseq = {a.i_seq: a.name.strip() for a in out_atoms}
  in_const = set()
  for dp in out_g.const_dihedral_proxies:
    for i in dp.i_seqs:
      nm = name_by_iseq.get(i, "")
      if nm in expected_h:
        in_const.add(nm)
  missing = set(expected_h) - in_const
  assert not missing, \
    "geostd EKB cif no longer routes %s through CONST torsions; test " \
    "no longer exercises the CONST pathway. Pin a local cif copy or " \
    "refresh the test." % sorted(missing)


# Two-residue ZZZ model used by the selection test. Residue 2 is residue 1
# translated by (+10, 0, 0); the internal geometry (and hence the CONST
# value_angle = 60 deg) is identical across copies. The pdb_interpretation
# pass yields two CONST proxies, one per residue.
pdb_str_two_residues = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
HETATM    1  X   ZZZ A   1       0.000   0.000   0.000  1.00 20.00           C
HETATM    2  Y   ZZZ A   1       1.520   0.000   0.000  1.00 20.00           C
HETATM    3  Z   ZZZ A   1       2.064   1.434   0.000  1.00 20.00           C
HETATM    4  HX  ZZZ A   1      -0.314   0.769   0.583  1.00 20.00           H
HETATM    5  X   ZZZ A   2      10.000   0.000   0.000  1.00 20.00           C
HETATM    6  Y   ZZZ A   2      11.520   0.000   0.000  1.00 20.00           C
HETATM    7  Z   ZZZ A   2      12.064   1.434   0.000  1.00 20.00           C
HETATM    8  HX  ZZZ A   2       9.686   0.769   0.583  1.00 20.00           H
END
"""


def test_const_dihedral_proxies_survive_selection():
  """Non-identity model.select() correctly filters and remaps
  const_dihedral_proxies. We build a 2-residue model (2 CONST proxies),
  then select residue 2 only — forcing the GRM to (a) drop residue 1's
  proxy and (b) remap residue 2's i_seqs from [4,5,6,7] down to [0,1,2,3].
  Then verify the selected model's riding-H manager still parameterizes
  HX via the CONST path."""
  model = _build_model_with_custom_cif(pdb_str_two_residues, custom_cif)
  g_full = model.get_restraints_manager().geometry
  assert len(g_full.const_dihedral_proxies) == 2, \
    "expected one CONST proxy per residue, got %d" \
    % len(g_full.const_dihedral_proxies)
  assert model.size() == 8, \
    "expected 8 atoms across 2 residues, got %d" % model.size()

  # Keep residue 2 only (atoms 4..7).
  sel = flex.bool(model.size(), False)
  for i in range(4, 8):
    sel[i] = True
  model_sel = model.select(sel)
  assert model_sel.size() == 4

  g_sel = model_sel.get_restraints_manager().geometry
  assert len(g_sel.const_dihedral_proxies) == 1, \
    "select() should drop residue-1's proxy and keep residue-2's; got %d" \
    % len(g_sel.const_dihedral_proxies)

  # Every i_seq referenced by the surviving proxy must be remapped into
  # the new [0, 4) range — proves the GRM ran proxy_select with the
  # selection's iselection rather than carrying stale full-model i_seqs.
  dp = g_sel.const_dihedral_proxies[0]
  for i in dp.i_seqs:
    assert 0 <= i < 4, \
      "proxy i_seq %d outside selected model's [0,4) range — i_seq " \
      "remapping failed" % i

  # The CONST value_angle should round-trip unchanged — residue 2 is just
  # a translation of residue 1, so the geometry (and the cif value 60 deg)
  # are identical.
  assert abs(abs(dp.angle_ideal) - 60.0) < 0.01, \
    "surviving proxy's angle_ideal should be +/-60 (cif value); got %s" \
    % dp.angle_ideal

  # And the riding-H manager built on the selected model must still
  # parameterize HX via the CONST path.
  model_sel.setup_riding_h_manager(use_ideal_dihedral=True)
  rhm = model_sel.get_riding_h_manager()
  assert rhm is not None, "riding-H manager should be built on selected model"
  atoms = model_sel.get_hierarchy().atoms()
  hx_iseq = next(a.i_seq for a in atoms if a.name.strip() == "HX")
  assert rhm.h_parameterization[hx_iseq] is not None, \
    "HX must be parameterized after model.select() to residue 2 only"


# Synthetic ligand with two H atoms driven by *different* third-neighbor
# sources: HA gets a Var torsion (HA-W-X-Y), HD gets a CONST torsion
# (HD-Z-Y-X). Designed to exercise both branches of find_third_neighbors's
# dispatcher in a single pass, on disjoint H atoms — distinct from the
# override test (which puts both torsions on the same H).
custom_cif_mixed = """\
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 ZZM ZZM test_mixed non-polymer 6 4 .

data_comp_ZZM
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
 ZZM W  C  C    0.00
 ZZM X  C  C    0.00
 ZZM Y  C  C    0.00
 ZZM Z  C  C    0.00
 ZZM HA H  H    0.00
 ZZM HD H  H    0.00

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 ZZM W  X   single 1.520 0.020
 ZZM X  Y   single 1.520 0.020
 ZZM Y  Z   single 1.520 0.020
 ZZM W  HA  single 0.970 0.020
 ZZM Z  HD  single 0.970 0.020

loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
 ZZM HA W  X   109.500 1.500
 ZZM W  X  Y   109.500 1.500
 ZZM X  Y  Z   109.500 1.500
 ZZM Y  Z  HD  109.500 1.500

loop_
_chem_comp_tor.comp_id
_chem_comp_tor.id
_chem_comp_tor.atom_id_1
_chem_comp_tor.atom_id_2
_chem_comp_tor.atom_id_3
_chem_comp_tor.atom_id_4
_chem_comp_tor.value_angle
_chem_comp_tor.value_angle_esd
_chem_comp_tor.period
 ZZM Var_01   HA W X Y 180.000 10.0   1
 ZZM CONST_01 HD Z Y X  60.000 0.000  0
"""

pdb_str_zzm = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
HETATM    1  W   ZZM A   1       0.000   0.000   0.000  1.00 20.00           C
HETATM    2  X   ZZM A   1       1.520   0.000   0.000  1.00 20.00           C
HETATM    3  Y   ZZM A   1       2.064   1.434   0.000  1.00 20.00           C
HETATM    4  Z   ZZM A   1       3.564   1.434   0.000  1.00 20.00           C
HETATM    5  HA  ZZM A   1      -0.314   0.769   0.583  1.00 20.00           H
HETATM    6  HD  ZZM A   1       3.878   2.203   0.583  1.00 20.00           H
END
"""


def test_const_and_var_on_different_h():
  """Both connectivity branches fire in a single find_third_neighbors pass
  on disjoint H atoms: HA's b1 comes from the Var helper, HD's b1 comes
  from the CONST helper. Verifies the dispatcher correctly routes each
  proxy without cross-contamination."""
  model = _build_model_with_custom_cif(pdb_str_zzm, custom_cif_mixed)
  g = model.get_restraints_manager().geometry
  assert len(g.const_dihedral_proxies) == 1, \
    "expected exactly one CONST proxy (got %d)" % len(g.const_dihedral_proxies)
  assert len(g.dihedral_proxies) == 1, \
    "expected exactly one Var proxy (got %d)" % len(g.dihedral_proxies)

  from mmtbx.hydrogens import connectivity
  cm = connectivity.determine_connectivity(
    pdb_hierarchy=model.get_hierarchy(),
    geometry_restraints=g)
  atoms = model.get_hierarchy().atoms()
  ha_iseq = next(a.i_seq for a in atoms if a.name.strip() == "HA")
  hd_iseq = next(a.i_seq for a in atoms if a.name.strip() == "HD")

  # HA: Var path -> dihedral derived from coords with periodicity-1 ideal
  # at 180.0. The result lands at ~+/-180 (period-1 has a unique minimum).
  ha_b1 = cm.h_connectivity[ha_iseq].b1
  assert "dihedral_ideal" in ha_b1, "HA b1 should carry dihedral_ideal"
  assert abs(abs(ha_b1["dihedral_ideal"]) - 180.0) < 5.0, \
    "HA expected Var-derived ~+/-180 deg; got %s" % ha_b1["dihedral_ideal"]

  # HD: CONST path -> b1.dihedral_ideal is the cif's value_angle (+/-60
  # depending on registry canonicalization sign-flip).
  hd_b1 = cm.h_connectivity[hd_iseq].b1
  assert "dihedral_ideal" in hd_b1, "HD b1 should carry dihedral_ideal"
  assert abs(abs(hd_b1["dihedral_ideal"]) - 60.0) < 0.01, \
    "HD expected CONST-derived +/-60 deg; got %s" % hd_b1["dihedral_ideal"]


def test_refinement_unaffected_by_const_proxies():
  """The refinement target and gradients must not consume
  const_dihedral_proxies. Clearing the list and recomputing must yield an
  identical target and gradient — proving const_dihedral_proxies is purely
  a connectivity-side annotation."""
  model = _build_model_with_custom_cif(pdb_str, custom_cif)
  sites_cart = model.get_sites_cart()
  g = model.get_restraints_manager().geometry
  assert len(g.const_dihedral_proxies) >= 1

  e_with = g.energies_sites(
    sites_cart=sites_cart, compute_gradients=True)

  saved_const = g.const_dihedral_proxies
  g.const_dihedral_proxies = geometry_restraints.shared_dihedral_proxy()
  e_without = g.energies_sites(
    sites_cart=sites_cart, compute_gradients=True)
  g.const_dihedral_proxies = saved_const

  assert abs(e_with.target - e_without.target) < 1e-9, \
    "refinement target changed when const_dihedral_proxies cleared: %s vs %s" \
    % (e_with.target, e_without.target)
  grad_diff = (e_with.gradients - e_without.gradients).norm()
  assert grad_diff < 1e-9, \
    "refinement gradients changed when const_dihedral_proxies cleared: %s" \
    % grad_diff


ZZN_CIF = """\
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
ZZN ZZN "synthetic no-torsion fixture" non-polymer 5 4 .
data_comp_ZZN
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
ZZN X  C C 0.0
ZZN Y  C C 0.0
ZZN Z  C C 0.0
ZZN D  C C 0.0
ZZN HX H H 0.0
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
ZZN X HX single 1.090 0.020
ZZN X Y  single 1.500 0.020
ZZN Y Z  single 1.500 0.020
ZZN Y D  single 1.500 0.020
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
ZZN HX X Y 109.470 3.000
ZZN X  Y Z 109.470 3.000
ZZN X  Y D 109.470 3.000
ZZN Z  Y D 109.470 3.000
"""

ZZN_PDB = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
HETATM    1  X   ZZN A   1       6.000   6.000   6.000  1.00 20.00           C
HETATM    2  Y   ZZN A   1       7.500   6.000   6.000  1.00 20.00           C
HETATM    3  Z   ZZN A   1       8.250   7.299   6.000  1.00 20.00           C
HETATM    4  D   ZZN A   1       8.250   4.701   6.000  1.00 20.00           C
HETATM    5  HX  ZZN A   1       5.470   6.918   6.000  1.00 20.00           H
END
"""

def test_solo_h_no_torsion_uses_staggered_default():
    """Solo H whose cif lacks a torsion is parameterized as alg1b with
    b1.dihedral_ideal == 180.0 (staggered default)."""
    model = _build_model_with_custom_cif(ZZN_PDB, ZZN_CIF)
    model.setup_riding_h_manager(use_ideal_dihedral=True)
    rh = model.get_riding_h_manager()
    assert rh is not None

    # Find HX
    atoms = model.get_hierarchy().atoms()
    ih = None
    for atom in atoms:
        if atom.name.strip() == 'HX':
            ih = atom.i_seq
            break
    assert ih is not None, 'HX not found in model'

    # Re-derive connectivity to inspect b1 directly
    from mmtbx.hydrogens import connectivity
    cm = connectivity.determine_connectivity(
        pdb_hierarchy       = model.get_hierarchy(),
        geometry_restraints = model.get_restraints_manager().geometry,
        mon_lib_srv         = model.get_mon_lib_srv())
    b1 = cm.h_connectivity[ih].b1
    assert 'iseq' in b1, 'b1.iseq not populated for HX'
    assert 'dihedral_ideal' in b1, 'b1.dihedral_ideal not populated for HX'
    assert abs(b1['dihedral_ideal'] - 180.0) < 1e-6, \
        'expected staggered default 180.0, got %s' % b1['dihedral_ideal']

    # End-to-end: HX is parameterized (h_parameterization entry is non-None)
    h_para = rh.h_parameterization
    assert h_para[ih] is not None, 'HX not parameterized'


ZZD_CIF = """\
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
ZZD ZZD "synthetic missing-atom fixture" non-polymer 5 4 .
data_comp_ZZD
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
ZZD A  C C 0.0
ZZD B  C C 0.0
ZZD C  C C 0.0
ZZD D  C C 0.0
ZZD HB H H 0.0
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
ZZD A B  single 1.500 0.020
ZZD B C  single 1.500 0.020
ZZD C D  single 1.500 0.020
ZZD B HB single 1.090 0.020
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
ZZD A  B C 109.470 3.000
ZZD A  B HB 109.470 3.000
ZZD HB B C 109.470 3.000
ZZD B  C D 109.470 3.000
"""

# No _chem_comp_tor block: HB has no torsion, so it enters the no-dihedral
# branch in process_a0_angles_and_third_neighbors_without_dihedral, which is
# what we want to exercise.

# PDB omits atom A — only B, C, D, HB are present.
ZZD_PDB_MISSING_A = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
HETATM    1  B   ZZD A   1       7.500   6.000   6.000  1.00 20.00           C
HETATM    2  C   ZZD A   1       8.250   7.299   6.000  1.00 20.00           C
HETATM    3  D   ZZD A   1       9.750   7.299   6.000  1.00 20.00           C
HETATM    4  HB  ZZD A   1       6.970   6.918   6.000  1.00 20.00           H
END
"""

def test_missing_heavy_atom_refuses_h_placement():
    """When parent (B) is missing an expected heavy neighbor (A) in the
    model, HB is wiped — number_non_h_neighbors == 0 — so the H is not
    parameterized."""
    model = _build_model_with_custom_cif(ZZD_PDB_MISSING_A, ZZD_CIF)

    atoms = model.get_hierarchy().atoms()
    ih = None
    for atom in atoms:
        if atom.name.strip() == 'HB':
            ih = atom.i_seq
            break
    assert ih is not None, 'HB not found in model'

    from mmtbx.hydrogens import connectivity
    cm = connectivity.determine_connectivity(
        pdb_hierarchy       = model.get_hierarchy(),
        geometry_restraints = model.get_restraints_manager().geometry,
        mon_lib_srv         = model.get_mon_lib_srv())
    no = cm.h_connectivity[ih]
    assert no.number_non_h_neighbors == 0, \
        'expected HB wiped (number_non_h_neighbors=0), got %d' % \
        no.number_non_h_neighbors


if __name__ == "__main__":
  test_const_proxy_populates_b1()
  test_var_overrides_const()
  test_const_and_var_on_different_h()
  test_ekb_const_h_placed()
  test_const_dihedral_proxies_survive_selection()
  test_refinement_unaffected_by_const_proxies()
  test_solo_h_no_torsion_uses_staggered_default()
  test_missing_heavy_atom_refuses_h_placement()
  print("OK")
