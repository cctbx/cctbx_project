
from __future__ import division
from mmtbx.secondary_structure import base_pairing
import iotbx.pdb.hierarchy
from libtbx.test_utils import approx_equal
import libtbx.load_env
from libtbx import Auto
import libtbx.phil
import os

def exercise () :
  db = base_pairing.db
  gu_classes = db.get_pair_saenger_classes("G-U")
  assert (gu_classes[0][0] == "XXVII")
  if libtbx.env.has_module("probe") and libtbx.env.has_module("reduce"):
    assert (db.get_atoms("A-U", "WWT", True) == [('H61', 'O4'), ('N1', 'H3')])
    assert (db.get_atoms("DA-DT", "WWT", False) == [('N6', 'O4'), ('N1', 'N3')])
    assert (db.get_atoms("DT-DA", "WWT", False) == [('O4', 'N6'), ('N3', 'N1')])
    assert (db.get_atoms("G-C", "WWT", False) == [('O6', 'N4'), ('N1', 'N3'),
      ('N2', 'O2')])
    assert (db.get_atoms("C-G", "WWT", False) == [('N4', 'O6'), ('N3', 'N1'),
      ('O2', 'N2')])
    assert db.get_pair_type("A-U", [('H61', 'O4'), ('N1', 'H3')], True) == "XX"
    assert db.get_pair_type("A-U", [('N1', 'H3'), ('H61', 'O4')], True) == "XX"
    assert db.get_pair_type("C-G", [('N4', 'O6'), ('N3', 'N1'), ('O2', 'N2')], False) == "XIX"
  else:
    print "Skipping: probe and/or reduce not available"
  if libtbx.env.has_module("phenix_regression") :
    pdb_file = libtbx.env.find_in_repositories(
      relative_path="phenix_regression/pdb/1u8d.pdb",
      test=os.path.isfile)
    from iotbx.file_reader import any_file
    pdb_in = any_file(pdb_file).file_object
    hierarchy = pdb_in.hierarchy
    hierarchy.atoms().reset_i_seq()
    bp_phil = libtbx.phil.parse(base_pairing.dna_rna_params_str)
    params = bp_phil.fetch(source=libtbx.phil.parse("""
      base_pair {
        base1 = chain A and resseq 18
        base2 = chain A and resseq 78
        leontis_westhof_class = *Auto
      }
      base_pair {
        base1 = chain A and resseq 21
        base2 = chain A and resseq 75
        leontis_westhof_class = *Auto
      }
      base_pair {
        base1 = chain A and resseq 33
        base2 = chain A and resseq 66
        leontis_westhof_class = *Auto
      }""")).extract()
    base_pairing.identify_base_pairs(
      pdb_hierarchy=hierarchy,
      base_pairs=params.base_pair,
      use_hydrogens=False,
      distance_ideal=3.0)
    classes = [ bp.saenger_class for bp in params.base_pair ]
    assert (classes == ["XIX", "XX", "V"])
    # and now in reverse order...
    params = bp_phil.fetch(source=libtbx.phil.parse("""
      base_pair {
        base1 = chain A and resseq 78
        base2 = chain A and resseq 18
        leontis_westhof_class = *Auto
      }
      base_pair {
        base1 = chain A and resseq 75
        base2 = chain A and resseq 21
        leontis_westhof_class = *Auto
      }
      base_pair {
        base1 = chain A and resseq 66
        base2 = chain A and resseq 33
        leontis_westhof_class = *Auto
      }""")).extract()
    base_pairing.identify_base_pairs(
      pdb_hierarchy=hierarchy,
      base_pairs=params.base_pair,
      use_hydrogens=False,
      distance_ideal=3.0)
    classes = [ bp.saenger_class for bp in params.base_pair ]
    assert (classes == ["XIX", "XX", "V"])
    # TODO: test h-bond extraction in either order???
  # rotateable base
  pdb_in = iotbx.pdb.hierarchy.input(pdb_string="""\
ATOM   3988  P    DA G  10       2.095 -23.407  14.671  1.00 24.76           P
ATOM   3989  OP1  DA G  10       2.702 -22.768  13.479  1.00 24.54           O
ATOM   3990  OP2  DA G  10       1.098 -22.688  15.497  1.00 26.02           O
ATOM   3991  O5'  DA G  10       3.272 -23.890  15.629  1.00 22.57           O
ATOM   3992  C5'  DA G  10       2.981 -24.494  16.871  1.00 21.71           C
ATOM   3993  C4'  DA G  10       4.289 -24.903  17.501  1.00 21.29           C
ATOM   3994  O4'  DA G  10       4.850 -25.979  16.721  1.00 20.23           O
ATOM   3995  C3'  DA G  10       5.340 -23.803  17.483  1.00 20.46           C
ATOM   3996  O3'  DA G  10       5.492 -23.313  18.807  1.00 20.49           O
ATOM   3997  C2'  DA G  10       6.602 -24.471  16.925  1.00 18.76           C
ATOM   3998  C1'  DA G  10       6.243 -25.951  16.910  1.00 17.49           C
ATOM   3999  N9   DA G  10       6.839 -26.749  15.842  1.00 14.88           N
ATOM   4000  C8   DA G  10       6.713 -26.565  14.492  1.00 15.01           C
ATOM   4001  N7   DA G  10       7.363 -27.452  13.774  1.00 15.05           N
ATOM   4002  C5   DA G  10       7.945 -28.283  14.718  1.00 13.74           C
ATOM   4003  C6   DA G  10       8.766 -29.426  14.610  1.00 12.94           C
ATOM   4004  N6   DA G  10       9.149 -29.944  13.443  1.00 13.28           N
ATOM   4005  N1   DA G  10       9.181 -30.024  15.746  1.00 12.16           N
ATOM   4006  C2   DA G  10       8.797 -29.503  16.921  1.00 14.06           C
ATOM   4007  N3   DA G  10       8.029 -28.435  17.152  1.00 13.57           N
ATOM   4008  C4   DA G  10       7.629 -27.865  15.998  1.00 13.67           C
""")
  #open("tmp1.pdb", "w").write(pdb_in.hierarchy.as_pdb_string())
  base = pdb_in.hierarchy.only_atom_group()
  xyz = base.atoms().extract_xyz().deep_copy()
  base_pairing.flip_base(base, angle=180)
  xyz_new = base.atoms().extract_xyz()
  assert approx_equal(xyz_new.rms_difference(xyz), 2.45942)
  #open("tmp2.pdb", "w").write(pdb_in.hierarchy.as_pdb_string())
  # TODO other bases?
  print "OK"

if __name__ == "__main__" :
  exercise()
