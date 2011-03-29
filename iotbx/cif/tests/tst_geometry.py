from cctbx import crystal, xray
from cctbx.array_family import flex
from libtbx.test_utils import show_diff
from iotbx.cif import geometry
from cStringIO import StringIO

def exercise_cif_from_cctbx():
  quartz = xray.structure(
    crystal_symmetry=crystal.symmetry(
      (5.01,5.01,5.47,90,90,120), "P6222"),
    scatterers=flex.xray_scatterer([
      xray.scatterer("Si", (1/2.,1/2.,1/3.)),
      xray.scatterer("O", (0.197,-0.197,0.83333))]))
  s = StringIO()
  loop = geometry.distances_as_cif_loop(
    quartz.pair_asu_table(distance_cutoff=2),
    site_labels=quartz.scatterers().extract_labels(),
    sites_frac=quartz.sites_frac()).loop
  print >> s, loop
  assert not show_diff(s.getvalue(), """\
loop_
  _geom_bond_atom_site_label_1
  _geom_bond_atom_site_label_2
  _geom_bond_distance
  _geom_bond_site_symmetry_2
   Si O 1.6160 4_554
   Si O 1.6160 2_554
   Si O 1.6160 3_664
   Si O 1.6160 5_664

""")
  s = StringIO()
  loop = geometry.angles_as_cif_loop(
    quartz.pair_asu_table(distance_cutoff=2),
    site_labels=quartz.scatterers().extract_labels(),
    sites_frac=quartz.sites_frac()).loop
  print >> s, loop
  assert not show_diff(s.getvalue(), """\
loop_
  _geom_angle_atom_site_label_1
  _geom_angle_atom_site_label_2
  _geom_angle_atom_site_label_3
  _geom_angle
  _geom_angle_site_symmetry_1
  _geom_angle_site_symmetry_3
   O Si O 101.3 2_554 4_554
   O Si O 111.3 3_664 4_554
   O Si O 116.1 3_664 2_554
   O Si O 116.1 5_664 4_554
   O Si O 111.3 5_664 2_554
   O Si O 101.3 5_664 3_664
   Si O Si 146.9 3 5

""")

if __name__ == '__main__':
  exercise_cif_from_cctbx()
  print "OK"
