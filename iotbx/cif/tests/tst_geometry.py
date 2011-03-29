from cctbx import crystal, sgtbx, xray
from cctbx.eltbx import covalent_radii
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

def exercise_hbond_as_cif_loop():
  xs = sucrose()
  radii = [
    covalent_radii.table(elt).radius() for elt in
    xs.scattering_type_registry().type_index_pairs_as_dict() ]
  asu_mappings = xs.asu_mappings(
    buffer_thickness=2*max(radii) + 0.5)
  pair_asu_table = crystal.pair_asu_table(asu_mappings)
  pair_asu_table.add_covalent_pairs(
    xs.scattering_types(),
    tolerance=0.5)
  hbonds = [
    geometry.hbond(1,5, sgtbx.rt_mx('-X,0.5+Y,2-Z')),
    geometry.hbond(5,14, sgtbx.rt_mx('-X,-0.5+Y,1-Z')),
    geometry.hbond(7,10, sgtbx.rt_mx('1+X,+Y,+Z')),
    geometry.hbond(10,0),
    geometry.hbond(12,14, sgtbx.rt_mx('-1-X,0.5+Y,1-Z')),
    geometry.hbond(14,12, sgtbx.rt_mx('-1-X,-0.5+Y,1-Z')),
    geometry.hbond(16,7)
  ]
  loop = geometry.hbonds_as_cif_loop(
    hbonds, pair_asu_table, xs.scatterers().extract_labels(),
    sites_frac=xs.sites_frac()).loop
  s = StringIO()
  print >> s, loop
  assert not show_diff(s.getvalue(), """\
loop_
  _geom_hbond_atom_site_label_D
  _geom_hbond_atom_site_label_H
  _geom_hbond_atom_site_label_A
  _geom_hbond_distance_DH
  _geom_hbond_distance_HA
  _geom_hbond_distance_DA
  _geom_hbond_angle_DHA
  _geom_hbond_site_symmetry_A
   O2 H2 O4 0.8200 2.0636 2.8635 165.0 2_557
   O4 H4 O9 0.8200 2.0559 2.8736 174.9 2_546
   O5 H5 O7 0.8200 2.0496 2.8589 169.0 1_655
   O7 H7 O1 0.8200 2.0573 2.8617 166.8 .
   O8 H8 O9 0.8200 2.1407 2.8943 152.8 2_456
   O9 H9 O8 0.8200 2.1031 2.8943 162.1 2_446
   O10 H10 O5 0.8200 2.0167 2.7979 159.1 .

""")


def sucrose():
  return xray.structure(
    crystal_symmetry=crystal.symmetry(
      unit_cell=(7.783, 8.7364, 10.9002, 90, 102.984, 90),
      space_group_symbol='hall:  P 2yb'),
    scatterers=flex.xray_scatterer((
      xray.scatterer( #0
                      label='O1',
                      site=(-0.131694, 0.935444, 0.877178),
                      u=(0.000411, 0.000231, 0.000177,
                         -0.000017, 0.000112, -0.000003)),
      xray.scatterer( #1
                      label='O2',
                      site=(-0.214093, 0.788106, 1.081133),
                      u=(0.000809, 0.000424, 0.000219,
                       0.000105, 0.000177, 0.000051)),
      xray.scatterer( #2
                      label='H2',
                      site=(-0.227111, 0.876745, 1.102068),
                      u=0.050326),
      xray.scatterer( #3
                      label='O3',
                      site=(-0.145099, 0.520224, 0.848374),
                      u=(0.000860, 0.000240, 0.000587,
                         -0.000149, 0.000303, -0.000093)),
      xray.scatterer( #4
                      label='H3',
                      site=(-0.073357, 0.454114, 0.840907),
                      u=0.064297),
      xray.scatterer( #5
                      label='O4',
                      site=(0.202446, 0.586695, 0.808946),
                      u=(0.000848, 0.000352, 0.000315,
                         0.000294, 0.000256, 0.000111)),
      xray.scatterer( #6
                      label='H4',
                      site=(0.234537, 0.569702, 0.743559),
                      u=0.052973),
      xray.scatterer( #7
                      label='O5',
                      site=(0.247688, 0.897720, 0.729153),
                      u=(0.000379, 0.000386, 0.000270,
                         -0.000019, 0.000106, 0.000001)),
      xray.scatterer( #8
                      label='H5',
                      site=(0.323942, 0.957695, 0.764491),
                      u=0.040241),
      xray.scatterer( #9
                      label='O6',
                      site=(-0.183682, 1.239489, 0.712253),
                      u=(0.000240, 0.000270, 0.000208,
                         0.000035, 0.000002, -0.000057)),
      xray.scatterer( #10
                      label='O7',
                      site=(-0.460944, 1.095427, 0.826725),
                      u=(0.000454, 0.000438, 0.000330,
                         0.000029, 0.000103, 0.000091)),
      xray.scatterer( #11
                      label='H7',
                      site=(-0.360598, 1.061356, 0.849066),
                      u=0.048134),
      xray.scatterer( #12
                      label='O8',
                      site=(-0.590098, 1.234857, 0.477504),
                      u=(0.000341, 0.000376, 0.000283,
                         0.000054, -0.000022, 0.000051)),
      xray.scatterer( #13
                      label='H8',
                      site=(-0.596179, 1.326096, 0.493669),
                      u=0.041878),
      xray.scatterer( #14
                      label='O9',
                      site=(-0.294970, 1.016028, 0.425505),
                      u=(0.000509, 0.000266, 0.000191,
                         -0.000031, 0.000056, -0.000050)),
      xray.scatterer( #15
                      label='H9',
                      site=(-0.305666, 0.937198, 0.463939),
                      u=0.035859),
      xray.scatterer( #16
                      label='O10',
                      site=(0.121231, 1.098881, 0.529558),
                      u=(0.000428, 0.000442, 0.000264,
                         0.000037, 0.000149, 0.000061)),
      xray.scatterer( #17
                      label='H10',
                      site=(0.164296, 1.026179, 0.573427),
                      u=0.042667),
      xray.scatterer( #18
                      label='O11',
                      site=(-0.108540, 0.986958, 0.671415),
                      u=(0.000301, 0.000160, 0.000170,
                         -0.000018, 0.000028, 0.000003)),
      xray.scatterer( #19
                      label='C1',
                      site=(-0.205698, 0.782086, 0.859435),
                      u=(0.000378, 0.000249, 0.000211,
                         -0.000046, 0.000020, 0.000007)),
      xray.scatterer( #20
                      label='H1',
                      site=(-0.282033, 0.773748, 0.774955),
                      u=0.033128),
      xray.scatterer( #21
                      label='C2',
                      site=(-0.315620, 0.763019, 0.956870),
                      u=(0.000478, 0.000387, 0.000268,
                         0.000006, 0.000115, 0.000085)),
      xray.scatterer( #22
                      label='H2B',
                      site=(-0.413052, 0.834909, 0.939327),
                      u=0.042991),
      xray.scatterer( #23
                      label='C3',
                      site=(-0.057676, 0.663433, 0.874279),
                      u=(0.000481, 0.000203, 0.000211,
                         -0.000039, 0.000097, -0.000026)),
      xray.scatterer( #24
                      label='H3B',
                      site=(0.010861, 0.664198, 0.961490),
                      u=0.033010),
      xray.scatterer( #25
                      label='C4',
                      site=(0.064486, 0.697929, 0.785789),
                      u=(0.000516, 0.000251, 0.000171,
                         0.000108, 0.000066, 0.000006)),
      xray.scatterer( #26
                      label='H2A',
                      site=(-0.364431, 0.660423, 0.951038),
                      u=0.042991),
      xray.scatterer( #27
                      label='H4B',
                    site=(-0.001643, 0.690588, 0.698163),
                    u=0.034151),
      xray.scatterer( #28
                      label='C5',
                      site=(0.134552, 0.859023, 0.812706),
                      u=(0.000346, 0.000319, 0.000114,
                         0.000016, 0.000027, 0.000034)),
      xray.scatterer( #29
                      label='H5B',
                      site=(0.204847, 0.862433, 0.899373),
                      u=0.028874),
      xray.scatterer( #30
                      label='C6',
                      site=(-0.013996, 0.975525, 0.800238),
                      u=(0.000321, 0.000194, 0.000125,
                         -0.000031, 0.000037, -0.000006)),
      xray.scatterer( #31
                      label='H6',
                      site=(0.037624, 1.075652, 0.827312),
                      u=0.023850),
      xray.scatterer( #32
                      label='C7',
                      site=(-0.130206, 1.141642, 0.624140),
                      u=(0.000313, 0.000143, 0.000184,
                         -0.000000, 0.000045, -0.000001)),
      xray.scatterer( #33
                      label='C8',
                      site=(0.043311, 1.203061, 0.603290),
                      u=(0.000354, 0.000269, 0.000247,
                         -0.000020, 0.000085, 0.000024)),
      xray.scatterer( #34
                      label='H8A',
                      site=(0.023319, 1.300964, 0.560480),
                      u=0.034035),
      xray.scatterer( #35
                      label='H8B',
                      site=(0.123967, 1.219261, 0.684061),
                      u=0.034035),
      xray.scatterer( #36
                      label='C9',
                      site=(-0.285452, 1.143037, 0.507357),
                      u=(0.000294, 0.000203, 0.000181,
                         0.000009, 0.000062, 0.000042)),
      xray.scatterer( #37
                      label='H9B',
                      site=(-0.273079, 1.235040, 0.458631),
                      u=0.026214),
      xray.scatterer( #38
                      label='C10',
                      site=(-0.444628, 1.167026, 0.565148),
                      u=(0.000289, 0.000234, 0.000210,
                         0.000018, 0.000036, 0.000040)),
      xray.scatterer( #39
                      label='H10B',
                      site=(-0.480820, 1.069195, 0.595409),
                      u=0.029493),
      xray.scatterer( #40
                      label='C11',
                      site=(-0.371995, 1.272625, 0.676806),
                      u=(0.000324, 0.000199, 0.000264,
                         0.000054, 0.000089, 0.000009)),
      xray.scatterer( #41
                      label='H11',
                      site=(-0.388144, 1.379127, 0.648319),
                      u=0.031434),
      xray.scatterer( #42
                      label='C12',
                      site=(-0.453277, 1.252159, 0.788677),
                      u=(0.000396, 0.000422, 0.000263,
                         0.000043, 0.000125, -0.000061)),
      xray.scatterer( #43
                      label='H12B',
                      site=(-0.571879, 1.293694, 0.768457),
                      u=0.041343),
      xray.scatterer( #44
                      label='H12A',
                      site=(-0.385488, 1.310416, 0.858834),
                      u=0.041343)
    )))

if __name__ == '__main__':
  exercise_hbond_as_cif_loop()
  exercise_cif_from_cctbx()
  print "OK"
