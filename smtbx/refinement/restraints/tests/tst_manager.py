from __future__ import absolute_import, division, print_function

from smtbx.refinement import restraints
from smtbx.refinement.restraints import adp_restraints
import smtbx.development
from cctbx import crystal
from cctbx import sgtbx
from cctbx.array_family import flex
import cctbx.geometry_restraints
from six.moves import cStringIO as StringIO
import sys

def exercise_manager():
  xray_structure = smtbx.development.sucrose()
  xray_structure.scatterers()[10].set_use_u_iso_only()
  asu_mappings = xray_structure.asu_mappings(buffer_thickness=3.5)
  pair_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
  scattering_types = xray_structure.scatterers().extract_scattering_types()
  pair_asu_table.add_covalent_pairs(
    scattering_types, exclude_scattering_types=flex.std_string(("H","D")))
  pair_sym_table = pair_asu_table.extract_pair_sym_table()
  # setup adp restraint proxies
  adp_similarity_proxies = \
    adp_restraints.adp_similarity_restraints(
      pair_sym_table=pair_sym_table).proxies
  rigid_bond_proxies = \
    adp_restraints.rigid_bond_restraints(
      pair_sym_table=pair_sym_table).proxies
  rigu_proxies = \
    adp_restraints.rigu_restraints(
      pair_sym_table=pair_sym_table).proxies
  isotropic_adp_proxies = \
    adp_restraints.isotropic_adp_restraints(
      xray_structure=xray_structure,
      pair_sym_table=pair_sym_table).proxies
  bond_proxies = cctbx.geometry_restraints.shared_bond_simple_proxy()
  bond_proxies.append(
    cctbx.geometry_restraints.bond_simple_proxy(
      i_seqs=(3,23), distance_ideal=1.44, weight=2))
  bond_proxies.append(
    cctbx.geometry_restraints.bond_simple_proxy(
      i_seqs=(5,25), distance_ideal=1.44, weight=2))
  bond_proxies.append(
    cctbx.geometry_restraints.bond_simple_proxy(
      i_seqs=(1,21), distance_ideal=1.44, weight=2))
  angle_proxies = cctbx.geometry_restraints.shared_angle_proxy()
  angle_proxies.append(
    cctbx.geometry_restraints.angle_proxy(
      i_seqs=(25,28,30),angle_ideal=110,weight=2))
  angle_proxies.append(
    cctbx.geometry_restraints.angle_proxy(
      i_seqs=(23,25,28),angle_ideal=110,weight=2))
  angle_proxies.append(
    cctbx.geometry_restraints.angle_proxy(
      i_seqs=(19,23,25),angle_ideal=110,weight=2))
  bond_similarity_proxies=cctbx.geometry_restraints.shared_bond_similarity_proxy()
  bond_similarity_proxies.append(
    cctbx.geometry_restraints.bond_similarity_proxy(
      i_seqs=((14,36),(12,38)),
      weights=(10,10),
      sym_ops=(sgtbx.rt_mx(),sgtbx.rt_mx())))
  chirality_proxies=cctbx.geometry_restraints.shared_chirality_proxy()
  chirality_proxies.append(
    cctbx.geometry_restraints.chirality_proxy(
      i_seqs=(14,36,12,38),
      weight=10**4,
      volume_ideal=1.2,
      both_signs=False))
  # setup restraints manager
  manager = restraints.manager(
    bond_proxies=bond_proxies,
    angle_proxies=angle_proxies,
    bond_similarity_proxies=bond_similarity_proxies,
    adp_similarity_proxies=adp_similarity_proxies,
    rigid_bond_proxies=rigid_bond_proxies,
    rigu_proxies=rigu_proxies,
    isotropic_adp_proxies=isotropic_adp_proxies,
    chirality_proxies=chirality_proxies)
  sio = StringIO()
  manager.show_sorted(xray_structure, max_items=1, f=sio)
  if sys.platform.startswith("win") and sys.version_info[:2] < (2,6):
    # This appears to be a windows-specific bug with string formatting
    # for python versions prior to 2.6, where the exponent is printed
    # with 3 digits rather than 2.
    pass
  else:
    assert sio.getvalue() == """\
Bond restraints: 3
Sorted by residual:
bond O3
     C3
  ideal  model  delta    sigma   weight residual
  1.440  1.422  0.018 7.07e-01 2.00e+00 6.58e-04
... (remaining 2 not shown)

Bond angle restraints: 3
Sorted by residual:
angle C3
      C4
      C5
    ideal   model   delta    sigma   weight residual
   110.00  108.00    2.00 7.07e-01 2.00e+00 8.03e+00
... (remaining 2 not shown)

Chirality restraints: 1
Sorted by residual:
chirality O9
          C9
          O8
          C10
  both_signs  ideal   model   delta    sigma   weight residual
    False      1.20    2.52   -1.32 1.00e-02 1.00e+04 1.75e+04

Bond similarity restraints: 1
Sorted by residual:
               delta    sigma   weight rms_deltas residual sym.op.
bond O9-C9    -0.010 3.16e-01 1.00e+01   9.93e-03 9.87e-05
     O8-C10    0.010 3.16e-01 1.00e+01

ADP similarity restraints: 24
Sorted by residual:
scatterers O7
           C12
          delta    sigma   weight rms_deltas residual
 U11  -1.02e+00 8.00e-02 1.56e+02   5.93e-01 4.95e+02
 U22  -1.03e+00 8.00e-02 1.56e+02
 U33  -1.03e+00 8.00e-02 1.56e+02
 U12  -4.23e-03 8.00e-02 1.56e+02
 U13  -3.49e-03 8.00e-02 1.56e+02
 U23   5.66e-03 8.00e-02 1.56e+02
... (remaining 23 not shown)

Rigid bond restraints: 60
Sorted by residual:
scatterers O7
           C12
   delta_z    sigma   weight residual
 -6.42e-01 1.00e-02 1.00e+04 4.12e+03
... (remaining 59 not shown)

Rigu bond restraints: 60
Sorted by residual:
scatterers O2
           C2
   delta_z    sigma   weight residual
  1.27e-03 6.36e-03 2.47e+04 4.01e-02
 -9.08e-03 6.36e-03 2.47e+04 2.04e+00
 -9.08e-03 6.36e-03 2.47e+04 2.40e+00
... (remaining 59 not shown)

Isotropic ADP restraints: 22
Sorted by residual:
scatterer O3
         delta    sigma   weight rms_deltas residual
 U11  1.20e-03 2.00e-01 2.50e+01   1.34e-02 4.06e-02
 U22 -2.46e-02 2.00e-01 2.50e+01
 U33  2.34e-02 2.00e-01 2.50e+01
 U12 -8.14e-03 2.00e-01 2.50e+01
 U13  9.78e-03 2.00e-01 2.50e+01
 U23 -8.63e-03 2.00e-01 2.50e+01
... (remaining 21 not shown)

"""

if __name__ == '__main__':
  exercise_manager()
  print('OK')
