from libtbx.test_utils import approx_equal, show_diff
from smtbx.refinement import restraints
from smtbx.refinement.restraints import adp_restraints
from smtbx.refinement.restraints.tests import trial_structure
from cctbx import crystal
from cctbx import sgtbx
from cctbx.array_family import flex
import cctbx.geometry_restraints
import cctbx.adp_restraints.flags
from cStringIO import StringIO
import sys

def exercise_manager(verbose=0):
  xray_structure = trial_structure()
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
  isotropic_adp_proxies = \
    adp_restraints.isotropic_adp_restraints(
      xray_structure=xray_structure,
      pair_sym_table=pair_sym_table).proxies
  # setup geometry restraint proxies
  max_bond_distance = 3.5
  bond_params_table = cctbx.geometry_restraints.bond_params_table(
    xray_structure.scatterers().size())
  asu_mappings = xray_structure.asu_mappings(buffer_thickness=max_bond_distance*3)
  bond_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
  bond_simple_proxies= cctbx.geometry_restraints.shared_bond_simple_proxy()
  bond_simple_proxies.append(
    cctbx.geometry_restraints.bond_simple_proxy(
      i_seqs=(3,24), distance_ideal=1.44, weight=2))
  bond_simple_proxies.append(
    cctbx.geometry_restraints.bond_simple_proxy(
      i_seqs=(5,26), distance_ideal=1.44, weight=2))
  bond_simple_proxies.append(
    cctbx.geometry_restraints.bond_simple_proxy(
      i_seqs=(1,21), distance_ideal=1.44, weight=2))
  for proxy in bond_simple_proxies:
    i_seq, j_seq = proxy.i_seqs
    bond_params_table.update(
      i_seq=i_seq,
      j_seq=j_seq,
      params=proxy)
  angle_proxies = cctbx.geometry_restraints.shared_angle_proxy()
  angle_proxies.append(
    cctbx.geometry_restraints.angle_proxy(
      i_seqs=(26,28,30),angle_ideal=110,weight=2))
  angle_proxies.append(
    cctbx.geometry_restraints.angle_proxy(
      i_seqs=(24,26,28),angle_ideal=110,weight=2))
  angle_proxies.append(
    cctbx.geometry_restraints.angle_proxy(
      i_seqs=(19,24,26),angle_ideal=110,weight=2))
  bond_similarity_proxies=cctbx.geometry_restraints.shared_bond_similarity_proxy()
  bond_similarity_proxies.append(
    cctbx.geometry_restraints.bond_similarity_proxy(
      i_seqs=((14,36),(12,38)),
      weights=(10,10),
      sym_ops=(sgtbx.rt_mx(),sgtbx.rt_mx())))
  # setup restraints manager
  manager = restraints.manager(
    bond_params_table=bond_params_table,
    angle_proxies=angle_proxies,
    bond_similarity_proxies=bond_similarity_proxies,
    adp_similarity_proxies=adp_similarity_proxies,
    rigid_bond_proxies=rigid_bond_proxies,
    isotropic_adp_proxies=isotropic_adp_proxies)
  sio = StringIO()
  manager.show_sorted(
    sites_cart=xray_structure.sites_cart(),
    site_labels=xray_structure.scatterers().extract_labels(),
    u_cart=xray_structure.scatterers().extract_u_cart(
      unit_cell=xray_structure.unit_cell()),
    u_iso=xray_structure.scatterers().extract_u_iso(),
    use_u_aniso=xray_structure.use_u_aniso(),
    max_items=1,
    f=sio)
  energies_adps = manager.energies_adps(
    sites_cart=xray_structure.sites_cart(),
    u_cart=xray_structure.scatterers().extract_u_cart(
      unit_cell=xray_structure.unit_cell()),
    u_iso=xray_structure.scatterers().extract_u_iso(),
    use_u_aniso=xray_structure.use_u_aniso(),
    flags=None,
    compute_gradients=True,
    normalization=True)
  energies_sites = manager.energies_sites(
    sites_cart=xray_structure.sites_cart(),
    unit_cell=xray_structure.unit_cell(),
    flags=None,
    compute_gradients=True,
    normalization=True)
  energies_adps.show(f=sio)
  energies_sites.show(f=sio)
  if sys.platform.startswith("win") and sys.version_info[:2] < (2,6):
    # This appears to be a windows-specific bug with string formatting
    # for python versions prior to 2.6, where the exponent is printed
    # with 3 digits rather than 2.
    pass
  else:
    assert not show_diff(sio.getvalue(), """\
Bond restraints: 3
Sorted by residual:
bond O3
     C3
  ideal  model  delta    sigma   weight residual
  1.440  1.421  0.019 7.07e-01 2.00e+00 7.04e-04
... (remaining 2 not shown)

Bond angle restraints: 3
Sorted by residual:
angle C3
      C4
      C5
    ideal   model   delta    sigma   weight residual
   110.00  108.10    1.90 7.07e-01 2.00e+00 7.20e+00
... (remaining 2 not shown)

Bond similarity restraints: 1
Sorted by residual:
               delta    sigma   weight rms_deltas residual
bond O9-C9    -0.010 3.16e-01 1.00e+01   1.02e-02 1.04e-04
     O8-C10    0.010 3.16e-01 1.00e+01

ADP similarity restraints: 24
Sorted by residual:
scatterers O3
           C3
          delta    sigma   weight rms_deltas residual
 U11   1.70e-02 8.00e-02 1.56e+02   1.61e-02 3.65e-01
 U22   2.94e-03 8.00e-02 1.56e+02
 U33   4.24e-02 8.00e-02 1.56e+02
 U12  -6.75e-03 8.00e-02 1.56e+02
 U13   6.37e-03 8.00e-02 1.56e+02
 U23  -5.96e-03 8.00e-02 1.56e+02
... (remaining 23 not shown)

Rigid bond restraints: 60
Sorted by residual:
scatterers O7
           C12
   delta_z    sigma   weight residual
 -6.39e-01 1.00e-02 1.00e+04 4.09e+03
... (remaining 59 not shown)

Isotropic ADP restraints: 22
Sorted by residual:
scatterer O3
         delta    sigma   weight rms_deltas residual
 U11  1.12e-03 2.00e-01 2.50e+01   1.33e-02 4.01e-02
 U22 -2.44e-02 2.00e-01 2.50e+01
 U33  2.32e-02 2.00e-01 2.50e+01
 U12 -8.65e-03 2.00e-01 2.50e+01
 U13  9.38e-03 2.00e-01 2.50e+01
 U23 -8.44e-03 2.00e-01 2.50e+01
... (remaining 21 not shown)

target: 19.762
  adp_similarity_residual_sum (n=24): 2.56474
  rigid_bond_residual_sum (n=60): 4463.41
  isotropic_adp_residual_sum (n=22): 0.243628
target: 2.21912
  bond_residual_sum (n=3): 0.00152702
  angle_residual_sum (n=3): 15.5322
  bond_similarity_residual_sum (n=1): 0.000103781
  norm of gradients: 160.014
""")
  if (0 or verbose): print sio.getvalue()
  assert energies_sites.number_of_restraints == 7
  assert approx_equal(energies_sites.residual_sum, 15.533823616839785)
  assert approx_equal(energies_sites.target, 2.2191176595485405)
  assert energies_adps.number_of_restraints == 226
  assert approx_equal(energies_adps.residual_sum, 4466.2193728738066)
  assert approx_equal(energies_adps.target, 19.762032623335426)
  assert approx_equal(energies_adps.adp_similarity_deviation(),
    (0.0011807357231092134, 0.016113630010596704, 0.0062318217917095619))
  assert approx_equal(energies_adps.rigid_bond_deviation(),
    (3.3307350961686577e-05, 0.63945903723502284, 0.086249743961269706))
  assert approx_equal(energies_adps.isotropic_adp_deviation(),
    (0.0015253519360392781, 0.013343889050210288, 0.0052497588362602453))
  # test flags
  adp_flags = cctbx.adp_restraints.flags.flags(
    adp_similarity=True,
    rigid_bond=True)
  geometric_flags = cctbx.geometry_restraints.flags.flags(
    bond=True)
  energies_adps = manager.energies_adps(
    sites_cart=xray_structure.sites_cart(),
    u_cart=xray_structure.scatterers().extract_u_cart(
      unit_cell=xray_structure.unit_cell()),
    u_iso=xray_structure.scatterers().extract_u_iso(),
    use_u_aniso=xray_structure.use_u_aniso(),
    flags=adp_flags)
  assert energies_adps.isotropic_adp_proxies is None
  assert energies_adps.isotropic_adp_residual_sum == 0
  assert approx_equal(energies_adps.residual_sum, 4465.9757449271165)
  assert approx_equal(energies_adps.target, energies_adps.residual_sum)
  energies_sites = manager.energies_sites(
    sites_cart=xray_structure.sites_cart(),
    flags=geometric_flags)
  assert energies_sites.angle_proxies is None
  assert energies_sites.angle_residual_sum == 0
  assert approx_equal(energies_sites.residual_sum, 0.0015270208454704576)
  assert approx_equal(energies_sites.target, energies_sites.residual_sum)

def run(verbose=0):
  exercise_manager(verbose=verbose)
  print "OK"

if __name__ == "__main__":
  run(verbose=("--verbose" in sys.argv[1:]))
