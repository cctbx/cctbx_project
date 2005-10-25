from mmtbx import ncs
import mmtbx.ncs.restraints
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from cctbx import adptbx
from cctbx.array_family import flex
from libtbx.itertbx import count
from libtbx.utils import Sorry, format_cpu_times
from libtbx.test_utils import eps_eq, show_diff
import libtbx.load_env
from cStringIO import StringIO
import sys, os

def finite_difference_site_gradients(
      ncs_operators,
      sites_cart,
      sites_average,
      eps=1.e-3):
  sites_cart = sites_cart.deep_copy()
  gradients = []
  for i_site,site in enumerate(sites_cart):
    grad = []
    for i_x in xrange(3):
      t = []
      for signed_eps in [eps, -eps]:
        site_mod = list(site)
        site_mod[i_x] += signed_eps
        sites_cart[i_site] = site_mod
        energies_sites = ncs_operators.energies_sites(
          sites_cart=sites_cart,
          compute_gradients=False,
          sites_average=sites_average)
        t.append(energies_sites.residual_sum)
      grad.append((t[0]-t[1])/(2*eps))
    gradients.append(grad)
    sites_cart[i_site] = site
  return gradients

def exercise_two_models_with_holes(processed_pdb):
  selection_strings=["chain A", "chain B", "chain C", "chain D"]
  group = ncs.restraints.group(
    processed_pdb=processed_pdb,
    reference_selection_string=None,
    selection_strings=selection_strings,
    coordinate_sigma=0.05,
    b_factor_weight=0.4321)
  sites_cart = processed_pdb.all_chain_proxies.stage_1.get_sites_cart()
  ncs_operators = group.operators(sites_cart=None)
  out = StringIO()
  ncs_operators.show(sites_cart=None, out=out, prefix="{*")
  assert not show_diff(out.getvalue(), """\
{*NCS operator 1:
{*  Reference selection: "chain A"
{*      Other selection: "chain B"
{*  Number of atom pairs: 26
{*  Rotation={{-0.928964, 0.31534, -0.193873},
{*            {0.322256, 0.43122, -0.842734},
{*            {-0.182146, -0.845346, -0.502208}}
{*  Translation={{164.149}, {-12.231}, {44.3922}}
{*  Histogram of differences:
{*    0.120681 - 0.267436: 1
{*    0.267436 - 0.414191: 4
{*    0.414191 - 0.560947: 8
{*    0.560947 - 0.707702: 3
{*    0.707702 - 0.854458: 6
{*    0.854458 - 1.001213: 4
{*  RMS difference with respect to the reference: 0.641057
{*NCS operator 2:
{*  Reference selection: "chain A"
{*      Other selection: "chain C"
{*  Number of atom pairs: 32
{*  Rotation={{-0.988874, -0.139883, -0.0506023},
{*            {0.0139383, -0.425808, 0.904706},
{*            {-0.1481, 0.893935, 0.423021}}
{*  Translation={{177.315}, {18.2319}, {1.24026}}
{*  Histogram of differences:
{*    0.291817 - 0.418479: 2
{*    0.418479 - 0.545141: 9
{*    0.545141 - 0.671802: 5
{*    0.671802 - 0.798464: 10
{*    0.798464 - 0.925126: 3
{*    0.925126 - 1.051787: 3
{*  RMS difference with respect to the reference: 0.677436
{*NCS operator 3:
{*  Reference selection: "chain A"
{*      Other selection: "chain D"
{*  Number of atom pairs: 24
{*  Rotation={{0.950594, -0.191857, 0.244055},
{*            {-0.192933, -0.981014, -0.0197252},
{*            {0.243205, -0.0283355, -0.969561}}
{*  Translation={{6.77797}, {56.2476}, {-5.96327}}
{*  Histogram of differences:
{*    0.270982 - 0.414517: 3
{*    0.414517 - 0.558053: 3
{*    0.558053 - 0.701588: 7
{*    0.701588 - 0.845124: 6
{*    0.845124 - 0.988659: 2
{*    0.988659 - 1.132195: 3
{*  RMS difference with respect to the reference: 0.724248
""")
  energies_sites_no_gradients = ncs_operators.energies_sites(
    sites_cart=None, compute_gradients=False)
  assert energies_sites_no_gradients.number_of_restraints == 114
  assert eps_eq(energies_sites_no_gradients.residual_sum, 7173.85077391)
  assert eps_eq(energies_sites_no_gradients.target, 7173.85077391)
  assert energies_sites_no_gradients.gradients is None
  assert eps_eq(energies_sites_no_gradients.rms_with_respect_to_average,
    [0.41634839722020417, 0.36587351734294055,
     0.40045430793729697, 0.39662478265194723])
  energies_sites = ncs_operators.energies_sites(sites_cart=None)
  assert energies_sites_no_gradients.number_of_restraints \
      == energies_sites.number_of_restraints
  assert energies_sites_no_gradients.residual_sum \
      == energies_sites.residual_sum
  assert energies_sites_no_gradients.target \
      == energies_sites.target
  assert eps_eq(energies_sites.gradients.norm(), 3397.62130783)
  assert eps_eq(energies_sites.rms_with_respect_to_average,
   energies_sites_no_gradients.rms_with_respect_to_average)
  out = StringIO()
  energies_sites.show_distances_to_average(out=out, prefix="#^")
  assert not show_diff(out.getvalue(), """\
#^NCS selection: "chain A"
#^                     Distance to NCS average
#^  " N   GLN A   1 ":   0.4263
#^  " CA  GLN A   1 ":   0.2141
#^  " C   GLN A   1 ":   0.4052
...
#^  " CA  THR B   6 ":   0.4231
#^  " C   THR B   6 ":   0.6203
#^NCS selection: "chain C"
#^                     Distance to NCS average
#^  " N   GLN C   1 ":   0.4135
#^  " CA  GLN C   1 ":   0.5070
...
#^  " CA BSER D   4 ":   0.4004
#^  " C   SER D   4 ":   0.6999
#^  " N   THR D   6 ":   0.3685
#^  " CA  THR D   6 ":   0.4179
#^  " C   THR D   6 ":   0.4045
""", selections=[range(5),range(60,66),range(-5,0)])
  for ag,fg in zip(energies_sites.gradients,
                   finite_difference_site_gradients(
                     ncs_operators=ncs_operators,
                     sites_cart=sites_cart,
                     sites_average=energies_sites.sites_average)):
    assert eps_eq(ag, fg)
  #
  u_isos = processed_pdb.xray_structure().extract_u_iso_or_u_equiv()
  eng = group.energies_adp_iso(
    u_isos=u_isos, average_power=1, compute_gradients=False)
  energies_adp_iso_no_gradients = eng
  assert eng.number_of_restraints == 114
  assert eps_eq(eng.residual_sum, 1.13989071027)
  assert eps_eq(eng.target, eng.residual_sum)
  assert eng.gradients is None
  assert eps_eq(eng.rms_with_respect_to_average,
    [3.8343001243481694, 4.3047451177562364,
     3.7243446215734459, 4.1272033923027323])
  energies_adp_iso = group.energies_adp_iso(u_isos=u_isos, average_power=1)
  assert energies_adp_iso.number_of_restraints == eng.number_of_restraints
  assert eps_eq(energies_adp_iso.residual_sum, eng.residual_sum)
  assert eps_eq(energies_adp_iso.target, eng.target)
  assert eps_eq(energies_adp_iso.gradients.norm(), 4.56964511297)
  assert eps_eq(energies_adp_iso.rms_with_respect_to_average,
                eng.rms_with_respect_to_average)
  out = StringIO()
  energies_adp_iso.show_differences_to_average(out=out, prefix="Y$")
  assert not show_diff(out.getvalue().replace(" 7.66 = ", " 7.65 = "), """\
Y$NCS selection: "chain A"
Y$                       B-iso   NCS ave  Difference
Y$  " N   GLN A   1 ":   11.77 -   12.08 =  -0.3133
Y$  " CA  GLN A   1 ":    9.09 -    9.28 =  -0.1933
Y$  " C   GLN A   1 ":   15.40 -   12.05 =   3.3500
...
Y$  " CA  THR A   6 ":    3.98 -    7.65 =  -3.6750
Y$  " C   THR A   6 ":   14.64 -   10.99 =   3.6475
Y$NCS selection: "chain B"
Y$                       B-iso   NCS ave  Difference
Y$  " N   GLU B   2 ":   12.87 -   10.35 =   2.5225
...
Y$  " CA BSER D   4 ":    5.23 -    8.71 =  -3.4750
Y$  " C   SER D   4 ":    7.29 -   10.45 =  -3.1575
Y$  " N   THR D   6 ":    4.55 -    8.69 =  -4.1425
Y$  " CA  THR D   6 ":    8.78 -    7.65 =   1.1250
Y$  " C   THR D   6 ":   10.80 -   10.99 =  -0.1925
""", selections=[range(5),range(32,37),range(-5,0)])
  for average_power in [1, 0.69, 0.35]:
    finite_difference_gradients = flex.double()
    eps = 1.e-6
    for i_u_iso in xrange(u_isos.size()):
      rs = []
      for signed_eps in [eps, -eps]:
        u_isos_eps = u_isos.deep_copy()
        u_isos_eps[i_u_iso] += signed_eps
        energies = group.energies_adp_iso(
          u_isos=u_isos_eps,
          average_power=average_power,
          compute_gradients=False)
        rs.append(energies.residual_sum)
      finite_difference_gradients.append((rs[0]-rs[1])/(2*eps))
    energies_adp_iso = group.energies_adp_iso(
      u_isos=u_isos, average_power=average_power)
    assert eps_eq(energies_adp_iso.gradients, finite_difference_gradients)
  #
  groups = ncs.restraints.groups()
  groups.members.append(group)
  for coordinate_sigma,b_factor_weight in [
        (None,1.234),
        (0.1,None)]:
    groups.members.append(
      ncs.restraints.group(
        processed_pdb=processed_pdb,
        reference_selection_string=None,
        selection_strings=selection_strings,
        coordinate_sigma=coordinate_sigma,
        b_factor_weight=b_factor_weight))
  energies_adp_iso = groups.energies_adp_iso(u_isos=u_isos, average_power=1)
  assert energies_adp_iso.number_of_restraints == 228
  assert eps_eq(energies_adp_iso.residual_sum, 4.39521386804)
  assert eps_eq(energies_adp_iso.target, energies_adp_iso.residual_sum)
  assert eps_eq(energies_adp_iso.gradients.norm(), 17.6197309019)
  energies_adp_iso = groups.energies_adp_iso(
    u_isos=u_isos, average_power=1, normalization=True)
  assert energies_adp_iso.number_of_restraints == 228
  assert eps_eq(energies_adp_iso.residual_sum, 4.39521386804)
  assert eps_eq(energies_adp_iso.target, energies_adp_iso.residual_sum/228)
  assert eps_eq(energies_adp_iso.gradients.norm(), 17.6197309019/228)
  for rms in energies_adp_iso.rms_with_respect_to_averages:
    if (rms is not None):
      assert eps_eq(
        rms, energies_adp_iso_no_gradients.rms_with_respect_to_average)
  assert energies_adp_iso.rms_with_respect_to_averages[2] is None
  energies_sites = groups.energies_sites(sites_cart=sites_cart)
  assert energies_sites.number_of_restraints == 228
  assert eps_eq(energies_sites.residual_sum, 8967.31346739)
  assert eps_eq(energies_sites.target, energies_sites.residual_sum)
  assert eps_eq(energies_sites.gradients.norm(), 4247.02663479)
  energies_sites = groups.energies_sites(
    sites_cart=sites_cart, normalization=True)
  assert energies_sites.number_of_restraints == 228
  assert eps_eq(energies_sites.residual_sum, 8967.31346739)
  assert eps_eq(energies_sites.target, energies_sites.residual_sum/228)
  assert eps_eq(energies_sites.gradients.norm(), 4247.02663479/228)
  for rms in energies_sites.rms_with_respect_to_averages:
    if (rms is not None):
      assert eps_eq(
        rms, energies_sites_no_gradients.rms_with_respect_to_average)
  assert energies_sites.rms_with_respect_to_averages[1] is None
  out = StringIO()
  groups.show_adp_iso_differences_to_average(
    u_isos=u_isos, out=out, prefix="W@")
  assert not show_diff(out.getvalue(), """\
W@NCS restraint group 1:
W@  weight: 0.4321
W@  NCS selection: "chain A"
W@                         B-iso   NCS ave  Difference
W@    " N   GLN A   1 ":   11.77 -   12.08 =  -0.3133
...
W@    " C   THR D   6 ":   10.80 -   10.99 =  -0.1925
W@NCS restraint group 3:
W@  b_factor_weight: None  =>  restraints disabled
""", selections=[range(5),range(-3,0)])
  out = StringIO()
  groups.show_operators(sites_cart=None, out=out, prefix="K&")
  assert not show_diff(out.getvalue(), """\
K&NCS restraint group 1:
K&  NCS operator 1:
K&    Reference selection: "chain A"
...
K&      0.854458 - 1.001213: 4
K&    RMS difference with respect to the reference: 0.641057
K&  NCS operator 2:
K&    Reference selection: "chain A"
K&        Other selection: "chain C"
K&    Number of atom pairs: 32
...
K&      0.845124 - 0.988659: 2
K&      0.988659 - 1.132195: 3
K&    RMS difference with respect to the reference: 0.724248
""", selections=[range(3),range(15,21),range(-3,0)])
  out = StringIO()
  groups.show_sites_distances_to_average(sites_cart=None, out=out, prefix="[")
  assert not show_diff(out.getvalue(), """\
[NCS restraint group 1:
[  coordinate_sigma: 0.05
[  weight:  400
[  NCS selection: "chain A"
[                       Distance to NCS average
[    " N   GLN A   1 ":   0.4263
...
[    " C   THR D   6 ":   0.4045
[NCS restraint group 2:
[  coordinate_sigma: None  =>  restraints disabled
[NCS restraint group 3:
[  coordinate_sigma: 0.1
[  weight:  100
[  NCS selection: "chain A"
[                       Distance to NCS average
[    " N   GLN A   1 ":   0.4263
...
[    " C   THR D   6 ":   0.4045
""", selections=[range(6),range(124,133),range(-1,0)])
  #
  out = StringIO()
  groups.show_unrestrained_atoms(out=out, prefix="%&")
  assert not show_diff(out.getvalue(), """\
%&Atoms without NCS restraints:
%&  Model 1, PDB serial number: 1
%&    Conformer 1, PDB altLoc: "A"
%&      Chain 2, PDB chainID: "B", segID: "    "
%&        "ALA B   3 "
%&          " CA  "
%&          " C   "
%&          " N   "
%&      Chain 3, PDB chainID: "C", segID: "    "
%&        "ALA C   3 "
%&          " CA  "
%&          " C   "
%&          " N   "
%&      Chain 4, PDB chainID: "D", segID: "    "
%&        "ALA D   3 "
%&          " CA  "
%&          " C   "
%&          " N   "
%&    Conformer 3, PDB altLoc: "C"
%&      Chain 4, PDB chainID: "D", segID: "    "
%&        "SER D   4 "
%&          " CA C"
%&  Model 2, PDB serial number: 2
%&    Conformer 1, PDB altLoc: "A"
%&      Chain 2, PDB chainID: "B", segID: "    "
%&        "ALA B   3 "
%&          " CA  "
%&          " C   "
%&          " N   "
%&      Chain 3, PDB chainID: "C", segID: "    "
%&        "ALA C   3 "
%&          " CA  "
%&          " C   "
%&          " N   "
%&      Chain 4, PDB chainID: "D", segID: "    "
%&        "ALA D   3 "
%&          " CA  "
%&          " C   "
%&          " N   "
%&    Conformer 3, PDB altLoc: "C"
%&      Chain 4, PDB chainID: "D", segID: "    "
%&        "SER D   4 "
%&          " CA C"
""")

def exercise(args):
  verbose = "--verbose" in args
  if (verbose):
    log = sys.stdout
  else:
    log = None
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  ncs_dir = libtbx.env.find_in_repositories(
    relative_path="regression/ncs",
    test=os.path.isdir)
  if (ncs_dir is None):
    print "Skipping exercise(): input files not available"
  else:
    processed_pdb = monomer_library.pdb_interpretation.process(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      file_name=os.path.join(ncs_dir, "ambiguous_alignment.pdb"),
      log=log)
    try:
      ncs.restraints.group(
        processed_pdb=processed_pdb,
        reference_selection_string=None,
        selection_strings=["chain A", "chain B"],
        coordinate_sigma=None,
        b_factor_weight=None)
    except Sorry, e:
      assert not show_diff(str(e), """\
Ambiguous alignment of NCS restraints selections:
  Reference selection: "chain A"
      Other selection: "chain B"
    Number of possible alginments: 2
    Alignment 1:
      ---  GLU
      GLU  GLU
      ALA  ALA
      SER  SER
      SER  SER
      ---  THR
    Alignment 2:
      GLU  GLU
      ---  GLU
      ALA  ALA
      SER  SER
      SER  SER
      ---  THR""")
    else:
      raise RuntimeError("Exception expected.")
    #
    processed_pdb = monomer_library.pdb_interpretation.process(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      file_name=os.path.join(ncs_dir, "two_models_with_holes.pdb"),
      log=log)
    exercise_two_models_with_holes(processed_pdb=processed_pdb)
  print format_cpu_times()

if (__name__ == "__main__"):
  exercise(sys.argv[1:])
