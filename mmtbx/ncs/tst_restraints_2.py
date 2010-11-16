from mmtbx import ncs
import mmtbx.ncs.restraints
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from cctbx.array_family import flex
from libtbx.utils import Sorry, format_cpu_times
from libtbx.test_utils import Exception_expected, eps_eq, show_diff
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
  group = ncs.restraints.group.from_atom_selections(
    processed_pdb=processed_pdb,
    reference_selection_string=None,
    selection_strings=selection_strings,
    coordinate_sigma=0.05,
    b_factor_weight=0.4321,
    special_position_warnings_only=False)
  sites_cart = processed_pdb.all_chain_proxies.pdb_atoms.extract_xyz()
  ncs_operators = group.operators(sites_cart=sites_cart)
  out = StringIO()
  ncs_operators.show(sites_cart=sites_cart, out=out, prefix="{*")
  assert not show_diff(out.getvalue(), """\
{*NCS operator 1:
{*  Reference selection: "chain A"
{*      Other selection: "chain B"
{*  Number of atom pairs: 22
{*  Rotation={{-0.925533, 0.322815, -0.197938},
{*            {0.329616, 0.429511, -0.840758},
{*            {-0.186393, -0.843393, -0.503931}}
{*  Translation={{163.62}, {-13.0292}, {44.8533}}
{*  Histogram of differences:
{*    0.092573 - 0.238983: 1
{*    0.238983 - 0.385393: 2
{*    0.385393 - 0.531803: 7
{*    0.531803 - 0.678214: 3
{*    0.678214 - 0.824624: 4
{*    0.824624 - 0.971034: 5
{*  RMS difference with respect to the reference: 0.653687
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
    sites_cart=sites_cart, compute_gradients=False)
  assert energies_sites_no_gradients.number_of_restraints == 110
  assert eps_eq(energies_sites_no_gradients.residual_sum, 7014.03969257)
  assert eps_eq(energies_sites_no_gradients.target, 7014.03969257)
  assert energies_sites_no_gradients.gradients is None
  assert eps_eq(energies_sites_no_gradients.rms_with_respect_to_average,
    [0.41226641576521778, 0.38139080907663186,
     0.39748408968570492, 0.40001937328488651])
  energies_sites = ncs_operators.energies_sites(sites_cart=sites_cart)
  assert energies_sites_no_gradients.number_of_restraints \
      == energies_sites.number_of_restraints
  assert energies_sites_no_gradients.residual_sum \
      == energies_sites.residual_sum
  assert energies_sites_no_gradients.target \
      == energies_sites.target
  assert eps_eq(energies_sites.gradients.norm(), 3349.99455344)
  assert eps_eq(energies_sites.rms_with_respect_to_average,
   energies_sites_no_gradients.rms_with_respect_to_average)
  site_labels = [
    '"'+atom.pdb_label_columns()+'"' for atom in
      processed_pdb.all_chain_proxies.pdb_atoms]
  out = StringIO()
  energies_sites.show_distances_to_average(
    site_labels=site_labels, out=out, prefix="#^")
  assert not show_diff(out.getvalue(), """\
#^NCS selection: "chain A"
#^                     Distance to NCS average
#^  " N   GLN A   1 ":   0.4263
#^  " CA  GLN A   1 ":   0.2141
#^  " C   GLN A   1 ":   0.4052
...
#^  " CA  THR B   6 ":   0.4001
#^  " C   THR B   6 ":   0.6281
#^NCS selection: "chain C"
#^                     Distance to NCS average
#^  " N   GLN C   1 ":   0.4135
#^  " CA  GLN C   1 ":   0.5070
...
#^  " C   SER D   4 ":   0.6943
#^  " CA BSER D   4 ":   0.4444
#^  " N   THR D   6 ":   0.3724
#^  " CA  THR D   6 ":   0.4129
#^  " C   THR D   6 ":   0.4017
""", selections=[range(5),range(56,62),range(-5,0)])
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
  assert eng.number_of_restraints == 110
  assert eps_eq(eng.residual_sum, 1.11021057745)
  assert eps_eq(eng.target, eng.residual_sum)
  assert eng.gradients is None
  assert eps_eq(eng.rms_with_respect_to_average,
    [3.8233537528289001, 4.4894247897900934,
     3.71150443476373, 4.0839849076232442])
  energies_adp_iso = group.energies_adp_iso(u_isos=u_isos, average_power=1)
  assert energies_adp_iso.number_of_restraints == eng.number_of_restraints
  assert eps_eq(energies_adp_iso.residual_sum, eng.residual_sum)
  assert eps_eq(energies_adp_iso.target, eng.target)
  assert eps_eq(energies_adp_iso.gradients.norm(), 4.50764745473)
  assert eps_eq(energies_adp_iso.rms_with_respect_to_average,
                eng.rms_with_respect_to_average)
  out = StringIO()
  energies_adp_iso.show_differences_to_average(
    site_labels=site_labels, out=out, prefix="Y$")
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
Y$  " C   SER D   4 ":    7.29 -   10.45 =  -3.1575
Y$  " CA BSER D   4 ":    5.23 -    7.00 =  -1.7667
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
      ncs.restraints.group.from_atom_selections(
        processed_pdb=processed_pdb,
        reference_selection_string=None,
        selection_strings=selection_strings,
        coordinate_sigma=coordinate_sigma,
        b_factor_weight=b_factor_weight,
        special_position_warnings_only=False))
  energies_adp_iso = groups.energies_adp_iso(u_isos=u_isos, average_power=1)
  assert energies_adp_iso.number_of_restraints == 220
  assert eps_eq(energies_adp_iso.residual_sum, 4.2807726061)
  assert eps_eq(energies_adp_iso.target, energies_adp_iso.residual_sum)
  assert eps_eq(energies_adp_iso.gradients.norm(), 17.3806790658)
  energies_adp_iso = groups.energies_adp_iso(
    u_isos=u_isos, average_power=1, normalization=True)
  assert energies_adp_iso.number_of_restraints == 220
  assert eps_eq(energies_adp_iso.residual_sum, 4.2807726061)
  assert eps_eq(energies_adp_iso.target, energies_adp_iso.residual_sum/220)
  assert eps_eq(energies_adp_iso.gradients.norm(), 17.3806790658/220)
  for rms in energies_adp_iso.rms_with_respect_to_averages:
    if (rms is not None):
      assert eps_eq(
        rms, energies_adp_iso_no_gradients.rms_with_respect_to_average)
  assert energies_adp_iso.rms_with_respect_to_averages[2] is None
  energies_sites = groups.energies_sites(sites_cart=sites_cart)
  assert energies_sites.number_of_restraints == 220
  assert eps_eq(energies_sites.residual_sum, 8767.54961571)
  assert eps_eq(energies_sites.target, energies_sites.residual_sum)
  assert eps_eq(energies_sites.gradients.norm(), 4187.49319181)
  energies_sites = groups.energies_sites(
    sites_cart=sites_cart, normalization=True)
  assert energies_sites.number_of_restraints == 220
  assert eps_eq(energies_sites.residual_sum, 8767.54961571)
  assert eps_eq(energies_sites.target, energies_sites.residual_sum/220)
  assert eps_eq(energies_sites.gradients.norm(), 4187.49319181/220)
  for rms in energies_sites.rms_with_respect_to_averages:
    if (rms is not None):
      assert eps_eq(
        rms, energies_sites_no_gradients.rms_with_respect_to_average)
  assert energies_sites.rms_with_respect_to_averages[1] is None
  out = StringIO()
  groups.show_adp_iso_differences_to_average(
    u_isos=u_isos, site_labels=site_labels, out=out, prefix="W@")
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
  groups.show_operators(sites_cart=sites_cart, out=out, prefix="K&")
  assert not show_diff(out.getvalue(), """\
K&NCS restraint group 1:
K&  NCS operator 1:
K&    Reference selection: "chain A"
...
K&      0.824624 - 0.971034: 5
K&    RMS difference with respect to the reference: 0.653687
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
  groups.show_sites_distances_to_average(
    sites_cart=sites_cart, site_labels=site_labels, out=out, prefix="[")
  assert not show_diff(out.getvalue(), """\
[NCS restraint group 1:
[  coordinate_sigma: 0.05
[  weight:  400
[  NCS selection: "chain A"
[                       Distance to NCS average
[    " N   GLN A   1 ":   0.4263
...
[    " C   THR D   6 ":   0.4017
[NCS restraint group 2:
[  coordinate_sigma: None  =>  restraints disabled
[NCS restraint group 3:
[  coordinate_sigma: 0.1
[  weight:  100
[  NCS selection: "chain A"
[                       Distance to NCS average
[    " N   GLN A   1 ":   0.4263
...
[    " C   THR D   6 ":   0.4017
""", selections=[range(6),range(120,129),range(-1,0)])
  #
  selection = groups.selection_restrained()
  assert selection.size() == 132
  assert selection.count(True) == 110
  out = StringIO()
  processed_pdb.show_atoms_without_ncs_restraints(
    ncs_restraints_groups=groups, out=out, prefix="%&")
  assert not show_diff(out.getvalue(), """\
%&Atoms without NCS restraints:
%&MODEL        1
%&ATOM     20  N   ALA B   3     111.517   7.175  -8.669  1.00 10.73
%&ATOM     21  CA  ALA B   3     112.152   8.103  -8.026  1.00 16.28
%&ATOM     22  C   ALA B   3     111.702   8.243  -5.903  1.00  9.19
%&ATOM     24  CA  SER B   4     111.797   9.689  -4.364  1.00  8.91
%&TER
%&ATOM     38  N   ALA C   3     109.043  27.391  28.663  1.00 15.05
%&ATOM     39  CA  ALA C   3     109.073  26.531  28.433  1.00  3.14
%&ATOM     40  C   ALA C   3     108.930  26.867  26.637  1.00 15.22
%&TER
%&ATOM     57  N   ALA D   3      65.439   9.903 -12.471  1.00  7.59
%&ATOM     58  CA  ALA D   3      65.019   9.566 -11.201  1.00 15.62
%&ATOM     59  C   ALA D   3      65.679  11.045 -11.097  1.00  2.65
%&ATOM     61  CA CSER D   4      65.657  12.870  -8.333  1.00  5.84
%&TER
%&ENDMDL
%&MODEL        2
%&ATOM     86  N   ALA B   3     111.984   7.364  -8.288  1.00  3.16
%&ATOM     87  CA  ALA B   3     112.389   8.456  -7.544  1.00 12.73
%&ATOM     88  C   ALA B   3     111.615   8.267  -6.238  1.00  7.17
%&ATOM     90  CA  SER B   4     111.707   9.131  -4.317  1.00 13.83
%&TER
%&ATOM    104  N   ALA C   3     109.237  27.879  29.334  1.00  6.45
%&ATOM    105  CA  ALA C   3     109.043  26.399  28.156  1.00 15.21
%&ATOM    106  C   ALA C   3     108.983  26.760  27.178  1.00 11.73
%&TER
%&ATOM    123  N   ALA D   3      65.188  10.335 -12.959  1.00  6.94
%&ATOM    124  CA  ALA D   3      65.543   9.614 -11.242  1.00  5.18
%&ATOM    125  C   ALA D   3      65.064  11.402 -10.648  1.00 16.01
%&ATOM    127  CA CSER D   4      65.127  12.969  -9.218  1.00  9.03
%&TER
%&ENDMDL
""")
  #
  for group in groups.members:
    assert group.registry.number_of_additional_isolated_sites == 0
  groups.register_additional_isolated_sites(number=10)
  for group in groups.members:
    assert group.registry.number_of_additional_isolated_sites == 10
  groups.register_additional_isolated_sites(number=3)
  for group in groups.members:
    assert group.registry.number_of_additional_isolated_sites == 13

def exercise(args):
  verbose = "--verbose" in args
  if (verbose):
    log = sys.stdout
  else:
    log = None
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  ncs_dir = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/ncs",
    test=os.path.isdir)
  if (ncs_dir is None):
    print "Skipping exercise(): input files not available"
  else:
    for file_name in ["simple.pdb", "ambiguous_alignment.pdb"]:
      processed_pdb = monomer_library.pdb_interpretation.process(
        mon_lib_srv=mon_lib_srv,
        ener_lib=ener_lib,
        file_name=os.path.join(ncs_dir, file_name),
        log=log)
      group = ncs.restraints.group.from_atom_selections(
        processed_pdb=processed_pdb,
        reference_selection_string=None,
        selection_strings=["chain A", "chain B"],
        coordinate_sigma=None,
        b_factor_weight=None,
        special_position_warnings_only=False)
      assert len(group.selection_pairs) == 1
      assert list(group.selection_pairs[0][0]) == [0,1,2,3]
      if (file_name == "simple.pdb"):
        assert list(group.selection_pairs[0][1]) == [4,5,6,7]
      else:
        assert list(group.selection_pairs[0][1]) == [4,6,7,8]
    #
    processed_pdb = monomer_library.pdb_interpretation.process(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      file_name=os.path.join(ncs_dir, "no_match.pdb"),
      log=log)
    try:
      ncs.restraints.group.from_atom_selections(
        processed_pdb=processed_pdb,
        reference_selection_string=None,
        selection_strings=["chain A", "chain B"],
        coordinate_sigma=None,
        b_factor_weight=None,
        special_position_warnings_only=False)
    except Sorry, e:
      assert not show_diff(str(e), '''\
NCS restraints selections do not produce any pairs of matching atoms:
  Reference selection: "chain A"
      Other selection: "chain B"''')
    else: raise Exception_expected
    try:
      ncs.restraints.group.from_atom_selections(
        processed_pdb=processed_pdb,
        reference_selection_string=None,
        selection_strings=["chain A", "chain C"],
        coordinate_sigma=None,
        b_factor_weight=None,
        special_position_warnings_only=False)
    except Sorry, e:
      assert not show_diff(str(e), '''\
NCS restraints selections produce only one pair of matching atoms:
  Reference selection: "chain A"
      Other selection: "chain C"''')
    else: raise Exception_expected
    processed_pdb = monomer_library.pdb_interpretation.process(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      file_name=os.path.join(ncs_dir, "special_position.pdb"),
      log=log)
    try:
      ncs.restraints.group.from_atom_selections(
        processed_pdb=processed_pdb,
        reference_selection_string=None,
        selection_strings=["chain A", "chain D"],
        coordinate_sigma=None,
        b_factor_weight=None,
        special_position_warnings_only=False)
    except Sorry, e:
      assert not show_diff(str(e), '''\
NCS selection includes an atom on a special position:
  Selection: "chain A"
    "ATOM    456  NE  ARG A  60 .*.     N  "''')
    else: raise Exception_expected
    log = StringIO()
    group = ncs.restraints.group.from_atom_selections(
      processed_pdb=processed_pdb,
      reference_selection_string=None,
      selection_strings=["chain A", "chain D"],
      coordinate_sigma=None,
      b_factor_weight=None,
      special_position_warnings_only=True,
      log=log)
    assert not show_diff(log.getvalue(), '''\
WARNING: NCS selection includes an atom on a special position:
  Selection: "chain A"
    "ATOM    456  NE  ARG A  60 .*.     N  "
''')
    assert len(group.selection_pairs) == 1
    assert list(group.selection_pairs[0][0]) == [0,1,2,3,4,5,6,8,9,10]
    assert list(group.selection_pairs[0][1]) == [11,12,13,14,15,16,17,19,20,21]
    #
    processed_pdb = monomer_library.pdb_interpretation.process(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      file_name=os.path.join(ncs_dir, "two_models_with_holes.pdb"),
      log=log)
    try:
      ncs.restraints.group.from_atom_selections(
        processed_pdb=processed_pdb,
        reference_selection_string=None,
        selection_strings=["chain A", "chain B", "chain B and resname SER"],
        coordinate_sigma=None,
        b_factor_weight=None,
        special_position_warnings_only=False)
    except Sorry, e:
      assert not show_diff(str(e), '''\
Two different NCS operators applied to same pair of atoms:
       Reference selection: "chain A"
  Previous other selection: "chain B"
   Current other selection: "chain B and resname SER"
    "ATOM      7  N   SER A   4 .*.        "
    "ATOM     23  N   SER B   4 .*.        "''')
    exercise_two_models_with_holes(processed_pdb=processed_pdb)
  print format_cpu_times()

if (__name__ == "__main__"):
  exercise(sys.argv[1:])
