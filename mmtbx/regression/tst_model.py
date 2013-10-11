from __future__ import division
from cctbx.array_family import flex
import os
import mmtbx.model
import libtbx.load_env
from libtbx import easy_pickle
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from cStringIO import StringIO
from mmtbx import utils
from libtbx.utils import format_cpu_times, null_out
from libtbx.test_utils import approx_equal

def exercise():
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  pdb_file = libtbx.env.find_in_repositories(
                   relative_path="phenix_regression/pdb/enk.pdb", test=os.path.isfile)
  processed_pdb_file = monomer_library.pdb_interpretation.process(
                                       mon_lib_srv               = mon_lib_srv,
                                       ener_lib                  = ener_lib,
                                       file_name                 = pdb_file,
                                       raw_records               = None,
                                       force_symmetry            = True)
################
  geometry = processed_pdb_file.geometry_restraints_manager(
                                                    show_energies      = False,
                                                    plain_pairs_radius = 5.0)
  restraints_manager = mmtbx.restraints.manager(geometry      = geometry,
                                                normalization = False)
  xray_structure = processed_pdb_file.xray_structure()
  selection = flex.bool(xray_structure.scatterers().size(), True)
  mol = mmtbx.model.manager(
    restraints_manager = restraints_manager,
    xray_structure = xray_structure,
    pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy)
  mol.xray_structure.scattering_type_registry(table = "wk1995")
################



  mol_copy = mol.deep_copy()
  assert mol.number_of_ordered_solvent_molecules() == 9
  mol.write_pdb_file(out = open("test_model_out.pdb","w"))
  mol.write_pdb_file(out = open("test_model_out_nohydrogens.pdb","w"))
  mol = mol.remove_solvent()
  assert mol.number_of_ordered_solvent_molecules() == 0
  mol.write_pdb_file(out = open("test_model_out_nosolvent.pdb","w"))
  mol.write_pdb_file(out = open("test_model_out_noO.pdb","w"))
  mol.show_geometry_statistics(ignore_hd = True)
  mol.show_geometry_statistics(ignore_hd = True)

################
  geometry = processed_pdb_file.geometry_restraints_manager(
                                                    show_energies      = False,
                                                    plain_pairs_radius = 5.0)
  restraints_manager = mmtbx.restraints.manager(geometry      = geometry,
                                                normalization = False)
  xray_structure = processed_pdb_file.xray_structure()
  selection = flex.bool(xray_structure.scatterers().size(), True)
  mol_other = mmtbx.model.manager(
    restraints_manager = restraints_manager,
    xray_structure = xray_structure,
    pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy)
  mol_other.xray_structure.scattering_type_registry(table = "wk1995")
################

  mol_other.show_geometry_statistics(ignore_hd = True)
  print
  mol.show_geometry_statistics(ignore_hd = True)

#####
  class iso: pass
  iso.sphere_radius = 5.0
  iso.distance_power = 1.57
  iso.average_power = 0.58
  iso.wilson_b_weight_auto = False
  iso.wilson_b_weight = None
  iso.plain_pairs_radius = 5.0
  iso.refine_ap_and_dp = False
  iso.use_u_local_only = False
#####

  mol.show_adp_statistics()
  print
  mol.show_adp_statistics()

  rm = mol.restraints_manager

  mol_copy.write_pdb_file(out = open("XXX.pdb","w"))
  mol_copy.write_pdb_file(out = open("XXXr.pdb","w"))


def exercise_2():
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  pdb_file = libtbx.env.find_in_repositories(
                   relative_path="phenix_regression/pdb/adp_out_stat.pdb", test=os.path.isfile)
  params = monomer_library.pdb_interpretation.master_params.extract()
  params.nonbonded_weight = 16
  processed_pdb_file = monomer_library.pdb_interpretation.process(
                                       mon_lib_srv               = mon_lib_srv,
                                       ener_lib                  = ener_lib,
                                       params                    = params,
                                       file_name                 = pdb_file,
                                       raw_records               = None,
                                       force_symmetry            = True)
  geometry = processed_pdb_file.geometry_restraints_manager(
                                                    show_energies      = False,
                                                    plain_pairs_radius = 5.0)
  restraints_manager = mmtbx.restraints.manager(geometry      = geometry,
                                                normalization = False)
  xray_structure = processed_pdb_file.xray_structure()
  selection = flex.bool(xray_structure.scatterers().size(), True)
  mol = mmtbx.model.manager(
    restraints_manager = restraints_manager,
    xray_structure = xray_structure,
    pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy)
  mol.xray_structure.scattering_type_registry(table = "wk1995")

  out = StringIO()
  adp_stat = mol.show_adp_statistics(out = out)
  expected_result = \
  """|-ADP statistics-------------------------------------------------------|
| Atom    | Number of   | Isotropic or equivalent| Anisotropy lmin/max |
| type    |iso    aniso | min     max     mean   | min   max    mean   |
| - - - - |- - - - - - -| - - - - - - - - - - - -| - - - - - - - - - - |
| all     : 14     7      1.31    18.27   13.11    0.02  0.73   0.20   |
| all(noH): 9      5      1.45    18.27   15.14    0.03  0.73   0.18   |
| Sol.    : 1      1      1.45    2.00    1.73     0.73  0.73   0.73   |
| Mac.    : 8      4      15.00   18.27   17.38    0.03  0.05   0.04   |
| Hyd.    : 5      2      1.31    17.15   9.04     0.02  0.50   0.26   |
| - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  |
|    Distribution of isotropic (or equivalent) ADP for non-H atoms:    |
| Bin#      value range     #atoms | Bin#      value range     #atoms  |
|   0:     1.450 -   3.132:    2   |   5:     9.860 -  11.542:    0    |
|   1:     3.132 -   4.814:    0   |   6:    11.542 -  13.224:    0    |
|   2:     4.814 -   6.496:    0   |   7:    13.224 -  14.906:    0    |
|   3:     6.496 -   8.178:    0   |   8:    14.906 -  16.588:    1    |
|   4:     8.178 -   9.860:    0   |   9:    16.588 -  18.270:   11    |
|                            =>continue=>                              |
| - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  |
|                     Distribution of anisotropy:                      |
| Bin#      value range     #atoms | Bin#      value range     #atoms  |
|   0:     0.028 -   0.098:    4   |   5:     0.377 -   0.447:    0    |
|   1:     0.098 -   0.168:    0   |   6:     0.447 -   0.517:    0    |
|   2:     0.168 -   0.238:    0   |   7:     0.517 -   0.587:    0    |
|   3:     0.238 -   0.307:    0   |   8:     0.587 -   0.657:    0    |
|   4:     0.307 -   0.377:    0   |   9:     0.657 -   0.726:    1    |
|                            =>continue=>                              |
|----------------------------------------------------------------------|
"""
  assert out.getvalue() == expected_result
  # XXX phenix GUI support (see wxtbx.adp_statistics)
  stats = mol.adp_statistics()
  tables = stats.format_tables()
  stats_pkl = easy_pickle.dumps(stats)
  stats2 = easy_pickle.loads(stats_pkl)
  tables2 = stats.format_tables()
  assert (tables2 == tables)
  t1,t2,t3 = tables2
  assert (len(t1) == 5) and (len(t2) == 10) and (len(t3) == 10)
  assert (t1[0][1:3] == ["14","7"]) and (t1[-1][1:3] == ["5","2"])
  assert (t2[1][1] == "3.132 - 4.814")
  p1, p2 = stats.format_plots()
  y1, yrange1 = p1
  assert (len(y1) == 10)

def exercise_3():
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/arg_h_hohh1_sh.pdb",
    test=os.path.isfile)
  params = monomer_library.pdb_interpretation.master_params.extract()
  params.nonbonded_weight = 16
  processed_pdb_files_srv = utils.process_pdb_file_srv(
    pdb_interpretation_params=params)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  #
  geometry = processed_pdb_file.geometry_restraints_manager(
                                                    show_energies      = False,
                                                    plain_pairs_radius = 5.0)
  restraints_manager = mmtbx.restraints.manager(geometry      = geometry,
                                                normalization = False)
  xray_structure = processed_pdb_file.xray_structure()
  out = StringIO()
  mol = mmtbx.model.manager(
    restraints_manager = restraints_manager,
    xray_structure = xray_structure,
    pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy,
    log = out)
  #
  mol.idealize_h()
  assert out.getvalue().splitlines()[0] == \
  "X-H deviation from ideal before regularization (bond): mean= 0.201 max= 0.636"
  assert out.getvalue().splitlines()[1] == \
  "X-H deviation from ideal after  regularization (bond): mean= 0.000 max= 0.000"

def exercise_4():
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/lysozyme_noH.pdb",
    test=os.path.isfile)
  processed_pdb_files_srv = utils.process_pdb_file_srv()
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  #
  geometry = processed_pdb_file.geometry_restraints_manager(
                                                    show_energies      = False,
                                                    plain_pairs_radius = 5.0)
  restraints_manager = mmtbx.restraints.manager(geometry      = geometry,
                                                normalization = False)
  xray_structure = processed_pdb_file.xray_structure()
  out = StringIO()
  mol = mmtbx.model.manager(
    restraints_manager = restraints_manager,
    xray_structure = xray_structure,
    pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy,
    log = out)
  #
  result = mol.isolated_atoms_selection()
  solvent_sel = mol.solvent_selection()
  assert solvent_sel.count(True)+1 == result.count(True)

def exercise_convert_atom() :
  from iotbx.pdb import hierarchy
  from cctbx import crystal
  from cctbx.xray import anomalous_scatterer_group
  coords = [ (2.12, 0., 0.), (0., 2.12, 0.), (0, 0, 2.12), (0,0,0),
             (-2.12, 0, 0), (0, -2.12, 0), (0, 0, -2.12) ]
  root = hierarchy.root()
  model = hierarchy.model()
  chain = hierarchy.chain(id="S")
  root.append_model(model)
  model.append_chain(chain)
  for k, xyz in enumerate(coords) :
    rg = hierarchy.residue_group(resseq=str(k+1))
    ag = hierarchy.atom_group(resname="HOH")
    atom = hierarchy.atom()
    atom.serial = str(k+1)
    atom.name = " O  "
    atom.occ = 1.0
    atom.b = 20.0
    atom.xyz=  xyz
    atom.element = " O"
    ag.append_atom(atom)
    rg.append_atom_group(ag)
    chain.append_residue_group(rg)
  symm = crystal.symmetry(
    space_group_symbol="P1",
    unit_cell=(10,10,10,90,90,90))
  open("tmp_tst_model_5.pdb", "w").write(root.as_pdb_string(symm))
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=null_out())
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = ["tmp_tst_model_5.pdb"])
  geometry = processed_pdb_file.geometry_restraints_manager(
                                                    show_energies      = False,
                                                    plain_pairs_radius = 5.0)
  restraints_manager = mmtbx.restraints.manager(geometry      = geometry,
                                                normalization = False)
  xray_structure = processed_pdb_file.xray_structure()
  pdb_atoms = processed_pdb_file.all_chain_proxies.pdb_hierarchy.atoms()
  mol = mmtbx.model.manager(
    restraints_manager = restraints_manager,
    xray_structure = xray_structure,
    pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy,
    log = null_out())
  mol.convert_atom(
    i_seq=4,
    scattering_type="Mg2+",
    atom_name="MG",
    element="MG",
    charge=2,
    residue_name="MG",
    initial_occupancy=0.99,
    chain_id='X')
  mol.anomalous_scatterer_groups = [
    anomalous_scatterer_group(
      iselection=flex.size_t([4]),
      f_prime=0,
      f_double_prime=0,
      refine=["f_prime","f_double_prime"],
      selection_string="element MG",
      update_from_selection=True), ]
  mol.geometry_minimization(nonbonded=True)
  # if the nonbonded type is set correctly, the nonbonded restraints should
  # not push the
  for atom in mol.pdb_hierarchy(sync_with_xray_structure=True).atoms() :
    xyz_max = max([ abs(n) for n in atom.xyz])
    assert (xyz_max < 2.5)
  mol = mol.select(flex.size_t([1,2,3,4,5,6]))
  assert mol.update_anomalous_groups(out=null_out())
  isel = mol.anomalous_scatterer_groups[0].iselection
  assert list(isel) == [3]

pdb_file_exercise_h_counts="""
CRYST1    8.228   11.366   10.991  90.00  90.00  90.00 P 1
ATOM      1  N  AGLY B   1       3.102   3.878   3.794  0.70  7.85           N
ATOM      2  CA AGLY B   1       3.985   4.960   4.314  0.70  6.79           C
ATOM      3  C  AGLY B   1       5.454   4.600   4.162  0.70  5.59           C
ATOM      4  O  AGLY B   1       5.756   3.431   3.920  0.70  6.04           O
ATOM      5  H1 AGLY B   1       2.660   4.173   3.080  0.70  7.85           H
ATOM      6  H2 AGLY B   1       3.597   3.175   3.566  0.70  7.85           H
ATOM      7  H3 AGLY B   1       2.522   3.640   4.425  0.70  7.85           H
ATOM      8  HA2AGLY B   1       3.815   5.782   3.828  0.70  6.79           H
ATOM      9  HA3AGLY B   1       3.798   5.111   5.254  0.70  6.79           H
ATOM     10  N  BGLY B   1       4.731   2.369   3.426  0.30  7.85           N
ATOM     11  CA BGLY B   1       5.847   3.227   3.917  0.30  6.79           C
ATOM     12  C  BGLY B   1       5.408   4.661   4.138  0.30  5.59           C
ATOM     13  O  BGLY B   1       4.212   4.956   4.163  0.30  6.04           O
ATOM     14  H1 BGLY B   1       4.655   1.653   3.949  0.30  7.85           H
ATOM     15  H2 BGLY B   1       3.971   2.832   3.444  0.30  7.85           H
ATOM     16  H3 BGLY B   1       4.901   2.110   2.592  0.30  7.85           H
ATOM     17  HA2BGLY B   1       6.182   2.874   4.756  0.30  6.79           H
ATOM     18  HA3BGLY B   1       6.569   3.222   3.269  0.30  6.79           H
ATOM     19  N   CYS B   2       6.380   5.555   4.293  1.00  5.95           N
ATOM     20  CA  CYS B   2       6.109   6.982   4.531  1.00  5.17           C
ATOM     21  C   CYS B   2       5.169   7.274   5.709  1.00  4.74           C
ATOM     22  O   CYS B   2       5.516   7.016   6.861  1.00  4.51           O
ATOM     23  CB  CYS B   2       5.631   7.673   3.244  1.00  5.99           C
ATOM     24  SG  CYS B   2       5.163   9.405   3.466  1.00  5.51           S
ATOM     25  H   CYS B   2       7.222   5.382   4.256  1.00  5.95           H
ATOM     26  HA  CYS B   2       6.956   7.396   4.758  1.00  5.17           H
ATOM     27  HB2 CYS B   2       6.346   7.642   2.590  1.00  5.99           H
ATOM     28  HB3 CYS B   2       4.856   7.199   2.905  1.00  5.99           H
ATOM     29  HG  CYS B   2       4.807   9.852   2.411  1.00  5.51           H
ATOM     30  N   CYS B   3       3.986   7.810   5.419  1.00  4.88           N
ATOM     31  CA  CYS B   3       3.018   8.129   6.461  1.00  5.24           C
ATOM     32  C   CYS B   3       2.431   6.857   7.062  1.00  6.13           C
ATOM     33  O   CYS B   3       1.976   6.851   8.205  1.00  7.85           O
ATOM     34  CB  CYS B   3       1.896   9.003   5.897  1.00  4.92           C
ATOM     35  SG  CYS B   3       2.450  10.600   5.255  1.00  6.13           S
ATOM     36  H   CYS B   3       3.720   7.999   4.623  1.00  4.88           H
ATOM     37  HA  CYS B   3       3.461   8.623   7.169  1.00  5.24           H
ATOM     38  HB2 CYS B   3       1.466   8.526   5.170  0.30  4.92           H
ATOM     39  HB3 CYS B   3       1.252   9.176   6.602  0.10  4.92           H
ATOM     40  HG  CYS B   3       3.237  10.413   4.368  1.00  6.13           H
TER
ATOM     41  O   HOH C   1       1.194   0.871   7.026  1.00  7.85           O
ATOM     42  H1  HOH C   1       2.011   1.046   7.112  1.00  7.85           H
ATOM     43  H2  HOH C   1       0.786   1.605   7.030  1.00  7.85           H
ATOM     44  O   HOH C   2       2.387   2.996   6.915  1.00  6.79           O
ATOM     45  H1  HOH C   2       3.061   3.497   6.915  1.00  6.79           H
ATOM     46  H2  HOH C   2       1.712   3.496   6.915  1.00  6.79           H
ATOM     47  O   HOH C   3       6.169   4.601   8.563  1.00  5.59           O
ATOM     48  D1  HOH C   3       6.992   4.768   8.578  1.00  5.59           D
ATOM     49  D2  HOH C   3       5.784   5.301   8.302  1.00  5.59           D
TER
END
"""

def exercise_h_counts():
  of = open("exercise_h_counts.pdb", "w")
  print >> of, pdb_file_exercise_h_counts
  of.close()
  import mmtbx.utils
  model = mmtbx.utils.model_simple(pdb_file_names=["exercise_h_counts.pdb"],
    scattering_table="n_gaussian")
  hc = model.h_counts()
  assert approx_equal(hc.h_count                , 26   , 0.01)
  assert approx_equal(hc.h_occ_sum              , 19.4 , 0.01)
  assert approx_equal(hc.h_fraction_of_total    , 50.52, 0.01)
  assert approx_equal(hc.hrot_count             , 8    , 0.01)
  assert approx_equal(hc.hrot_occ_sum           , 5    , 0.01)
  assert approx_equal(hc.hrot_fraction_of_total , 13.02, 0.01)
  model.show_h_counts()

def run():
  exercise()
  exercise_2()
  exercise_3()
  exercise_4()
  exercise_convert_atom()
  exercise_h_counts()
  print format_cpu_times()

if (__name__ == "__main__"):
  run()
