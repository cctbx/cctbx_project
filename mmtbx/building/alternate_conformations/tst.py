
from __future__ import division
from mmtbx.building import alternate_conformations
import mmtbx.utils
import iotbx.pdb.hierarchy
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out

def exercise_utils () :
  #--- coord_stats_with_flips
  pdb_1 = iotbx.pdb.hierarchy.input(pdb_string="""\
ATOM   1639  N   PHE A 113      18.514  41.994  54.886  1.00 12.36           N
ATOM   1640  CA  PHE A 113      19.737  42.742  54.898  1.00 12.97           C
ATOM   1641  C   PHE A 113      19.514  44.206  55.111  1.00 11.88           C
ATOM   1642  O   PHE A 113      18.442  44.758  54.823  1.00 13.75           O
ATOM   1643  CB  PHE A 113      20.587  42.460  53.684  1.00 20.08           C
ATOM   1644  CG  PHE A 113      19.893  42.692  52.399  1.00 17.15           C
ATOM   1645  CD1 PHE A 113      19.182  41.689  51.798  1.00 20.34           C
ATOM   1646  CD2 PHE A 113      20.029  43.895  51.745  1.00 21.22           C
ATOM   1647  CE1 PHE A 113      18.577  41.887  50.601  1.00 23.49           C
ATOM   1648  CE2 PHE A 113      19.411  44.105  50.542  1.00 24.87           C
ATOM   1649  CZ  PHE A 113      18.675  43.097  49.976  1.00 21.81           C
ATOM   1650  H   PHE A 113      17.875  42.342  54.427  1.00 14.83           H
ATOM   1651  HA  PHE A 113      20.250  42.435  55.661  1.00 15.57           H
ATOM   1652  HB2 PHE A 113      21.367  43.036  53.708  1.00 24.10           H
ATOM   1653  HB3 PHE A 113      20.865  41.531  53.707  1.00 24.10           H
ATOM   1654  HD1 PHE A 113      19.099  40.866  52.224  1.00 24.40           H
ATOM   1655  HD2 PHE A 113      20.518  44.581  52.138  1.00 25.47           H
ATOM   1656  HE1 PHE A 113      18.077  41.204  50.215  1.00 28.18           H
ATOM   1657  HE2 PHE A 113      19.488  44.927  50.113  1.00 29.84           H
ATOM   1658  HZ  PHE A 113      18.261  43.230  49.154  1.00 26.17           H
""")
  pdb_2 = iotbx.pdb.hierarchy.input(pdb_string="""\
ATOM   1639  N   PHE A 113      18.514  41.994  54.886  1.00 12.36           N
ATOM   1640  CA  PHE A 113      19.737  42.742  54.898  1.00 12.97           C
ATOM   1641  C   PHE A 113      19.514  44.206  55.111  1.00 11.88           C
ATOM   1642  O   PHE A 113      18.442  44.758  54.823  1.00 13.75           O
ATOM   1643  CB  PHE A 113      20.587  42.460  53.684  1.00 20.08           C
ATOM   1644  CG  PHE A 113      19.893  42.692  52.399  1.00 17.15           C
ATOM   1645  CD1 PHE A 113      20.174  43.797  51.643  1.00 20.34           C
ATOM   1646  CD2 PHE A 113      18.888  41.845  51.991  1.00 21.22           C
ATOM   1647  CE1 PHE A 113      19.512  44.033  50.483  1.00 23.49           C
ATOM   1648  CE2 PHE A 113      18.224  42.071  50.816  1.00 24.87           C
ATOM   1649  CZ  PHE A 113      18.548  43.166  50.057  1.00 21.81           C
ATOM   1650  H   PHE A 113      17.875  42.342  54.427  1.00 14.83           H
ATOM   1651  HA  PHE A 113      20.250  42.435  55.661  1.00 15.57           H
ATOM   1652  HB2 PHE A 113      21.367  43.036  53.708  1.00 24.10           H
ATOM   1653  HB3 PHE A 113      20.865  41.531  53.707  1.00 24.10           H
ATOM   1654  HD1 PHE A 113      20.840  44.386  51.919  1.00 24.40           H
ATOM   1655  HD2 PHE A 113      18.681  41.096  52.501  1.00 25.47           H
ATOM   1656  HE1 PHE A 113      19.730  44.776  49.967  1.00 28.18           H
ATOM   1657  HE2 PHE A 113      17.560  41.484  50.533  1.00 29.84           H
ATOM   1658  HZ  PHE A 113      18.093  43.330  49.263  1.00 26.17           H
""")
  hierarchy_1 = pdb_1.hierarchy
  hierarchy_2 = pdb_2.hierarchy
  pdb_atoms_1 = hierarchy_1.atoms()
  sites_1 = pdb_atoms_1.extract_xyz()
  sites_2 = hierarchy_2.atoms().extract_xyz()
  result = alternate_conformations.coord_stats_with_flips(sites_1, sites_2, pdb_atoms_1)
  assert (approx_equal(result.rmsd, 0.2, eps=0.001))
  assert (approx_equal(result.max_dev, 0.452427, eps=0.00001))
  #--- get_selection_gap
  assert (alternate_conformations.get_selection_gap([1,2,3,4,5],[6,7,8,9]) == 0)
  assert (alternate_conformations.get_selection_gap([6,7], [1,2,3,4,5]) == 0)
  assert (alternate_conformations.get_selection_gap([1,2,3,4,5],[7,8,9]) == 1)
  assert (alternate_conformations.get_selection_gap([1,2,3,4,5],[4,5,7,8,9]) == -2)
  assert (alternate_conformations.get_selection_gap([4,5,7,8,9],[1,2,3,4]) == -1)
  #--- score_rotamers
  n_outliers = alternate_conformations.score_rotamers(hierarchy_2,
    flex.bool(sites_2.size(), True))
  assert (n_outliers == 0)
  #--- get_partial_omit_map
  xrs = pdb_1.input.xray_structure_simple()
  f_calc = abs(xrs.structure_factors(d_min=1.5).f_calc())
  f_calc.set_observation_type_xray_amplitude()
  flags = f_calc.generate_r_free_flags()
  fmodel = mmtbx.utils.fmodel_simple(
    xray_structures=[xrs],
    f_obs=f_calc,
    r_free_flags=flags,
    scattering_table="n_gaussian",
    bulk_solvent_and_scaling=False,
    skip_twin_detection=True)
  sel = pdb_1.hierarchy.atom_selection_cache().selection("name CZ")
  assert (sel.count(True) == 1)
  two_fofc_map, fofc_map = alternate_conformations.get_partial_omit_map(
    fmodel=fmodel,
    selection=sel)
  site = xrs.sites_frac()[sel.iselection()[0]] # CZ coordinate
  assert (fofc_map.tricubic_interpolation(site) > 20)

def exercise_rejoin () :
  from mmtbx.building.alternate_conformations.tst_build_simple import pdb_raw
  pdb_in = iotbx.pdb.hierarchy.input(pdb_string=pdb_raw)
  hierarchy = pdb_in.hierarchy
  params = alternate_conformations.rejoin_phil.extract()
  n_modified = alternate_conformations.rejoin_split_single_conformers(
    pdb_hierarchy=hierarchy,
    params=params,
    model_error_ml=0.5,
    log=null_out())
  assert (n_modified == 3), n_modified # Gln5
  # split residue 6 without changing coordinates, set occupancy very low
  chain = hierarchy.only_model().chains()[0]
  rg6 = chain.residue_groups()[5]
  ag = rg6.only_atom_group().detached_copy()
  for atom in ag.atoms() :
    atom.occ = 0.05
  for atom in rg6.atoms() :
    atom.occ = 0.95
  rg6.only_atom_group().altloc = 'A'
  ag.altloc = 'B'
  rg6.append_atom_group(ag)
  n_modified = alternate_conformations.rejoin_split_single_conformers(
    pdb_hierarchy=hierarchy.deep_copy(),
    params=params,
    model_error_ml=0.5,
    log=null_out())
  assert (n_modified == 1), n_modified
  # now with higher B-factors for all atoms
  for atom in hierarchy.atoms() :
    atom.b = atom.b * 10
  n_modified = alternate_conformations.rejoin_split_single_conformers(
    pdb_hierarchy=hierarchy,
    params=params,
    log=null_out())
  assert (n_modified == 1), n_modified
  # TODO more needed...

if (__name__ == "__main__") :
  exercise_utils()
  exercise_rejoin()
  print "OK"
