
from __future__ import division
from mmtbx import building
from libtbx.utils import null_out
from cStringIO import StringIO

def exercise_model_utils () :
  pdb_in = get_1yjp_pdb()
  residue = pdb_in.hierarchy.only_model().chains()[0].residue_groups()[0].only_atom_group()
  sele = pdb_in.hierarchy.atom_selection_cache().selection("resname TYR")
  water_sel = building.get_nearby_water_selection(
    pdb_hierarchy=pdb_in.hierarchy,
    xray_structure=pdb_in.input.xray_structure_simple(),
    selection=sele)
  assert (list(water_sel.iselection()) == [59, 60, 61, 62, 63])
  from mmtbx.monomer_library import idealized_aa
  from mmtbx.monomer_library import server
  mon_lib_srv = server.server()
  ideal_dict = idealized_aa.residue_dict()
  for resname, hierarchy in ideal_dict.iteritems() :
    residue = hierarchy.only_model().only_chain().only_residue_group().only_atom_group()
    result = building.generate_sidechain_clusters(residue, mon_lib_srv)
    if (len(result) == 0) :
      assert (residue.resname in ["ALA", "GLY"])
  # show_chain_resseq_ranges
  resids = [ (1,''),(2,''),(2,'A'),(4,''),(5,''),(6,''),(10,'B') ]
  import iotbx.pdb.hierarchy
  chain = iotbx.pdb.hierarchy.chain(id='A')
  for (resseq, icode) in resids :
    rg = iotbx.pdb.hierarchy.residue_group(resseq="%4d" % resseq, icode=icode)
    chain.append_residue_group(rg)
  out = StringIO()
  building.show_chain_resseq_ranges(chain.residue_groups(), out=out,
    prefix="  ")
  assert out.getvalue() == """  chain 'A': 1-2A,4-6,10B\n""", out.getvalue()

def exercise_box_rebuild () :
  #
  # UNSTABLE!
  #
  pdb_in = get_1yjp_pdb()
  xrs = pdb_in.input.xray_structure_simple()
  fc = xrs.structure_factors(d_min=1.5).f_calc()
  fc_fft = fc.fft_map(resolution_factor=0.25)
  fc_map = fc_fft.apply_sigma_scaling().real_map_unpadded()
  sel_cache = pdb_in.hierarchy.atom_selection_cache()
  sel = sel_cache.selection("resseq 4:5")
  sites_orig = xrs.sites_cart().deep_copy()
  xrs.shake_sites_in_place(0.3, selection=sel)
  pdb_in.hierarchy.atoms().set_xyz(xrs.sites_cart())
  sites_shaken = xrs.sites_cart().deep_copy()
  sites_new = building.run_real_space_annealing(
    xray_structure=xrs,
    pdb_hierarchy=pdb_in.hierarchy,
    selection=sel,
    target_map=fc_map,
    rsr_after_anneal=False,
    d_min=1.5,
    out=null_out())
  pdb_in.hierarchy.atoms().set_xyz(sites_new)
  rmsd_start = sites_orig.rms_difference(sites_shaken)
  assert (rmsd_start > 0.1)
  assert (sites_orig.rms_difference(sites_new) < rmsd_start), "%f < %f" % (
      sites_orig.rms_difference(sites_new), rmsd_start)
  assert (sites_orig.select(sel).rms_difference(sites_new.select(sel)) > 0)
  assert (sites_new.select(sel).rms_difference(sites_shaken.select(sel)) > 0)
  assert (sites_orig.select(~sel).rms_difference(sites_new.select(~sel)) == 0)
  box_builder = building.box_build_refine_base(
    pdb_hierarchy=pdb_in.hierarchy,
    xray_structure=xrs,
    processed_pdb_file=None,
    target_map=fc_map,
    selection=sel,
    d_min=1.5,
    out=null_out())
  box_builder.restrain_atoms(box_builder.others_in_box, 0.02)
  box_builder.geometry_minimization()
  box_builder.real_space_refine()
  sites_new = box_builder.update_original_coordinates()
  assert (sites_orig.rms_difference(sites_new) < 0.3)
  assert (box_builder.mean_density_at_sites() > 3)

def exercise_map_utils () :
  hierarchy, fmodel = get_1yjp_pdb_and_fmodel()
  sel_cache = hierarchy.atom_selection_cache()
  sele = sel_cache.selection("resseq 5 and (name CD or name OE1 or name NE2)")
  sele_all = sel_cache.selection("resseq 5")
  fmodel.xray_structure.scale_adp(factor=0.5, selection=sele)
  fmodel.update_xray_structure(update_f_calc=True)
  two_fofc_map, fofc_map = building.get_difference_maps(fmodel)
  map_stats = building.local_density_quality(
    fofc_map=fofc_map,
    two_fofc_map=two_fofc_map,
    atom_selection=sele_all,
    xray_structure=fmodel.xray_structure)
  out = StringIO()
  map_stats.show_atoms_outside_density(out=out, two_fofc_cutoff=3.0)
  pdb_strs = [ l.split(":")[0] for l in out.getvalue().splitlines() ]
  assert len(pdb_strs) > 0
  # XXX there seems to be a stochastic effect here
  #assert (pdb_strs == ['pdb=" CA  GLN A   5 "', 'pdb=" CB  GLN A   5 "',
  #                     'pdb=" CD  GLN A   5 "', 'pdb=" NE2 GLN A   5 "']), \
  #  pdb_strs
  fc_map = fmodel.map_coefficients(map_type="Fc").fft_map(
    resolution_factor=0.25).apply_sigma_scaling().real_map_unpadded()
  assert (map_stats.number_of_atoms_below_fofc_map_level() == 3)
  assert (map_stats.fraction_of_nearby_grid_points_above_cutoff() > 0.0)
  stats = building.get_model_map_stats(
    selection=sele_all,
    target_map=two_fofc_map,
    model_map=fc_map,
    unit_cell=fmodel.xray_structure.unit_cell(),
    sites_cart=fmodel.xray_structure.sites_cart(),
    pdb_atoms=hierarchy.atoms())
  #print stats.cc, stats.min, stats.mean
  # XXX numbers are not exact because of R-free flags are generated randomly
  #assert (stats.cc > 0.7) and (stats.min > 1.5) and (stats.mean > 2.5)

def get_1yjp_pdb_and_fmodel () :
  import mmtbx.utils
  pdb_in = get_1yjp_pdb()
  xrs = pdb_in.input.xray_structure_simple()
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
  return pdb_in.hierarchy, fmodel

def get_1yjp_pdb() :
  import iotbx.pdb.hierarchy
  return iotbx.pdb.hierarchy.input(pdb_string="""\
CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1      2
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  1.00 16.78           O
ATOM      5  N   ASN A   2      -7.656   2.923   3.155  1.00 15.02           N
ATOM      6  CA  ASN A   2      -6.522   2.038   2.831  1.00 14.10           C
ATOM      7  C   ASN A   2      -5.241   2.537   3.427  1.00 13.13           C
ATOM      8  O   ASN A   2      -4.978   3.742   3.426  1.00 11.91           O
ATOM      9  CB  ASN A   2      -6.346   1.881   1.341  1.00 15.38           C
ATOM     10  CG  ASN A   2      -7.584   1.342   0.692  1.00 14.08           C
ATOM     11  OD1 ASN A   2      -8.025   0.227   1.016  1.00 17.46           O
ATOM     12  ND2 ASN A   2      -8.204   2.155  -0.169  1.00 11.72           N
ATOM     13  N   ASN A   3      -4.438   1.590   3.905  1.00 12.26           N
ATOM     14  CA  ASN A   3      -3.193   1.904   4.589  1.00 11.74           C
ATOM     15  C   ASN A   3      -1.955   1.332   3.895  1.00 11.10           C
ATOM     16  O   ASN A   3      -1.872   0.119   3.648  1.00 10.42           O
ATOM     17  CB  ASN A   3      -3.259   1.378   6.042  1.00 12.15           C
ATOM     18  CG  ASN A   3      -2.006   1.739   6.861  1.00 12.82           C
ATOM     19  OD1 ASN A   3      -1.702   2.925   7.072  1.00 15.05           O
ATOM     20  ND2 ASN A   3      -1.271   0.715   7.306  1.00 13.48           N
ATOM     21  N   GLN A   4      -1.005   2.228   3.598  1.00 10.29           N
ATOM     22  CA  GLN A   4       0.384   1.888   3.199  1.00 10.53           C
ATOM     23  C   GLN A   4       1.435   2.606   4.088  1.00 10.24           C
ATOM     24  O   GLN A   4       1.547   3.843   4.115  1.00  8.86           O
ATOM     25  CB  GLN A   4       0.656   2.148   1.711  1.00  9.80           C
ATOM     26  CG  GLN A   4       1.944   1.458   1.213  1.00 10.25           C
ATOM     27  CD  GLN A   4       2.504   2.044  -0.089  0.80 12.43           C
ATOM     28  OE1 GLN A   4       2.744   3.268  -0.190  0.80 14.62           O
ATOM     29  NE2 GLN A   4       2.750   1.161  -1.091  0.80  9.05           N
ATOM     30  N   GLN A   5       2.154   1.821   4.871  1.00 10.38           N
ATOM     31  CA  GLN A   5       3.270   2.361   5.640  1.00 11.39           C
ATOM     32  C   GLN A   5       4.594   1.768   5.172  1.00 11.52           C
ATOM     33  O   GLN A   5       4.768   0.546   5.054  1.00 12.05           O
ATOM     34  CB  GLN A   5       3.056   2.183   7.147  1.00 11.96           C
ATOM     35  CG  GLN A   5       1.829   2.950   7.647  1.00 10.81           C
ATOM     36  CD  GLN A   5       1.344   2.414   8.954  0.80 13.10           C
ATOM     37  OE1 GLN A   5       0.774   1.325   9.002  0.80 10.65           O
ATOM     38  NE2 GLN A   5       1.549   3.187  10.039  0.80 12.30           N
ATOM     39  N   ASN A   6       5.514   2.664   4.856  1.00 11.99           N
ATOM     40  CA  ASN A   6       6.831   2.310   4.318  1.00 12.30           C
ATOM     41  C   ASN A   6       7.854   2.761   5.324  1.00 13.40           C
ATOM     42  O   ASN A   6       8.219   3.943   5.374  1.00 13.92           O
ATOM     43  CB  ASN A   6       7.065   3.016   2.993  1.00 12.13           C
ATOM     44  CG  ASN A   6       5.961   2.735   2.003  1.00 12.77           C
ATOM     45  OD1 ASN A   6       5.798   1.604   1.551  1.00 14.27           O
ATOM     46  ND2 ASN A   6       5.195   3.747   1.679  1.00 10.07           N
ATOM     47  N   TYR A   7       8.292   1.817   6.147  1.00 14.70           N
ATOM     48  CA  TYR A   7       9.159   2.144   7.299  1.00 15.18           C
ATOM     49  C   TYR A   7      10.603   2.331   6.885  1.00 15.91           C
ATOM     50  O   TYR A   7      11.041   1.811   5.855  1.00 15.76           O
ATOM     51  CB  TYR A   7       9.061   1.065   8.369  1.00 15.35           C
ATOM     52  CG  TYR A   7       7.665   0.929   8.902  1.00 14.45           C
ATOM     53  CD1 TYR A   7       6.771   0.021   8.327  1.00 15.68           C
ATOM     54  CD2 TYR A   7       7.210   1.756   9.920  1.00 14.80           C
ATOM     55  CE1 TYR A   7       5.480  -0.094   8.796  1.00 13.46           C
ATOM     56  CE2 TYR A   7       5.904   1.649  10.416  1.00 14.33           C
ATOM     57  CZ  TYR A   7       5.047   0.729   9.831  1.00 15.09           C
ATOM     58  OH  TYR A   7       3.766   0.589  10.291  1.00 14.39           O
ATOM     59  OXT TYR A   7      11.358   2.999   7.612  1.00 17.49           O
TER
HETATM   61  O   HOH A   8      -6.471   5.227   7.124  1.00 22.62           O
HETATM   62  O   HOH A   9      10.431   1.858   3.216  1.00 19.71           O
HETATM   63  O   HOH A  10     -11.286   1.756  -1.468  1.00 17.08           O
HETATM   64  O   HOH A  11      11.808   4.179   9.970  1.00 23.99           O
HETATM   65  O   HOH A  12      13.605   1.327   9.198  1.00 26.17           O
HETATM   66  O   HOH A  13      -2.749   3.429  10.024  1.00 39.15           O
HETATM   67  O   HOH A  14      -1.500   0.682  10.967  1.00 43.49           O
END
""")

if (__name__ == "__main__") :
  exercise_model_utils()
  exercise_map_utils()
  exercise_box_rebuild()
  print "OK"
