from __future__ import absolute_import, division, print_function
import time
import mmtbx.refinement.real_space
import iotbx.pdb
from libtbx.utils import null_out
from scitbx.matrix import rotate_point_around_axis
from mmtbx.rotamer.rotamer_eval import RotamerEval
from mmtbx.idealized_aa_residues import sample_rotamers
from cctbx import maptbx
from scitbx.array_family import flex

pdb_good = """\
CRYST1   21.207   17.255   22.844  90.00  90.00  90.00 P 1
ATOM      1  N   TRP L 148      14.260  10.871  17.844  1.00 22.69           N
ATOM      2  CA  TRP L 148      14.094  10.389  16.478  1.00 19.52           C
ATOM      3  C   TRP L 148      14.872   9.095  16.265  1.00 21.40           C
ATOM      4  O   TRP L 148      15.903   8.867  16.898  1.00 21.95           O
ATOM      5  CB  TRP L 148      14.553  11.451  15.476  1.00 18.93           C
ATOM      6  CG  TRP L 148      14.396  11.033  14.046  1.00 16.05           C
ATOM      7  CD1 TRP L 148      15.328  10.415  13.265  1.00 15.60           C
ATOM      8  CD2 TRP L 148      13.234  11.204  13.225  1.00 15.47           C
ATOM      9  NE1 TRP L 148      14.819  10.190  12.008  1.00 17.97           N
ATOM     10  CE2 TRP L 148      13.535  10.665  11.958  1.00 16.42           C
ATOM     11  CE3 TRP L 148      11.968  11.759  13.437  1.00 16.97           C
ATOM     12  CZ2 TRP L 148      12.618  10.666  10.909  1.00 15.78           C
ATOM     13  CZ3 TRP L 148      11.059  11.759  12.395  1.00 17.54           C
ATOM     14  CH2 TRP L 148      11.389  11.216  11.147  1.00 17.61           C
ATOM     15  HA  TRP L 148      13.155  10.207  16.317  1.00 19.52           H
ATOM     16  HB2 TRP L 148      14.028  12.255  15.610  1.00 18.93           H
ATOM     17  HB3 TRP L 148      15.492  11.640  15.629  1.00 18.93           H
ATOM     18  HD1 TRP L 148      16.184  10.180  13.542  1.00 15.60           H
ATOM     19  HE1 TRP L 148      15.239   9.813  11.359  1.00 17.97           H
ATOM     20  HE3 TRP L 148      11.742  12.122  14.264  1.00 16.97           H
ATOM     21  HZ2 TRP L 148      12.834  10.306  10.079  1.00 15.78           H
ATOM     22  HZ3 TRP L 148      10.215  12.126  12.525  1.00 17.54           H
ATOM     23  HH2 TRP L 148      10.757  11.230  10.464  1.00 17.61           H
ATOM     24  N   TYR L 192      11.041   5.312   9.916  1.00 20.96           N
ATOM     25  CA  TYR L 192      11.151   6.712  10.309  1.00 18.77           C
ATOM     26  C   TYR L 192       9.991   7.112  11.214  1.00 20.42           C
ATOM     27  O   TYR L 192       8.825   6.971  10.844  1.00 20.12           O
ATOM     28  CB  TYR L 192      11.189   7.614   9.074  1.00 15.79           C
ATOM     29  CG  TYR L 192      12.351   7.332   8.148  1.00 19.38           C
ATOM     30  CD1 TYR L 192      12.232   6.414   7.112  1.00 18.58           C
ATOM     31  CD2 TYR L 192      13.566   7.984   8.308  1.00 13.17           C
ATOM     32  CE1 TYR L 192      13.291   6.154   6.263  1.00 17.30           C
ATOM     33  CE2 TYR L 192      14.631   7.730   7.464  1.00 19.20           C
ATOM     34  CZ  TYR L 192      14.487   6.814   6.444  1.00 19.28           C
ATOM     35  OH  TYR L 192      15.544   6.558   5.601  1.00 22.20           O
ATOM     36  HA  TYR L 192      11.976   6.840  10.802  1.00 18.77           H
ATOM     37  HB2 TYR L 192      10.370   7.488   8.569  1.00 15.79           H
ATOM     38  HB3 TYR L 192      11.258   8.537   9.363  1.00 15.79           H
ATOM     39  HD1 TYR L 192      11.426   5.967   6.988  1.00 18.58           H
ATOM     40  HD2 TYR L 192      13.666   8.602   8.996  1.00 13.17           H
ATOM     41  HE1 TYR L 192      13.196   5.536   5.574  1.00 17.30           H
ATOM     42  HE2 TYR L 192      15.439   8.174   7.584  1.00 19.20           H
ATOM     43  HH  TYR L 192      16.207   7.024   5.821  1.00 22.20           H
ATOM     44  N   PHE L 209       6.405   6.669   9.777  1.00 21.36           N
ATOM     45  CA  PHE L 209       5.861   6.089   8.555  1.00 22.67           C
ATOM     46  C   PHE L 209       6.783   5.000   8.016  1.00 26.96           C
ATOM     47  O   PHE L 209       7.997   5.049   8.212  1.00 22.94           O
ATOM     48  CB  PHE L 209       5.650   7.172   7.495  1.00 24.63           C
ATOM     49  CG  PHE L 209       6.906   7.909   7.122  1.00 23.08           C
ATOM     50  CD1 PHE L 209       7.686   7.482   6.060  1.00 22.08           C
ATOM     51  CD2 PHE L 209       7.306   9.027   7.834  1.00 22.40           C
ATOM     52  CE1 PHE L 209       8.841   8.158   5.715  1.00 26.03           C
ATOM     53  CE2 PHE L 209       8.460   9.707   7.493  1.00 21.71           C
ATOM     54  CZ  PHE L 209       9.229   9.271   6.433  1.00 20.59           C
ATOM     55  HA  PHE L 209       5.000   5.686   8.750  1.00 22.67           H
ATOM     56  HB2 PHE L 209       5.299   6.759   6.691  1.00 24.63           H
ATOM     57  HB3 PHE L 209       5.015   7.822   7.835  1.00 24.63           H
ATOM     58  HD1 PHE L 209       7.430   6.732   5.574  1.00 22.08           H
ATOM     59  HD2 PHE L 209       6.792   9.325   8.549  1.00 22.40           H
ATOM     60  HE1 PHE L 209       9.357   7.862   5.000  1.00 26.03           H
ATOM     61  HE2 PHE L 209       8.719  10.457   7.979  1.00 21.71           H
ATOM     62  HZ  PHE L 209      10.006   9.727   6.202  1.00 20.59           H
TER
END
"""

pdb_poor = """\
CRYST1   21.207   17.255   22.844  90.00  90.00  90.00 P 1
ATOM      1  N   TRP L 148      14.260  10.871  17.844  1.00 22.69           N
ATOM      2  CA  TRP L 148      14.094  10.389  16.478  1.00 19.52           C
ATOM      3  C   TRP L 148      14.872   9.095  16.265  1.00 21.40           C
ATOM      4  O   TRP L 148      15.903   8.867  16.898  1.00 21.95           O
ATOM      5  CB  TRP L 148      14.553  11.451  15.476  1.00 18.93           C
ATOM      6  CG  TRP L 148      14.396  11.033  14.046  1.00 16.05           C
ATOM      7  CD1 TRP L 148      15.328  10.415  13.265  1.00 15.60           C
ATOM      8  CD2 TRP L 148      13.234  11.204  13.225  1.00 15.47           C
ATOM      9  NE1 TRP L 148      14.819  10.190  12.008  1.00 17.97           N
ATOM     10  CE2 TRP L 148      13.535  10.665  11.958  1.00 16.42           C
ATOM     11  CE3 TRP L 148      11.968  11.759  13.437  1.00 16.97           C
ATOM     12  CZ2 TRP L 148      12.618  10.666  10.909  1.00 15.78           C
ATOM     13  CZ3 TRP L 148      11.059  11.759  12.395  1.00 17.54           C
ATOM     14  CH2 TRP L 148      11.389  11.216  11.147  1.00 17.61           C
ATOM     15  HA  TRP L 148      13.155  10.207  16.317  1.00 19.52           H
ATOM     16  HB2 TRP L 148      14.028  12.255  15.610  1.00 18.93           H
ATOM     17  HB3 TRP L 148      15.492  11.640  15.629  1.00 18.93           H
ATOM     18  HD1 TRP L 148      16.184  10.180  13.542  1.00 15.60           H
ATOM     19  HE1 TRP L 148      15.239   9.813  11.359  1.00 17.97           H
ATOM     20  HE3 TRP L 148      11.742  12.122  14.264  1.00 16.97           H
ATOM     21  HZ2 TRP L 148      12.834  10.306  10.079  1.00 15.78           H
ATOM     22  HZ3 TRP L 148      10.215  12.126  12.525  1.00 17.54           H
ATOM     23  HH2 TRP L 148      10.757  11.230  10.464  1.00 17.61           H
ATOM     24  N   TYR L 192      11.041   5.312   9.916  1.00 20.96           N
ATOM     25  CA  TYR L 192      11.151   6.712  10.309  1.00 18.77           C
ATOM     26  C   TYR L 192       9.991   7.112  11.214  1.00 20.42           C
ATOM     27  O   TYR L 192       8.825   6.971  10.844  1.00 20.12           O
ATOM     28  CB  TYR L 192      11.189   7.614   9.074  1.00 15.79           C
ATOM     29  CG  TYR L 192      12.585   7.885   8.559  1.00 19.38           C
ATOM     30  CD1 TYR L 192      12.786   8.551   7.356  1.00 18.58           C
ATOM     31  CD2 TYR L 192      13.702   7.478   9.276  1.00 13.17           C
ATOM     32  CE1 TYR L 192      14.060   8.801   6.883  1.00 17.30           C
ATOM     33  CE2 TYR L 192      14.980   7.724   8.810  1.00 19.20           C
ATOM     34  CZ  TYR L 192      15.152   8.386   7.613  1.00 19.28           C
ATOM     35  OH  TYR L 192      16.422   8.633   7.145  1.00 22.20           O
ATOM     36  HA  TYR L 192      11.976   6.840  10.802  1.00 18.77           H
ATOM     37  HB2 TYR L 192      10.687   7.189   8.360  1.00 15.79           H
ATOM     38  HB3 TYR L 192      10.785   8.467   9.297  1.00 15.79           H
ATOM     39  HD1 TYR L 192      12.050   8.832   6.862  1.00 18.58           H
ATOM     40  HD2 TYR L 192      13.589   7.030  10.083  1.00 13.17           H
ATOM     41  HE1 TYR L 192      14.178   9.248   6.075  1.00 17.30           H
ATOM     42  HE2 TYR L 192      15.718   7.444   9.300  1.00 19.20           H
ATOM     43  HH  TYR L 192      16.991   8.329   7.682  1.00 22.20           H
ATOM     44  N   PHE L 209       6.405   6.669   9.777  1.00 21.36           N
ATOM     45  CA  PHE L 209       5.861   6.089   8.555  1.00 22.67           C
ATOM     46  C   PHE L 209       6.783   5.000   8.016  1.00 26.96           C
ATOM     47  O   PHE L 209       7.997   5.049   8.212  1.00 22.94           O
ATOM     48  CB  PHE L 209       5.650   7.172   7.495  1.00 24.63           C
ATOM     49  CG  PHE L 209       6.906   7.909   7.122  1.00 23.08           C
ATOM     50  CD1 PHE L 209       7.686   7.482   6.060  1.00 22.08           C
ATOM     51  CD2 PHE L 209       7.306   9.027   7.834  1.00 22.40           C
ATOM     52  CE1 PHE L 209       8.841   8.158   5.715  1.00 26.03           C
ATOM     53  CE2 PHE L 209       8.460   9.707   7.493  1.00 21.71           C
ATOM     54  CZ  PHE L 209       9.229   9.271   6.433  1.00 20.59           C
ATOM     55  HA  PHE L 209       5.000   5.686   8.750  1.00 22.67           H
ATOM     56  HB2 PHE L 209       5.299   6.759   6.691  1.00 24.63           H
ATOM     57  HB3 PHE L 209       5.015   7.822   7.835  1.00 24.63           H
ATOM     58  HD1 PHE L 209       7.430   6.732   5.574  1.00 22.08           H
ATOM     59  HD2 PHE L 209       6.792   9.325   8.549  1.00 22.40           H
ATOM     60  HE1 PHE L 209       9.357   7.862   5.000  1.00 26.03           H
ATOM     61  HE2 PHE L 209       8.719  10.457   7.979  1.00 21.71           H
ATOM     62  HZ  PHE L 209      10.006   9.727   6.202  1.00 20.59           H
TER      63      PHE L 209
END
"""

def main(use_mask):
  #
  # SETUP INPUTS: read good model, compute map from it
  #
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_good)
  model_good = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  mon_lib_srv = model_good.get_mon_lib_srv()
  with open("model_good.pdb", "w") as fo:
    fo.write(model_good.model_as_pdb())
  xrs_good = model_good.get_xray_structure()
  f_calc = xrs_good.structure_factors(d_min = 2).f_calc()
  crystal_gridding = maptbx.crystal_gridding(
    unit_cell        = xrs_good.unit_cell(),
    space_group_info = xrs_good.space_group_info(),
    symmetry_flags   = maptbx.use_space_group_symmetry,
    step             = 0.5)
  fft_map = f_calc.fft_map(crystal_gridding=crystal_gridding)
  fft_map.apply_sigma_scaling()
  target_map = fft_map.real_map_unpadded()
  #
  # Read poor model (same as good model except one residue side chain)
  #
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_poor)
  model_poor = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  with open("model_poor.pdb", "w") as fo:
    fo.write(model_poor.model_as_pdb())
  h_poor = model_poor.get_hierarchy()
  sites_cart_all = model_poor.get_sites_cart()
  uc = model_poor.crystal_symmetry().unit_cell()
  #
  # Track states
  #
  states = model_poor.get_states_collector()
  #
  # This is the residue (number 192) we will be working on (sampling rotamers)
  #
  for m in h_poor.models():
    for c in m.chains():
      for r in c.residues():
        if r.resseq_as_int()==192:
          residue = r
          selection = r.atoms().extract_i_seq()
          break
  #
  # Mask neighjbors. If any trial rotamer atom hits mask value 0 it will be
  # rejected.
  #
  if(use_mask):
    sel_to_mask = ~flex.bool(model_poor.size(), selection)
    mask = maptbx.mask(
      xray_structure = model_poor.get_xray_structure().select(sel_to_mask),
      n_real         = crystal_gridding.n_real(),
      mask_value_inside_molecule = 0,
      mask_value_outside_molecule = 1,
      solvent_radius = 0,
      atom_radius = None)
  #
  # Sample all rotamers
  #
  clusters = mmtbx.refinement.real_space.aa_residue_axes_and_clusters(
    residue         = residue,
    mon_lib_srv     = mon_lib_srv,
    backbone_sample = False).clusters
  rotamer_eval = RotamerEval()
  nested_loop = sample_rotamers.get_nested_loop(
    n=len(clusters), fine=False, start=0, end=360)
  sites_cart = residue.atoms().extract_xyz()
  score_best = -1.e9
  sites_cart_best = None
  for angles in nested_loop:
    sites_cart_moved = sites_cart.deep_copy()
    score=0
    for i, angle in enumerate(angles):
      cl = clusters[i]
      for atom_to_rotate in cl.atoms_to_rotate:
        new_site_cart = rotate_point_around_axis(
          axis_point_1 = sites_cart_moved[cl.axis[0]],
          axis_point_2 = sites_cart_moved[cl.axis[1]],
          point        = sites_cart_moved[atom_to_rotate],
          angle        = angle,
          deg          = True)
        sites_cart_moved[atom_to_rotate] = new_site_cart
    residue.atoms().set_xyz(sites_cart_moved)
    sites_frac_moved = uc.fractionalize(sites_cart_moved)
    # Check if trial rotamer bumps into neighor atoms
    hit_neighbour = False
    if(use_mask):
      for site_frac in sites_frac_moved:
        mv = mask.value_at_closest_grid_point(site_frac)
        if abs(mv)<1.e-6: hit_neighbour = True
        score += target_map.eight_point_interpolation(site_frac)
    # Evaluate rotemeric state
    fl = str(rotamer_eval.evaluate_residue_2(residue = residue)).strip().upper()
    # Accept rotamer if it is not OUTLIER and does not clash with neighbours
    if(fl in ["ALLOWED", "FAVORED"] and not hit_neighbour):
      tmp = sites_cart_all.deep_copy()
      tmp = tmp.set_selected(selection, sites_cart_moved)
      states.add(sites_cart=tmp)
      if(score > score_best):
        score_best = score
        sites_cart_best = sites_cart_moved.deep_copy()
  #
  # Write multi-model file that contains all accepted rotamers
  #
  states.write(file_name="all_states_UseMask%s.pdb"%str(use_mask),
    crystal_symmetry=model_poor.crystal_symmetry())
  #
  # Set best fitting rotamer and write final model
  #
  residue.atoms().set_xyz(sites_cart_best)
  with open("model_fitted.pdb", "w") as fo:
    fo.write(model_poor.model_as_pdb())
  #
  # Assert fitted model matches the answer
  #
  if use_mask:
    s1 = model_good.get_sites_cart()
    s2 = model_poor.get_sites_cart()
    max_dist = flex.max(flex.sqrt((s1 - s2).dot()))
    assert max_dist < 0.1

if(__name__ == "__main__"):
  t0 = time.time()
  for use_mask in [True, False]:
    main(use_mask = use_mask)
  print(time.time()-t0)
