from __future__ import division

import iotbx.pdb
from cctbx.array_family import flex
from cctbx import adp_restraints # import dependency
import random
import mmtbx.command_line.real_space_refine as rs
from libtbx import group_args
from mmtbx.refinement import real_space
from mmtbx.rotamer.rotamer_eval import RotamerEval
from cctbx import miller
from cctbx import maptbx

if(1):
  random.seed(0)
  flex.set_random_seed(0)

pdb_str = """\
CRYST1   16.960   19.455   19.841  90.00  90.00  90.00 P 1
SCALE1      0.058962  0.000000  0.000000        0.00000
SCALE2      0.000000  0.051401  0.000000        0.00000
SCALE3      0.000000  0.000000  0.050401        0.00000
ATOM      1  N   ASP A  18      14.043  11.263  12.015  1.00100.00           N
ATOM      2  CA  ASP A  18      12.990  12.237  12.161  1.00100.00           C
ATOM      3  C   ASP A  18      12.453  12.618  10.780  1.00100.00           C
ATOM      4  O   ASP A  18      13.209  12.935   9.857  1.00100.00           O
ATOM      5  CB  ASP A  18      13.509  13.475  12.888  1.00100.00           C
ATOM      6  CG  ASP A  18      12.403  14.532  13.085  1.00100.00           C
ATOM      7  OD1 ASP A  18      11.437  14.212  13.836  1.00100.00           O
ATOM      8  OD2 ASP A  18      12.519  15.608  12.543  1.00100.00           O
ATOM      9  N   ASN A  19      11.135  12.599  10.637  1.00100.00           N
ATOM     10  CA  ASN A  19      10.462  12.955   9.415  1.00100.00           C
ATOM     11  C   ASN A  19      10.932  12.139   8.190  1.00100.00           C
ATOM     12  O   ASN A  19      11.745  12.583   7.411  1.00100.00           O
ATOM     13  CB  ASN A  19      10.574  14.466   9.101  1.00100.00           C
ATOM     14  CG  ASN A  19       9.667  14.910   8.017  1.00100.00           C
ATOM     15  OD1 ASN A  19       8.567  14.372   7.811  1.00100.00           O
ATOM     16  ND2 ASN A  19      10.066  15.976   7.273  1.00100.00           N
ATOM     17  N   TYR A  20      10.418  10.903   8.137  1.00100.00           N
ATOM     18  CA  TYR A  20      10.810   9.981   7.034  1.00100.00           C
ATOM     19  C   TYR A  20       9.703   8.988   6.760  1.00100.00           C
ATOM     20  O   TYR A  20       9.364   8.158   7.667  1.00100.00           O
ATOM     21  CB  TYR A  20      12.115   9.280   7.395  1.00100.00           C
ATOM     22  CG  TYR A  20      13.316  10.200   7.509  1.00100.00           C
ATOM     23  CD1 TYR A  20      14.072  10.516   6.392  1.00100.00           C
ATOM     24  CD2 TYR A  20      13.705  10.721   8.703  1.00100.00           C
ATOM     25  CE1 TYR A  20      15.189  11.316   6.470  1.00100.00           C
ATOM     26  CE2 TYR A  20      14.836  11.519   8.815  1.00100.00           C
ATOM     27  CZ  TYR A  20      15.542  11.836   7.693  1.00100.00           C
ATOM     28  OH  TYR A  20      16.673  12.662   7.779  1.00100.00           O
ATOM     29  N   ARG A  21       9.160   8.989   5.577  1.00100.00           N
ATOM     30  CA  ARG A  21       8.073   8.140   5.126  1.00100.00           C
ATOM     31  C   ARG A  21       6.863   8.252   6.048  1.00100.00           C
ATOM     32  O   ARG A  21       6.208   7.218   6.352  1.00100.00           O
ATOM     33  CB  ARG A  21       8.542   6.674   5.053  1.00100.00           C
ATOM     34  CG  ARG A  21       8.909   6.063   6.409  1.00100.00           C
ATOM     35  CD  ARG A  21       9.409   4.633   6.259  1.00100.00           C
ATOM     36  NE  ARG A  21       9.737   4.045   7.579  1.00100.00           N
ATOM     37  CZ  ARG A  21       8.851   3.418   8.335  1.00100.00           C
ATOM     38  NH1 ARG A  21       7.585   3.273   7.930  1.00100.00           N
ATOM     39  NH2 ARG A  21       9.220   2.935   9.504  1.00100.00           N
ATOM     40  N   GLY A  22       6.554   9.436   6.463  1.00100.00           N
ATOM     41  CA  GLY A  22       5.410   9.707   7.327  1.00100.00           C
ATOM     42  C   GLY A  22       5.748   9.622   8.821  1.00100.00           C
ATOM     43  O   GLY A  22       5.479  10.518   9.577  1.00100.00           O
ATOM     44  N   TYR A  23       6.330   8.487   9.201  1.00100.00           N
ATOM     45  CA  TYR A  23       6.695   8.251  10.596  1.00100.00           C
ATOM     46  C   TYR A  23       7.852   9.136  11.044  1.00100.00           C
ATOM     47  O   TYR A  23       8.641   9.570  10.203  1.00100.00           O
ATOM     48  CB  TYR A  23       7.052   6.786  10.831  1.00100.00           C
ATOM     49  CG  TYR A  23       7.531   6.486  12.226  1.00100.00           C
ATOM     50  CD1 TYR A  23       6.627   6.386  13.294  1.00100.00           C
ATOM     51  CD2 TYR A  23       8.870   6.238  12.480  1.00100.00           C
ATOM     52  CE1 TYR A  23       7.067   6.121  14.552  1.00100.00           C
ATOM     53  CE2 TYR A  23       9.341   5.961  13.770  1.00100.00           C
ATOM     54  CZ  TYR A  23       8.418   5.898  14.784  1.00100.00           C
ATOM     55  OH  TYR A  23       8.841   5.617  16.089  1.00100.00           O
ATOM     56  N   SER A  24       7.948   9.408  12.310  1.00100.00           N
ATOM     57  CA  SER A  24       9.013  10.264  12.855  1.00100.00           C
ATOM     58  C   SER A  24       9.607   9.636  14.109  1.00100.00           C
ATOM     59  O   SER A  24       8.916   9.382  15.096  1.00100.00           O
ATOM     60  CB  SER A  24       8.482  11.655  13.176  1.00100.00           C
ATOM     61  OG  SER A  24       8.072  12.331  11.990  1.00100.00           O
ATOM     62  N   LEU A  25      10.933   9.449  14.093  1.00100.00           N
ATOM     63  CA  LEU A  25      11.651   8.890  15.233  1.00100.00           C
ATOM     64  C   LEU A  25      12.390   9.961  16.015  1.00100.00           C
ATOM     65  O   LEU A  25      11.825  10.650  16.815  1.00100.00           O
ATOM     66  CB  LEU A  25      12.616   7.785  14.748  1.00100.00           C
ATOM     67  CG  LEU A  25      12.075   6.386  14.572  1.00100.00           C
ATOM     68  CD1 LEU A  25      11.031   6.304  13.433  1.00100.00           C
ATOM     69  CD2 LEU A  25      13.212   5.385  14.266  1.00100.00           C
TER      70      LEU A  25
END
"""

def get_pdb_inputs(pdb_str):
  raw_records = flex.std_string(pdb_str.splitlines())
  processed_pdb_file = rs.get_processed_pdb_object(raw_records=raw_records,
    rama_potential=None, log = None)
  xrs = processed_pdb_file.xray_structure(show_summary = False)
  geometry_restraints_manager = rs.get_geometry_restraints_manager(
    processed_pdb_file = processed_pdb_file,
    xray_structure     = xrs)
  pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  return group_args(
    ph  = pdb_hierarchy,
    grm = geometry_restraints_manager,
    xrs = xrs)

def get_tmo(inp, d_min):
  sel = inp.ph.atom_selection_cache().selection(
    string = "name CA or name CB or name N or name O or name C")
  f_calc = inp.xrs.select(
    selection=sel).structure_factors(d_min=d_min).f_calc()
  fft_map = f_calc.fft_map(resolution_factor = 0.25)
  fft_map.apply_sigma_scaling()
  target_map = fft_map.real_map_unpadded()
  return group_args(
    data             = target_map,
    miller_array     = f_calc,
    crystal_gridding = fft_map)

def shake_sites(xrs, random, shift, grm=None):
  from mmtbx.dynamics import cartesian_dynamics
  if(random):
    xrs.shake_sites_in_place(mean_distance = shift)
  else:
    grad_calc = cartesian_dynamics.gradients_calculator_geometry_restraints(
      restraints_manager = grm)
    cartesian_dynamics.run(
      xray_structure       = xrs,
      gradients_calculator = grad_calc,
      temperature          = 1000,
      n_steps              = 100000,
      time_step            = 0.0005,
      stop_cm_motion       = True,
      stop_at_diff         = shift)
  return xrs

def compute_map(target_map, xray_structure):
  mc = target_map.miller_array.structure_factors_from_scatterers(
    xray_structure = xray_structure).f_calc()
  fft_map = miller.fft_map(
    crystal_gridding     = target_map.crystal_gridding,
    fourier_coefficients = mc)
  fft_map.apply_sigma_scaling()
  return fft_map.real_map_unpadded()

def show(pdb_hierarchy, tm, xrs, grm, prefix):
  map = compute_map(target_map=tm, xray_structure=xrs)
  cc = flex.linear_correlation(
    x=map.as_1d(),
    y=tm.data.as_1d()).coefficient()
  es = grm.energies_sites(sites_cart = xrs.sites_cart())
  rmsd_a = es.angle_deviations()[2]
  rmsd_b = es.bond_deviations()[2]
  print "%s: overall CC: %6.4f rmsd_bonds=%6.3f rmsd_angles=%6.3f"%(
    prefix, cc, rmsd_b, rmsd_a)
  pdb_hierarchy.adopt_xray_structure(xrs)
  rotamer_manager = RotamerEval()
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for residue in chain.residues():
        sites_cart = residue.atoms().extract_xyz()
        sel = maptbx.grid_indices_around_sites(
          unit_cell  = xrs.unit_cell(),
          fft_n_real = map.focus(),
          fft_m_real = map.all(),
          sites_cart = sites_cart,
          site_radii = flex.double(sites_cart.size(), 2))
        ccr = flex.linear_correlation(
          x=map.select(sel).as_1d(),
          y=tm.data.select(sel).as_1d()).coefficient()
        fmt = "%s: %4s %10s CC: %6.4f"
        print fmt%(prefix, residue.resname, rotamer_manager.evaluate_residue(residue),ccr)

def exercise(d_min=5, random_seed=1111111):
  inp = get_pdb_inputs(pdb_str=pdb_str)
  xrs_good = inp.xrs.deep_copy_scatterers()
  target_map = get_tmo(inp=inp, d_min = d_min)
  if(1):
    inp.ph.write_pdb_file(file_name="strat.pdb")
  xrs_poor = shake_sites(xrs=xrs_good.deep_copy_scatterers(), random=False,
    shift=2, grm=inp.grm)
  if(1):
    inp.ph.adopt_xray_structure(xrs_poor)
    inp.ph.write_pdb_file(file_name="poor.pdb")
  #
  show(prefix="GOOD",pdb_hierarchy = inp.ph, tm=target_map, xrs=xrs_good, grm=inp.grm.geometry)
  #
  for use_reference_torsion in [False, True]:
    inp.ph.adopt_xray_structure(xrs_good)
    random.seed(random_seed)
    flex.set_random_seed(random_seed)
    print "*"*79
    print "use_reference_torsion:", use_reference_torsion
    print "*"*79
    show(prefix="START",pdb_hierarchy = inp.ph, tm=target_map, xrs=xrs_poor, grm=inp.grm.geometry)
    #
    if(use_reference_torsion):
      for model in inp.ph.models():
        for chain in model.chains():
          for residue in chain.residues():
            inp.grm.geometry.generic_restraints_manager.reference_manager.add_torsion_restraints(
              pdb_hierarchy = inp.ph,
              sites_cart    = residue.atoms().extract_xyz(),
              selection     = residue.atoms().extract_i_seq(),
              sigma         = 1.0)
    #
    tmp = xrs_poor.deep_copy_scatterers()
    rsr_simple_refiner = real_space.simple(
      target_map                  = target_map.data,
      selection                   = flex.bool(tmp.scatterers().size(), True),
      real_space_gradients_delta  = d_min/4,
      max_iterations              = 500,
      geometry_restraints_manager = inp.grm.geometry)
    refined = real_space.refinery(
      refiner          = rsr_simple_refiner,
      optimize_weight  = True,
      xray_structure   = tmp,
      start_trial_weight_value = 50,
      rms_bonds_limit  = 0.02,
      rms_angles_limit = 2.0)
    assert refined.sites_cart_result is not None
    tmp = tmp.replace_sites_cart(refined.sites_cart_result)
    show(prefix="FINAL",pdb_hierarchy = inp.ph, tm=target_map, xrs=tmp, grm=inp.grm.geometry)
    if(1):
      inp.ph.adopt_xray_structure(tmp)
      inp.ph.write_pdb_file(file_name="final_%s.pdb"%str(use_reference_torsion))


if (__name__ == "__main__") :
  exercise()
  print "OK"
