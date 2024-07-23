from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
import iotbx.pdb
from libtbx import group_args
from libtbx.utils import user_plus_sys_time
import mmtbx.f_model
from mmtbx.dynamics import simulated_annealing as sa
import random
from mmtbx.dynamics import cartesian_dynamics
from mmtbx.refinement import real_space
import mmtbx.utils
import mmtbx.refinement.real_space.individual_sites
from libtbx.utils import null_out
import mmtbx.model

pdb_str_1 = """\
CRYST1   26.960   29.455   29.841  90.00  90.00  90.00 P 21 21 21
ATOM      1  N   ASP A  18      14.467  12.258  12.105  1.00 21.41           N
ATOM      2  CA  ASP A  18      13.225  13.003  12.273  1.00 34.70           C
ATOM      3  C   ASP A  18      12.437  13.061  10.968  1.00 35.50           C
ATOM      4  O   ASP A  18      12.962  13.483   9.937  1.00 38.48           O
ATOM      5  CB  ASP A  18      13.514  14.418  12.777  1.00 38.15           C
ATOM      6  CG  ASP A  18      12.251  15.231  12.986  1.00 27.39           C
ATOM      7  OD1 ASP A  18      11.675  15.164  14.092  1.00 35.78           O
ATOM      8  OD2 ASP A  18      11.836  15.939  12.044  1.00 30.67           O
ATOM      9  N   ASN A  19      11.178  12.632  11.029  1.00 22.94           N
ATOM     10  CA  ASN A  19      10.286  12.621   9.871  1.00 22.82           C
ATOM     11  C   ASN A  19      10.847  11.832   8.689  1.00 20.25           C
ATOM     12  O   ASN A  19      11.419  12.404   7.762  1.00 33.24           O
ATOM     13  CB  ASN A  19       9.930  14.049   9.443  1.00 28.83           C
ATOM     14  CG  ASN A  19       8.884  14.085   8.345  1.00 26.44           C
ATOM     15  OD1 ASN A  19       8.101  13.149   8.185  1.00 22.64           O
ATOM     16  ND2 ASN A  19       8.867  15.171   7.581  1.00 28.91           N
ATOM     17  N   TYR A  20      10.677  10.514   8.730  1.00 23.81           N
ATOM     18  CA  TYR A  20      11.164   9.645   7.666  1.00 35.69           C
ATOM     19  C   TYR A  20      10.078   8.675   7.210  1.00 28.94           C
ATOM     20  O   TYR A  20       9.800   7.684   7.886  1.00 25.26           O
ATOM     21  CB  TYR A  20      12.403   8.874   8.129  1.00 33.37           C
ATOM     22  CG  TYR A  20      13.017   7.996   7.061  1.00 31.48           C
ATOM     23  CD1 TYR A  20      13.880   8.527   6.111  1.00 36.91           C
ATOM     24  CD2 TYR A  20      12.740   6.636   7.007  1.00 21.36           C
ATOM     25  CE1 TYR A  20      14.446   7.729   5.134  1.00 36.56           C
ATOM     26  CE2 TYR A  20      13.301   5.831   6.034  1.00 29.56           C
ATOM     27  CZ  TYR A  20      14.153   6.382   5.100  1.00 35.08           C
ATOM     28  OH  TYR A  20      14.714   5.584   4.130  1.00 38.73           O
ATOM     29  N   ARG A  21       9.475   8.975   6.062  1.00 38.95           N
ATOM     30  CA  ARG A  21       8.416   8.153   5.478  1.00 38.77           C
ATOM     31  C   ARG A  21       7.232   7.966   6.426  1.00 27.69           C
ATOM     32  O   ARG A  21       7.166   6.988   7.171  1.00 22.82           O
ATOM     33  CB  ARG A  21       8.964   6.795   5.023  1.00 20.00           C
ATOM     34  CG  ARG A  21       7.967   5.947   4.248  1.00 20.00           C
ATOM     35  CD  ARG A  21       8.578   4.619   3.832  1.00 20.00           C
ATOM     36  NE  ARG A  21       7.634   3.793   3.085  1.00 20.00           N
ATOM     37  CZ  ARG A  21       7.914   2.585   2.606  1.00 20.00           C
ATOM     38  NH1 ARG A  21       9.116   2.058   2.795  1.00 20.00           N
ATOM     39  NH2 ARG A  21       6.993   1.904   1.938  1.00 20.00           N
ATOM     40  N   GLY A  22       6.299   8.912   6.390  1.00 24.85           N
ATOM     41  CA  GLY A  22       5.121   8.854   7.235  1.00 29.53           C
ATOM     42  C   GLY A  22       5.429   9.164   8.687  1.00 33.22           C
ATOM     43  O   GLY A  22       5.303  10.306   9.128  1.00 30.06           O
ATOM     44  N   TYR A  23       5.835   8.140   9.432  1.00 27.25           N
ATOM     45  CA  TYR A  23       6.161   8.302  10.844  1.00 34.16           C
ATOM     46  C   TYR A  23       7.511   8.988  11.024  1.00 23.48           C
ATOM     47  O   TYR A  23       8.309   9.064  10.089  1.00 39.30           O
ATOM     48  CB  TYR A  23       6.161   6.946  11.553  1.00 29.65           C
ATOM     49  CG  TYR A  23       4.829   6.231  11.512  1.00 34.88           C
ATOM     50  CD1 TYR A  23       3.865   6.460  12.485  1.00 30.77           C
ATOM     51  CD2 TYR A  23       4.536   5.326  10.500  1.00 32.29           C
ATOM     52  CE1 TYR A  23       2.646   5.809  12.452  1.00 39.91           C
ATOM     53  CE2 TYR A  23       3.320   4.670  10.458  1.00 30.45           C
ATOM     54  CZ  TYR A  23       2.379   4.915  11.436  1.00 37.13           C
ATOM     55  OH  TYR A  23       1.167   4.264  11.398  1.00 38.59           O
ATOM     56  N   SER A  24       7.760   9.485  12.231  1.00 36.25           N
ATOM     57  CA  SER A  24       9.014  10.164  12.536  1.00 26.44           C
ATOM     58  C   SER A  24       9.914   9.295  13.408  1.00 20.40           C
ATOM     59  O   SER A  24       9.451   8.660  14.355  1.00 39.27           O
ATOM     60  CB  SER A  24       8.743  11.501  13.229  1.00 21.07           C
ATOM     61  OG  SER A  24       7.969  12.355  12.404  1.00 34.95           O
ATOM     62  N   LEU A  25      11.202   9.273  13.082  1.00 22.77           N
ATOM     63  CA  LEU A  25      12.170   8.483  13.834  1.00 36.60           C
ATOM     64  C   LEU A  25      13.164   9.380  14.565  1.00 28.89           C
ATOM     65  O   LEU A  25      12.775  10.223  15.373  1.00 38.59           O
ATOM     66  CB  LEU A  25      12.917   7.510  12.914  1.00 27.50           C
ATOM     67  CG  LEU A  25      12.178   6.262  12.416  1.00 38.78           C
ATOM     68  CD1 LEU A  25      11.174   6.594  11.319  1.00 22.10           C
ATOM     69  CD2 LEU A  25      13.171   5.212  11.938  1.00 26.92           C
TER
END
"""

def shake_sites(xrs, random, shift, grm=None):
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

def get_pdb_inputs(pdb_str):
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  model = mmtbx.model.manager(model_input = pdb_inp, log = null_out())
  model.process(make_restraints=True)
  return group_args(
    ph  = model.get_hierarchy(),
    grm = model.get_restraints_manager(),
    xrs = model.get_xray_structure())

def exercise_1():
  random.seed(0)
  flex.set_random_seed(0)
  pi = get_pdb_inputs(pdb_str=pdb_str_1)
  f_obs = abs(pi.xrs.structure_factors(d_min = 2.5).f_calc())
  r_free_flags = f_obs.generate_r_free_flags(use_lattice_symmetry=False)
  if(0):
    pi.ph.adopt_xray_structure(pi.xrs)
    pi.ph.write_pdb_file(file_name="start.pdb",
      crystal_symmetry = pi.xrs.crystal_symmetry())
  xrs_poor = shake_sites(xrs = pi.xrs.deep_copy_scatterers(), random=False,
   shift = 1.5, grm=pi.grm)
  if(0):
    pi.ph.adopt_xray_structure(xrs_poor)
    pi.ph.write_pdb_file(file_name="poor.pdb",
      crystal_symmetry = xrs_poor.crystal_symmetry())

  fmodel = mmtbx.f_model.manager(
    f_obs          = f_obs,
    r_free_flags   = r_free_flags,
    xray_structure = xrs_poor)
  print("start r_work:", fmodel.r_work())
  #
  params = sa.master_params().extract()
  params.start_temperature=3000
  params.final_temperature=0
  params.cool_rate = 100
  params.number_of_steps = 100
  params.update_grads_shift = 0.
  #
  sa.run(
    params = params,
    fmodel = fmodel,
    restraints_manager = pi.grm,
    wx                 = 20,
    wc                 = 1,
    verbose            = True)
  #
  r = fmodel.r_work()
  print("final r_work:", r)
  assert r < 0.03, r
  dist = flex.mean(flex.sqrt((pi.xrs.sites_cart() -
          fmodel.xray_structure.sites_cart()).dot()))
  print("Distance(refined, answer): %6.4f"%dist)
  assert dist < 0.25, dist
  if(0):
    pi.ph.adopt_xray_structure(fmodel.xray_structure)
    pi.ph.write_pdb_file(file_name="refined.pdb",
      crystal_symmetry = fmodel.xray_structure.crystal_symmetry())

def exercise_2(d_min = 1.5):
  random.seed(2679941)
  flex.set_random_seed(2679941)
  for shake in [True, False]:
    pi = get_pdb_inputs(pdb_str=pdb_str_1)
    f_obs = abs(pi.xrs.structure_factors(d_min = d_min).f_calc())
    r_free_flags = f_obs.generate_r_free_flags(use_lattice_symmetry=False)
    xrs_poor = pi.xrs.deep_copy_scatterers()
    if(shake):
      xrs_poor = shake_sites(xrs = pi.xrs.deep_copy_scatterers(), random=False,
       shift = 2.0, grm=pi.grm)
    fmodel = mmtbx.f_model.manager(
      f_obs          = f_obs,
      r_free_flags   = r_free_flags,
      xray_structure = xrs_poor)
    print("start r_work:", fmodel.r_work())
    #
    f_calc = pi.xrs.structure_factors(d_min = d_min).f_calc()
    fft_map = f_calc.fft_map(resolution_factor=0.25)
    fft_map.apply_sigma_scaling()
    target_map = fft_map.real_map_unpadded()
    # find optimal weight
    rsr_simple_refiner = mmtbx.refinement.real_space.individual_sites.simple(
      target_map                  = target_map,
      selection                   = flex.bool(pi.xrs.scatterers().size(), True),
      real_space_gradients_delta  = d_min/4,
      max_iterations              = 150,
      geometry_restraints_manager = pi.grm.geometry)
    refined = mmtbx.refinement.real_space.individual_sites.refinery(
      refiner                  = rsr_simple_refiner,
      xray_structure           = xrs_poor.deep_copy_scatterers(),
      start_trial_weight_value = 1,
      rms_bonds_limit          = 0.02,
      rms_angles_limit         = 2)
    print(refined.weight_final, refined.rms_bonds_final, refined.rms_angles_final)
    #
    params = sa.master_params().extract()
    params.start_temperature=5000
    params.final_temperature=0
    params.cool_rate = 100
    params.number_of_steps = 100
    params.update_grads_shift = 0. # does not change runtime visibly
    #
    sa.run(
      params             = params,
      fmodel             = fmodel,
      real_space         = True,
      target_map         = target_map,
      restraints_manager = pi.grm,
      wx                 = refined.weight_final,
      wc                 = 1.,
      verbose            = True)
    #
    r = fmodel.r_work()
    print("final r_work:", r)
    if(shake):
      assert r < 0.07, r
    else:
      assert r < 0.13, r
    dist = flex.mean(flex.sqrt((pi.xrs.sites_cart() -
            fmodel.xray_structure.sites_cart()).dot()))
    print("Distance(refined, answer): %6.4f"%dist)
    if(shake):
      assert dist < 0.35, r
    else:
      assert dist < 0.14, dist
    if(0):
      pi.ph.adopt_xray_structure(fmodel.xray_structure)
      pi.ph.write_pdb_file(file_name="refined.pdb",
        crystal_symmetry = fmodel.xray_structure.crystal_symmetry())

def exercise_3():
  pi = get_pdb_inputs(pdb_str=pdb_str_1)
  xrs = pi.xrs.deep_copy_scatterers()
  sites_cart_start = xrs.sites_cart()
  states_collector = mmtbx.utils.states(
    pdb_hierarchy  = pi.ph)
  #
  params = sa.master_params().extract()
  params.start_temperature=5000
  params.final_temperature=0
  params.cool_rate = 100
  params.number_of_steps = 100
  params.update_grads_shift = 0.
  params.time_step = 0.0005
  params.interleave_minimization=True
  #
  sa.run(
    params = params,
    xray_structure     = xrs,
    restraints_manager = pi.grm,
    states_collector   = states_collector)
  states_collector.write(file_name = "all.pdb")


if(__name__ == "__main__"):
  timer = user_plus_sys_time()
  exercise_1()
  exercise_2()
  exercise_3()
  print("Time: %6.2f" % timer.elapsed())
