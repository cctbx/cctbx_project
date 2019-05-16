from __future__ import absolute_import, division, print_function
from iotbx.pdb.tst_pdb import dump_pdb
from cctbx import geometry_restraints
from iotbx.pymol import pml_stick, pml_write
from cctbx.array_family import flex
from scitbx.matrix import col, sqr
from scitbx.math import euler_angles_as_matrix
from libtbx.test_utils import approx_equal, not_approx_equal
import random
import sys
from six.moves import range
from six.moves import zip

if (1): # fixed random seed to avoid rare failures
  random.seed(0)
  flex.set_random_seed(0)

class residual_functor(object):

  def __init__(self, restraint_type, **keyword_arguments):
    assert "sites" not in keyword_arguments
    self.restraint_type = restraint_type
    self.keyword_arguments = keyword_arguments

  def __call__(self, sites):
    self.keyword_arguments["sites"] = sites
    return self.restraint_type(**self.keyword_arguments).residual()
    del self.keyword_arguments["sites"]

def finite_differences(sites, residual_obj, eps=1.e-6):
  sites_mod = sites[:]
  gradients = []
  for i_site,site in enumerate(sites):
    grad = []
    for i_x in range(3):
      t = []
      for signed_eps in [eps, -eps]:
        add = [0]*3
        add[i_x] = signed_eps
        sites_mod[i_site] = site + col(add)
        r = residual_obj(sites_mod)
        t.append(r)
        sites_mod[i_site] = site
      grad.append((t[0]-t[1])/(2*eps))
    gradients.append(grad)
  return gradients

def exercise_bond():
  for i_trial in range(5):
    sites = random_sites(n_sites=2)
    b = geometry_restraints.bond(sites=sites, distance_ideal=0, weight=0)
    distance_ideal = b.distance_model
    sigma = distance_ideal * 0.01
    weight = 1 / sigma**2
    residual_obj = residual_functor(
      restraint_type=geometry_restraints.bond,
      distance_ideal=distance_ideal,
      weight=weight)
    for i_pert in range(5):
      if (i_pert == 0):
        sites_mod = sites
      else:
        sites_mod = []
        for site in sites:
          shift = col([random.uniform(-3,3)*sigma for i in range(3)])
          sites_mod.append(site+shift)
      b = geometry_restraints.bond(
        sites=sites_mod, distance_ideal=distance_ideal, weight=weight)
      ag = b.gradients()
      fg = finite_differences(sites_mod, residual_obj)
      for analytical,finite in zip(ag,fg):
        assert approx_equal(analytical, finite)

def exercise_nonbonded(nonbonded_type, repulsion_function):
  for i_trial in range(5):
    sites = random_sites(n_sites=2)
    r = nonbonded_type(
      sites=sites,
      vdw_distance=1,
      function=repulsion_function)
    vdw_distance = r.delta
    residual_obj = residual_functor(
      restraint_type=nonbonded_type,
      vdw_distance=vdw_distance,
      function=repulsion_function)
    for i_pert in range(5):
      if (i_pert == 0):
        sites_mod = sites
      else:
        sites_mod = []
        for site in sites:
          shift = col([random.uniform(-1,1)*vdw_distance*0.05
            for i in range(3)])
          sites_mod.append(site+shift)
      r = nonbonded_type(
        sites=sites_mod,
        vdw_distance=vdw_distance,
        function=repulsion_function)
      if (   not hasattr(repulsion_function, "irexp")
          or repulsion_function.irexp <= 2
          or not_approx_equal(r.residual(), 0)):
        fg = flex.vec3_double(finite_differences(sites_mod, residual_obj)) \
          .as_double()
        ag = flex.vec3_double(r.gradients()).as_double()
        scale = max(1, flex.mean(flex.abs(fg)))
        assert approx_equal(ag/scale, fg/scale, eps=1.e-2)
  sites = [[0,0,0], [0,0,0]]
  vdw_distance = 1.8
  for i_step in range(1,100):
    delta = i_step/50.
    sites[1][0] = delta
    r = nonbonded_type(
      sites=[col(site) for site in sites],
      vdw_distance=vdw_distance,
      function=repulsion_function)
    assert approx_equal(delta, r.delta)

def exercise_angle():
  eps = 1
  n_trials = 20
  for i_trial in range(n_trials):
    if (i_trial < -99):
      if (i_trial == n_trials-1):
        eps = 0
      else:
        eps *= 0.1
      sites = [col(site) for site in [(-1,0,0),(0,0,0),(1,eps,0)]]
      angle_ideal = 170
    else:
      while 1:
        sites = random_sites(n_sites=3)
        angle_ideal = geometry_restraints.angle(
          sites, angle_ideal=0, weight=0).angle_model
        if (angle_ideal > 10):
          break
    sigma = angle_ideal * 0.1
    weight = 1 / sigma**2
    residual_obj = residual_functor(
      restraint_type=geometry_restraints.angle,
      angle_ideal=angle_ideal,
      weight=weight)
    for i_pert in range(5):
      if (i_pert == 0):
        sites_mod = sites
      else:
        sites_mod = []
        for site in sites:
          shift = col([random.uniform(3,3)*sigma for i in range(3)])
          sites_mod.append(site+shift)
      a = geometry_restraints.angle(
        sites=sites_mod, angle_ideal=angle_ideal, weight=weight)
      if (eps != 0):
        ag = a.gradients()
      fg = finite_differences(sites_mod, residual_obj)
      if (eps == 0):
        for finite in fg:
          print(finite)
          print()
      else:
        for analytical,finite in zip(ag,fg):
          if (0):
            print(analytical)
            print(finite)
          if (0):
            for x,y in zip(analytical, finite):
              if (y == 0): print(None, end=' ')
              else: print(x/y, end=' ')
            print()
          if (0):
            print()
          assert approx_equal(analytical, finite,
                              eps=max(1.e-6,max(analytical)*1.e-6))

def random_site(max_coordinate):
  return col([random.uniform(-max_coordinate, max_coordinate)
    for i in range(3)])

def random_sites(n_sites, max_coordinate=10, min_distance=0.5,
                 max_trials=1000):
  sites = []
  run_away_counter = 0
  while 1:
    run_away_counter += 1
    assert run_away_counter < max_trials
    new_site = random_site(max_coordinate)
    ok = True
    for site in sites:
      if ((new_site-site).norm_sq() < min_distance**2):
        ok = False
        break
    if (ok):
      sites.append(new_site)
      run_away_counter = 0
      if (len(sites) == n_sites):
        break
  return sites

def exercise_dihedral_core(sites, angle_ideal, angle_esd, periodicity,
                                  angle_model):
  assert periodicity > 0
  for signed_periodicity in [periodicity, -periodicity]:
    dih = geometry_restraints.dihedral(
      sites=sites,
      angle_ideal=angle_ideal,
      weight=1./angle_esd**2,
      periodicity=signed_periodicity)
    if (angle_model is not None):
      assert approx_equal(angle_model, dih.angle_model)
    ag = flex.vec3_double(dih.gradients())
    residual_obj = residual_functor(
      restraint_type=geometry_restraints.dihedral,
      angle_ideal=dih.angle_ideal,
      weight=dih.weight,
      periodicity=signed_periodicity)
    fg = flex.vec3_double(finite_differences(sites, residual_obj)).as_double()
    ag = ag.as_double()
    scale = max(1, flex.mean(flex.abs(fg)))
    assert approx_equal(ag/scale, fg/scale)

def exercise_dihedral():
  sites = [col(site) for site in [
    (1,0,0), (0,0,0), (0,1,0), (1,0,1)]]
  angle_ideal = -45
  for flip in range(2):
    if (flip != 0):
      sites = [sites[i] for i in [0,2,1,3]]
      angle_ideal *= -1
    dih = geometry_restraints.dihedral(
      sites=sites,
      angle_ideal=angle_ideal,
      weight=1)
    assert approx_equal(dih.angle_ideal, angle_ideal)
    assert approx_equal(dih.angle_model, angle_ideal)
    assert approx_equal(dih.delta, 0)
    for i_trial in range(20):
      sites_mod = []
      for site in sites:
        shift = col([random.uniform(-.1,.1) for i in range(3)])
        sites_mod.append(site+shift)
      exercise_dihedral_core(
        sites_mod, angle_ideal, angle_esd=1, periodicity=1, angle_model=None)
  for sites,angle_ideal,angle_esd,period,angle_model in dihedral_test_data:
    sites = [col(site) for site in sites]
    exercise_dihedral_core(
      sites, angle_ideal, angle_esd, period, angle_model)

def improper_permutation(sites):
  return [sites[0], sites[1], sites[3], sites[2]]

def exercise_chirality(verbose=0):
  # monomer library: CA N CB C
  sites = [col(site) for site in [
    (27.660, 9.903, 2.078),
    (28.049, 9.675, 3.486),
    (28.183, 11.269, 1.625),
    (28.085, 8.759, 1.165)]]
  chir = geometry_restraints.chirality(
    sites=sites,
    volume_ideal=-2.48,
    both_signs=False,
    weight=1)
  assert approx_equal(chir.volume_model, -2.411548478)
  assert approx_equal(chir.residual(), 0.00468561086412)
  dih = geometry_restraints.dihedral(
    sites=improper_permutation(sites),
    angle_ideal=35.26439,
    weight=1)
  assert approx_equal(dih.angle_model, 32.2587249641)
  assert approx_equal(dih.residual(), 9.03402230782)
  if (verbose):
    dump_pdb("sites.pdb", sites)
    print("volume_ideal:", chir.volume_ideal)
    print("volume_model:", chir.volume_model)
    print("angle_ideal:", dih.angle_ideal)
    print("angle model:", dih.angle_model)
  for i_trial in range(50):
    volume_ideal = chir.volume_ideal
    for both_signs in [False, True]:
      sites_mod = []
      for site in sites:
        shift = col([random.uniform(-.5,.5) for i in range(3)])
        sites_mod.append(site+shift)
      if (both_signs): volume_ideal = abs(volume_ideal)
      c = geometry_restraints.chirality(
        sites=sites_mod,
        volume_ideal=volume_ideal,
        both_signs=both_signs,
        weight=500*1/0.2**2)
      residual_obj = residual_functor(
        restraint_type=geometry_restraints.chirality,
        volume_ideal=c.volume_ideal,
        both_signs=c.both_signs,
        weight=c.weight)
      fg = flex.vec3_double(
        finite_differences(sites_mod, residual_obj)).as_double()
      ag = flex.vec3_double(
        c.gradients()).as_double()
      scale = max(1, flex.mean(flex.abs(fg)))
      assert approx_equal(ag/scale, fg/scale)
      d = geometry_restraints.dihedral(
        sites=improper_permutation(sites_mod),
        angle_ideal=dih.angle_ideal,
        weight=750)
      ag_dih = dih.gradients()
      if (verbose and i_trial == 0 and not both_signs):
        dump_pdb("sites_mod.pdb", sites_mod)
        max_g_len = 0
        for g in ag_dih:
          max_g_len = max(max_g_len, abs(col(g)))
        for g in flex.vec3_double(fg):
          max_g_len = max(max_g_len, abs(col(g)))
        sticks = []
        for site,g in zip(improper_permutation(sites_mod),ag_dih):
          sticks.append(
            pml_stick(
              begin=site,
              end=site+col(g)/max_g_len,
              colors=[[1,0,0]]*2,
              width=0.01))
        pml_write(f=open("dih.pml", "w"), label="dih", sticks=sticks)
        sticks = []
        for site,g in zip(sites_mod,flex.vec3_double(fg)):
          sticks.append(
            pml_stick(
              begin=site,
              end=site+col(g)/max_g_len,
              colors=[[0,1,0]]*2,
              width=0.01))
        pml_write(f=open("chir.pml", "w"), label="chir", sticks=sticks)

def exercise_planarity():
  weights = flex.double([1,2,3,4])
  for points,norm in [([(1,1,0), (-1,-1,0), (-1,1,0), (1,-1,0)], (0,0,1)),
                      ([(0,1,1), (0,-1,-1), (0,-1,1), (0,1,-1)], (1,0,0)),
                      ([(1,0,1), (-1,0,-1), (-1,0,1), (1,0,-1)], (0,1,0)),
                      ([(1,1,-1), (-1,-1,1), (-1,1,-1), (1,-1,1)], (0,1,1))]:
    norm = col(norm)
    norm /= abs(norm)
    for i_trial in range(5):
      for i_pert in range(5):
        if (i_trial == 0 and i_pert == 0):
          sites = [col(v) for v in points]
          rot_norm = norm
        else:
          shift = col([random.uniform(-10,10) for i in range(3)])
          rot = sqr(euler_angles_as_matrix(
            [random.uniform(0,360) for i in range(3)]))
          rot_norm = rot * norm
          if (i_pert == 0):
            sites = [rot*col(v)+shift for v in points]
          else:
            sites = []
            for v in points:
              f = 0.1
              pert = col([(random.random()-0.5)*f for i in range(3)])
              sites.append(rot*col(v)+shift+pert)
        pl = geometry_restraints.planarity(sites=sites, weights=weights)
        n = col(pl.normal())
        gradients_analytical = pl.gradients()
        if (i_pert == 0):
          assert abs(abs(n.dot(rot_norm))-1) < 1.e-5
          assert approx_equal(pl.residual(), 0)
          for grad in gradients_analytical:
            assert approx_equal(grad, [0,0,0])
        assert approx_equal(pl.residual(), pl.lambda_min())
        residual_obj = residual_functor(
          restraint_type=geometry_restraints.planarity,
          weights=weights)
        gradients_finite = finite_differences(sites, residual_obj)
        assert approx_equal(gradients_finite, gradients_analytical)

def exercise():
  exercise_bond()
  for irexp in [1,2,3,4,5]:
    for rexp in [3,4]:
      exercise_nonbonded(
        nonbonded_type=geometry_restraints.nonbonded_prolsq,
        repulsion_function=geometry_restraints.prolsq_repulsion_function(
          irexp=irexp,
          rexp=rexp))
  for irexp in [1,2,3,4,5]:
    exercise_nonbonded(
      nonbonded_type=geometry_restraints.nonbonded_inverse_power,
      repulsion_function=geometry_restraints.inverse_power_repulsion_function(
        nonbonded_distance_cutoff=1.e20,
        irexp=irexp))
  for exponent in [1,2,3]:
    exercise_nonbonded(
      nonbonded_type=geometry_restraints.nonbonded_cos,
      repulsion_function=geometry_restraints.cos_repulsion_function(
        max_residual=13,
        exponent=exponent))
  for norm_height_at_vdw_distance in [0.1,0.2,0.3]:
    exercise_nonbonded(
      nonbonded_type=geometry_restraints.nonbonded_gaussian,
      repulsion_function=geometry_restraints.gaussian_repulsion_function(
        max_residual=12,
        norm_height_at_vdw_distance=norm_height_at_vdw_distance))
  exercise_angle()
  exercise_dihedral()
  exercise_chirality(verbose="--verbose" in sys.argv[1:])
  exercise_planarity()
  print("OK")

dihedral_test_data = [
[[ [69.141,9.390,2.567],
[70.651,9.597,2.216],
[70.898,11.039,2.558],
[69.552,11.796,2.401],
], 35,15,3,-27.623156 ],
[[ [69.094,9.368,7.421],
[70.388,10.106,7.552],
[70.526,11.172,6.519],
[71.741,11.347,6.212],
], 0,30,2,-149.24817 ],
[[ [66.200,0.138,5.659],
[67.055,-1.079,5.550],
[68.451,-0.802,4.924],
[68.274,-0.562,3.439],
], 180,15,3,-73.186553 ],
[[ [67.330,-1.591,9.213],
[67.723,-0.520,10.245],
[69.155,0.070,10.088],
[70.080,-1.137,10.323],
], 180,15,3,-63.28745 ],
[[ [64.929,-1.607,9.124],
[63.628,-2.233,9.453],
[62.488,-1.199,9.105],
[61.125,-2.031,9.143],
], 180,15,3,166.85678 ],
[[ [70.128,2.016,15.960],
[70.973,2.866,16.769],
[71.244,4.032,15.846],
[70.108,4.001,14.818],
], -25,15,3,21.478587 ],
[[ [51.224,5.278,19.162],
[49.792,5.189,19.476],
[49.560,5.913,20.813],
[50.073,5.054,21.945],
], 180,15,3,-74.662005 ],
[[ [47.467,9.582,20.922],
[46.748,10.422,21.963],
[47.382,10.282,23.330],
[46.730,9.309,24.166],
], 180,15,3,-96.181112 ],
[[ [46.748,10.422,21.963],
[47.382,10.282,23.330],
[46.730,9.309,24.166],
[45.949,8.736,25.060],
], 0,15,4,-95.23452 ],
[[ [54.391,1.288,14.581],
[55.651,1.273,15.387],
[56.467,-0.004,15.273],
[55.681,-1.123,15.821],
], 180,15,3,-63.913454 ],
[[ [67.779,9.293,11.519],
[68.887,10.169,12.058],
[70.267,9.371,12.083],
[71.437,10.417,12.273],
], 180,15,3,164.95936 ],
[[ [74.272,17.569,20.340],
[74.302,18.830,21.056],
[73.048,19.608,20.974],
[72.642,20.465,21.916],
], 0,15,4,152.18051 ],
[[ [67.020,12.001,18.091],
[66.093,10.976,18.685],
[66.532,9.757,19.391],
[65.579,8.603,19.413],
], 60,15,3,-155.18126 ],
[[ [59.525,6.469,10.069],
[60.307,5.499,9.152],
[60.931,6.381,7.983],
[61.840,5.442,7.183],
], 180,15,3,173.24322 ],
[[ [54.325,-2.812,2.710],
[54.755,-1.497,2.220],
[54.490,-1.281,0.739],
[54.841,-2.246,-0.385],
], 180,15,3,-48.495407 ],
[[ [54.449,2.584,5.084],
[55.475,3.162,5.951],
[56.953,2.788,5.575],
[57.159,1.313,5.585],
], 180,15,3,-58.143716 ],
[[ [55.273,5.394,6.776],
[55.070,6.795,6.699],
[53.864,7.253,7.664],
[52.594,6.446,7.313],
], 180,15,3,-55.174556 ],
[[ [55.070,6.795,6.699],
[53.864,7.253,7.664],
[52.594,6.446,7.313],
[52.327,5.272,8.025],
], 90,20,2,93.499752 ],
[[ [59.901,2.524,-2.158],
[60.037,1.248,-1.326],
[60.057,1.643,0.158],
[60.110,0.416,1.003],
], 180,15,3,176.77171 ],
[[ [34.183,3.032,3.282],
[34.441,2.116,4.414],
[35.660,1.339,3.940],
[35.767,1.558,2.468],
], -25,15,3,-16.279611 ],
[[ [35.660,1.339,3.940],
[35.767,1.558,2.468],
[35.251,2.993,2.262],
[34.183,3.032,3.282],
], -30,15,3,-36.045859 ],
[[ [18.701,2.935,18.700],
[17.376,2.652,19.264],
[17.387,1.554,20.356],
[17.769,0.234,19.673],
], 180,15,3,-66.444499 ],
[[ [17.350,2.078,24.829],
[18.086,1.833,26.121],
[17.850,2.399,27.443],
[17.154,2.723,28.532],
], 0,15,4,53.638935 ],
[[ [15.882,6.548,24.854],
[15.228,7.545,25.766],
[16.213,8.309,26.673],
[17.715,8.339,26.304],
], 180,15,3,18.478984 ],
[[ [14.761,12.218,19.123],
[15.765,11.584,18.129],
[17.272,11.907,18.285],
[17.609,13.327,17.950],
], 180,15,3,-69.045664 ],
[[ [15.999,13.559,4.210],
[17.156,14.185,3.544],
[16.653,14.774,2.121],
[15.890,13.732,1.302],
], 180,15,3,-50.961962 ],
[[ [19.166,15.041,4.775],
[19.811,16.046,5.668],
[19.418,17.518,5.404],
[19.549,18.202,4.065],
], 180,15,3,-54.972852 ],
[[ [19.418,17.518,5.404],
[19.549,18.202,4.065],
[20.987,18.306,3.555],
[21.515,19.426,3.487],
], 0,30,2,-111.91728 ],
[[ [21.126,16.443,21.482],
[21.298,16.312,20.045],
[20.406,15.288,19.376],
[20.706,13.843,19.714],
], 180,15,3,-68.764347 ],
[[ [23.143,21.668,10.296],
[22.933,22.521,9.149],
[22.565,23.890,9.595],
[22.302,24.797,8.392],
], 180,15,3,-178.47172 ],
]

if (__name__ == "__main__"):
  exercise()
