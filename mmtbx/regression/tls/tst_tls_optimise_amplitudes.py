from __future__ import absolute_import, division, print_function

import numpy
from mmtbx.tls.utils import TLSMatrices
from mmtbx.tls.optimise_amplitudes import \
  MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator, \
  OptimiseAmplitudes, OptimisationWeights
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from pytest import raises
from six.moves import zip
from six.moves import range

rran = numpy.random.random
iran = numpy.random.choice

numpy.random.seed(123232)

TLSMatrices.set_precision(12)

T_ONLY = [1.,1.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
L_ONLY = [0.,0.,0.,0.,0.,0.,1.,1.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]

TLS_MATRICES = [
  [0.88875254,0.266694873,0.493991604,-0.280846445,0.018132653,0.128202717,0.109792777,0.143780355,0.111757325,0.00745827,0.044408043,0.004683883,-0.031880059,0.087486387,-0.035440044,0.053723107,-0.093285664,-0.126741771,0.078707102,0.047063985,0.125165723],
  [0.244377462,0.459148081,0.710013128,-0.008898371,0.07279876,0.036277241,3.755718081,2.739489028,2.134579946,-0.200385262,0.074588032,-0.023382499,0.00048379,0.229469716,-0.187382678,-0.05705213,0.308356068,0.219549925,0.260721724,0.137599684,-0.308839858],
  [0.29790555,0.765985455,0.362144247,0.188051876,-0.07468928,0.176222813,1.157009238,1.603540537,1.146982202,0.076470473,0.076337332,-0.036514342,-0.057413883,-0.4930283,-0.141063091,0.199094828,0.180380061,0.564268854,0.271176005,0.192103145,-0.122966178],
  [0.550613174,0.642391951,0.512539873,0.186990637,-0.018948302,0.182865203,0.696794749,1.168473012,1.459255982,0.122241174,-0.070965716,-0.127405648,-0.126672963,-0.057326201,0.307248593,0.23402513,-0.081293525,-0.457364662,-0.270995819,0.131106455,0.207966488],
  [0.801839264,0.111567149,0.473172234,0.111448705,0.01322435,0.053023611,0.119600352,0.158201942,0.501511288,-0.00097325,0.038852731,-0.020179123,0.094859837,0.027409079,0.027958886,-0.036977684,0.115684103,-0.023249574,-0.167599437,0.028274946,-0.21054394],
  [0.476189761,0.185513316,0.494838237,0.079787304,0.186996542,0.089705044,6.468559896,0.03302114,1.360918838,-0.092774147,-0.470033825,0.142680478,0.007373318,0.322513481,-0.087204939,-0.102689013,0.085596047,-0.037679444,-0.183875771,-0.048267033,-0.092969365],
  [0.271946335,0.353828607,0.011539801,-0.012010653,-0.03669219,-0.006760052,23.89487928,3.564598018,0.086574587,-1.913464478,0.89104049,-0.011635907,0.00083053,-0.00120248,-0.06596398,0.059314995,0.042902408,0.001580523,-0.013344186,0.103940739,-0.043732938],
  [0.436887477,0.558232952,0.172244124,0.292168093,-0.130328728,-0.074615062,3.137593554,0.025205525,0.920908537,-0.108005337,0.00472811,0.015510299,-0.000822308,0.118599046,-0.151446932,0.029991796,0.042430554,-0.009368142,0.000888338,0.153578145,-0.041608246],
  [0.507487457,0.513473662,0.465150239,-0.024656406,0.052695005,0.030032721,9.358808433,1.113449297,12.326819619,-0.988256782,1.235061392,-1.371684512,0.454760772,-0.522429055,0.566566556,0.34528893,-0.057886461,-0.28449733,-0.171754982,-0.103504915,-0.396874311],
  [0.708808312,0.55187043,0.429938391,0.197257684,-3.9774e-05,-0.054427743,0.178345211,0.297789551,0.087920317,-0.003702672,-0.002424311,0.000192857,-0.029460632,0.037160795,0.061827225,0.276788373,0.000147768,-0.235898673,0.018602792,-0.125798933,0.029312864],
  [0.151113427,0.654574237,0.240882778,-0.018212974,-0.014025963,0.100180189,0.285310135,2.881490187,11.011691132,-0.437653779,-0.698453455,0.157185441,-0.1905967,-0.166253585,-0.151567002,0.150109019,0.444896847,-0.230918972,0.078079958,0.50626905,-0.254300147],
]

def get_optimisation_test_set(
  n_grp,
  n_dst,
  n_atm,
  atomic_amplitude=1.0,
  random_dataset_order=False,
  ):
  """Generate a test set of random uijs for optimisation"""

  assert n_grp > 0
  assert n_grp <= len(TLS_MATRICES)
  assert n_dst > 0
  assert n_atm > 3

  target_uijs = numpy.zeros((n_dst, n_atm, 6))

  # Default to equal weights
  target_weights = flex.double(flex.grid((n_dst, n_atm)), 1.0)

  # Create atomic base
  atomic_values_numpy = rran(n_atm*6).reshape((n_atm,6))
  real_atomic_amps = atomic_amplitude * rran(n_atm)
  atomic_base = flex.sym_mat3_double(atomic_values_numpy)

  # Create atomic total and add to target
  atomic_uijs = real_atomic_amps.reshape((n_atm,1)) * atomic_values_numpy
  for i_dst in range(n_dst):
    target_uijs[i_dst] = atomic_uijs

  # Generate a random set of amplitudes
  real_group_amps = rran(n_grp*n_dst).reshape((n_grp,n_dst))

  # Output lists of base uijs and atoms they are associated with
  base_uijs = []
  base_sels = []

  # Mapping of elements to datasets
  if random_dataset_order is True:
    dataset_hash = flex.size_t(numpy.concatenate([flex.random_permutation(n_dst) for _ in range(n_grp)]))
  else:
    dataset_hash = flex.size_t(list(range(n_dst)) * n_grp)

  # base counter
  i_cnt = 0
  for i_grp in range(n_grp):
    # Number of atoms in this group -- at least 2
    if n_grp == 1:
      n_atm_this = n_atm
    else:
      n_atm_this = max(2,iran(n_atm))
    # Which atoms are covered by this group
    i_sel = iran(n_atm, size=n_atm_this, replace=False)
    # Generate some uijs for this
    tls_m = TLSMatrices(TLS_MATRICES[i_grp])
    for _ in range(n_dst):
      # Which dataset?
      i_dst = dataset_hash[i_cnt]
      i_cnt += 1
      # Generate random coordinates
      x_ = 50.0*(-0.5+rran(n_atm_this*3).reshape((n_atm_this,3)))
      x_ = flex.vec3_double(x_)
      # Generate the base elements
      u_ = tls_m.uijs(x_, origin=(0.,0.,0.))
      base_uijs.append(flex.sym_mat3_double(u_))
      base_sels.append(flex.size_t(i_sel))
      # Generate the amplitude-multiplied equivalent
      u_m = (real_group_amps[i_grp,i_dst]*tls_m).uijs(x_, origin=(0.,0.,0.))
      target_uijs[i_dst, i_sel, :] += numpy.array(u_m)
  # Reshape
  target_uijs = flex.sym_mat3_double(target_uijs.reshape((n_dst*n_atm,6)))
  target_uijs.reshape(flex.grid((n_dst,n_atm)))

  return (
    target_uijs, target_weights, \
    base_uijs, base_sels, dataset_hash, \
    atomic_base, \
    real_group_amps.flatten(), real_atomic_amps,
  )

def resort_amplitudes_by_dataset_hash(
  n_grp,
  n_dst,
  dataset_hash,
  real_group_amps,
  ):
  """Sorted the real amplitudes into the order that they will be returned by the optimiser"""
  sorted_group_amps = []
  i_cnt = 0
  for i_grp in range(n_grp):
    for _ in range(n_dst):
      i_dst = dataset_hash[i_cnt]
      i_cnt += 1
      sorted_group_amps.append(real_group_amps[i_grp*n_dst+i_dst])
  return sorted_group_amps

def calculate_expected_f_and_g_least_squares(
  n_grp,
  n_dst,
  n_atm,
  target_uijs,
  target_weights,
  base_uijs,
  base_sels,
  dataset_hash,
  atomic_base,
  current_amplitudes_base,
  atomic_mask = None,
  ):
  """Calculate what we expect from the current functional and gradients"""

  # Atomics always start with unitary amplitudes
  current_amplitudes_atomic = flex.double(n_atm, 1.0)

  if atomic_mask is None:
    atomic_mask = flex.bool(n_dst, True)

  functional = 0.0
  gradients = numpy.zeros(n_grp*n_dst+n_atm)

  target_uijs_numpy = numpy.array(target_uijs).reshape((n_dst,n_atm,6))
  target_weights_numpy = numpy.array(target_weights).reshape((n_dst,n_atm))
  atomic_uijs_numpy = numpy.array(current_amplitudes_atomic).reshape((n_atm,1)) * numpy.array(atomic_base)

  for i_dst in range(n_dst):
    b_sel = [i_b for i_b, i_d in enumerate(dataset_hash) if i_d==i_dst]
    total_uij = numpy.zeros((n_atm, 6))
    # Add up contributions from different components
    for i_b in b_sel:
      amp = current_amplitudes_base[i_b]
      uijs = base_uijs[i_b]
      atms = base_sels[i_b]
      assert len(uijs) == len(atms)
      for i_atm, u in zip(atms, uijs):
        total_uij[i_atm] += amp*numpy.array(u)
    total_uij += atomic_uijs_numpy
    delta_uijs = (target_uijs_numpy[i_dst] - total_uij)

    # Calculate functional (least squares)
    for u, w in zip(delta_uijs, target_weights_numpy[i_dst]):
      assert len(u) == 6
      functional += numpy.sum(w * u * u)

    # Calculate gradients (for each base element)
    for i_b in b_sel:
      uijs = base_uijs[i_b]
      atms = base_sels[i_b]
      assert len(uijs) == len(atms)
      for i_atm, ug in zip(atms, uijs):
        ud = delta_uijs[i_atm]
        w = target_weights_numpy[i_dst, i_atm]
        gradients[i_b] += numpy.sum(-2.0 * w * ud * ug)

    if atomic_mask[i_dst]:
      for i_atm in range(n_atm):
        ur = atomic_base[i_atm]
        ud = delta_uijs[i_atm]
        w = target_weights_numpy[i_dst, i_atm]
        w2 = float(n_dst) / sum(atomic_mask)
        gradients[-n_atm+i_atm] += numpy.sum(-2.0 * w * ud * ur * w2)

  return functional, gradients

def calculate_expected_f_and_g_sum_of_amplitudes(
  n_grp, n_dst, n_atm,
  target_weights,
  dataset_hash,
  current_amplitudes,
  atomic_base,
  overall_weight,
  atomic_mask = None,
  ):

  mult = overall_weight

  dataset_weights = target_weights.as_numpy_array().mean(axis=1)
  assert len(dataset_weights) == n_dst

  amplitude_sums = numpy.zeros(n_dst)

  assert len(current_amplitudes) == n_grp*n_dst + n_atm

  for i, a in enumerate(current_amplitudes[:n_grp*n_dst]):
    amplitude_sums[dataset_hash[i]] += a

  atomic_amps = (numpy.array(current_amplitudes[-n_atm:]).reshape((n_atm,1)) * numpy.array(atomic_base))
  atomic_amps = atomic_amps[:,0:3].mean(axis=1)
  assert len(atomic_amps) == n_atm

  # Add to all datasets
  amplitude_sums[:] += atomic_amps.sum()

  f = mult * (amplitude_sums * dataset_weights).sum()
  g = flex.double(
    [mult * dataset_weights[i_dst] for i_dst in dataset_hash] +
    [mult * dataset_weights.sum()]*n_atm
  )

  return f,g

def calculate_expected_f_and_g_sum_of_amplitudes_squared(
  n_grp,
  n_dst,
  n_atm,
  target_weights,
  dataset_hash,
  current_amplitudes,
  atomic_base,
  overall_weight,
  atomic_mask = None,
  ):

  mult = overall_weight

  dataset_weights = target_weights.as_numpy_array().mean(axis=1)
  assert len(dataset_weights) == n_dst

  amplitude_sums = numpy.zeros(n_dst)

  assert len(current_amplitudes) == n_grp*n_dst + n_atm

  for i, a in enumerate(current_amplitudes[:n_grp*n_dst]):
    amplitude_sums[dataset_hash[i]] += a

  atomic_amps = (numpy.array(current_amplitudes[-n_atm:]).reshape((n_atm,1)) * numpy.array(atomic_base))
  atomic_amps = atomic_amps[:,0:3].mean(axis=1)
  assert len(atomic_amps) == n_atm

  # Add to all datasets
  amplitude_sums[:] += atomic_amps.sum()

  f = mult * (amplitude_sums * amplitude_sums * dataset_weights).sum()
  g = flex.double(
    [2.0 * mult * amplitude_sums[i_dst] * dataset_weights[i_dst] for i_dst in dataset_hash] +
    [2.0 * mult * (amplitude_sums * dataset_weights).sum()]*n_atm
  )

  return f,g

def calculate_expected_f_and_g_sum_of_squared_amplitudes(
  n_grp,
  n_dst,
  n_atm,
  target_weights,
  dataset_hash,
  current_amplitudes,
  atomic_base,
  overall_weight,
  atomic_mask=None,
  ):

  mult = overall_weight

  dataset_weights = target_weights.as_numpy_array().mean(axis=1)
  assert len(dataset_weights) == n_dst

  dataset_weights_base = numpy.array([dataset_weights[i_dst] for i_dst in dataset_hash])

  assert len(current_amplitudes) == n_grp*n_dst + n_atm

  group_amplitudes = current_amplitudes[:n_grp*n_dst]

  atomic_amps = (numpy.array(current_amplitudes[-n_atm:]).reshape((n_atm,1)) * numpy.array(atomic_base))
  atomic_amps = atomic_amps[:,0:3].mean(axis=1)
  assert len(atomic_amps) == n_atm

  f = mult * (
    (group_amplitudes * group_amplitudes * dataset_weights_base).sum() +
    (atomic_amps * atomic_amps * dataset_weights.sum()).sum()
  )
  g = flex.double(
    list(2.0 * mult * group_amplitudes * dataset_weights_base) +
    list(2.0 * mult * atomic_amps * dataset_weights.sum())
  )

  return f,g

def tst_optimisation_weights():

  wgts_dbl = rran(3)
  wgts_dict = dict(
    sum_of_amplitudes = wgts_dbl[0],
    sum_of_squared_amplitudes = wgts_dbl[1],
    sum_of_amplitudes_squared = wgts_dbl[2],
  )

  # Initialise from dict
  weights = OptimisationWeights(**wgts_dict)
  assert approx_equal(weights.sum_of_amplitudes, wgts_dbl[0])
  assert approx_equal(weights.sum_of_squared_amplitudes, wgts_dbl[1])
  assert approx_equal(weights.sum_of_amplitudes_squared, wgts_dbl[2])

  print('OK')

def tst_optimise_amplitudes_single_group():
  """Check that very basic optimisation reaches expected values regardless of starting values"""

  n_grp = 1
  n_atm = 10

  for n_dst in [1,5]:
    for g_amp_start in [0.0, 1.0, None]:
      for r_amp_mult in [0.0, 1.0]:

        # Get a random set of values to test
        target_uijs, target_weights, \
          base_uijs, base_sels, dataset_hash, \
          atomic_base, \
          real_group_amps, real_atomic_amps \
            = get_optimisation_test_set(n_grp, n_dst, n_atm, atomic_amplitude=r_amp_mult)

        if r_amp_mult == 0.0:
          assert not real_atomic_amps.any()

        # Starting values
        if g_amp_start is not None:
          base_amplitudes_start = flex.double(n_grp*n_dst, g_amp_start)
        else:
          base_amplitudes_start = flex.double(rran(n_grp*n_dst))

        # Since not given, amplitudes start at 1.0 for atomic
        atomic_amplitudes_start = flex.double(n_atm, 1.0)

        # Optimise
        opt = OptimiseAmplitudes(
          target_uijs = target_uijs,
          target_weights = target_weights,
          base_amplitudes = base_amplitudes_start,
          base_uijs = base_uijs,
          base_atom_indices = base_sels,
          base_dataset_hash = dataset_hash,
          atomic_uijs = atomic_base,
          atomic_amplitudes = None,
          atomic_optimisation_mask = None,
          optimisation_weights = None,
          convergence_tolerance=1e-08,
        ).run()

        assert approx_equal(list(opt.initial), list(base_amplitudes_start) + list(atomic_amplitudes_start), 1e-6)
        assert approx_equal(list(opt.result), list(real_group_amps) + list(real_atomic_amps), 1e-6)

  print('OK')

def tst_optimise_amplitudes_multiple_groups_with_atomic():
  """Check that multiple partial groups optimise correctly with atomic level"""

  n_grp = 3
  n_dst = 5
  n_atm = 10

  # Get a random set of values to test
  target_uijs, target_weights, \
    base_uijs, base_sels, dataset_hash, \
    atomic_base, \
    real_group_amps, real_atomic_amps \
      = get_optimisation_test_set(n_grp, n_dst, n_atm)

  # Starting values
  base_amplitudes_start = flex.double(n_grp*n_dst, 0.0)
  atomic_amplitudes_start = flex.double(n_atm, 0.0)

  opt = OptimiseAmplitudes(
    target_uijs = target_uijs,
    target_weights = target_weights,
    base_amplitudes = base_amplitudes_start,
    base_uijs = base_uijs,
    base_atom_indices = base_sels,
    base_dataset_hash = dataset_hash,
    atomic_uijs = atomic_base,
    atomic_amplitudes = atomic_amplitudes_start,
    atomic_optimisation_mask = None,
    optimisation_weights = None,
    convergence_tolerance=1e-08,
  ).run()

  assert approx_equal(list(opt.initial), list(base_amplitudes_start) + list(atomic_amplitudes_start), 1e-6)
  assert approx_equal(list(opt.result), list(real_group_amps) + list(real_atomic_amps), 1e-6)

  ########################################################
  # Check that shuffling the input lists has the expected effect
  ########################################################

  n_base = n_grp * n_dst
  i_perm = iran(n_base, size=n_base, replace=False)
  # Reorder the base elements by random permutation
  base_uijs = [base_uijs[i] for i in i_perm]
  base_sels = [base_sels[i] for i in i_perm]
  dataset_hash = flex.size_t([dataset_hash[i] for i in i_perm])
  base_amplitudes_start = flex.double([base_amplitudes_start[i] for i in i_perm])

  opt = OptimiseAmplitudes(
    target_uijs = target_uijs,
    target_weights = target_weights,
    base_amplitudes = base_amplitudes_start,
    base_uijs = base_uijs,
    base_atom_indices = base_sels,
    base_dataset_hash = dataset_hash,
    atomic_uijs = atomic_base,
    atomic_amplitudes = atomic_amplitudes_start,
    atomic_optimisation_mask = None,
    optimisation_weights = None,
    convergence_tolerance=1e-08,
  ).run()

  assert approx_equal(list(opt.initial), list(base_amplitudes_start) + list(atomic_amplitudes_start), 1e-6)
  assert approx_equal(list(opt.result), list([real_group_amps[i] for i in i_perm]) + list(real_atomic_amps), 1e-6)

  print('OK')

def tst_optimise_amplitudes_multiple_groups_permuted_dataset_order():

  n_grp = 3
  n_dst = 5
  n_atm = 10

  # Get a random set of values to test
  target_uijs, target_weights, \
    base_uijs, base_sels, dataset_hash, \
    atomic_base, \
    real_group_amps, real_atomic_amps \
      = get_optimisation_test_set(n_grp, n_dst, n_atm, random_dataset_order=True)

  sorted_group_amps = resort_amplitudes_by_dataset_hash(
    n_grp=n_grp, n_dst=n_dst,
    dataset_hash=dataset_hash,
    real_group_amps=real_group_amps,
  )

  # Starting values
  base_amplitudes_start = flex.double(n_grp*n_dst, 0.0)
  atomic_amplitudes_start = flex.double(n_atm, 0.0)

  opt = OptimiseAmplitudes(
    target_uijs = target_uijs,
    target_weights = target_weights,
    base_amplitudes = base_amplitudes_start,
    base_uijs = base_uijs,
    base_atom_indices = base_sels,
    base_dataset_hash = dataset_hash,
    atomic_uijs = atomic_base,
    atomic_amplitudes = atomic_amplitudes_start,
    atomic_optimisation_mask = None,
    optimisation_weights = None,
    convergence_tolerance=1e-08,
  ).run()

  assert approx_equal(list(opt.initial), list(base_amplitudes_start) + list(atomic_amplitudes_start), 1e-6)
  assert approx_equal(list(opt.result), list(sorted_group_amps) + list(real_atomic_amps), 1e-6)

  print('OK')

def tst_functional_gradient_calculator():
  """Check that the functional and gradients are calculated as expected"""

  n_grp = 3
  n_dst = 5
  n_atm = 10

  for atomic_mask_size in [None, 1, 2, 3]:
    if atomic_mask_size is not None:
      atomic_mask = flex.bool(n_dst, False)
      n_sel = max(n_dst, atomic_mask_size)
      atomic_mask.set_selected(flex.size_t(iran(n_dst, size=n_sel, replace=False)), True)
      assert sum(atomic_mask) == n_sel
    else:
      atomic_mask = None

    target_uijs, target_weights, \
      base_uijs, base_sels, dataset_hash, \
      atomic_base, \
      real_group_amps, real_atomic_amps \
        = get_optimisation_test_set(n_grp, n_dst, n_atm)

    # Starting values - randomise
    base_amplitudes_start = flex.random_double(n_grp*n_dst)

    ########################################################
    # Check least-squares f/g
    ########################################################

    functional, gradients = calculate_expected_f_and_g_least_squares(
      n_grp = n_grp,
      n_dst = n_dst,
      n_atm = n_atm,
      target_uijs = target_uijs,
      target_weights = target_weights,
      base_uijs = base_uijs,
      base_sels = base_sels,
      dataset_hash = dataset_hash,
      atomic_base = atomic_base,
      current_amplitudes_base = base_amplitudes_start,
      atomic_mask = atomic_mask,
    )

    f_g_calculator = MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
      target_uijs = target_uijs,
      target_weights = target_weights,
      base_amplitudes = base_amplitudes_start,
      base_uijs = base_uijs,
      base_atom_indices = base_sels,
      base_dataset_hash = dataset_hash,
      atomic_uijs = atomic_base,
      weight_sum_of_amplitudes = 0.0,
      weight_sum_of_amplitudes_squared = 0.0,
      weight_sum_of_squared_amplitudes = 0.0,
    )
    if atomic_mask:
      f_g_calculator.set_atomic_optimisation_mask(atomic_mask)

    f, g = f_g_calculator.compute_functional_and_gradients()

    assert approx_equal(f, functional)
    assert approx_equal(list(g), list(gradients))

    ########################################################
    # check that correct solution gives no functional and gradient
    ########################################################

    correct_amplitudes = flex.double(list(real_group_amps) + list(real_atomic_amps))

    # Manually set amplitudes
    f_g_calculator.set_current_amplitudes(correct_amplitudes)

    f, g = f_g_calculator.compute_functional_and_gradients()

    assert approx_equal(f, 0.0)
    assert g.all_approx_equal(0.0)

    ########################################################
    # Calculate gradients for weight_sum_of_amplitudes
    ########################################################

    # !!! Still using the correct solution so that those functional and gradient components are zero

    wgt = rran()

    f_g_calculator = MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
      target_uijs = target_uijs,
      target_weights = target_weights,
      base_amplitudes = base_amplitudes_start,
      base_uijs = base_uijs,
      base_atom_indices = base_sels,
      base_dataset_hash = dataset_hash,
      atomic_uijs = atomic_base,
      weight_sum_of_amplitudes = wgt,
      weight_sum_of_amplitudes_squared = 0.0,
      weight_sum_of_squared_amplitudes = 0.0,
    )
    if atomic_mask:
      f_g_calculator.set_atomic_optimisation_mask(atomic_mask)

    f_g_calculator.set_current_amplitudes(correct_amplitudes)

    functional, gradients = calculate_expected_f_and_g_sum_of_amplitudes(
      n_grp = n_grp,
      n_dst = n_dst,
      n_atm = n_atm,
      target_weights = target_weights,
      dataset_hash = dataset_hash,
      current_amplitudes = f_g_calculator.get_current_amplitudes(),
      atomic_base = atomic_base,
      overall_weight = wgt,
    )

    f, g = f_g_calculator.compute_functional_and_gradients()

    assert approx_equal(f, functional)
    assert approx_equal(list(g), list(gradients))

    ########################################################
    # Calculate gradients for weight_sum_of_amplitudes_squared
    ########################################################

    # !!! Still using the correct solution so that those functional and gradient components are zero

    wgt = rran()

    f_g_calculator = MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
      target_uijs = target_uijs,
      target_weights = target_weights,
      base_amplitudes = base_amplitudes_start,
      base_uijs = base_uijs,
      base_atom_indices = base_sels,
      base_dataset_hash = dataset_hash,
      atomic_uijs = atomic_base,
      weight_sum_of_amplitudes = 0.0,
      weight_sum_of_amplitudes_squared = wgt,
      weight_sum_of_squared_amplitudes = 0.0,
    )
    if atomic_mask:
      f_g_calculator.set_atomic_optimisation_mask(atomic_mask)

    f_g_calculator.set_current_amplitudes(correct_amplitudes)

    functional, gradients = calculate_expected_f_and_g_sum_of_amplitudes_squared(
      n_grp = n_grp,
      n_dst = n_dst,
      n_atm = n_atm,
      target_weights = target_weights,
      dataset_hash = dataset_hash,
      current_amplitudes = f_g_calculator.get_current_amplitudes(),
      atomic_base = atomic_base,
      overall_weight = wgt,
    )

    f, g = f_g_calculator.compute_functional_and_gradients()

    assert approx_equal(f, functional)
    assert approx_equal(list(g), list(gradients))

    ########################################################
    # Calculate gradients for weight_sum_of_squared_amplitudes
    ########################################################

    # !!! Still using the correct solution so that those functional and gradient components are zero

    wgt = rran()

    f_g_calculator = MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
      target_uijs = target_uijs,
      target_weights = target_weights,
      base_amplitudes = base_amplitudes_start,
      base_uijs = base_uijs,
      base_atom_indices = base_sels,
      base_dataset_hash = dataset_hash,
      atomic_uijs = atomic_base,
      weight_sum_of_amplitudes = 0.0,
      weight_sum_of_amplitudes_squared = 0.0,
      weight_sum_of_squared_amplitudes = wgt,
    )
    if atomic_mask:
      f_g_calculator.set_atomic_optimisation_mask(atomic_mask)

    f_g_calculator.set_current_amplitudes(correct_amplitudes)

    functional, gradients = calculate_expected_f_and_g_sum_of_squared_amplitudes(
      n_grp = n_grp,
      n_dst = n_dst,
      n_atm = n_atm,
      target_weights = target_weights,
      dataset_hash = dataset_hash,
      current_amplitudes = f_g_calculator.get_current_amplitudes(),
      atomic_base = atomic_base,
      overall_weight = wgt,
    )

    f, g = f_g_calculator.compute_functional_and_gradients()

    assert approx_equal(f, functional)
    assert approx_equal(list(g), list(gradients))

    ########################################################
    # Turn on all weights and non-optimal solution
    ########################################################

    wgts = rran(3)

    f_g_calculator = MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
      target_uijs = target_uijs,
      target_weights = target_weights,
      base_amplitudes = base_amplitudes_start,
      base_uijs = base_uijs,
      base_atom_indices = base_sels,
      base_dataset_hash = dataset_hash,
      atomic_uijs = atomic_base,
      weight_sum_of_amplitudes = wgts[0],
      weight_sum_of_amplitudes_squared = wgts[1],
      weight_sum_of_squared_amplitudes = wgts[2],
    )
    if atomic_mask:
      f_g_calculator.set_atomic_optimisation_mask(atomic_mask)

    flsq, glsq = calculate_expected_f_and_g_least_squares(
      n_grp = n_grp,
      n_dst = n_dst,
      n_atm = n_atm,
      target_uijs = target_uijs,
      target_weights = target_weights,
      base_uijs = base_uijs,
      base_sels = base_sels,
      dataset_hash = dataset_hash,
      atomic_base = atomic_base,
      current_amplitudes_base = base_amplitudes_start,
      atomic_mask = atomic_mask,
    )

    fw1, gw1 = calculate_expected_f_and_g_sum_of_amplitudes(
      n_grp = n_grp,
      n_dst = n_dst,
      n_atm = n_atm,
      target_weights = target_weights,
      dataset_hash = dataset_hash,
      current_amplitudes = f_g_calculator.get_current_amplitudes(),
      atomic_base = atomic_base,
      overall_weight = wgts[0],
    )

    fw2, gw2 = calculate_expected_f_and_g_sum_of_amplitudes_squared(
      n_grp = n_grp,
      n_dst = n_dst,
      n_atm = n_atm,
      target_weights = target_weights,
      dataset_hash = dataset_hash,
      current_amplitudes = f_g_calculator.get_current_amplitudes(),
      atomic_base = atomic_base,
      overall_weight = wgts[1],
    )

    fw3, gw3 = calculate_expected_f_and_g_sum_of_squared_amplitudes(
      n_grp = n_grp,
      n_dst = n_dst,
      n_atm = n_atm,
      target_weights = target_weights,
      dataset_hash = dataset_hash,
      current_amplitudes = f_g_calculator.get_current_amplitudes(),
      atomic_base = atomic_base,
      overall_weight = wgts[2],
    )

    f, g = f_g_calculator.compute_functional_and_gradients()

    assert approx_equal(f, flsq+fw1+fw2+fw3)
    assert approx_equal(list(g), list(glsq+gw1+gw2+gw3))

    ########################################################
    # Make some amplitudes negative and recalculate - negative amplitudes should be set to zero
    ########################################################

    # Starting values
    base_amplitudes_start = flex.double(n_grp*n_dst, 0.0)

    functional, gradients = calculate_expected_f_and_g_least_squares(
      n_grp = n_grp,
      n_dst = n_dst,
      n_atm = n_atm,
      target_uijs = target_uijs,
      target_weights = target_weights,
      base_uijs = base_uijs,
      base_sels = base_sels,
      dataset_hash = dataset_hash,
      atomic_base = atomic_base,
      current_amplitudes_base = base_amplitudes_start,
      atomic_mask = atomic_mask,
    )

    f_g_calculator = MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
      target_uijs = target_uijs,
      target_weights = target_weights,
      base_amplitudes = base_amplitudes_start,
      base_uijs = base_uijs,
      base_atom_indices = base_sels,
      base_dataset_hash = dataset_hash,
      atomic_uijs = atomic_base,
      weight_sum_of_amplitudes = 0.0,
      weight_sum_of_amplitudes_squared = 0.0,
      weight_sum_of_squared_amplitudes = 0.0,
    )
    if atomic_mask:
      f_g_calculator.set_atomic_optimisation_mask(atomic_mask)

    f, g = f_g_calculator.compute_functional_and_gradients()

    assert approx_equal(f, functional)
    assert approx_equal(list(g), list(gradients))

    # Now make some amplitudes negative
    base_amplitudes_start.set_selected(flex.size_t(iran(n_grp*n_dst)), -0.5)
    base_amplitudes_start.set_selected(flex.size_t(iran(n_grp*n_dst)), -200.)

    f_g_calculator = MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
      target_uijs = target_uijs,
      target_weights = target_weights,
      base_amplitudes = base_amplitudes_start,
      base_uijs = base_uijs,
      base_atom_indices = base_sels,
      base_dataset_hash = dataset_hash,
      atomic_uijs = atomic_base,
      weight_sum_of_amplitudes = 0.0,
      weight_sum_of_amplitudes_squared = 0.0,
      weight_sum_of_squared_amplitudes = 0.0,
    )
    if atomic_mask:
      f_g_calculator.set_atomic_optimisation_mask(atomic_mask)

    f, g = f_g_calculator.compute_functional_and_gradients()

    assert approx_equal(f, functional)
    assert approx_equal(list(g), list(gradients))

  print('OK')

def tst_functional_gradient_calculator_permuted_dataset_order():

  n_grp = 3
  n_dst = 5
  n_atm = 10

  for atomic_mask_size in [None, 1, 2, 3]:
    if atomic_mask_size is not None:
      atomic_mask = flex.bool(n_dst, False)
      n_sel = max(n_dst, atomic_mask_size)
      atomic_mask.set_selected(flex.size_t(iran(n_dst, size=n_sel, replace=False)), True)
      assert sum(atomic_mask) == n_sel
    else:
      atomic_mask = None

    target_uijs, target_weights, \
      base_uijs, base_sels, dataset_hash, \
      atomic_base, \
      real_group_amps, real_atomic_amps \
        = get_optimisation_test_set(n_grp, n_dst, n_atm, random_dataset_order=True)

    # Starting values - randomise
    base_amplitudes_start = flex.random_double(n_grp*n_dst)

    wgts = rran(3)

    f_g_calculator = MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
      target_uijs = target_uijs,
      target_weights = target_weights,
      base_amplitudes = base_amplitudes_start,
      base_uijs = base_uijs,
      base_atom_indices = base_sels,
      base_dataset_hash = dataset_hash,
      atomic_uijs = atomic_base,
      weight_sum_of_amplitudes = wgts[0],
      weight_sum_of_amplitudes_squared = wgts[1],
      weight_sum_of_squared_amplitudes = wgts[2],
    )
    if atomic_mask:
      f_g_calculator.set_atomic_optimisation_mask(atomic_mask)

    flsq, glsq = calculate_expected_f_and_g_least_squares(
      n_grp = n_grp,
      n_dst = n_dst,
      n_atm = n_atm,
      target_uijs = target_uijs,
      target_weights = target_weights,
      base_uijs = base_uijs,
      base_sels = base_sels,
      dataset_hash = dataset_hash,
      atomic_base = atomic_base,
      current_amplitudes_base = base_amplitudes_start,
      atomic_mask = atomic_mask,
    )

    fw1, gw1 = calculate_expected_f_and_g_sum_of_amplitudes(
      n_grp = n_grp,
      n_dst = n_dst,
      n_atm = n_atm,
      target_weights = target_weights,
      dataset_hash = dataset_hash,
      current_amplitudes = f_g_calculator.get_current_amplitudes(),
      atomic_base = atomic_base,
      overall_weight = wgts[0],
    )

    fw2, gw2 = calculate_expected_f_and_g_sum_of_amplitudes_squared(
      n_grp = n_grp,
      n_dst = n_dst,
      n_atm = n_atm,
      target_weights = target_weights,
      dataset_hash = dataset_hash,
      current_amplitudes = f_g_calculator.get_current_amplitudes(),
      atomic_base = atomic_base,
      overall_weight = wgts[1],
    )

    fw3, gw3 = calculate_expected_f_and_g_sum_of_squared_amplitudes(
      n_grp = n_grp,
      n_dst = n_dst,
      n_atm = n_atm,
      target_weights = target_weights,
      dataset_hash = dataset_hash,
      current_amplitudes = f_g_calculator.get_current_amplitudes(),
      atomic_base = atomic_base,
      overall_weight = wgts[2],
    )

    f, g = f_g_calculator.compute_functional_and_gradients()

    assert approx_equal(f, flsq+fw1+fw2+fw3)
    assert approx_equal(list(g), list(glsq+gw1+gw2+gw3))

  print('OK')

def tst_functional_gradient_calculator_invalid_arguments():
  """Check errors are raised as expected"""

  n_grp = 3
  n_dst = 5
  n_atm = 10

  target_uijs, target_weights, \
    base_uijs, base_sels, dataset_hash, \
    atomic_base, \
    real_group_amps, real_atomic_amps \
      = get_optimisation_test_set(n_grp, n_dst, n_atm)

  ########################################################
  # Check expected error messages are raised
  ########################################################

  # Starting values
  base_amplitudes_start = flex.double(n_grp*n_dst, 1.0)

  wgt_kw_args = dict(
    weight_sum_of_amplitudes = 0.0,
    weight_sum_of_amplitudes_squared = 0.0,
    weight_sum_of_squared_amplitudes = 0.0,
  )

  # Should not error
  f_g_calculator = MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
    target_uijs = target_uijs,
    target_weights = target_weights,
    base_amplitudes = base_amplitudes_start,
    base_uijs = base_uijs,
    base_atom_indices = base_sels,
    base_dataset_hash = dataset_hash,
    atomic_uijs = atomic_base,
    **wgt_kw_args
  )
  f, g = f_g_calculator.compute_functional_and_gradients()

  # target_uijs

  msg = "invalid target_uijs: must be 2-dimensional flex array (currently 3)"
  with raises(Exception) as e:
    f_g_calculator = MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
      target_uijs = flex.sym_mat3_double(flex.grid((n_dst-1,n_atm,5)), (1.,1.,1.,0.,0.,0.)),
      target_weights = target_weights,
      base_amplitudes = base_amplitudes_start,
      base_uijs = base_uijs,
      base_atom_indices = base_sels,
      base_dataset_hash = dataset_hash,
      atomic_uijs = atomic_base,
      **wgt_kw_args
    )
  assert msg == str(e.value), '"{}" does not match "{}"'.format(msg, str(e.value))

  # target_weights

  msg = "invalid dimension of target_weights (dimension 3): must be same dimension as target_uijs (dimension 2)"
  with raises(Exception) as e:
    f_g_calculator = MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
      target_uijs = target_uijs,
      target_weights = flex.double(flex.grid((n_dst,n_atm,5)), 1.0),
      base_amplitudes = base_amplitudes_start,
      base_uijs = base_uijs,
      base_atom_indices = base_sels,
      base_dataset_hash = dataset_hash,
      atomic_uijs = atomic_base,
      **wgt_kw_args
    )
  assert msg == str(e.value)

  msg = "incompatible dimension of target_weights (axis 0): must be same size as target_uijs ({} != {})".format(n_dst, n_dst-1)
  with raises(Exception) as e:
    f_g_calculator = MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
      target_uijs = flex.sym_mat3_double(flex.grid((n_dst-1,n_atm)), (1.,1.,1.,0.,0.,0.)),
      target_weights = target_weights,
      base_amplitudes = base_amplitudes_start,
      base_uijs = base_uijs,
      base_atom_indices = base_sels,
      base_dataset_hash = dataset_hash,
      atomic_uijs = atomic_base,
      **wgt_kw_args
    )
  assert msg == str(e.value)

  msg = "incompatible dimension of target_weights (axis 1): must be same size as target_uijs ({} != {})".format(n_atm, n_atm+1)
  with raises(Exception) as e:
    f_g_calculator = MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
      target_uijs = flex.sym_mat3_double(flex.grid((n_dst,n_atm+1)), (1.,1.,1.,0.,0.,0.)),
      target_weights = target_weights,
      base_amplitudes = base_amplitudes_start,
      base_uijs = base_uijs,
      base_atom_indices = base_sels,
      base_dataset_hash = dataset_hash,
      atomic_uijs = atomic_base,
      **wgt_kw_args
    )
  assert msg == str(e.value)

  # base components

  msg = "invalid input base components. base_amplitudes (length {}), base_uijs (length {}) and base_atom_indices (length {}) must all be the same length"
  with raises(Exception) as e:
    f_g_calculator = MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
      target_uijs = target_uijs,
      target_weights = target_weights,
      base_amplitudes = flex.double(n_grp*n_dst-1, 1.0),
      base_uijs = base_uijs,
      base_atom_indices = base_sels,
      base_dataset_hash = dataset_hash,
      atomic_uijs = atomic_base,
      **wgt_kw_args
    )
  assert msg.format(n_grp*n_dst-1, n_grp*n_dst, n_grp*n_dst) == str(e.value)
  with raises(Exception) as e:
    f_g_calculator = MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
      target_uijs = target_uijs,
      target_weights = target_weights,
      base_amplitudes = base_amplitudes_start,
      base_uijs = base_uijs[:-1],
      base_atom_indices = base_sels[:-1],
      base_dataset_hash = dataset_hash,
      atomic_uijs = atomic_base,
      **wgt_kw_args
    )
  assert msg.format(n_grp*n_dst, n_grp*n_dst-1, n_grp*n_dst-1) == str(e.value)
  with raises(Exception) as e:
    f_g_calculator = MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
      target_uijs = target_uijs,
      target_weights = target_weights,
      base_amplitudes = base_amplitudes_start,
      base_uijs = base_uijs[:-1],
      base_atom_indices = base_sels,
      base_dataset_hash = dataset_hash,
      atomic_uijs = atomic_base,
      **wgt_kw_args
    )
  assert msg.format(n_grp*n_dst, n_grp*n_dst-1, n_grp*n_dst) == str(e.value)
  with raises(Exception) as e:
    f_g_calculator = MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
      target_uijs = target_uijs,
      target_weights = target_weights,
      base_amplitudes = base_amplitudes_start,
      base_uijs = base_uijs,
      base_atom_indices = base_sels[:-1],
      base_dataset_hash = dataset_hash,
      atomic_uijs = atomic_base,
      **wgt_kw_args
    )
  assert msg.format(n_grp*n_dst, n_grp*n_dst, n_grp*n_dst-1) == str(e.value)

  msg = "incompatible pair (element 2) in base_uijs/base_atom_indices: pairwise elements must be the same length ({} and {})".format(n_atm+1, len(base_uijs[2]))
  with raises(Exception) as e:
    f_g_calculator = MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
      target_uijs = target_uijs,
      target_weights = target_weights,
      base_amplitudes = base_amplitudes_start,
      base_uijs = base_uijs[:2] + [flex.sym_mat3_double(n_atm+1)] + base_uijs[3:],
      base_atom_indices = base_sels,
      base_dataset_hash = dataset_hash,
      atomic_uijs = atomic_base,
      **wgt_kw_args
    )
  assert msg == str(e.value)

  msg = "incompatible pair (element 2) in base_uijs/base_atom_indices: pairwise elements must be the same length ({} and {})".format(len(base_sels[2]), n_atm+1)
  with raises(Exception) as e:
    f_g_calculator = MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
      target_uijs = target_uijs,
      target_weights = target_weights,
      base_amplitudes = base_amplitudes_start,
      base_uijs = base_uijs,
      base_atom_indices = base_sels[:2] + [flex.size_t_range(n_atm+1)] + base_sels[3:],
      base_dataset_hash = dataset_hash,
      atomic_uijs = atomic_base,
      **wgt_kw_args
    )
  assert msg == str(e.value)

  msg = "invalid selection in base_atom_indices ({}): attempting to select atom outside of array (size {})".format(n_atm, n_atm)
  with raises(Exception) as e:
    f_g_calculator = MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
      target_uijs = target_uijs,
      target_weights = target_weights,
      base_amplitudes = base_amplitudes_start,
      base_uijs = base_uijs[:2] + [flex.sym_mat3_double(n_atm+1)] + base_uijs[3:],
      base_atom_indices = base_sels[:2] + [flex.size_t_range(n_atm+1)] + base_sels[3:],
      base_dataset_hash = dataset_hash,
      atomic_uijs = atomic_base,
      **wgt_kw_args
    )
  assert msg == str(e.value)

  # dataset_hash

  msg = "invalid base_dataset_hash (length {}): must be same length as base_amplitudes, base_uijs & base_atom_indices (length {})".format(len(dataset_hash)-1, len(dataset_hash))
  with raises(Exception) as e:
    f_g_calculator = MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
      target_uijs = target_uijs,
      target_weights = target_weights,
      base_amplitudes = base_amplitudes_start,
      base_uijs = base_uijs,
      base_atom_indices = base_sels,
      base_dataset_hash = flex.size_t(list(dataset_hash)[:-1]),
      atomic_uijs = atomic_base,
      **wgt_kw_args
    )
  assert msg == str(e.value), (msg, str(e.value))

  msg = "invalid value in base_dataset_hash ({}): attempts to select element outside range of target_uijs (size {})".format(n_dst-1, n_dst-1)
  with raises(Exception) as e:
    f_g_calculator = MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
      target_uijs = flex.sym_mat3_double(flex.grid((n_dst-1,n_atm)), (1.,1.,1.,0.,0.,0.)),
      target_weights = flex.double(flex.grid((n_dst-1,n_atm)), 1.0), # need to resize this also
      base_amplitudes = base_amplitudes_start,
      base_uijs = base_uijs,
      base_atom_indices = base_sels,
      base_dataset_hash = dataset_hash,
      atomic_uijs = atomic_base,
      **wgt_kw_args
    )
  assert msg == str(e.value), (msg, str(e.value))

  msg = "Dataset index {} is not present in base_dataset_hash -- this dataset has no base elements associated with it.".format(n_dst)
  with raises(Exception) as e:
    f_g_calculator = MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
      target_uijs = flex.sym_mat3_double(flex.grid((n_dst+1,n_atm)), (1.,1.,1.,0.,0.,0.)),
      target_weights = flex.double(flex.grid((n_dst+1,n_atm)), 1.0), # need to resize this also
      base_amplitudes = base_amplitudes_start,
      base_uijs = base_uijs,
      base_atom_indices = base_sels,
      base_dataset_hash = dataset_hash,
      atomic_uijs = atomic_base,
      **wgt_kw_args
    )
  assert msg == str(e.value), (msg, str(e.value))

  # atomic uijs

  msg = "invalid size of atomic_uijs ({}): must match 2nd dimension of target_uijs ({})".format(n_atm-1, n_atm)
  with raises(Exception) as e:
    f_g_calculator = MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
      target_uijs = target_uijs,
      target_weights = target_weights,
      base_amplitudes = base_amplitudes_start,
      base_uijs = base_uijs,
      base_atom_indices = base_sels,
      base_dataset_hash = dataset_hash,
      atomic_uijs = flex.sym_mat3_double(n_atm-1),
      **wgt_kw_args
    )
  assert msg == str(e.value), (msg, str(e.value))

  # atomic mask

  f_g_calculator = MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
    target_uijs = target_uijs,
    target_weights = target_weights,
    base_amplitudes = base_amplitudes_start,
    base_uijs = base_uijs,
    base_atom_indices = base_sels,
    base_dataset_hash = dataset_hash,
    atomic_uijs = atomic_base,
    **wgt_kw_args
  )
  # should not error
  f_g_calculator.set_atomic_optimisation_mask(flex.bool(n_dst, True))
  f_g_calculator.set_atomic_optimisation_mask(flex.bool(n_dst, False))
  # should error
  msg = "Input array (size {}) must be the same length as number of datasets ({})".format(n_dst-1, n_dst)
  with raises(Exception) as e:
    f_g_calculator.set_atomic_optimisation_mask(flex.bool(n_dst-1, True))
  assert msg == str(e.value)
  msg = "Input array (size {}) must be the same length as number of datasets ({})".format(n_dst+1, n_dst)
  with raises(Exception) as e:
    f_g_calculator.set_atomic_optimisation_mask(flex.bool(n_dst+1, True))
  assert msg == str(e.value)

  # setting amplitudes

  # should not error
  f_g_calculator.set_current_amplitudes(flex.double(n_grp*n_dst+n_atm))
  # should error
  msg = "Input array (size {}) must be the same length as current_amplitudes (size {})".format(n_grp*n_dst+n_atm-1, n_grp*n_dst+n_atm)
  with raises(Exception) as e:
    f_g_calculator.set_current_amplitudes(flex.double(n_grp*n_dst+n_atm-1))
  assert msg == str(e.value)

  print('OK')

if __name__ == "__main__":
  tst_optimisation_weights()
  tst_optimise_amplitudes_single_group()
  tst_optimise_amplitudes_multiple_groups_with_atomic()
  tst_optimise_amplitudes_multiple_groups_permuted_dataset_order()
  tst_functional_gradient_calculator()
  tst_functional_gradient_calculator_permuted_dataset_order()
  tst_functional_gradient_calculator_invalid_arguments()
