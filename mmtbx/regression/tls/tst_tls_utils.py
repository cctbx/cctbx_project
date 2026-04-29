from __future__ import absolute_import, division, print_function

import math
import numpy
from mmtbx.tls.utils import *
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal, not_approx_equal
from libtbx.test_utils import raises
from six.moves import zip
from six.moves import range

rran = numpy.random.random
iran = numpy.random.choice

numpy.random.seed(123232)

TLSMatrices.set_precision(12)
TLSAmplitudes.set_precision(12)

# Maximum tolerable rounding error
MAT_RND_TOL = (0.51)*(10**(-TLSMatrices.get_precision()))
AMP_RND_TOL = (0.51)*(10**(-TLSAmplitudes.get_precision()))
# For errors in the last decimal place
LAST_DP_TOL = 10.0 * max(10**(-TLSMatrices.get_precision()),
              10**(-TLSAmplitudes.get_precision()))

TEST_LENGTHS = [1,2,3]
TEST_TARGETS = [0.25,1.0,3.3]

TLSMATRICES = [
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

# Matrices that are valid if tolerance > 1e-6
LS = (180. / math.pi)**2
INVALID_TLS_MATRICES = [
    [1.,1.,-1e-6,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
    [0.,0.,0.,0.,0.,0.,LS,LS,-LS*1e-6,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
  ]

def get_TLS_values():
  vals = rran(21)
  obj = TLSMatrices(values=vals)
  assert approx_equal(obj.get("TLS"), vals, MAT_RND_TOL)
  return obj, vals

def get_TLS_amplitudes(n):
  vals = rran(n)
  obj = TLSAmplitudes(values=vals)
  assert approx_equal(obj.values, vals, AMP_RND_TOL)
  return obj, vals

def tst_set_precision():
  for cl in [TLSMatrices, TLSAmplitudes]:
    initial = cl.get_precision()
    for i in range(5):
      precision = i
      cl.set_precision(precision)
      assert precision == cl.get_precision()
    # Check is class method
    precision = 2
    a = cl(list(range(21)))
    b = cl(list(range(21)))
    a.set_precision(precision)
    assert b.get_precision() == precision
    # Set back to original
    cl.set_precision(initial)
    assert cl.get_precision() == initial

  print('OK')

def tst_set_tolerance():
  for cl in [TLSMatrices, TLSAmplitudes]:
    initial = cl.get_tolerance()
    for i in range(6):
      tolerance = 10.0**(-i)
      cl.set_tolerance(tolerance)
      assert tolerance == cl.get_tolerance()
    # Check is class method
    tolerance = 0.1
    a = cl(list(range(21)))
    b = cl(list(range(21)))
    a.set_tolerance(tolerance)
    assert b.get_tolerance() == tolerance
    # Set back to original
    cl.set_tolerance(initial)
    assert cl.get_tolerance() == initial

  print('OK')

def tst_TLSMatrices():
  """Check classes operating as expected (copying where should be copying rather than referencing...)"""

  for _ in range(5):

    # A
    a, a_vals = get_TLS_values()

    # A - with zero trace of S
    a2_vals = a_vals.copy()
    a2_vals[20] = -(
        round(a2_vals[12], a.get_precision()) +
        round(a2_vals[16], a.get_precision())
        )

    # B
    b, b_vals = get_TLS_values()

    # B2 - Test initialisation from matrices
    b2 = TLSMatrices(T=b_vals[0:6], L=b_vals[6:12], S=b_vals[12:21])

    # Addition of another instance
    c = a + b
    c_vals = (a_vals.round(a.get_precision()) + b_vals.round(b.get_precision()))

    # Multiplication by scalar
    d = (5.0*a)
    d_vals = (5.0*a_vals.round(a.get_precision()))

    # Check vals stored and extracted properly
    assert approx_equal(a.T+a.L+a.S, a_vals, MAT_RND_TOL) # matrices are returned as tuples
    assert approx_equal(a.get("TLS"), a_vals, MAT_RND_TOL)
    assert approx_equal(a.get("T"), a_vals[0:6], MAT_RND_TOL)
    assert approx_equal(a.get("L"), a_vals[6:12], MAT_RND_TOL)
    assert approx_equal(a.get("S"), a_vals[12:21], MAT_RND_TOL)
    assert approx_equal(b.get("TLS"), b_vals, MAT_RND_TOL)
    assert approx_equal(b.get("T"), b_vals[0:6], MAT_RND_TOL)
    assert approx_equal(b.get("L"), b_vals[6:12], MAT_RND_TOL)
    assert approx_equal(b.get("S"), b_vals[12:21], MAT_RND_TOL)
    # Check initialisation from matrices
    assert approx_equal(b2.get("TLS"), b_vals, MAT_RND_TOL)
    # Checks addition (needs greater tolerance)
    assert approx_equal(c.get("TLS"), c_vals, MAT_RND_TOL)
    assert approx_equal(c.get("T"), c_vals[0:6], MAT_RND_TOL)
    assert approx_equal(c.get("L"), c_vals[6:12], MAT_RND_TOL)
    assert approx_equal(c.get("S"), c_vals[12:21], MAT_RND_TOL)
    # Checks multiplication by scalar (need greater tolerance)
    assert approx_equal(d.get("TLS"), d_vals, MAT_RND_TOL)
    assert approx_equal(d.get("T"), d_vals[0:6], MAT_RND_TOL)
    assert approx_equal(d.get("L"), d_vals[6:12], MAT_RND_TOL)
    assert approx_equal(d.get("S"), d_vals[12:21], MAT_RND_TOL)

    # Check values are not being added/multiplied in place
    for _ in range(5):
      # Check addition
      assert approx_equal((a+a).get(), (a+a).get(), 0.0)
      # Check left- and right-multiplication
      assert approx_equal((5.0*a).get(), (5.0*a).get(), 0.0)
      assert approx_equal((a*5.0).get(), (a*5.0).get(), 0.0)
      # Check commutativity
      assert approx_equal((5.0*a).get(), (a*5.0).get(), 0.0)
    assert approx_equal(a.get("TLS"), a_vals, MAT_RND_TOL)

    # Check copies - reset
    a_copy = a.copy()
    a_copy.reset()
    assert approx_equal(a.get(), a_vals, MAT_RND_TOL)
    # Check copies - set
    a_copy = a.copy()
    a_copy.set(list(range(21)))
    assert approx_equal(a.get(), a_vals, MAT_RND_TOL)
    # Check initialisation from other does not affect original
    a_copy = TLSMatrices(a)
    assert approx_equal(a_copy.get(), a_vals, MAT_RND_TOL)
    a_copy.reset()
    assert approx_equal(a.get(), a_vals, MAT_RND_TOL)

    # Check addition in place
    a.add(b)
    assert approx_equal(a.get(), c_vals, MAT_RND_TOL)
    assert approx_equal(b.get(), b_vals, MAT_RND_TOL) # Check b unchanged
    # Check multiplication in place
    a.multiply(10.0)
    assert approx_equal(a.get(), 10.0*c_vals, MAT_RND_TOL)
    # set back to original values
    a.set(a_vals, "TLS")
    assert approx_equal(a.get(), a_vals, MAT_RND_TOL)

    # Order of the letters should not matter (always returns T-L-S)
    assert a.get("TL") == a.get("LT")
    assert a.get("LS") == a.get("SL")
    assert a.get("TS") == a.get("ST")
    assert a.get("TLS") == a.get("STL")
    assert a.get("TLS") == a.get("LST")
    assert a.get("TLS") == a.get("TSL")
    assert a.get("TLS") == a.get("LTS")
    assert a.get("TLS") == a.get("SLT")

    # Check values can be replaced correctly
    assert approx_equal(a.get(), a_vals, MAT_RND_TOL) # Check values are correct at start
    # T
    new_t = rran(6).tolist()
    a.set(values=new_t, component_string="T")
    assert approx_equal(a.get("T"), new_t, MAT_RND_TOL)
    # L
    new_l = rran(6).tolist()
    a.set(values=new_l, component_string="L")
    assert approx_equal(a.get("L"), new_l, MAT_RND_TOL)
    # S but not Szz
    new_s = rran(8).tolist()
    a.set(new_s, "S", include_szz=False)
    new_s_all = new_s + [-(round(new_s[0], a.get_precision()) + round(new_s[4], a.get_precision()))]
    assert len(new_s_all) == 9
    assert approx_equal(a.get("S"), new_s_all, MAT_RND_TOL)
    assert approx_equal(a.get("S")[8], new_s_all[8]) # This number is pre-rounded so should be "exact""
    # S
    new_s = rran(9).tolist()
    a.set(new_s, "S", include_szz=True)
    assert approx_equal(a.get("S"), new_s, MAT_RND_TOL)
    # all
    assert approx_equal(a.get(), new_t+new_l+new_s, MAT_RND_TOL)
    assert approx_equal(a.get("TLS"), new_t+new_l+new_s, MAT_RND_TOL)
    # set all values except Szz
    a.set(a_vals[:-1], "TLS", include_szz=False)
    assert approx_equal(a.get(), a2_vals, MAT_RND_TOL)
    assert approx_equal(a.get()[20], a2_vals[20]) # This number is pre-rounded in a2_vals
    # set back to original values
    a.set(a_vals, "TLS", include_szz=True)
    assert approx_equal(a.get(), a_vals, MAT_RND_TOL)

  print('OK')

def tst_TLSMatrices_counts():
  """Test that number of non-zero parameters are identifed correctly, etc"""

  letter_combinations = ["T","L","S","TL","LS","TS","TLS"]

  for _ in range(5):

    # Get values & object
    a, a_vals = get_TLS_values()
    assert (a_vals!=0.0).all()
    # Test reset function
    a.reset()
    assert approx_equal(a.get(), [0.0]*21, 0.0)

    # Test setting values to zero
    for letts in [""]+letter_combinations:
      # Get a fresh set of matrices each time
      a, a_vals = get_TLS_values()
      assert (a_vals!=0.0).all()
      # Set selected matrices to zero
      for l in letts:
        v = a.get(l)
        a.set([0.0]*len(v), l)
      # Test ANY
      for t_l in letter_combinations:
        non_zero_l = ''.join(set(t_l).difference(letts))
        ret = a.any(t_l)
        assert ret == bool(non_zero_l)
        # High tolerance -- always false
        assert not a.any(t_l, 1e3)
      # Test NPARAMS (number of non-zero values)
      if letts == "":
        n_non_zero = 21
      else:
        n_non_zero = 21 - len(a.get(letts))
      assert a.n_params(free=False, non_zero=True) == n_non_zero
      assert a.n_params(free=True, non_zero=True) == n_non_zero - int("S" not in letts)

    # Test setting values to non-zero
    for letts in [""]+letter_combinations:
      # Get a fresh set of matrices each time
      a = TLSMatrices()

      # Set random values that should be rounded to zero
      for l in letts:
        v = list(a.get(l))
        i = numpy.random.randint(len(v))
        v[i] = 0.49*(10.**(-a.get_precision()))
        a.set(v, l)
      # Test ANY
      for t_l in letter_combinations:
        ret = a.any(t_l)
        assert ret == False # ALL values should have rounded to zero
        # High tolerance -- always false
        assert not a.any(t_l, 1e3)
      # Test NPARAMS
      assert a.n_params(free=False, non_zero=True) == 0
      assert a.n_params(free=True, non_zero=True) == 0
      assert a.n_params(free=False, non_zero=False) == 21
      assert a.n_params(free=True, non_zero=False) == 20

      # Set selected matrices to non-zero
      for l in letts:
        v = list(a.get(l))
        i = numpy.random.randint(len(v))
        # Set values just abovd the rounding/tolerance limit
        v[i] = 2.0*max(10.**(-a.get_precision()), a.get_tolerance())
        a.set(v, l)
      # Test ANY
      for t_l in letter_combinations:
        non_zero_l = ''.join(set(t_l).intersection(letts))
        ret = a.any(t_l)
        assert ret == bool(non_zero_l)
        # High tolerance -- always false
        assert not a.any(t_l, 1e3)
      # Test NPARAMS (number of non-zero values)
      if letts == "":
        n_non_zero = 0
      else:
        n_non_zero = len(a.get(letts))
      assert a.n_params(free=False, non_zero=True) == n_non_zero
      assert a.n_params(free=True, non_zero=True) == n_non_zero - int("S" in letts)
      assert a.n_params(free=False, non_zero=False) == 21
      assert a.n_params(free=True, non_zero=False) == 20

  print('OK')

def tst_TLSMatrices_fails():
  """Check that exceptions are raised where expected"""

  # Tolerances & Precision
  with raises(Exception) as e:
    TLSMatrices.set_precision(2.3) # must be int
  with raises(Exception) as e:
    TLSMatrices.set_tolerance(-0.1) # must be greater than 0.0
  with raises(Exception) as e:
    TLSMatrices.set_tolerance(0.0) # must be greater than 0.0

  # Check doesn't error with the correct length
  a = TLSMatrices(T=list(range(6)), L=list(range(6)), S=list(range(9)))

  # Check length of input matrices
  for tt,ll,ss in [
      (5, 6, 9),
      (7, 6, 9),
      (6, 5, 9),
      (6, 7, 9),
      (6, 6, 8),
      (6, 6, 10),
      ]:
    with raises(Exception) as e:
      a = TLSMatrices(rran(tt),rran(ll),rran(ss))
    assert "Python argument types in" in str(e.value)

  # Check setting values with incorrect size array
  msg = "Mismatch between the length of the selected matrices and the length of the input array"
  a = TLSMatrices()
  for l, sel in [
      (6,  "T"),
      (6,  "L"),
      (9,  "S"),
      (12, "TL"),
      (15, "TS"),
      (15, "LS"),
      (21, "TLS"),
      ]:
    a.set(rran(l), sel)
    with raises(Exception) as e:
      a.set(rran(l-1), sel)
    assert msg == str(e.value)
    with raises(Exception) as e:
      a.set(rran(l+1), sel)
    assert msg == str(e.value)
    if "S" in sel:
      a.set(rran(l-1), sel, include_szz=False)
      with raises(Exception) as e:
        a.set(rran(l-2), sel, include_szz=False)
      assert msg == str(e.value)
      with raises(Exception) as e:
        a.set(rran(l), sel, include_szz=False)
      assert msg == str(e.value)

  # Check length of input array
  for i in [20,22]:
    with raises(Exception) as e:
      a = TLSMatrices(list(range(i)))
    assert "Input values must have length 21" == str(e.value)

  # Invalid letters in the selection strings
  a = TLSMatrices(TLSMATRICES[0])
  for s in ["", "F", "TLSD", "2T"]:
    with raises(ValueError) as e:
      a.any(s)
    if not s:
      assert "Empty string provided: '{}'".format(s) == str(e.value)
    else:
      assert "Invalid letters in string (not T, L or S): '{}'".format(s) == str(e.value)

    with raises(ValueError) as e:
      a.get(s)
    if not s:
      assert "Empty string provided: '{}'".format(s) == str(e.value)
    else:
      assert "Invalid letters in string (not T, L or S): '{}'".format(s) == str(e.value)

  # Multipliers must be positive
  with raises(Exception) as e:
    a.multiply(-0.1) # must be greater than 0.0
  assert "Multiplier must be positive" == str(e.value)

  a.multiply(0.0) # can be equal to zero!

  # Negative tolerances!
  msg = "Tolerance provided must either be positive or -1"
  # Allowed to be -1
  a.any("TLS", -1)
  # Any
  with raises(Exception) as e:
    a.any("TLS", -0.1)
  assert msg == str(e.value)
  # Decompose
  with raises(Exception) as e:
    a.decompose(-0.1)
  assert msg == str(e.value)
  # Valid
  with raises(Exception) as e:
    a.is_valid(-0.1)
  assert msg == str(e.value)
  # Normalise
  with raises(Exception) as e:
    a.normalise([(1,2,3)],(0,0,0),1.0,-0.1)
  assert msg == str(e.value)

  # Negative/zero-value target values!
  msg = "target must be positive"
  with raises(Exception) as e:
    a.normalise([(1,2,3)],(0,0,0),0.0)
  assert msg == str(e.value)
  with raises(Exception) as e:
    a.normalise([(1,2,3)],(0,0,0),-1.0)
  assert msg == str(e.value)

  print('OK')

def tst_TLSMatrices_uijs():
  """Check that uijs are generated from coordinates as expected"""

  from mmtbx.tls import tlso, uaniso_from_tls_one_group
  from scitbx.linalg.eigensystem import real_symmetric

  for vals in TLSMATRICES:

    # Initialise
    a = TLSMatrices(vals)
    # Test sets of matrices are valid
    assert a.is_valid()

    # Get a set of random coordinates
    n_atm = 5
    coords = rran(n_atm*3).reshape(n_atm,3)
    origin = rran(3)

    # Reference Uijs
    tls_obj = tlso(t=a.T, l=a.L, s=a.S, origin=(0.0,0.0,0.0)) # Set origin to zero and subtract from coords as manual check
    uij_ref = uaniso_from_tls_one_group(tls_obj, coords-origin, zeroize_trace=False)
    uij_ref_np = numpy.array(uij_ref)

    # Uijs from the class being tested
    uij_test = a.uijs(coords, origin)
    uij_test_np = numpy.array(uij_test)

    # Compare (should be identical because class calls the uaniso function)
    assert approx_equal(uij_ref_np.flatten(), uij_test_np.flatten(), 0.0)

    # Calculate average eigenvalue of output uijs
    # TODO: verify that real_symmetric returns non-dict, I assume its an eigensystem-y instance
    eigs = [real_symmetric(u).values() for u in uij_test]
    orig_mean_eig = numpy.mean(eigs)

    for target in TEST_TARGETS:
      # Extract copy to test normalisation
      b = a.copy()
      # Normalise
      mult = b.normalise(coords, origin, target=target)
      assert approx_equal(orig_mean_eig, target*mult, LAST_DP_TOL)
      # Extract normalised uijs
      uij_norm = b.uijs(coords, origin)
      uij_norm_np = numpy.array(uij_norm)
      # Check average eigenvalue
      # TODO: see above todo
      eigs = [real_symmetric(u).values() for u in uij_norm]
      mean_eig_norm = numpy.mean(eigs)
      assert approx_equal(mean_eig_norm, target, LAST_DP_TOL)
      # Check that normalised values related to original values by multiplier
      uij_comp_np = mult*uij_norm_np
      assert approx_equal(uij_comp_np.flatten(), uij_test_np.flatten(), LAST_DP_TOL)

  print('OK')

def tst_TLSAmplitudes():
  """Exercise the TLSAmplitudes class"""

  for length in TEST_LENGTHS:

    # Initialise from array
    a, a_vals = get_TLS_amplitudes(length)

    # Initialise from length
    b = TLSAmplitudes(length)
    assert b.size() == length
    # Check values intialised to 1
    assert approx_equal(b.values, [1.0]*length, 0.0)
    # Set values using set function
    b_vals = rran(length)
    b.set(b_vals)
    assert approx_equal(b.values, b_vals, AMP_RND_TOL)

    # Addition of another instance
    c = a + b
    c_vals = (a_vals.round(a.get_precision()) + b_vals.round(b.get_precision()))

    # Multiplication by scalar
    d = (5.0*a)
    d_vals = (5.0*a_vals.round(a.get_precision()))

    # Intialise from other
    e = TLSAmplitudes(a)

    # Check initial values
    assert approx_equal(a.values, a_vals, AMP_RND_TOL)
    assert approx_equal(b.values, b_vals, AMP_RND_TOL)
    assert approx_equal(c.values, c_vals, AMP_RND_TOL)
    assert approx_equal(d.values, d_vals, AMP_RND_TOL)
    assert approx_equal(e.values, a_vals, AMP_RND_TOL)
    assert approx_equal(e.values, a.values, 0.0)

    assert a.any()
    assert b.any()
    assert c.any()
    assert d.any()
    assert e.any()

    assert not a.any(1e3)
    assert not b.any(1e3)
    assert not c.any(1e3)
    assert not d.any(1e3)
    assert not e.any(1e3)

    assert a.size() == length
    assert b.size() == length
    assert c.size() == length
    assert d.size() == length
    assert e.size() == length

    # Check get
    assert approx_equal(a.values, a.get(list(range(a.size()))), 0.0)
    assert approx_equal(b.values, b.get(list(range(b.size()))), 0.0)
    assert approx_equal(c.values, c.get(list(range(c.size()))), 0.0)
    assert approx_equal(d.values, d.get(list(range(d.size()))), 0.0)
    assert approx_equal(e.values, e.get(list(range(e.size()))), 0.0)

    # Check reset works
    e.reset()
    assert approx_equal(e.values, [1.0]*e.size(), 0.0)
    assert e.any()
    # Check that a is unchanged by e.reset())
    assert approx_equal(a.values, a_vals, AMP_RND_TOL)

    # Check values are not being added/multiplied in place
    for _ in range(5):
      # Check addition
      assert approx_equal((a+a).values, (a+a).values, 0.0)
      # Check left- and right-multiplication
      assert approx_equal((5.0*a).values, (5.0*a).values, 0.0)
      assert approx_equal((a*5.0).values, (a*5.0).values, 0.0)
      # Check commutativity
      assert approx_equal((5.0*a).values, (a*5.0).values, 0.0)
    # Final check against initial values
    assert approx_equal(a.values, a_vals, AMP_RND_TOL)

    # Check indexing
    for i in range(a.size()):
      assert approx_equal(a[i], a_vals[i], AMP_RND_TOL)

    # Check addition in place
    a.add(b)
    assert approx_equal(a.values, c_vals, AMP_RND_TOL)
    assert approx_equal(b.values, b_vals, AMP_RND_TOL) # Check b unchanged
    # Check multiplication in place
    mult = 10.0
    a.multiply(mult)
    assert approx_equal(a.values, mult*c_vals, AMP_RND_TOL)
    # set back to original values
    a.set(a_vals)
    assert approx_equal(a.values, a_vals, AMP_RND_TOL)
    assert approx_equal(b.values, b_vals, AMP_RND_TOL)

    # Test setting subset of values
    selection = iran(a.size(), size=min(3,a.size()), replace=False)
    new_vals = rran(len(selection)).tolist()
    chg_a_vals = a_vals.copy()
    assert approx_equal(a_vals, chg_a_vals, 0.0)
    for new_idx, chg_idx in enumerate(selection):
      assert chg_a_vals[chg_idx] != new_vals[new_idx] # Check that new values are different!
      chg_a_vals[chg_idx] = new_vals[new_idx]
    assert approx_equal(chg_a_vals[selection], new_vals, AMP_RND_TOL) # Check that values are set correctly
    a_chg = a.copy()
    assert approx_equal(a_chg.values, a_vals, AMP_RND_TOL)
    a_chg.set(new_vals, list(selection.astype(numpy.uint)))
    assert approx_equal(a_chg.values, chg_a_vals, AMP_RND_TOL)

    # Test reset and n_params
    for _ in range(5):

      # All non-zero
      assert a_vals.all()
      assert a.any()
      assert approx_equal(a.values, a_vals, AMP_RND_TOL)
      assert a.n_params(non_zero=True) == a.size()
      assert a.n_params(non_zero=False) == a.size()

      # Set all values to one
      a.reset()
      assert a.any()
      assert approx_equal(a.values, [1.0]*a.size(), 0.0)
      assert a.n_params(non_zero=True) == a.size()
      assert a.n_params(non_zero=False) == a.size()
      # Set some values to non-one
      selection = iran(a.size(), size=min(3,a.size()), replace=False).tolist()
      new_vals = rran(len(selection)).tolist()
      a.set(new_vals, selection)
      for i, v in enumerate(a.values):
        if i in selection:
          assert approx_equal(v, new_vals[selection.index(i)], AMP_RND_TOL)
        else:
          assert v == 1.0
      assert a.n_params(non_zero=True) == a.size() - new_vals.count(0.0)
      assert a.n_params(non_zero=False) == a.size()

      # Set all values to zero
      a.zero_values()
      assert not a.any()
      assert approx_equal(a.values, [0.0]*a.size(), 0.0)
      assert a.n_params(non_zero=True) == 0
      assert a.n_params(non_zero=False) == a.size()
      # Set some values to non-zero
      selection = iran(a.size(), size=min(3,a.size()), replace=False).tolist()
      new_vals = rran(len(selection)).tolist()
      a.set(new_vals, selection)
      for i, v in enumerate(a.values):
        if i in selection:
          assert approx_equal(v, new_vals[selection.index(i)], AMP_RND_TOL)
        else:
          assert v == 0.0
      assert a.n_params(non_zero=True) == len(new_vals) - new_vals.count(0.0)
      assert a.n_params(non_zero=False) == a.size()
      assert a.any()

      # Set all values to original
      a.set(a_vals.tolist(), list(range(a.size()))) # set all values by selection, just to mix things up
      assert approx_equal(a.values, a_vals, AMP_RND_TOL)
      assert a.n_params(non_zero=True) == a.size()
      assert a.n_params(non_zero=False) == a.size()
      # Set some values to zero
      new_vals = [0.0, 100.0, 0.0, 34.5, 1e-4, -1e-4][:a.size()] # make sure values not longer than a
      selection = iran(a.size(), size=len(new_vals), replace=False).tolist()
      assert len(selection) == len(new_vals)
      a.set(new_vals, selection)
      for i, v in enumerate(a.values):
        if i in selection:
          assert approx_equal(v, new_vals[selection.index(i)], AMP_RND_TOL)
        else:
          assert approx_equal(v, a_vals[i], AMP_RND_TOL)
      assert a.n_params(non_zero=True) == a.size() - numpy.round(new_vals, a.get_precision()).tolist().count(0.0)
      assert a.n_params(non_zero=False) == a.size()

      # Set all values to original
      a.set(a_vals.tolist(), list(range(a.size()))) # set all values by selection, just to mix things up
      assert approx_equal(a.values, a_vals, AMP_RND_TOL)
      assert a.n_params(non_zero=True) == a.size()
      assert a.n_params(non_zero=False) == a.size()
      # Set some values to negative
      selection = iran(a.size(), size=min(4,a.size()), replace=False).tolist()
      new_vals = [-1.0*v for v in a.get(selection)]
      a.set(new_vals, selection)
      for i, v in enumerate(a.values):
        if i in selection:
          assert approx_equal(v, -a_vals[i], AMP_RND_TOL)
        else:
          assert approx_equal(v, a_vals[i], AMP_RND_TOL)
      assert a.n_params(non_zero=True) == a.size()
      assert a.n_params(non_zero=False) == a.size()
      # Zero negative values and test they are zero
      a.zero_negative_values()
      for i, v in enumerate(a.values):
        if i in selection:
          assert v == 0.0
        else:
          assert approx_equal(v, a_vals[i], AMP_RND_TOL)
      assert a.n_params(non_zero=True) == a.size() - len(selection)
      assert a.n_params(non_zero=False) == a.size()

      # Reset for next loop, etc
      a.set(a_vals.tolist(), list(range(a.size()))) # set all values by selection, just to mix things up
      assert approx_equal(a.values, a_vals, AMP_RND_TOL)
      assert a.n_params(non_zero=True) == a.size()
      assert a.n_params(non_zero=False) == a.size()

  print('OK')

def tst_TLSAmplitudes_fails():
  """Check that exceptions are raised where expected"""

  # Tolerances & Precision
  with raises(Exception) as e:
    TLSAmplitudes.set_precision(2.3) # must be int
  with raises(Exception) as e:
    TLSAmplitudes.set_tolerance(-0.1) # must be greater than 0.0
  with raises(Exception) as e:
    TLSAmplitudes.set_tolerance(0.0) # must be greater than 0.0

  with raises(Exception) as e:
    TLSAmplitudes(0)

  a = TLSAmplitudes(10)
  b = TLSAmplitudes(11)

  # Adding incompatible lengths
  with raises(Exception) as e:
    a+b
  assert "TLSAmplitudes must have the same length" == str(e.value)
  with raises(Exception) as e:
    a.add(b)
  assert "TLSAmplitudes must have the same length" == str(e.value)

  # Negative tolerances!
  msg = "Tolerance provided must either be positive or -1"
  # Allowed to be -1
  a.any(-1)
  # Any
  with raises(Exception) as e:
    a.any(-0.1)
  assert msg == str(e.value)

  # Invalid selections!
  with raises(Exception) as e:
    a.get([])
  assert "No indices given for selection" == str(e.value)
  with raises(Exception) as e:
    a.get([a.size()])
  assert "Selection indices out of range of TLSAmplitudes" == str(e.value)
  with raises(Exception) as e:
    a.get(list(range(a.size()))+[0])
  assert "Selection indices cannot be longer than TLSAmplitudes" == str(e.value)

  # Valid, but odd, selections -- maybe change later
  a.get([1,1,1,1])

  # Invalid selections -- setting values
  with raises(Exception) as e:
    a.set(list(range(a.size()-1)))
  assert "Input array must be the same length as TLSAmplitudes" == str(e.value)
  with raises(Exception) as e:
    a.set(list(range(a.size()+1)))
  assert "Input array must be the same length as TLSAmplitudes" == str(e.value)
  # No selection
  with raises(Exception) as e:
    a.set([],[])
  assert "No indices given for selection" == str(e.value)
  # Negative indices!
  with raises(Exception) as e:
    a.set([],[])
  assert "No indices given for selection" == str(e.value)
  # Mismatching size
  for n in range(1,a.size()):
    with raises(Exception) as e:
      a.set(rran(n-1), iran(a.size(), size=n, replace=False).astype(numpy.uint))
    assert "Input values must be the same length as input selection" == str(e.value)
    with raises(Exception) as e:
      a.set(rran(n+1), iran(a.size(), size=n, replace=False).astype(numpy.uint))
    assert "Input values must be the same length as input selection" == str(e.value)
  # Negative indices
  with raises(OverflowError) as e:
    a.set([1], [-1])
  assert ("can't convert negative value to unsigned" == str(e.value) or
          "negative overflow" in str(e.value) or
          "can't convert negative value to unsigned" in str(e.value))

  # Negative/zero-value target values!
  msg = "target must be positive"
  with raises(Exception) as e:
    a.normalise(0.0)
  assert msg == str(e.value)
  with raises(Exception) as e:
    a.normalise(-1.0)
  assert msg == str(e.value)

  print('OK')

def tst_TLSMatricesAndAmplitudes():
  """Exercise the TLSMatricesAndAmplitudes class"""

  for length in TEST_LENGTHS:

    # Iterate through test matrices
    for a_mats in TLSMATRICES:
      a_amps = rran(length)
      # Initialise from existing classes (objects are passed by reference!)
      a_m = TLSMatrices(a_mats)
      a_a = TLSAmplitudes(a_amps)
      a = TLSMatricesAndAmplitudes(a_m, a_a)

      # Check values passed through
      assert approx_equal(a.matrices.get("TLS"), a_mats, MAT_RND_TOL)
      assert approx_equal(a.amplitudes.values     , a_amps, AMP_RND_TOL)
      # Double check that values are the same
      assert approx_equal(a_m.T, a.matrices.T, 0.0)
      assert approx_equal(a_m.L, a.matrices.L, 0.0)
      assert approx_equal(a_m.S, a.matrices.S, 0.0)
      assert approx_equal(a_a.values, a.amplitudes.values, 0.0)

      # Initialise from other
      b = TLSMatricesAndAmplitudes(a)
      assert approx_equal(b.matrices.get("TLS"), a.matrices.get("TLS"), 0.0)
      assert approx_equal(b.amplitudes.values  , a.amplitudes.values, 0.0)

      # Initialise from copy
      c = a.copy()
      assert approx_equal(c.matrices.get("TLS"), a.matrices.get("TLS"), 0.0)
      assert approx_equal(c.amplitudes.values  , a.amplitudes.values, 0.0)

      # Initialise from length
      d = TLSMatricesAndAmplitudes(length)

      # Initialise straight from values
      e = TLSMatricesAndAmplitudes(a_mats, a_amps)

      # Check that a is using objects by reference, and that b & c are initialised by values from a
      a.reset()
      assert approx_equal(a.matrices.get("TLS"), [0.0]*21, 0.0)
      assert approx_equal(a.amplitudes.values  , [1.0]*length, 0.0)
      # Original objects should also be affected (same object)
      assert approx_equal(a_m.get("TLS")     , [0.0]*21, 0.0)
      assert approx_equal(a_a.values          , [1.0]*length, 0.0)
      # Check b & c are not reset (copied by reference)
      assert approx_equal(b.matrices.get("TLS"), a_mats, MAT_RND_TOL)
      assert approx_equal(b.amplitudes.values  , a_amps, AMP_RND_TOL)
      assert approx_equal(c.matrices.get("TLS"), a_mats, MAT_RND_TOL)
      assert approx_equal(c.amplitudes.values  , a_amps, AMP_RND_TOL)
      # Reset b and check c unaffected
      b.reset()
      assert approx_equal(b.matrices.get("TLS"), [0.0]*21, 0.0)
      assert approx_equal(b.amplitudes.values  , [1.0]*length, 0.0)
      assert approx_equal(c.matrices.get("TLS"), a_mats, MAT_RND_TOL)
      assert approx_equal(c.amplitudes.values  , a_amps, AMP_RND_TOL)

      # Check d that defaults are intialised as expected
      assert approx_equal(d.matrices.get("TLS"), [0.0]*21, 0.0)
      assert approx_equal(d.amplitudes.values  , [1.0]*length, 0.0)

      # Check e is initialised correctly
      assert approx_equal(e.matrices.get("TLS"), a_mats, MAT_RND_TOL)
      assert approx_equal(e.amplitudes.values  , a_amps, AMP_RND_TOL)

  print('OK')

def tst_TLSMatricesAndAmplitudes_counts():
  """Test that number of non-zero parameters are identifed correctly, etc"""

  for length in TEST_LENGTHS:

    # Iterate through test matrices
    for a_mats in TLSMATRICES:
      a_amps = rran(length)
      a = TLSMatricesAndAmplitudes(a_mats, a_amps.tolist())
      assert approx_equal(a.matrices.get("TLS"), a_mats, MAT_RND_TOL)
      assert approx_equal(a.amplitudes.values     , a_amps, AMP_RND_TOL)

      # Check tolerances
      assert not a.is_null()
      assert not a.is_null(-1, -1)
      assert not a.is_null(0.0, 0.0)
      assert a.is_null(1e3, -1)
      assert a.is_null(-1, 1e3)
      assert a.is_null(1e3, 1e3)

      # Check null function + n_params
      assert not a.is_null()
      # All should be non-zero
      assert a.n_params(free=True, non_zero=True)   == 20 + length
      assert a.n_params(free=True, non_zero=False)  == 20 + length
      assert a.n_params(free=False, non_zero=True)  == 21 + length
      assert a.n_params(free=False, non_zero=False) == 21 + length
      # Set matrices to zero
      a.matrices.reset()
      assert a.is_null()
      # Matrices are null, amplitudes not
      assert a.n_params(free=True, non_zero=True)   == length
      assert a.n_params(free=True, non_zero=False)  == 20 + length
      assert a.n_params(free=False, non_zero=True)  == length
      assert a.n_params(free=False, non_zero=False) == 21 + length
      # Reset back to original values
      a.matrices.set(a_mats)
      assert not a.is_null()
      # All should be non-zero
      assert a.n_params(free=True, non_zero=True)   == 20 + length
      assert a.n_params(free=True, non_zero=False)  == 20 + length
      assert a.n_params(free=False, non_zero=True)  == 21 + length
      assert a.n_params(free=False, non_zero=False) == 21 + length
      # Set amplitudes to one -- NOT NULL
      a.amplitudes.reset()
      assert not a.is_null()
      # Amplitudes are all one
      assert a.n_params(free=True, non_zero=True)   == 20 + length
      assert a.n_params(free=True, non_zero=False)  == 20 + length
      assert a.n_params(free=False, non_zero=True)  == 21 + length
      assert a.n_params(free=False, non_zero=False) == 21 + length
      # Set amplitudes to zero -- NULL
      a.amplitudes.zero_values()
      assert a.is_null()
      # Amplitudes are all one
      assert a.n_params(free=True, non_zero=True)   == 20
      assert a.n_params(free=True, non_zero=False)  == 20 + length
      assert a.n_params(free=False, non_zero=True)  == 21
      assert a.n_params(free=False, non_zero=False) == 21 + length

      # Reset to original values
      a.amplitudes.set(a_amps)
      assert approx_equal(a.matrices.get("TLS"), a_mats, MAT_RND_TOL)
      assert approx_equal(a.amplitudes.values     , a_amps, AMP_RND_TOL)
      assert not a.is_null()

  print('OK')

def tst_TLSMatricesAndAmplitudes_fails():
  """Check that exceptions are raised where expected"""

  m=TLSMatrices()
  a=TLSAmplitudes(10)
  with raises(Exception) as e:
    TLSMatricesAndAmplitudes(m,m)
  with raises(Exception) as e:
    TLSMatricesAndAmplitudes(a,a)
  with raises(Exception) as e:
    TLSMatricesAndAmplitudes(a,m)

  with raises(Exception) as e:
    TLSMatricesAndAmplitudes(rran(20), rran(10))
  assert "Matrix values must have length 21" == str(e.value)

  with raises(Exception) as e:
    TLSMatricesAndAmplitudes(rran(20), rran(10))
  assert "Matrix values must have length 21" == str(e.value)

  with raises(Exception) as e:
    TLSMatricesAndAmplitudes(rran(10), rran(21))
  assert "Matrix values must have length 21" == str(e.value)

  with raises(Exception) as e:
    TLSMatricesAndAmplitudes(rran(21), [])
  assert "Amplitude values must have length greater than 0" == str(e.value)

  # For use throughout
  n_dst = 3
  n_atm = 2
  assert n_dst>2 # needed for later
  assert n_atm>1 # needed for later
  assert n_dst != n_atm # need different otherwise cannot test for swaps of n_dst and n_atm
  ma = TLSMatricesAndAmplitudes(n_dst)

  # Tolerances
  ma.is_null(0.0, 0.0) # values can be zero
  ma.is_null(-1, -1) # values can be -1
  msg = "Tolerance provided must either be positive or -1"
  with raises(Exception) as e:
    ma.is_null(-0.1, -1)
  assert msg == str(e.value)
  with raises(Exception) as e:
    ma.reset_if_null(-0.1, -1)
  assert msg == str(e.value)
  with raises(Exception) as e:
    ma.is_null(-1, -0.1)
  assert msg == str(e.value)
  with raises(Exception) as e:
    ma.reset_if_null(-1, -0.1)
  assert msg == str(e.value)

  # Coordinates sets for following tests
  sites = flex.vec3_double(rran(n_dst*n_atm*3).reshape(n_dst*n_atm, 3).tolist())
  sites.reshape(flex.grid(n_dst,n_atm))
  origins = rran(n_dst*3).reshape(n_dst,3)

  # Check this doesn't break it
  ma.copy().normalise_by_amplitudes(1.0)
  ma.copy().normalise_by_matrices(sites, origins, 1.0)
  # Normalise (target values)
  msg = "target must be positive"
  with raises(Exception) as e:
    ma.normalise_by_amplitudes(-0.1)
  assert msg == str(e.value)
  with raises(Exception) as e:
    ma.normalise_by_matrices(sites, origins, -0.1)
  assert msg == str(e.value)

  # Compatibility of sites/origins arrays (normalise and uijs)
  # Array has swapped axes
  msg = "Mismatch between the size of origins and first dimension of sites_carts"
  new_sites = flex.vec3_double(rran(n_dst*n_atm*3).reshape(n_dst*n_atm, 3).tolist())
  new_sites.reshape(flex.grid(n_atm,n_dst)) # swap axes
  with raises(Exception) as e:
    ma.normalise_by_matrices(new_sites, origins)
  assert msg == str(e.value)
  with raises(Exception) as e:
    ma.uijs(new_sites, origins)
  assert msg == str(e.value)
  # sites is too short/long
  msg = "Mismatch between the size of origins and first dimension of sites_carts"
  for n_tmp in [n_dst-1, n_dst+1]:
    new_sites = flex.vec3_double(rran(n_tmp*n_atm*3).reshape(n_tmp*n_atm, 3).tolist())
    new_sites.reshape(flex.grid(n_tmp,n_atm))
    with raises(Exception) as e:
      ma.normalise_by_matrices(new_sites, origins)
    assert msg == str(e.value)
    with raises(Exception) as e:
      ma.uijs(new_sites, origins)
    assert msg == str(e.value)
    # Sites/origins compatible but not same length as amplitudes
    new_origins = rran(n_tmp*3).reshape(n_tmp, 3).tolist()
    ma.copy().normalise_by_matrices(new_sites, new_origins) # Should not error -- does not use amplitudes
    with raises(Exception) as e:
      ma.uijs(new_sites, new_origins)
    assert "Mismatch between the size of TLSAmplitudes and the input arrays" == str(e.value)

  # Check dimension of sites_carts
  msg = "sites_carts must be 2-dimensional array of size (n_dst, n_atm)"
  # Make 3-d array of sites
  n_fake = 2 # length of new dimension
  new_sites = flex.vec3_double(rran(n_fake*n_dst*n_atm*3).reshape(n_fake*n_dst*n_atm, 3).tolist())
  new_sites.reshape(flex.grid(n_fake,n_atm,n_dst))
  with raises(Exception) as e:
    ma.normalise_by_matrices(new_sites, origins)
  assert msg == str(e.value)
  with raises(Exception) as e:
    ma.uijs(new_sites, origins)
  assert msg == str(e.value)
  # Make 1-d array of sites
  new_sites = flex.vec3_double(rran(n_atm*3).reshape(n_atm, 3).tolist())
  new_sites.reshape(flex.grid(n_atm))
  with raises(Exception) as e:
    ma.normalise_by_matrices(new_sites, origins)
  assert msg == str(e.value)
  with raises(Exception) as e:
    ma.uijs(new_sites, origins)
  assert msg == str(e.value)
  # Make other 1-d array
  new_sites = flex.vec3_double(rran(n_atm*3).reshape(n_atm, 3))
  with raises(Exception) as e:
    ma.normalise_by_matrices(new_sites, origins)
  assert msg == str(e.value)
  with raises(Exception) as e:
    ma.uijs(new_sites, origins)
  assert msg == str(e.value)

  # Check error when origins is 2-dimensional
  n_fake = 2 # length of new dimension
  new_origins = flex.vec3_double(rran(n_fake*n_atm*3).reshape(n_fake*n_atm, 3).tolist())
  new_origins.reshape(flex.grid(n_fake,n_atm))
  with raises(Exception) as e:
    ma.normalise_by_matrices(sites, new_origins)
  assert "Python argument types in" in str(e.value)

  # Make new sites/origins for subset of ampls for use below
  t_n_dst = n_dst - 1
  new_sites = flex.vec3_double(rran(t_n_dst*n_atm*3).reshape(t_n_dst*n_atm, 3).tolist())
  new_sites.reshape(flex.grid(t_n_dst,n_atm))
  new_origins = rran(t_n_dst*3).reshape(t_n_dst,3)

  # Selection compatibility with sites and origins
  msg = "Mismatch between the size of selection and the input arrays"
  for l in [t_n_dst-1, t_n_dst+1]:
    # Selection too short
    sel = iran(n_dst, size=l, replace=False)
    with raises(Exception) as e:
      ma.uijs(new_sites, new_origins, sel.astype(numpy.uint))
    assert msg == str(e.value)
  # Invalid selection
  sel = iran(n_dst, size=t_n_dst, replace=False)
  sel[-1] = n_dst # larger than allowed
  with raises(Exception) as e:
    ma.uijs(new_sites, new_origins, sel.astype(numpy.uint))
  assert "Selection indices out of range of TLSAmplitudes" == str(e.value)

  print('OK')

def tst_TLSMatricesAndAmplitudes_valid():
  """Check that whole object is valid/invalid based on components"""

  length = 10
  indices = [3,5]
  tol = 1e-6

  for fail_mat in INVALID_TLS_MATRICES:
    a_amps = rran(length)
    assert (a_amps<1.0).all()
    a = TLSMatricesAndAmplitudes(fail_mat, a_amps.tolist())
    for i in range(length):
      assert a.expand()[i].is_valid(tol)
      assert a.is_valid(flex.size_t([i]), tol)
    assert a.is_valid(tol)
    assert a.is_valid(flex.size_t(indices), tol)
    a.amplitudes.set([12., 1.1], flex.size_t(indices))
    assert not a.is_valid(tol)
    assert not a.is_valid(flex.size_t(indices), tol)
    assert a.is_valid(flex.size_t([i for i in range(length) if i not in indices]), tol)
    for i in range(length):
      if i in indices:
        assert not a.expand()[i].is_valid(tol)
        assert not a.is_valid(flex.size_t([i]), tol)
      else:
        assert a.expand()[i].is_valid(tol)
        assert a.is_valid(flex.size_t([i]), tol)

  print('OK')

def tst_TLSMatricesAndAmplitudes_uijs():
  """Test uijs generated correctly"""

  for length in TEST_LENGTHS:

    # Iterate through test matrices
    for a_mats in TLSMATRICES:
      a_amps = rran(length)
      a = TLSMatricesAndAmplitudes(a_mats, a_amps.tolist())
      assert approx_equal(a.matrices.get("TLS"), a_mats, MAT_RND_TOL)
      assert approx_equal(a.amplitudes.values, a_amps, AMP_RND_TOL)

      # Get a set of random coordinates
      n_atm = 5
      coords = rran(length*n_atm*3).reshape(length,n_atm,3)
      origns = rran(length*3).reshape(length,3)

      # Format coordinates into flex array
      sites = flex.vec3_double(coords.reshape(length*n_atm, 3).tolist())
      sites.reshape(flex.grid(length,n_atm))
      # Uijs from the testing class
      all_uijs = a.uijs(sites, origns)
      # Reformat to numpy array for slicing
      all_uijs_np = numpy.array(all_uijs).reshape(length,n_atm,6)

      # Check expand
      exp = a.expand() # these also used later
      for i, m in enumerate(exp):

        # Check expanded matrices have correct values
        exp_amp = round(a_amps[i], a.amplitudes.get_precision())
        exp_mat = numpy.round(a_mats, a.matrices.get_precision())
        assert approx_equal(m.get("TLS"), exp_amp*exp_mat, LAST_DP_TOL)

        # Check that uij from this class equals that expected from the container class
        this_uij = m.uijs(coords[i].tolist(), origns[i].tolist())
        assert approx_equal(numpy.array(this_uij).flatten(), all_uijs_np[i].flatten())

      # Reshape to 1d and round for convenience
      all_uijs_np = all_uijs_np.reshape(length*n_atm,6)

      # Test normalise by amplitudes
      for target in TEST_TARGETS:

        # Create copy to operate on
        new_a = a.copy()

        # Apply normalisation
        mult = new_a.normalise_by_amplitudes(target=target)

        # Average amplitude should now be target
        mean_amp = numpy.mean(new_a.amplitudes.values)
        assert approx_equal(mean_amp, target, LAST_DP_TOL)

        # Check expanded matrices are the same
        new_exp = new_a.expand()
        for new, orig in zip(new_exp, exp):
          assert approx_equal(new.get("TLS"), orig.get("TLS"))

      # Test normalise by matrices
      from scitbx.linalg.eigensystem import real_symmetric
      for target in TEST_TARGETS:

        # Create copy to operate on
        new_a = a.copy()

        # Apply normalisation
        new_a.normalise_by_matrices(sites, origns, target)

        # Average uij eigenvalue from Matrices object should now be target
        uijs = [new_a.matrices.uijs(coords[i], origns[i]) for i in range(length)]
        # TODO: see above todo
        eigenvalues = [[list(real_symmetric(u).values()) for u in uijs[i]] for i in range(length)]
        mean_eig = numpy.mean(eigenvalues)
        assert approx_equal(mean_eig, target, LAST_DP_TOL)

        # Check expanded matrices are the same
        new_exp = new_a.expand()
        for new, orig in zip(new_exp, exp):
          assert approx_equal(new.get("TLS"), orig.get("TLS"))

        # Output uijs should still be the same
        new_uijs = new_a.uijs(sites, origns)
        new_uijs_np = numpy.array(new_uijs)
        assert approx_equal(new_uijs_np.flatten(), all_uijs_np.flatten(), LAST_DP_TOL)

  print('OK')

def tst_TLSMatricesAndAmplitudesList():
  """Exercise the TLSMatricesAndAmplitudesList class"""

  for length in [1,2,3]:
    for n_dst in [1,2,3]:
      mal = TLSMatricesAndAmplitudesList(length=length, n_amplitudes=n_dst)
      assert mal.is_null()
      assert mal.size() == length
      assert mal.n_params(free=False, non_zero=False) == (21 + n_dst) * length
      assert mal.n_params(free=True, non_zero=False) == (20 + n_dst) * length
      assert mal.n_params(free=False, non_zero=True) == (n_dst) * length # initialised with all matrices zero
      assert mal.n_params(free=True, non_zero=True) == (n_dst) * length # initialised with all matrices zero
      for i, ma in enumerate(mal):
        assert ma.amplitudes.size() == n_dst
        ma.matrices.set(rran(21))
        assert not mal.is_null()
        ma.amplitudes.set(rran(n_dst))
        assert not mal.is_null()
        for j in range(length):
          if i==j: continue
          ma2 = mal.get(j)
          assert not_approx_equal(list(ma.amplitudes.values), list(ma2.amplitudes.values))
          assert not_approx_equal(list(ma.matrices.get("TLS")), list(ma2.matrices.get("TLS")))
        assert mal.n_params(free=False, non_zero=False) == (21 + n_dst) * length
        assert mal.n_params(free=True, non_zero=False) == (20 + n_dst) * length
        assert mal.n_params(free=False, non_zero=True) == 21*(i+1) + (n_dst) * length
        assert mal.n_params(free=True, non_zero=True) == 20*(i+1) + (n_dst) * length
      assert mal.n_params(free=False, non_zero=False) == (21 + n_dst) * length
      assert mal.n_params(free=True, non_zero=False) == (20 + n_dst) * length
      assert mal.n_params(free=False, non_zero=True) == (21 + n_dst) * length
      assert mal.n_params(free=True, non_zero=True) == (20 + n_dst) * length

      # Reset all matrices
      mal.reset_matrices()
      assert mal.is_null()
      for ma in mal:
        assert not ma.matrices.any("TLS")
        assert ma.amplitudes.any()
        assert ma.is_null()
      assert mal.n_params(free=False, non_zero=False) == (21 + n_dst) * length
      assert mal.n_params(free=True, non_zero=False) == (20 + n_dst) * length
      assert mal.n_params(free=False, non_zero=True) == (n_dst) * length
      assert mal.n_params(free=True, non_zero=True) == (n_dst) * length
      # This should mean all are null
      mal.reset_null_modes()
      assert mal.is_null()
      for ma in mal:
        assert not ma.matrices.any("TLS")
        assert approx_equal(ma.amplitudes.values, [1.0]*n_dst)
        assert ma.is_null()
      assert mal.n_params(free=False, non_zero=False) == (21 + n_dst) * length
      assert mal.n_params(free=True, non_zero=False) == (20 + n_dst) * length
      assert mal.n_params(free=False, non_zero=True) == (n_dst) * length
      assert mal.n_params(free=True, non_zero=True) == (n_dst) * length
      # Reset all to non-zero
      for i, ma in enumerate(mal):
        ma.matrices.set(rran(21))
        ma.amplitudes.set(rran(n_dst))
        assert not ma.is_null()
        assert not mal.is_null()
      assert not mal.is_null()
      assert mal.n_params(free=False, non_zero=False) == (21 + n_dst) * length
      assert mal.n_params(free=True, non_zero=False) == (20 + n_dst) * length
      assert mal.n_params(free=False, non_zero=True) == (21 + n_dst) * length
      assert mal.n_params(free=True, non_zero=True) == (20 + n_dst) * length

      # Zero selection of amplitudes
      sel = iran(length, size=max(1,length-1), replace=False)
      mal.zero_amplitudes(sel.astype(numpy.uint))
      if length == 1:
        assert mal.is_null()
      else:
        assert not mal.is_null()
      for i, ma in enumerate(mal):
        if i in sel:
          exp = True
        else:
          exp = False
        assert True == ma.matrices.any("TLS")
        assert exp == (not ma.amplitudes.any())
        assert exp == ma.is_null()
      assert mal.n_params(free=False, non_zero=False) == (21 + n_dst) * length
      assert mal.n_params(free=True, non_zero=False) == (20 + n_dst) * length
      assert mal.n_params(free=False, non_zero=True) == (21 + n_dst) * length - n_dst * len(sel)
      assert mal.n_params(free=True, non_zero=True) == (20 + n_dst) * length - n_dst * len(sel)
      # Only selection should be null
      mal.reset_null_modes()
      for i, ma in enumerate(mal):
        if i in sel:
          exp = True
          assert approx_equal(ma.amplitudes.values, [1.0]*n_dst)
        else:
          exp = False
          assert not_approx_equal(ma.amplitudes.values, [1.0]*n_dst)
        assert exp == (not ma.matrices.any("TLS"))
        assert exp == ma.is_null()
      assert mal.n_params(free=False, non_zero=False) == (21 + n_dst) * length
      assert mal.n_params(free=True, non_zero=False) == (20 + n_dst) * length
      assert mal.n_params(free=False, non_zero=True) == 21 * (length - len(sel)) + n_dst * length
      assert mal.n_params(free=True, non_zero=True) == 20 * (length - len(sel)) + n_dst * length
      # Reset all to non-zero
      for i, ma in enumerate(mal):
        ma.matrices.set(rran(21))
        ma.amplitudes.set(rran(n_dst))
        assert not ma.is_null()
        assert not mal.is_null()
      assert mal.n_params(free=False, non_zero=False) == (21 + n_dst) * length
      assert mal.n_params(free=True, non_zero=False) == (20 + n_dst) * length
      assert mal.n_params(free=False, non_zero=True) == (21 + n_dst) * length
      assert mal.n_params(free=True, non_zero=True) == (20 + n_dst) * length

      # Check modes with zero matrices are reset
      sel = iran(length, size=max(1,length-1), replace=False)
      for i, ma in enumerate(mal):
        if i in sel:
          ma.matrices.reset()
      mal.reset_null_modes()
      # Only selection should be null
      for i, ma in enumerate(mal):
        if i in sel:
          exp = True
          assert approx_equal(ma.amplitudes.values, [1.0]*n_dst)
        else:
          exp = False
          assert not_approx_equal(ma.amplitudes.values, [1.0]*n_dst)
        assert exp == (not ma.matrices.any("TLS"))
        assert exp == ma.is_null()
      assert mal.n_params(free=False, non_zero=False) == (21 + n_dst) * length
      assert mal.n_params(free=True, non_zero=False) == (20 + n_dst) * length
      assert mal.n_params(free=False, non_zero=True) == (21 + n_dst) * length - 21 * len(sel)
      assert mal.n_params(free=True, non_zero=True) == (20 + n_dst) * length - 20 * len(sel)
      # Reset all to non-zero
      for i, ma in enumerate(mal):
        ma.matrices.set(rran(21))
        ma.amplitudes.set(rran(n_dst))
        assert not ma.is_null()

      # Set random amplitudes to negative and then check they are zero-d
      isel = iran(length)
      asel = iran(n_dst, max(1,n_dst-1), replace=False)
      ma = mal.get(isel)
      ma.amplitudes.set(-rran(len(asel)), asel.astype(numpy.uint))
      prev_v = ma.amplitudes.values
      for i, v in enumerate(prev_v):
        if i in asel:
          assert v<0.0
        else:
          assert v>0.0
      mal.zero_negative_amplitudes()
      for i_ma, ma in enumerate(mal):
        if i_ma == isel:
          for i_v, v in enumerate(ma.amplitudes.values):
            if i_v in asel:
              assert v==0.0
            else:
              assert v==prev_v[i_v]
        else:
          assert ma.amplitudes.n_params(non_zero=True) == n_dst
      assert mal.n_params(free=False, non_zero=False) == (21 + n_dst) * length
      assert mal.n_params(free=True, non_zero=False) == (20 + n_dst) * length
      assert mal.n_params(free=False, non_zero=True) == (21 + n_dst) * length - len(asel)
      assert mal.n_params(free=True, non_zero=True) == (20 + n_dst) * length - len(asel)

  # Test initialisation from arrays
  mats = flex.double(rran(3*21).tolist())
  mats.reshape(flex.grid((3,21)))
  amps = flex.double(rran(3*10).tolist())
  amps.reshape(flex.grid((3,10)))
  mal = TLSMatricesAndAmplitudesList(mats, amps)
  assert not mal.is_null()
  for i, m in enumerate(mal):
    assert approx_equal(m.matrices.get(), mats[21*i:21*(i+1)])
    assert approx_equal(m.amplitudes.get(), amps[10*i:10*(i+1)])
  mats.set_selected(flex.bool(flex.grid(mats.all()), True), 0.0)
  amps.set_selected(flex.bool(flex.grid(amps.all()), True), 0.0)
  assert mats.all_eq(0.0)
  assert amps.all_eq(0.0)
  assert not mal.is_null()

  # Test reset
  mal.reset()
  assert mal.is_null()
  for ma in mal:
    assert ma.is_null()

  # Test uijs
  n_atm = 2
  n_dst = 3
  n_mod = 4
  mal = TLSMatricesAndAmplitudesList(length=n_mod, n_amplitudes=n_dst)
  for i, mats in enumerate(TLSMATRICES[:n_mod]):
    ma = mal[i]
    ma.matrices.set(mats)
    ma.amplitudes.set(rran(3))

  # For all modes
  coord_scale = 20.0 # set suitable scale for coordinates in angstrom
  rand_coords = (coord_scale*rran(n_dst*n_atm*3)).reshape(n_dst*n_atm, 3).tolist()
  sites = flex.vec3_double(rand_coords)
  sites.reshape(flex.grid(n_dst,n_atm))
  origins = rran(n_dst*3).reshape(n_dst,3)
  all_uijs = mal.uijs(sites, origins)
  sum_uijs = flex.sym_mat3_double(flex.grid((n_dst, n_atm)))
  assert all_uijs.all() == sum_uijs.all()
  for ma in mal:
    new_uijs = ma.uijs(sites, origins)
    assert new_uijs.all() == (n_dst, n_atm)
    sum_uijs = sum_uijs + new_uijs
  assert approx_equal(numpy.array(all_uijs).flatten(), numpy.array(sum_uijs).flatten())
  # And from an amplitudes selection
  n_tmp = max(1,n_dst-2)
  sel = iran(n_dst, size=n_tmp, replace=False)
  sites = flex.vec3_double(rran(n_tmp*n_atm*3).reshape(n_tmp*n_atm, 3).tolist())
  sites.reshape(flex.grid(n_tmp,n_atm))
  origins = rran(n_tmp*3).reshape(n_tmp,3)
  all_uijs = mal.uijs(sites, origins, sel.astype(numpy.uint))
  sum_uijs = flex.sym_mat3_double(flex.grid((n_tmp, n_atm)))
  assert all_uijs.all() == sum_uijs.all()
  for i, ma in enumerate(mal):
    new_uijs = ma.uijs(sites, origins, sel.astype(numpy.uint))
    assert new_uijs.all() == (n_tmp, n_atm)
    sum_uijs = sum_uijs + new_uijs
  assert approx_equal(numpy.array(all_uijs).flatten(), numpy.array(sum_uijs).flatten())

  # Check copy creates separate object
  mal = TLSMatricesAndAmplitudesList(length=n_mod, n_amplitudes=n_dst)
  for ma in mal:
    ma.matrices.set(rran(21))
    ma.amplitudes.set(rran(n_dst))
  mal_c = mal.copy()
  for i in range(n_mod):
    assert approx_equal(list(mal[i].matrices.get("TLS")), list(mal_c[i].matrices.get("TLS")))
    assert approx_equal(list(mal[i].amplitudes.values), list(mal_c[i].amplitudes.values))
  mal_c.reset_matrices()
  for i in range(n_mod):
    assert mal[i].matrices.any()
    assert not mal_c[i].matrices.any()
    assert approx_equal(list(mal[i].amplitudes.values), list(mal_c[i].amplitudes.values))

  print('OK')

def tst_TLSMatricesAndAmplitudesList_fails():
  """Check that exceptions are raised where expected"""

  n_mod = 4
  n_dst = 3
  n_atm = 2

  mal = TLSMatricesAndAmplitudesList(length=n_mod, n_amplitudes=n_dst)

  # Check indexing errors
  with raises(Exception) as e:
    mal[n_mod+1]
  assert "index out of range of TLSMatricesAndAmplitudesList" == str(e.value)

  # Tolerances
  mal.copy().reset_null_modes(0.0, 0.0) # values can be zero
  mal.copy().reset_null_modes(-1, -1) # values can be -1
  msg = "Tolerance provided must either be positive or -1"
  with raises(Exception) as e:
    mal.reset_null_modes(-0.1, -1)
  assert msg == str(e.value)
  with raises(Exception) as e:
    mal.reset_null_modes(-1, -0.1)
  assert msg == str(e.value)

  # Set some values
  for i, ma in enumerate(mal):
    ma.matrices.set(TLSMATRICES[i])
    ma.amplitudes.set(rran(n_dst))

  # Coordinates sets for following tests
  sites = flex.vec3_double(rran(n_dst*n_atm*3).reshape(n_dst*n_atm, 3).tolist())
  sites.reshape(flex.grid(n_dst,n_atm))
  origins = rran(n_dst*3).reshape(n_dst,3)

  # Compatibility of sites/origins arrays (normalise and uijs)
  # Array has swapped axes
  msg = "Mismatch between the size of origins and first dimension of sites_carts"
  new_sites = flex.vec3_double(rran(n_dst*n_atm*3).reshape(n_dst*n_atm, 3).tolist())
  new_sites.reshape(flex.grid(n_atm,n_dst)) # swap axes
  with raises(Exception) as e:
    mal.uijs(new_sites, origins)
  assert msg == str(e.value)
  # sites is too short/long
  msg = "Mismatch between the size of origins and first dimension of sites_carts"
  for n_tmp in [n_dst-1, n_dst+1]:
    new_sites = flex.vec3_double(rran(n_tmp*n_atm*3).reshape(n_tmp*n_atm, 3).tolist())
    new_sites.reshape(flex.grid(n_tmp,n_atm))
    with raises(Exception) as e:
      mal.uijs(new_sites, origins)
    assert msg == str(e.value)
    # Sites/origins compatible but not same length as amplitudes
    new_origins = rran(n_tmp*3).reshape(n_tmp, 3).tolist()
    with raises(Exception) as e:
      mal.uijs(new_sites, new_origins)
    assert "Mismatch between the size of TLSAmplitudes and the input arrays" == str(e.value)

  # Check dimension of sites_carts
  msg = "sites_carts must be 2-dimensional array of size (n_dst, n_atm)"
  # Make 3-d array of sites
  n_fake = 2 # length of new dimension
  new_sites = flex.vec3_double(rran(n_fake*n_dst*n_atm*3).reshape(n_fake*n_dst*n_atm, 3).tolist())
  new_sites.reshape(flex.grid(n_fake,n_atm,n_dst))
  with raises(Exception) as e:
    mal.uijs(new_sites, origins)
  assert msg == str(e.value)
  # Make 1-d array of sites
  new_sites = flex.vec3_double(rran(n_atm*3).reshape(n_atm, 3).tolist())
  new_sites.reshape(flex.grid(n_atm))
  with raises(Exception) as e:
    mal.uijs(new_sites, origins)
  assert msg == str(e.value)
  # Make other 1-d array
  new_sites = flex.vec3_double(rran(n_atm*3).reshape(n_atm, 3))
  with raises(Exception) as e:
    mal.uijs(new_sites, origins)
  assert msg == str(e.value)

  # Make new sites/origins for subset of ampls for use below
  t_n_dst = n_dst - 1
  new_sites = flex.vec3_double(rran(t_n_dst*n_atm*3).reshape(t_n_dst*n_atm, 3).tolist())
  new_sites.reshape(flex.grid(t_n_dst,n_atm))
  new_origins = rran(t_n_dst*3).reshape(t_n_dst,3)

  # Selection compatibility with sites and origins
  msg = "Mismatch between the size of selection and the input arrays"
  for l in [t_n_dst-1, t_n_dst+1]:
    # Selection too short
    sel = iran(n_dst, size=l, replace=False)
    with raises(Exception) as e:
      mal.uijs(new_sites, new_origins, sel.astype(numpy.uint))
    assert msg == str(e.value)
  # Invalid selection
  sel = iran(n_dst, size=t_n_dst, replace=False)
  sel[-1] = n_dst # larger than allowed
  with raises(Exception) as e:
    mal.uijs(new_sites, new_origins, sel.astype(numpy.uint))
  assert "Selection indices out of range of TLSAmplitudes" == str(e.value)

  print('OK')

if __name__ == "__main__":
  tst_set_precision()
  tst_set_tolerance()
  tst_TLSMatrices()
  tst_TLSMatrices_counts()
  tst_TLSMatrices_fails()
  tst_TLSMatrices_uijs()
  tst_TLSAmplitudes()
  tst_TLSAmplitudes_fails()
  tst_TLSMatricesAndAmplitudes()
  tst_TLSMatricesAndAmplitudes_counts()
  tst_TLSMatricesAndAmplitudes_fails()
  tst_TLSMatricesAndAmplitudes_valid()
  tst_TLSMatricesAndAmplitudes_uijs()
  tst_TLSMatricesAndAmplitudesList()
  tst_TLSMatricesAndAmplitudesList_fails()
