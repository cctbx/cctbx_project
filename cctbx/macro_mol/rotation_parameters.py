"""Module for rotation parameter conversions.

Example usage:

import rotation_parameters
a = rotation_parameters.amore_alpha_beta_gamma(params=(30,40,50))
print a.params
print a.matrix
c = rotation_parameters.cns_theta1_theta2_theta3(matrix=a.matrix)
print c.params
print c.matrix

To see a list of all available converters:

import rotation_parameters
for conv in rotation_parameters.get_converters():
  print conv.__doc__

With kind permission of A.G. Urzhumtsev, some functions are
based on the FORTRAN program CONVROT:
  Urzhumtseva, L.M. & Urzhumtsev, A.G. (1997).
  Tcl/Tk-based programs. II.
  CONVROT: a program to recalculate different rotation descriptions.
  J. Appl. Cryst. 30, 402-410.

Revision history:
  2002 Jan: Created (Ralf W. Grosse-Kunstleve)
"""
from __future__ import absolute_import, division, print_function

import math, types
from six.moves import range

class matrix33(object):
  "Minimal class for the handling of (3x3) matrices."

  def __init__(self, elems = None):
    if (elems is not None):
      self.elems = elems[:]
    else:
      self.elems = [1,0,0,
                    0,1,0,
                    0,0,1]

  def index1d(self, key):
    return key[0] * 3 + key[1]

  def __call__(self, i, j):
    return self.elems[self.index1d((i,j))]

  def __repr__(self):
    return (" %.6f %.6f %.6f\n" * 3) % tuple(self.elems)

  def trace(self):
    m = self.elems
    return m[0] + m[4] + m[8]

  def det(self):
    m = self.elems
    return (  m[0] * (m[4] * m[8] - m[5] * m[7])
            - m[1] * (m[3] * m[8] - m[5] * m[6])
            + m[2] * (m[3] * m[7] - m[4] * m[6]))

def degree_as_radians(angles):
  if (type(angles) in (tuple, list)):
    return [a * math.pi / 180 for a in angles]
  return angles * math.pi / 180

def degree_from_radians(angles):
  if (type(angles) in (tuple, list)):
    return [a / math.pi * 180 for a in angles]
  return angles / math.pi * 180

def adjust_cosine(c):
  if (c >  1.0): return  1.
  if (c < -1.0): return -1.
  return c

def safe_atan2(y, x):
  if (y == 0. and x == 0.): return 0.
  return math.atan2(y, x)

def fmod_positive(phi, period = 360.):
  phi = math.fmod(phi, period)
  if (phi < 0.): phi = phi + period
  return phi

def preprocess_angles(params):
  a_rad = degree_as_radians(params)
  c = [math.cos(a) for a in a_rad]
  s = [math.sin(a) for a in a_rad]
  return c, s

class converter_base(object):
  """Base class for conversions of rotation parameters.
    m = self.matrix is a (3x3) matrix that transforms
    cartesian coordinates of the search body according to:

    (x')   (m(0,0),m(0,1),m(0,2)) (x)
    (y') = (m(1,0),m(1,1),m(1,2)) (y)
    (z')   (m(2,0),m(2,1),m(2,2)) (z)"""

  def __init__(self, params = None, matrix = None, tolerance = 1.e-6):
    if ((params is None) == (matrix is None)):
      raise ArgumentError("Either parameters or matrix required.")
    if (params is not None):
      self.params = params
      self.matrix = self.params_to_matrix(params, tolerance)
    else:
      if (abs(matrix.det()-1) > tolerance):
        raise ValueError("Determinant of matrix is not equal to 1.")
      self.matrix = matrix
      self.params = self.matrix_to_params(matrix, tolerance)
      self.params = self.normalize()

def amore_alpha_beta_gamma_as_matrix(params, tolerance = 1.e-6):
  "Core routine for conversion of Euler angles."
  # Based on subroutine angro7 in CONVROT.
  c, s = preprocess_angles(params)
  return matrix33((
     c[0]*c[1]*c[2]-s[0]*s[2],
    -c[0]*c[1]*s[2]-s[0]*c[2],
     c[0]*s[1],
     s[0]*c[1]*c[2]+c[0]*s[2],
    -s[0]*c[1]*s[2]+c[0]*c[2],
     s[0]*s[1],
    -s[1]*c[2],
     s[1]*s[2],
     c[1]))

def amore_alpha_beta_gamma_from_matrix(matrix, tolerance = 1.e-6):
  "Core routine for conversion of Euler angles."
  # Based on subroutine angro7 in CONVROT.
  c = [0,0,0]
  a = [0,0,0]
  if (abs(matrix(2,2)) > 1.0 + tolerance):
    raise ValueError("Corrupt matrix.")
  c[1] = adjust_cosine(matrix(2,2))
  a[1]=math.acos(c[1])
  if (c[1] == 1.0):
    a[0]=safe_atan2(-matrix(0,1),matrix(1,1))
    a[2]=0.0
  elif (c[1] == -1.0):
    a[0]=safe_atan2(-matrix(0,1),matrix(1,1))
    a[2]=0.0
  else:
    a[0]=safe_atan2(matrix(1,2), matrix(0,2))
    a[2]=safe_atan2(matrix(2,1),-matrix(2,0))
  return degree_from_radians(a)

def amore_alpha_beta_gamma_normalize(params, alternative):
  a, b, g = params
  g = fmod_positive(g)
  if ((g >= 180.) == (alternative == 0)):
    return [fmod_positive(a+180.), fmod_positive(-b), g-180.]
  return [fmod_positive(a), fmod_positive(b), g]

def amore_kappa_l_m_n_as_matrix(params, tolerance = 1.e-6):
  "Core routine for conversion of direction cosines."
  # Based on subroutine angro17 in CONVROT.
  k = degree_as_radians(params[0])
  ck = math.cos(k)
  sk = math.sin(k)
  c = params[1:4]
  sc2=c[0]*c[0]+c[1]*c[1]+c[2]*c[2]
  if (abs(sc2-1.0) > tolerance):
    raise ValueError("Corrupt direction cosines: the sum of their squares is not equal to 1.")
  c = [x / math.sqrt(sc2) for x in c]
  return matrix33((
      (1.-ck)*c[0]*c[0]+ck,
      (1.-ck)*c[0]*c[1]-sk*c[2],
      (1.-ck)*c[0]*c[2]+sk*c[1],
      (1.-ck)*c[0]*c[1]+sk*c[2],
      (1.-ck)*c[1]*c[1]+ck,
      (1.-ck)*c[1]*c[2]-sk*c[0],
      (1.-ck)*c[0]*c[2]-sk*c[1],
      (1.-ck)*c[1]*c[2]+sk*c[0],
      (1.-ck)*c[2]*c[2]+ck))

def amore_kappa_l_m_n_from_matrix(matrix, tolerance = 1.e-6):
  "Core routine for conversion of direction cosines."
  ck = (matrix.trace() - 1.) / 2.
  if (-1.-tolerance > ck > 1.+tolerance):
    raise ValueError("Corrupt matrix.")
  if (ck >= 1.):
    return [0.,0.,0.,1.]
  if (ck < -1.): ck = -1.
  kappa = math.acos(ck)
  sk = math.sin(kappa)
  c = [(matrix(i,i) - ck) / (1.-ck) for i in range(3)]
  for i in range(3):
    if (c[i] < -tolerance):
      raise ValueError("Corrupt matrix.")
    if (c[i] <= 0.): c[i] = 0.
    else: c[i] = math.sqrt(c[i])
  # At this point the absolute values of the direction cosines
  # are given by c. The signs are determined by trial and error.
  # This is slower than a strictly deterministic algorithm, but
  # requires less attention to special cases. Therefore the
  # source code is more compact. Since there are only eight
  # trials, the runtime penalty is very modest.
  best_fc = None
  best_sum_d2 = None
  f = [0,0,0]
  for f[0] in (1,-1):
    for f[1] in (1,-1):
      for f[2] in (1,-1):
        fc = [f[i] * c[i] for i in range(3)]
        ok = 1
        sum_d2 = 0
        for i,j,k in ((0,1,2), (1,2,0), (2,0,1)):
          dij = abs(matrix(i,j) - ((1.-ck)*fc[i]*fc[j]-sk*fc[k]))
          dji = abs(matrix(j,i) - ((1.-ck)*fc[i]*fc[j]+sk*fc[k]))
          if (max(dij, dji) > tolerance):
            ok = 0
            break
          sum_d2 += dij*dij + dji*dji
        if (ok and (best_sum_d2 is None or best_sum_d2 > sum_d2)):
          best_fc = fc
          best_sum_d2 = sum_d2
  if (best_fc is None):
    raise ValueError("Corrupt matrix.")
  return [degree_from_radians(kappa)] + best_fc

def amore_kappa_l_m_n_normalize(kappa, l, m, n, alternative):
  kappa = fmod_positive(kappa)
  if ((kappa > 180.) == (alternative == 0)):
    return [360.-kappa, -l, -m, -n]
  return [kappa, l, m, n]

def generic_polar_normalize(a1, a2, kappa, alternative):
  kappa = fmod_positive(kappa)
  if ((kappa > 180.) == (alternative == 0)):
    return [fmod_positive(a1+180.), fmod_positive(-a2+180.), 360.-kappa]
  return [fmod_positive(a1), fmod_positive(a2), kappa]

class amore_alpha_beta_gamma(converter_base):
  "AMoRe alpha, beta, gamma (Crowther, 1972)"

  def params_to_matrix(self, params, tolerance = 1.e-6):
    return amore_alpha_beta_gamma_as_matrix(params, tolerance)

  def matrix_to_params(self, matrix, tolerance = 1.e-6):
    return amore_alpha_beta_gamma_from_matrix(matrix, tolerance)

  def normalize(self, alternative = 0):
    return amore_alpha_beta_gamma_normalize(self.params, alternative)

class amore_kappa_l_m_n(converter_base):
  "AMoRe kappa, l, m, n (Diamond, 1993; Rossmann, 1993)"

  def params_to_matrix(self, params, tolerance = 1.e6):
    return amore_kappa_l_m_n_as_matrix(params, tolerance)

  def matrix_to_params(self, matrix, tolerance = 1.e-6):
    return amore_kappa_l_m_n_from_matrix(matrix, tolerance)

  def normalize(self, alternative = 0):
    kappa, l, m, n = self.params
    return amore_kappa_l_m_n_normalize(kappa, l, m, n, alternative)

def cns_theta_123_as_p2m(params):
  return [params[0]+params[2], params[1], params[0]-params[2]]

def cns_theta_123_from_p2m(params):
  return [(params[0]+params[2])/2., params[1], (params[0]-params[2])/2.]

class cns_theta1_theta2_theta3(converter_base):
  "CNS theta1, theta2, theta3 (Rossmann, 1962)"

  def params_to_matrix(self, params, tolerance = 1.e-6):
    return amore_alpha_beta_gamma_as_matrix(
      [270.-params[2], -params[1], 90.-params[0]], tolerance)

  def matrix_to_params(self, matrix, tolerance = 1.e-6):
    a, b, g = amore_alpha_beta_gamma_from_matrix(matrix, tolerance)
    return [90.-g, -b, 270.-a]

  def normalize(self, alternative = 0):
    return amore_alpha_beta_gamma_normalize(self.params, alternative)

class cns_theta_plus_theta2_theta_minus(converter_base):
  "CNS theta+, theta2, theta- (Lattman, 1972)"

  def params_to_matrix(self, params, tolerance = 1.e-6):
    t1, t2, t3 = cns_theta_123_from_p2m(self.params)
    return amore_alpha_beta_gamma_as_matrix(
      [270.-t3, -t2, 90.-t1], tolerance)

  def matrix_to_params(self, matrix, tolerance = 1.e-6):
    a, b, g = amore_alpha_beta_gamma_from_matrix(matrix, tolerance)
    return cns_theta_123_as_p2m([90.-g, -b, 270.-a])

  def normalize(self, alternative = 0):
    tp = fmod_positive(self.params[0], 720.)
    t2 = fmod_positive(self.params[1])
    tm = fmod_positive(self.params[2]+360.,720.)-360.
    if ((tp >= 360.) == (alternative == 0)):
      return [tp-360., fmod_positive(-t2), tm]
    return [tp, t2, tm]

class cns_psi_phi_kappa(converter_base):
  "CNS psi, phi, kappa (Rossmann, 1962)"

  def params_to_matrix(self, params, tolerance = 1.e-6):
    c, s = preprocess_angles((180.-params[0], params[1]+180.))
    return amore_kappa_l_m_n_as_matrix(
      [params[2], s[0]*c[1], c[0], -s[0]*s[1]], tolerance)

  def matrix_to_params(self, matrix, tolerance = 1.e-6):
    kappa, l, m, n = amore_kappa_l_m_n_from_matrix(matrix, tolerance)
    psi = math.acos(m)
    phi = safe_atan2(-n, l)
    return degree_from_radians([math.pi-psi, math.pi+phi]) + [kappa]

  def normalize(self, alternative = 0):
    psi, phi, kappa = self.params
    phi, psi, kappa = generic_polar_normalize(phi, psi, kappa, alternative)
    return [psi, phi, kappa]

class cns_axis_x_axis_y_axis_z_axis_kappa(converter_base):
  "CNS axis_x, axis_y, axis_z, axis_kappa (Brunger, 1992)"

  def params_to_matrix(self, params, tolerance = 1.e6):
    return amore_kappa_l_m_n_as_matrix(
      [params[3], -params[0], -params[1], -params[2]], tolerance)

  def matrix_to_params(self, matrix, tolerance = 1.e-6):
    kappa, l, m, n = amore_kappa_l_m_n_from_matrix(matrix, tolerance)
    return [-l, -m, -n, kappa]

  def normalize(self, alternative = 0):
    ax, ay, az, kappa = self.params
    kappa, ax, ay, az = amore_kappa_l_m_n_normalize(
      kappa, ax, ay, az, alternative)
    return [ax, ay, az, kappa]

class ccp4_alpha_beta_gamma_fobs_fcalc(converter_base):
  "CCP4 alpha, beta, gamma, Fobs/Fcalc (Crowther, 1972)"

  def params_to_matrix(self, params, tolerance = 1.e-6):
    return amore_alpha_beta_gamma_as_matrix(params, tolerance)

  def matrix_to_params(self, matrix, tolerance = 1.e-6):
    return amore_alpha_beta_gamma_from_matrix(matrix, tolerance)

  def normalize(self, alternative = 0):
    return amore_alpha_beta_gamma_normalize(self.params, alternative)

class ccp4_alpha_beta_gamma_fcalc_fobs(converter_base):
  "CCP4 alpha, beta, gamma, Fcalc/Fobs (Crowther, 1972)"

  def params_to_matrix(self, params, tolerance = 1.e-6):
    return amore_alpha_beta_gamma_as_matrix(
      [-params[2], -params[1], -params[0]], tolerance)

  def matrix_to_params(self, matrix, tolerance = 1.e-6):
    a, b, g = amore_alpha_beta_gamma_from_matrix(matrix, tolerance)
    return [-g, -b, -a]

  def normalize(self, alternative = 0):
    return amore_alpha_beta_gamma_normalize(self.params, alternative)

def ccp4_phi_omega_kappa_fcalc_fobs_params_to_matrix(
      params, tolerance = 1.e-6):
  c, s = preprocess_angles(params[:2])
  return amore_kappa_l_m_n_as_matrix(
    [params[2], s[1] * c[0], s[1] * s[0], c[1]], tolerance)

def ccp4_phi_omega_kappa_fobs_fcalc_matrix_to_params(
      matrix, tolerance = 1.e-6):
  kappa, l, m, n = amore_kappa_l_m_n_from_matrix(matrix, tolerance)
  omega = math.acos(n)
  phi = safe_atan2(m, l)
  return degree_from_radians([phi, omega]) + [kappa]

class ccp4_phi_omega_kappa_fobs_fcalc(converter_base):
  "CCP4 phi, omega, kappa, Fobs/Fcalc (Crowther, 1972)"

  def params_to_matrix(self, params, tolerance = 1.e-6):
    return ccp4_phi_omega_kappa_fcalc_fobs_params_to_matrix(
      params, tolerance)

  def matrix_to_params(self, matrix, tolerance = 1.e-6):
    return ccp4_phi_omega_kappa_fobs_fcalc_matrix_to_params(
      matrix, tolerance)

  def normalize(self, alternative = 0):
    phi, omega, kappa = self.params
    return generic_polar_normalize(phi, omega, kappa, alternative)

class ccp4_phi_omega_kappa_fcalc_fobs(converter_base):
  "CCP4 phi, omega, kappa, Fcalc/Fobs (Crowther, 1972)"

  def params_to_matrix(self, params, tolerance = 1.e-6):
    return ccp4_phi_omega_kappa_fcalc_fobs_params_to_matrix(
      [params[0]-180., 180.-params[1], params[2]], tolerance)

  def matrix_to_params(self, matrix, tolerance = 1.e-6):
    phi, omega, kappa = ccp4_phi_omega_kappa_fobs_fcalc_matrix_to_params(
      matrix, tolerance)
    return [180.+phi, 180.-omega, kappa]

  def normalize(self, alternative = 0):
    phi, omega, kappa = self.params
    return generic_polar_normalize(phi, omega, kappa, alternative)

def get_converters():
  converters = (
    amore_alpha_beta_gamma,
    amore_kappa_l_m_n,
    cns_theta1_theta2_theta3,
    cns_theta_plus_theta2_theta_minus,
    cns_psi_phi_kappa,
    cns_axis_x_axis_y_axis_z_axis_kappa,
    ccp4_alpha_beta_gamma_fobs_fcalc,
    ccp4_alpha_beta_gamma_fcalc_fobs,
    ccp4_phi_omega_kappa_fobs_fcalc,
    ccp4_phi_omega_kappa_fcalc_fobs)
  # Check if the list is up to date.
  for global_attr in globals().values():
    if (hasattr(global_attr, "__bases__")):
      if (converter_base in global_attr.__bases__):
        assert global_attr in converters
  return converters

def get_converter_by_docstring(docstring):
  for conv in get_converters():
    if (conv.__doc__ == docstring): return conv
  return None

if (__name__ == "__main__"):
  # Exhaustive regression tests.
  # Usages:
  #
  #   python rotation_parameters.py
  #     Exhaustive self-consistency test for all conversions.
  #     Expected output: a list of the conversion types.
  #
  #   python rotation_parameters.py --CNS > cns.inp
  #   cns_solve < cns.inp > cns.out
  #   grep Distance cns.out | uniq
  #     All distances shown in cns.out should be zero.

  import sys

  def cns_input(matrix):
    print("rotman")
    print("matrix=" + (("(%.10g %.10g %.10g)" * 3) % matrix.elems))
    print("copy")
    for label, conv in (("eule", cns_theta1_theta2_theta3),
                        ("latt", cns_theta_plus_theta2_theta_minus),
                        ("sphe", cns_psi_phi_kappa),
                        ("axis", cns_axis_x_axis_y_axis_z_axis_kappa)):
      if (label == "axis"): fmt = "%s=(%.10g %.10g %.10g) %.10g"
      else:                 fmt = "%s=(%.10g %.10g %.10g)"
      print(fmt % tuple([label] + conv(matrix=matrix).params))
      print("? distance")
      print(fmt % tuple([label] + conv(matrix=matrix).normalize(1)))
      print("? distance")
    print("end")

  def compare_lists(list1, list2):
    for i in range(len(list1)):
      if (abs(list1[i] - list2[i]) > 1.e-5):
        print((" %.5f %.5f %.5f\n" * 3) % tuple(list1))
        print((" %.5f %.5f %.5f\n" * 3) % tuple(list2))
        raise RuntimeError("Matrix mismatch.")

  def check_conversion(conv, matrix):
    if ("--CNS" in sys.argv[1:]):
      cns_input(matrix)
      return
    rot0 = conv(matrix = matrix)
    rot1 = conv(params = rot0.params)
    rot2 = conv(params = rot0.normalize(1))
    if ("-v" in sys.argv[1:]):
      print("rot1", rot1.params)
      print("rot2", rot2.params)
    compare_lists(rot0.matrix.elems, rot1.matrix.elems)
    compare_lists(rot0.matrix.elems, rot2.matrix.elems)

  def check_45deg_incr(conv):
    for a1 in range(0, 361, 45):
      for a2 in range(0, 361, 45):
        for a3 in range(0, 361, 45):
          if ("-v" in sys.argv[1:]):
            print("Angles:", a1, a2, a3)
          matrix = cns_theta1_theta2_theta3(params = (a1, a2, a3)).matrix
          try:
            check_conversion(conv, matrix)
          except RuntimeError as e:
            print(e)
            print()
            return

  def random_angles():
    import random
    rng = random.random
    return [rng() * 180 for i in range(3)]

  def run():
    for conv in get_converters():
      if (not "--CNS" in sys.argv[1:]):
        print(conv.__doc__)
      for trial in range(10):
        matrix = cns_theta1_theta2_theta3(params = random_angles()).matrix
        check_conversion(conv, matrix)
      if (not "--RandomOnly" in sys.argv[1:]):
        check_45deg_incr(conv)
      if ("--CNS" in sys.argv[1:]):
        print("stop")
        break

  run()
