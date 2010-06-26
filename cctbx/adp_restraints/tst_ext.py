from libtbx.test_utils import approx_equal, not_approx_equal, eps_eq, show_diff
from iotbx.shelx import from_ins
from cctbx.array_family import flex
from cctbx import adptbx
from cctbx import uctbx
from cctbx import adp_restraints
from scitbx import matrix
import libtbx.load_env
import math, os, sys
from cStringIO import StringIO

def finite_difference_gradients(restraint_type,
                                proxy,
                                sites_cart=None,
                                u_cart=None,
                                u_iso=None,
                                use_u_aniso=None,
                                eps=1.e-8):
  def residual(restraint_type, proxy, sites_cart=None,
               u_cart=None, u_iso=None, use_u_aniso=None):
    if sites_cart is not None:
      return restraint_type(
        sites_cart=sites_cart,
        u_cart=u_cart,
        proxy=proxy).residual()
    elif u_iso is None:
      return restraint_type(u_cart=u_cart,proxy=proxy).residual()
    else:
      assert use_u_aniso is not None
      return restraint_type(
        u_cart=u_cart,
        u_iso=u_iso,
        use_u_aniso=use_u_aniso,
        proxy=proxy).residual()
  result_aniso = [(0,0,0,0,0,0)]*len(u_cart)
  result_iso = [0] * len(u_cart)
  if sites_cart is not None:
    assert len(sites_cart) == len(u_cart)
  for i in xrange(len(u_cart)):
    if u_iso is None:
      result_aniso_i = []
      for j in xrange(6):
        h = [0,0,0,0,0,0]
        h[j] = eps
        h = matrix.sym(sym_mat3=h)
        u_cart[i]=list((matrix.sym(sym_mat3=u_cart[i]) + h).as_sym_mat3())
        qp = residual(restraint_type, proxy,
                      sites_cart=sites_cart, u_cart=u_cart)
        u_cart[i]=list((matrix.sym(sym_mat3=u_cart[i]) - 2*h).as_sym_mat3())
        qm = residual(restraint_type, proxy,
                      sites_cart=sites_cart, u_cart=u_cart)
        dq = (qp-qm)/2
        result_aniso_i.append(dq/(eps))
      result_aniso[i] = result_aniso_i
    else:
      if use_u_aniso[i]:
        result_aniso_i = []
        for j in xrange(6):
          h = [0,0,0,0,0,0]
          h[j] = eps
          h = matrix.sym(sym_mat3=h)
          u_cart[i]=list((matrix.sym(sym_mat3=u_cart[i]) + h).as_sym_mat3())
          qp = residual(restraint_type, proxy,
                        u_cart=u_cart, u_iso=u_iso, use_u_aniso=use_u_aniso)
          u_cart[i]=list((matrix.sym(sym_mat3=u_cart[i]) - 2*h).as_sym_mat3())
          qm = residual(restraint_type, proxy,
                        u_cart=u_cart, u_iso=u_iso, use_u_aniso=use_u_aniso)
          dq = (qp-qm)/2
          result_aniso_i.append(dq/(eps))
        result_aniso[i] = result_aniso_i
      else:
        u_iso[i] += eps
        qp = residual(restraint_type, proxy,
                      u_cart=u_cart, u_iso=u_iso, use_u_aniso=use_u_aniso)
        u_iso[i] -= 2*eps
        qm = residual(restraint_type, proxy,
                      u_cart=u_cart, u_iso=u_iso, use_u_aniso=use_u_aniso)
        dq = (qp-qm)/2
        result_iso[i] = dq/(eps)
  return result_aniso, result_iso

result = [
  ['C1',   'C2',   0.0039,  0.0162,  0.0123],
  ['C1',   'N1',   0.0002,  0.0129,  0.0131],
  ['C2',   'C1',   0.0039,  0.0123,  0.0162],
  ['C3',   'C4',   0.0001,  0.0147,  0.0146],
  ['C3',   'C8',   0.0024,  0.0078,  0.0102],
  ['C4',   'C3',   0.0001,  0.0146,  0.0147],
  ['C4',   'C5',   0.0013,  0.0156,  0.0144],
  ['C5',   'C4',   0.0013,  0.0144,  0.0156],
  ['C5',   'C6',   0.0012,  0.0109,  0.0121],
  ['C6',   'C5',   0.0012,  0.0121,  0.0109],
  ['C6',   'C7',   0.0002,  0.0171,  0.0169],
  ['C6',   'O1',   0.0008,  0.0132,  0.0140],
  ['C7',   'C6',   0.0002,  0.0169,  0.0171],
  ['C7',   'C8',   0.0004,  0.0165,  0.0161],
  ['C8',   'C3',   0.0024,  0.0102,  0.0078],
  ['C8',   'C7',   0.0004,  0.0161,  0.0165],
  ['C9',   'O2',   0.0017,  0.0106,  0.0123],
  ['C11',  'O3',   0.0007,  0.0151,  0.0145],
  ['C11',  'N3',   0.0009,  0.0207,  0.0198],
  ['C12',  'C13',  0.0006,  0.0114,  0.0119],
  ['C12',  'N3',   0.0040,  0.0193,  0.0153],
  ['C13',  'C12',  0.0006,  0.0119,  0.0114],
  ['C13',  'O4',   0.0001,  0.0128,  0.0130],
  ['C13',  'N4',   0.0009,  0.0110,  0.0119],
  ['C14',  'N4',   0.0006,  0.0090,  0.0096],
  ['C16',  'C17',  0.0017,  0.0168,  0.0186],
  ['C16',  'C21',  0.0023,  0.0205,  0.0183],
  ['C17',  'C16',  0.0017,  0.0186,  0.0168],
  ['C17',  'C18',  0.0063,  0.0178,  0.0241],
  ['C18',  'C17',  0.0063,  0.0241,  0.0178],
  ['C18',  'C19',  0.0049,  0.0358,  0.0309],
  ['C19',  'C18',  0.0049,  0.0309,  0.0358],
  ['C19',  'C20',  0.0012,  0.0207,  0.0196],
  ['C20',  'C19',  0.0012,  0.0196,  0.0207],
  ['C20',  'C21',  0.0006,  0.0163,  0.0157],
  ['C21',  'C16',  0.0023,  0.0183,  0.0205],
  ['C21',  'C20',  0.0006,  0.0157,  0.0163],
  ['C22',  'N5',   0.0015,  0.0098,  0.0083],
  ['C23',  'C24',  0.0002,  0.0072,  0.0073],
  ['C24',  'C23',  0.0002,  0.0073,  0.0072],
  ['C25',  'C27',  0.0001,  0.0075,  0.0076],
  ['C27',  'C25',  0.0001,  0.0076,  0.0075],
  ['C28',  'O6',   0.0023,  0.0192,  0.0169],
  ['C28',  'O7',   0.0001,  0.0120,  0.0119],
  ['O1',   'C6',   0.0008,  0.0140,  0.0132],
  ['O2',   'C9',   0.0017,  0.0123,  0.0106],
  ['O3',   'C11',  0.0007,  0.0145,  0.0151],
  ['O4',   'C13',  0.0001,  0.0130,  0.0128],
  ['O6',   'C28',  0.0023,  0.0169,  0.0192],
  ['O7',   'C28',  0.0001,  0.0119,  0.0120],
  ['N1',   'C1',   0.0002,  0.0131,  0.0129],
  ['N3',   'C11',  0.0009,  0.0198,  0.0207],
  ['N3',   'C12',  0.0040,  0.0153,  0.0193],
  ['N4',   'C13',  0.0009,  0.0119,  0.0110],
  ['N4',   'C14',  0.0006,  0.0096,  0.0090],
  ['N5',   'C22',  0.0015,  0.0083,  0.0098]]

def exercise_rigid_bond_test():
  """
  Results compared with THMA11 (Ver. 20-04-91) - TLS Thermal Motion
  Analysis used as a part of WinGX (WinGX - Crystallographic Program
  System for Windows)
  """
  ins_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/enk_11i.res", test=os.path.isfile)
  if (ins_file is None):
    print "Skipping exercise_rigid_bond_test(): input file not available"
    return
  ins_xray_structure = from_ins.from_ins(file_name = ins_file)
  sites_frac = ins_xray_structure.sites_frac()
  sites_cart = ins_xray_structure.sites_cart()
  ustars = ins_xray_structure.scatterers().extract_u_star()
  scatterers = ins_xray_structure.scatterers()
  j = 0
  for site_cart_1,site_frac_1,ustar_1,scat_1 in zip(sites_cart,sites_frac,ustars,scatterers):
    for site_cart_2,site_frac_2,ustar_2, scat_2 in zip(sites_cart,sites_frac,ustars,scatterers):
      d = math.sqrt(flex.sum(flex.pow2(flex.double(site_cart_1)-\
                                       flex.double(site_cart_2))))
      if(d > 1.1 and d < 1.55):
        p = adp_restraints.rigid_bond_pair(site_frac_1,
                                           site_frac_2,
                                           ustar_1,
                                           ustar_2,
                                           ins_xray_structure.unit_cell())
        if(0):
          print "%4s %4s %7.4f %7.4f %7.4f" % \
                (scat_1.label,scat_2.label,p.delta_z(),p.z_12(),p.z_21())
        r = result[j]
        assert r[0] == scat_1.label
        assert r[1] == scat_2.label
        assert approx_equal(r[2], p.delta_z(), 1.e-4)
        assert approx_equal(r[3], p.z_12(), 1.e-4)
        assert approx_equal(r[4], p.z_21(), 1.e-4)
        j += 1
  assert j == 56

def exercise_rigid_bond():
  i_seqs = (1,2)
  weight = 1
  p = adp_restraints.rigid_bond_proxy(i_seqs=i_seqs,weight=weight)
  assert p.i_seqs == i_seqs
  assert p.weight == weight
  sites = ((1,2,3),(2,3,4))
  u_cart = ((1,2,3,4,5,6), (3,4,5,6,7,8))
  expected_gradients = ((-4, -4, -4, -8, -8, -8), (4, 4, 4, 8, 8, 8))
  r = adp_restraints.rigid_bond(sites=sites, u_cart=u_cart, weight=weight)
  assert r.weight == weight
  assert approx_equal(r.delta_z(), 6)
  assert approx_equal(r.residual(), 36)
  assert approx_equal(r.gradients(), expected_gradients)
  assert approx_equal(r.sites, sites)
  assert approx_equal(r.u_cart, u_cart)
  sites_cart = flex.vec3_double(((1,2,3),(2,5,4),(3,4,5)))
  u_cart = flex.sym_mat3_double(((1,2,3,4,5,6),
                                 (2,3,3,5,7,7),
                                 (3,4,5,3,7,8)))
  r = adp_restraints.rigid_bond(sites_cart=sites_cart,
                                u_cart=u_cart,
                                proxy=p)
  assert approx_equal(r.sites, sites_cart[1:3])
  assert approx_equal(r.u_cart, u_cart[1:3])
  assert approx_equal(r.weight, weight)
  unit_cell = uctbx.unit_cell([15,25,30,90,90,90])
  sites_frac = unit_cell.fractionalize(sites_cart=sites_cart)
  u_star = flex.sym_mat3_double([
    adptbx.u_cart_as_u_star(unit_cell, u_cart_i)
    for u_cart_i in u_cart])
  pair = adp_restraints.rigid_bond_pair(sites_frac[1],
                                     sites_frac[2],
                                     u_star[1],
                                     u_star[2],
                                     unit_cell)
  assert approx_equal(pair.delta_z(), r.delta_z())
  assert approx_equal(pair.z_12(), r.z_12())
  assert approx_equal(pair.z_21(), r.z_21())
  #
  gradients_aniso_cart = flex.sym_mat3_double(sites_cart.size(), (0,0,0,0,0,0))
  gradients_iso = flex.double(sites_cart.size(), 0)
  proxies = adp_restraints.shared_rigid_bond_proxy([p,p])
  residuals = adp_restraints.rigid_bond_residuals(
    sites_cart=sites_cart, u_cart=u_cart, proxies=proxies)
  assert approx_equal(residuals, (r.residual(),r.residual()))
  deltas = adp_restraints.rigid_bond_deltas(
    sites_cart=sites_cart, u_cart=u_cart, proxies=proxies)
  assert approx_equal(deltas, (r.delta_z(),r.delta_z()))
  residual_sum = adp_restraints.rigid_bond_residual_sum(
    sites_cart=sites_cart,
    u_cart=u_cart,
    proxies=proxies,
    gradients_aniso_cart=gradients_aniso_cart)
  assert approx_equal(residual_sum, 2 * r.residual())
  for g,e in zip(gradients_aniso_cart[1:3], r.gradients()):
    assert approx_equal(g, matrix.col(e)*2)
  fd_grads_aniso, fd_grads_iso = finite_difference_gradients(
    restraint_type=adp_restraints.rigid_bond,
    proxy=p,
    sites_cart=sites_cart,
    u_cart=u_cart)
  for g,e in zip(gradients_aniso_cart, fd_grads_aniso):
    assert approx_equal(g, matrix.col(e)*2)
  #
  # check frame invariance of residual
  #
  u_cart_1 = matrix.sym(sym_mat3=(0.1,0.2,0.05,0.03,0.02,0.01))
  u_cart_2 = matrix.sym(sym_mat3=(0.21,0.32,0.11,0.02,0.02,0.07))
  u_cart = (u_cart_1.as_sym_mat3(),u_cart_2.as_sym_mat3())
  site_cart_1 = matrix.col((1,2,3))
  site_cart_2 = matrix.col((3,1,4.2))
  sites = (tuple(site_cart_1),tuple(site_cart_2))
  a = adp_restraints.rigid_bond(sites=sites, u_cart=u_cart, weight=1)
  expected_residual = a.residual()
  gen = flex.mersenne_twister()
  for i in range(20):
    R = matrix.rec(gen.random_double_r3_rotation_matrix(),(3,3))
    u_cart_1_rot = R * u_cart_1 * R.transpose()
    u_cart_2_rot = R * u_cart_2 * R.transpose()
    u_cart = (u_cart_1_rot.as_sym_mat3(),u_cart_2_rot.as_sym_mat3())
    site_cart_1_rot = R * site_cart_1
    site_cart_2_rot = R * site_cart_2
    sites = (tuple(site_cart_1_rot),tuple(site_cart_2_rot))
    a = adp_restraints.rigid_bond(sites=sites, u_cart=u_cart, weight=1)
    assert approx_equal(a.residual(), expected_residual)

def exercise_adp_similarity():
  u_cart = ((1,3,2,4,3,6),(2,4,2,6,5,1))
  u_iso = (-1,-1)
  use_u_aniso = (True, True)
  weight = 1
  a = adp_restraints.adp_similarity(
    u_cart=u_cart, u_iso=u_iso, use_u_aniso=use_u_aniso, weight=weight)
  assert approx_equal(a.u_cart, u_cart)
  assert approx_equal(a.u_iso, u_iso)
  assert approx_equal(a.use_u_aniso, use_u_aniso)
  assert a.weight == weight
  assert approx_equal(a.residual(), 68)
  assert approx_equal(a.gradients(),
    ((-2.0, -2.0, 0.0, -8.0, -8.0, 20.0), (2.0, 2.0, -0.0, 8.0, 8.0, -20.0)))
  assert approx_equal(a.deltas(), (-1.0, -1.0, 0.0, -2.0, -2.0, 5.0))
  assert approx_equal(a.rms_deltas(), 2.7487370837451071)
  #
  u_cart = ((1,3,2,4,3,6),(-1,-1,-1,-1,-1,-1))
  u_iso = (-1,2)
  use_u_aniso = (True, False)
  a = adp_restraints.adp_similarity(
    u_cart=u_cart, u_iso=u_iso, use_u_aniso=use_u_aniso, weight=weight)
  assert approx_equal(a.u_cart, u_cart)
  assert approx_equal(a.u_iso, u_iso)
  assert approx_equal(a.use_u_aniso, use_u_aniso)
  assert a.weight == weight
  assert approx_equal(a.residual(), 124)
  assert approx_equal(a.gradients(),
    ((-2, 2, 0, 16, 12, 24), (2, -2, 0, -16, -12, -24)))
  assert approx_equal(a.deltas(), (-1, 1, 0, 4, 3, 6))
  assert approx_equal(a.rms_deltas(), 3.711842908553348)
  #
  i_seqs_aa = (1,2) # () - ()
  i_seqs_ai = (1,0) # () - o
  i_seqs_ia = (3,2) #  o - ()
  i_seqs_ii = (0,3) #  o - o
  p_aa = adp_restraints.adp_similarity_proxy(i_seqs=i_seqs_aa,weight=weight)
  p_ai = adp_restraints.adp_similarity_proxy(i_seqs=i_seqs_ai,weight=weight)
  p_ia = adp_restraints.adp_similarity_proxy(i_seqs=i_seqs_ia,weight=weight)
  p_ii = adp_restraints.adp_similarity_proxy(i_seqs=i_seqs_ii,weight=weight)
  assert p_aa.i_seqs == i_seqs_aa
  assert p_aa.weight == weight
  u_cart = flex.sym_mat3_double(((-1,-1,-1,-1,-1,-1),
                                 (1,2,2,4,3,6),
                                 (2,4,2,6,5,1),
                                 (-1,-1,-1,-1,-1,-1)))
  u_iso = flex.double((1,-1,-1,2))
  use_u_aniso = flex.bool((False, True,True,False))
  for p in (p_aa,p_ai,p_ia,p_ii):
    a = adp_restraints.adp_similarity(u_cart=u_cart,
                                      u_iso=u_iso,
                                      use_u_aniso=use_u_aniso,
                                      proxy=p)
    assert approx_equal(
      a.u_cart, (u_cart[p.i_seqs[0]],u_cart[p.i_seqs[1]]))
    assert approx_equal(a.weight, weight)
    #
    gradients_aniso_cart = flex.sym_mat3_double(u_cart.size(), (0,0,0,0,0,0))
    gradients_iso = flex.double(u_cart.size(), 0)
    proxies = adp_restraints.shared_adp_similarity_proxy([p,p])
    residuals = adp_restraints.adp_similarity_residuals(
      u_cart=u_cart, u_iso=u_iso, use_u_aniso=use_u_aniso, proxies=proxies)
    assert approx_equal(residuals, (a.residual(),a.residual()))
    deltas_rms = adp_restraints.adp_similarity_deltas_rms(
      u_cart=u_cart, u_iso=u_iso, use_u_aniso=use_u_aniso, proxies=proxies)
    assert approx_equal(deltas_rms, (a.rms_deltas(),a.rms_deltas()))
    residual_sum = adp_restraints.adp_similarity_residual_sum(
      u_cart=u_cart,
      u_iso=u_iso,
      use_u_aniso=use_u_aniso,
      proxies=proxies,
      gradients_aniso_cart=gradients_aniso_cart,
      gradients_iso=gradients_iso)
    assert approx_equal(residual_sum, 2 * a.residual())
    fd_grads_aniso, fd_grads_iso = finite_difference_gradients(
      restraint_type=adp_restraints.adp_similarity,
      proxy=p,
      u_cart=u_cart,
      u_iso=u_iso,
      use_u_aniso=use_u_aniso)
    for g,e in zip(gradients_aniso_cart, fd_grads_aniso):
      assert approx_equal(g,  matrix.col(e)*2)
    for g,e in zip(gradients_iso, fd_grads_iso):
      assert approx_equal(g, e*2)
  #
  # check frame invariance of residual
  #
  u_cart_1 = matrix.sym(sym_mat3=(0.1,0.2,0.05,0.03,0.02,0.01))
  u_cart_2 = matrix.sym(sym_mat3=(0.21,0.32,0.11,0.02,0.02,0.07))
  u_cart = (u_cart_1.as_sym_mat3(),u_cart_2.as_sym_mat3())
  u_iso = (-1, -1)
  use_u_aniso = (True, True)
  a = adp_restraints.adp_similarity(
    u_cart=u_cart, u_iso=u_iso, use_u_aniso=use_u_aniso, weight=1)
  expected_residual = a.residual()
  gen = flex.mersenne_twister()
  for i in range(20):
    R = matrix.rec(gen.random_double_r3_rotation_matrix(),(3,3))
    u_cart_1_rot = R * u_cart_1 * R.transpose()
    u_cart_2_rot = R * u_cart_2 * R.transpose()
    u_cart = (u_cart_1_rot.as_sym_mat3(),u_cart_2_rot.as_sym_mat3())
    a = adp_restraints.adp_similarity(
      u_cart=u_cart, u_iso=u_iso, use_u_aniso=use_u_aniso, weight=1)
    assert approx_equal(a.residual(), expected_residual)

def exercise_isotropic_adp():
  i_seq = 0
  weight = 2
  u_cart = (1,2,3,5,2,8)
  p = adp_restraints.isotropic_adp_proxy(
    i_seq=i_seq,
    weight=weight)
  assert p.i_seq == i_seq
  assert approx_equal(p.weight, weight)
  i = adp_restraints.isotropic_adp(
    u_cart=u_cart,
    weight=weight)
  expected_deltas = (-1, 0, 1, 5, 2, 8)
  expected_gradients = (-4, 0, 4, 40, 16, 64)
  assert approx_equal(i.u_cart, u_cart)
  assert approx_equal(i.weight, weight)
  assert approx_equal(i.deltas(), expected_deltas)
  assert approx_equal(i.rms_deltas(), 4.5704364002673632)
  assert approx_equal(i.residual(), 376.0)
  assert approx_equal(i.gradients(), expected_gradients)
  gradients_aniso_cart = flex.sym_mat3_double(1, (0,0,0,0,0,0))
  proxies = adp_restraints.shared_isotropic_adp_proxy([p,p])
  u_cart = flex.sym_mat3_double((u_cart,))
  residuals = adp_restraints.isotropic_adp_residuals(
    u_cart=u_cart, proxies=proxies)
  assert approx_equal(residuals, (i.residual(),i.residual()))
  deltas_rms = adp_restraints.isotropic_adp_deltas_rms(
    u_cart=u_cart, proxies=proxies)
  assert approx_equal(deltas_rms, (i.rms_deltas(),i.rms_deltas()))
  residual_sum = adp_restraints.isotropic_adp_residual_sum(
    u_cart=u_cart,
    proxies=proxies,
    gradients_aniso_cart=gradients_aniso_cart)
  assert approx_equal(residual_sum, 752.0)
  fd_grads_aniso, fd_grads_iso = finite_difference_gradients(
    restraint_type=adp_restraints.isotropic_adp,
    proxy=p,
    u_cart=u_cart)
  for g,e in zip(gradients_aniso_cart, fd_grads_aniso):
    assert approx_equal(g, matrix.col(e)*2)
  #
  # check frame invariance of residual
  #
  u_cart = matrix.sym(sym_mat3=(0.1,0.2,0.05,0.03,0.02,0.01))
  a = adp_restraints.isotropic_adp(u_cart=u_cart.as_sym_mat3(), weight=1)
  expected_residual = a.residual()
  gen = flex.mersenne_twister()
  for i in range(20):
    R = matrix.rec(gen.random_double_r3_rotation_matrix(),(3,3))
    u_cart_rot = R * u_cart * R.transpose()
    a = adp_restraints.isotropic_adp(u_cart=u_cart_rot.as_sym_mat3(), weight=1)
    assert approx_equal(a.residual(), expected_residual)

def exercise_proxy_show():
  if sys.platform.startswith("win") and sys.version_info[:2] < (2,6):
    # This appears to be a windows-specific bug with string formatting
    # for python versions prior to 2.6, where the exponent is printed
    # with 3 digits rather than 2.
    print "Skipping exercise_proxy_show()"
    return
  sites_cart = flex.vec3_double((
    (-3.1739,10.8317,7.5653),(-2.5419,9.7567,6.6306),
    (-3.3369,8.8794,4.5191),(-3.4640,9.9882,5.3896)))
  site_labels = ("C1", "C2", "O16", "N8")
  u_cart = flex.sym_mat3_double((
    (0.0153,0.0206,0.0234,0.0035,-0.0052,-0.0051),
    (0.0185,0.0109,0.0206,0.0005,-0.0010,0.0002),
    (0.0295,0.0203,0.0218,-0.0010,-0.0003,-0.0044),
    (0.0159,0.0154,0.0206,-0.0003,0.0004,0.0036)))
  u_iso = flex.double((-1,-1,-1,-1))
  use_u_aniso = flex.bool((True,True,True,True))
  #
  proxies = adp_restraints.shared_adp_similarity_proxy()
  sio = StringIO()
  proxies.show_sorted(
    by_value="residual",
    u_cart=flex.sym_mat3_double(),
    u_iso=flex.double(),
    use_u_aniso=flex.bool(),
    f=sio)
  assert not show_diff(sio.getvalue(), """\
ADP similarity restraints: 0
""")
  proxies = adp_restraints.shared_adp_similarity_proxy([
    adp_restraints.adp_similarity_proxy(i_seqs=[0,1],weight=25),
    adp_restraints.adp_similarity_proxy(i_seqs=[2,3],weight=0.3)])
  sio = StringIO()
  proxies.show_sorted(
    by_value="residual",
    u_cart=u_cart,
    u_iso=u_iso,
    use_u_aniso=use_u_aniso,
    f=sio,
    prefix=":")
  assert not show_diff(sio.getvalue(), """\
:ADP similarity restraints: 2
:Sorted by residual:
:scatterers 0
:           1
:         delta    sigma   weight rms_deltas residual
: U11 -3.20e-03 2.00e-01 2.50e+01   4.96e-03 5.54e-03
: U22  9.70e-03 2.00e-01 2.50e+01
: U33  2.80e-03 2.00e-01 2.50e+01
: U12  3.00e-03 2.00e-01 2.50e+01
: U13 -4.20e-03 2.00e-01 2.50e+01
: U23 -5.30e-03 2.00e-01 2.50e+01
:scatterers 2
:           3
:         delta    sigma   weight rms_deltas residual
: U11  1.36e-02 1.83e+00 3.00e-01   6.15e-03 1.02e-04
: U22  4.90e-03 1.83e+00 3.00e-01
: U33  1.20e-03 1.83e+00 3.00e-01
: U12 -7.00e-04 1.83e+00 3.00e-01
: U13 -7.00e-04 1.83e+00 3.00e-01
: U23 -8.00e-03 1.83e+00 3.00e-01
""")
  sio = StringIO()
  proxies.show_sorted(
    by_value="rms_deltas",
    site_labels=site_labels,
    u_cart=u_cart,
    u_iso=u_iso,
    use_u_aniso=use_u_aniso,
    f=sio,
    prefix="=")
  assert not show_diff(sio.getvalue(), """\
=ADP similarity restraints: 2
=Sorted by rms_deltas:
=scatterers O16
=           N8
=         delta    sigma   weight rms_deltas residual
= U11  1.36e-02 1.83e+00 3.00e-01   6.15e-03 1.02e-04
= U22  4.90e-03 1.83e+00 3.00e-01
= U33  1.20e-03 1.83e+00 3.00e-01
= U12 -7.00e-04 1.83e+00 3.00e-01
= U13 -7.00e-04 1.83e+00 3.00e-01
= U23 -8.00e-03 1.83e+00 3.00e-01
=scatterers C1
=           C2
=         delta    sigma   weight rms_deltas residual
= U11 -3.20e-03 2.00e-01 2.50e+01   4.96e-03 5.54e-03
= U22  9.70e-03 2.00e-01 2.50e+01
= U33  2.80e-03 2.00e-01 2.50e+01
= U12  3.00e-03 2.00e-01 2.50e+01
= U13 -4.20e-03 2.00e-01 2.50e+01
= U23 -5.30e-03 2.00e-01 2.50e+01
""")
  #
  proxies = adp_restraints.shared_isotropic_adp_proxy()
  sio = StringIO()
  proxies.show_sorted(
    by_value="residual",
    u_cart=flex.sym_mat3_double(),
    f=sio)
  assert not show_diff(sio.getvalue(), """\
Isotropic ADP restraints: 0
""")
  proxies = adp_restraints.shared_isotropic_adp_proxy([
    adp_restraints.isotropic_adp_proxy(i_seq=0,weight=25),
    adp_restraints.isotropic_adp_proxy(i_seq=2,weight=0.3)])
  sio = StringIO()
  proxies.show_sorted(
    by_value="residual",
    site_labels=site_labels,
    u_cart=u_cart,
    f=sio,
    prefix=" ")
  assert not show_diff(sio.getvalue(), """\
 Isotropic ADP restraints: 2
 Sorted by residual:
 scatterer C1
          delta    sigma   weight rms_deltas residual
  U11 -4.47e-03 2.00e-01 2.50e+01   4.27e-03 4.11e-03
  U22  8.33e-04 2.00e-01 2.50e+01
  U33  3.63e-03 2.00e-01 2.50e+01
  U12  3.50e-03 2.00e-01 2.50e+01
  U13 -5.20e-03 2.00e-01 2.50e+01
  U23 -5.10e-03 2.00e-01 2.50e+01
 scatterer O16
          delta    sigma   weight rms_deltas residual
  U11  5.63e-03 1.83e+00 3.00e-01   3.16e-03 2.69e-05
  U22 -3.57e-03 1.83e+00 3.00e-01
  U33 -2.07e-03 1.83e+00 3.00e-01
  U12 -1.00e-03 1.83e+00 3.00e-01
  U13 -3.00e-04 1.83e+00 3.00e-01
  U23 -4.40e-03 1.83e+00 3.00e-01
""")
  sio = StringIO()
  proxies.show_sorted(
    by_value="rms_deltas",
    u_cart=u_cart,
    f=sio,
    prefix="$")
  assert not show_diff(sio.getvalue(), """\
$Isotropic ADP restraints: 2
$Sorted by rms_deltas:
$scatterer 0
$         delta    sigma   weight rms_deltas residual
$ U11 -4.47e-03 2.00e-01 2.50e+01   4.27e-03 4.11e-03
$ U22  8.33e-04 2.00e-01 2.50e+01
$ U33  3.63e-03 2.00e-01 2.50e+01
$ U12  3.50e-03 2.00e-01 2.50e+01
$ U13 -5.20e-03 2.00e-01 2.50e+01
$ U23 -5.10e-03 2.00e-01 2.50e+01
$scatterer 2
$         delta    sigma   weight rms_deltas residual
$ U11  5.63e-03 1.83e+00 3.00e-01   3.16e-03 2.69e-05
$ U22 -3.57e-03 1.83e+00 3.00e-01
$ U33 -2.07e-03 1.83e+00 3.00e-01
$ U12 -1.00e-03 1.83e+00 3.00e-01
$ U13 -3.00e-04 1.83e+00 3.00e-01
$ U23 -4.40e-03 1.83e+00 3.00e-01
""")
  #
  proxies = adp_restraints.shared_rigid_bond_proxy()
  sio = StringIO()
  proxies.show_sorted(
    by_value="residual",
    sites_cart=flex.vec3_double(),
    u_cart=flex.sym_mat3_double(),
    f=sio)
  assert not show_diff(sio.getvalue(), """\
Rigid bond restraints: 0
""")
  proxies = adp_restraints.shared_rigid_bond_proxy([
    adp_restraints.rigid_bond_proxy(i_seqs=(0,1),weight=25),
    adp_restraints.rigid_bond_proxy(i_seqs=(0,2),weight=15),
    adp_restraints.rigid_bond_proxy(i_seqs=(2,3),weight=25),
    adp_restraints.rigid_bond_proxy(i_seqs=(3,1),weight=30)])
  sio = StringIO()
  proxies.show_sorted(
    by_value="residual",
    sites_cart=sites_cart,
    site_labels=site_labels,
    u_cart=u_cart,
    f=sio,
    prefix="*")
  assert not show_diff(sio.getvalue(), """\
*Rigid bond restraints: 4
*Sorted by residual:
*scatterers O16
*           N8
*   delta_z    sigma   weight residual
*  3.96e-03 2.00e-01 2.50e+01 3.92e-04
*scatterers C1
*           C2
*   delta_z    sigma   weight residual
*  1.08e-03 2.00e-01 2.50e+01 2.89e-05
*scatterers C1
*           O16
*   delta_z    sigma   weight residual
*  4.03e-04 2.58e-01 1.50e+01 2.44e-06
*scatterers N8
*           C2
*   delta_z    sigma   weight residual
*  1.54e-04 1.83e-01 3.00e+01 7.16e-07
""")
  sio = StringIO()
  proxies.show_sorted(
    by_value="delta",
    sites_cart=sites_cart,
    u_cart=u_cart,
    f=sio,
    prefix="||",
    max_items=2)
  assert not show_diff(sio.getvalue(), """\
||Rigid bond restraints: 4
||Sorted by delta:
||scatterers 2
||           3
||   delta_z    sigma   weight residual
||  3.96e-03 2.00e-01 2.50e+01 3.92e-04
||scatterers 0
||           1
||   delta_z    sigma   weight residual
||  1.08e-03 2.00e-01 2.50e+01 2.89e-05
||... (remaining 2 not shown)
""")

def exercise():
  exercise_proxy_show()
  exercise_adp_similarity()
  exercise_isotropic_adp()
  exercise_rigid_bond()
  exercise_rigid_bond_test()
  print "OK"

if (__name__ == "__main__"):
  exercise()
