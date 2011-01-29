
from __future__ import division
from math import sqrt, cos
from libtbx import adopt_init_args

sigma_base = sqrt(2.0) / 3.0

master_phil = libtbx.phil.parse("""
  restraints_weight = 1.0
    .type = float
  implicit
    .short_caption = Implicit hydrogens
    .help = Based on H-bond potential for CNS by Chapman lab
  {
    distance_ideal = 2.9
      .type = float
    distance_cut = 3.5
      .type = float
    theta_high = 155
      .type = float
    theta_low = 115
      .type = float
  }
  explicit
    .short_caption = Explicit hydrogens
    .help = Similar to Rosetta H-bond energy (Kortemme & Baker)
  {
    distance_ideal = 1.975
      .type = float
    distance_cut = None
      .type = float
    distance_sigma = 0.05
      .type = float
    theta_ideal = 180
      .type = float
    theta_sigma = 5
      .type = float
    psi_ideal = 155
      .type = float
    psi_sigma = 5
      .type = float
    relative_weights = 1.0 1.0 1.0
      .type = floats(size=3)
  }
  simple
    .short_caption = Simple distance-based potential
    .help = Pseudo-bond restraints
  {
    distance_ideal_h_o = 1.975
      .type = float
    distance_cut_h_o = 2.5
      .type = float
    distance_ideal_n_o = 2.9
      .type = float
    distance_cut_n_o = 3.5
      .type = float
    sigma = 0.05
      .type = float
    slack = 0.0
      .type = float
  }
""")

class implicit_proxy (object) :
  def __init__ (self,
                i_seqs, # donor, acceptor, acceptor base
                distance_ideal,
                distance_cut,
                theta_low,
                theta_high,
                weight=1.0) :
    assert (len(i_seqs) == 3))
    assert (distance_cut is None) or (distance_cut > distance_ideal)
    adopt_init_args(self, locals())

class explicit_proxy (object) :
  def __init__ (self,
                i_seqs, # donor, H, acceptor, acceptor base
                distance_ideal,
                distance_cut,
                theta_ideal,
                psi_ideal,
                weight=1.0,
                relative_weights=(1.0,1.0,1.0)) :
    assert (len(relative_weights) == 3)
    assert (len(i_seqs) == 4)
    assert (distance_cut is None) or (distance_cut > distance_ideal)
    adopt_init_args(self, locals())

class distance_proxy (object) :
  def __init__ (self,
                i_seqs,
                distance_ideal,
                distance_cut,
                sigma,
                slack,
                weight=1.0) :
    assert (len(i_seqs) == 2)
    assert (distance_cut is None) or (distance_cut > distance_ideal)
    adopt_init_args(self, locals())
    import cctbx.geometry_restraints
    self.cctbx_proxy = cctbx.geometry_restraints.bond_simple_proxy(
      i_seqs=i_seqs,
      distance_ideal=distance_ideal,
      weight=1/(sigma**2))

def hbond_target_and_gradients (proxies,
                                sites_cart,
                                gradient_array=None) :
  from scitbx.array_family import flex
  if (gradient_array is None) :
    gradient_array = flex.vec3_double(sites_cart.size(), (0,0,0))
  sum = 0.0
  for proxy in proxies :
    if isinstance(proxy, implicit_proxy) :
      sum += _implicit_target_and_gradients_fd(
        proxy=proxy,
        sites_cart=sites_cart)
    elif isinstance(proxy, explicit_proxy) :
      pass
    elif isinstance(proxy, distance_proxy) :
      pass
    else :
      assert False, "Not a recognized H-bond proxy type"
  return sum

# Fabiola et al. (2002) Protein Sci. 11:1415-23
# http://www.ncbi.nlm.nih.gov/pubmed/12021440
def _implicit_target_and_gradients_fd (proxy,
                                       sites_cart,
                                       gradients_array) :
  from cctbx.geometry import angle
  weight = proxy.weight
  i_seqs = proxy.i_seqs
  sites = (sites_cart[i_seqs[0]], sites_cart[i_seqs[1]], sites_cart[i_seqs[2]])
  hb_angle = angle(sites)
  #d_theta_d_xyz = hb_angle.d_angle_d_sites()
  delta_high = proxy.theta_high - hb_angle.angle_model
  delta_low = proxy.theta_low - hb_angle.angle_model
  if (abs(delta_high) < abs(delta_low)) :
    angle_ideal = proxy.theta_high
    delta_theta = delta_high
  else :
    angle_ideal = proxy.theta_low
    delta_theta = delta_low
  residual = _eval_energy_implicit(
    sites=sites,
    distance_ideal=distance_ideal,
    weight=proxy.weight,
    delta_theta=delta_theta)
  for j in range(3) :
    i_seq = i_seqs[j]
    grads_j = gradients_array[i_seq]
    for k in range(3) :
      sites[j][k] -= epsilon
      hb_angle = angle(sites)
      delta_theta = angle_ideal - hb_angle.angle_model
      e1 = _eval_energy_implicit(
        sites=sites,
        distance_ideal=distance_ideal,
        weight=proxy.weight,
        delta_theta=delta_theta)
      sites[j][k] += 2 * epsilon
      hb_angle = angle(sites)
      delta_theta = angle_ideal - hb_angle.angle_model
      e2 = _eval_energy_implicit(
        sites=sites,
        distance_ideal=distance_ideal,
        weight=proxy.weight,
        delta_theta=delta_theta)
      grad_j[k] += (e2 - e1) / epsilon
      sites[j][k] -= epsilon
    gradients_array[i_seq] = grads_j
  return residual

def _eval_energy_implicit (sites,
                           distance_ideal,
                           distance_cut,
                           weight,
                           delta_theta) :
  if (distance_ideal > distance_cut) :
    return 0.0
  from scitbx.matrix import rec
  v1 = rec(sites[0], (1,3))
  v2 = rec(sites[1], (1,3))
  bond_dist = abs(v2 - v1)
  sigma = sigma_base * distance_ideal
  energy = weight * ((sigma / bond_dist)**6 - (sigma / bond_dist)**4) * \
           cos(delta_theta)**4
  if ((distance_cut - distance_ideal) < 0.05) :
    energy *= (distance_cut - distance_ideal) / 0.05
  return energy

def _expicit_target_and_gradients_fd (proxy,
                                      sites_cart,
                                      gradients_array) :
  assert isinstance(proxy, explicit_proxy)
