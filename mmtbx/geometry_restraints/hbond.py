
from __future__ import division
from math import sqrt, cos

sigma_base = sqrt(2.0) / 3.0

class implicit_proxy (object) :
  def __init__ (self,
                i_seqs, # donor, acceptor, acceptor base
                distance_ideal,
                distance_cut,
                theta_low,
                theta_high,
                weight=1.0) :
    self.i_seqs = i_seqs
    self.distance_ideal = distance_ideal
    self.distance_cut = distance_cut
    self.theta_low = theta_low
    self.theta_high = theta_high
    self.weight = weight

def hbond_target_and_gradients (proxies,
                                sites_cart,
                                gradient_array=None) :
  from scitbx.array_family import flex
  if (gradient_array is None) :
    gradient_array = flex.vec3_double(sites_cart.size(), (0,0,0))
  sum = 0.0
  for proxy in proxies :
    if isinstance(proxy, implicit_proxy) :
      sum += _implicit_residue_and_gradients_fd(
        proxy=proxy,
        sites_cart=sites_cart)
    else :
      pass
  return sum

# Fabiola et al. (2002) Protein Sci. 11:1415-23
# http://www.ncbi.nlm.nih.gov/pubmed/12021440
def _implicit_residue_and_gradients_fd (proxy,
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
           cos(delta_low)**4
  if ((distance_cut - distance_ideal) < 0.05) :
    energy *= (distance_cut - distance_ideal) / 0.05
  return energy
