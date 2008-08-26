from scitbx.array_family import flex
from scitbx import matrix
import math
import sys

angle_scale = math.pi / 2

def euler_xyz_matrix(ea):
  """
  Mathematica code:
    rx = {{1, 0, 0}, {0, cx, -sx}, {0, sx, cx}}
    ry = {{cy, 0, sy}, {0, 1, 0}, {-sy, 0, cy}}
    rz = {{cz, -sz, 0}, {sz, cz, 0}, {0, 0, 1}}
    rx.ry.rz
  """
  sin, cos = math.sin, math.cos
  cx = cos(ea[0] * angle_scale)
  sx = sin(ea[0] * angle_scale)
  cy = cos(ea[1] * angle_scale)
  sy = sin(ea[1] * angle_scale)
  cz = cos(ea[2] * angle_scale)
  sz = sin(ea[2] * angle_scale)
  return (
              cy*cz,         -cy*sz,     sy,
     cz*sx*sy+cx*sz, cx*cz-sx*sy*sz, -cy*sx,
    -cx*cz*sy+sx*sz, cz*sx+cx*sy*sz,  cx*cy)

def euler_xyz_ea_d_as_omega_fixed_frame_matrix(ea):
  "Goldstein (A14.xyz) with sinus sign reversed"
  sin, cos = math.sin, math.cos
  cx = cos(ea[0] * angle_scale)
  sx = sin(ea[0] * angle_scale)
  cy = cos(ea[1] * angle_scale)
  sy = sin(ea[1] * angle_scale)
  return (
    1,  0,     sy,
    0, cx, -cy*sx,
    0, sx,  cx*cy)

def newton_euler_f(sites, pivot, d_potential_energy_d_site):
  "Schwieters & Clore (2001) equations 24"
  sum_grads = matrix.col((0,0,0))
  sum_moments = matrix.col((0,0,0))
  for site,grad in zip(sites, d_potential_energy_d_site):
    grad = -matrix.col(grad)
    sum_grads += grad
    sum_moments += (matrix.col(site) - pivot).cross(grad)
  return matrix.col(sum_moments.elems + sum_grads.elems)

class rigid_body(object):

  def __init__(self, sites):
    self.sites_orig = sites
    self.center_of_mass_orig = matrix.col(sites.mean())
    self.lt = matrix.col((0,0,0))
    self.ea = matrix.col((0,0,0))

  def rotation_matrix(self):
    return euler_xyz_matrix(ea=self.ea)

  def center_of_mass_moved(self):
    return self.center_of_mass_orig + self.lt

  def sites_moved(self):
    return \
      self.rotation_matrix() \
      * (self.sites_orig - self.center_of_mass_orig) \
      + self.center_of_mass_moved()

  def ea_gradients(self, energy_cart_function):
    sites_moved = self.sites_moved()
    energy_cart = energy_cart_function(
      nodes=sites_moved, homes=self.sites_orig)
    ne_f = newton_euler_f(
      sites=sites_moved,
      pivot=self.center_of_mass_moved(),
      d_potential_energy_d_site=energy_cart.gradients())
    f = list(-ne_f)
    c = matrix.sqr(euler_xyz_ea_d_as_omega_fixed_frame_matrix(
      ea=self.ea)).transpose()
    return list(angle_scale * c * matrix.col(f[:3])) + f[-3:]

def exercise(args):
  assert len(args) == 0
  sites = flex.vec3_double([
    (10.949, 12.815, 15.189),
    (10.405, 13.954, 15.917),
    (10.779, 15.262, 15.227)])

  class energy_cart(object):

    def __init__(self, nodes, homes):
      assert nodes.size() == homes.size()
      self.nodes = nodes
      self.homes = homes

    def functional(self):
      return flex.sum((self.nodes-self.homes).dot())

    def gradients(self):
      return 2*(self.nodes-self.homes)

  def incr_position(rb, i, delta):
    assert 0 <= i < 6
    if (i < 3):
      v = list(rb.ea)
      v[i] += delta
      rb.ea = matrix.col(v)
    else:
      v = list(rb.lt)
      v[i-3] += delta
      rb.lt = matrix.col(v)

  def ea_gradients_fd(rb, energy_cart_function, eps=1.e-6):
    result = []
    for i in xrange(6):
      fs = []
      incr_position(rb=rb, i=i, delta=eps)
      fs.append(energy_cart_function(
        nodes=rb.sites_moved(), homes=rb.sites_orig).functional())
      incr_position(rb=rb, i=i, delta=-eps)
      incr_position(rb=rb, i=i, delta=-eps)
      fs.append(energy_cart_function(
        nodes=rb.sites_moved(), homes=rb.sites_orig).functional())
      incr_position(rb=rb, i=i, delta=eps)
      result.append((fs[0]-fs[1])/(2*eps))
    return result

  def show_gradients(rb):
    an = rb.ea_gradients(energy_cart_function=energy_cart)
    fd = ea_gradients_fd(rb=rb, energy_cart_function=energy_cart)
    print "an ea:", an[:3]
    print "fd ea:", fd[:3]
    print "an lt:", an[3:]
    print "fd lt:", fd[3:]
    print

  rb = rigid_body(sites=sites)
  mt = flex.mersenne_twister()
  n_trials = 4
  for i in xrange(n_trials):
    rb.ea = matrix.col(mt.random_double_point_on_sphere()) * i
    rb.lt = matrix.col(mt.random_double_point_on_sphere()) * i
    show_gradients(rb=rb)

  print "OK"

if (__name__ == "__main__"):
  exercise(sys.argv[1:])
