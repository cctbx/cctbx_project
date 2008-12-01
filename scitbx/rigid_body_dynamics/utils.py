from scitbx.rigid_body_dynamics import featherstone
from scitbx import matrix

def center_of_mass_from_sites(sites):
  assert len(sites) != 0
  result = matrix.col((0,0,0))
  for site in sites:
    result += site
  result /= len(sites)
  return result

def inertia_from_sites(sites, pivot):
  m = [0] * 9
  for site in sites:
    x,y,z = site - pivot
    m[0] += y*y+z*z
    m[4] += x*x+z*z
    m[8] += x*x+y*y
    m[1] -= x*y
    m[2] -= x*z
    m[5] -= y*z
  m[3] = m[1]
  m[6] = m[2]
  m[7] = m[5]
  return matrix.sqr(m)

def spatial_inertia_from_sites(
      sites,
      mass=None,
      center_of_mass=None,
      alignment_T=None):
  if (mass is None):
    mass = len(sites)
  if (center_of_mass is None):
    center_of_mass = center_of_mass_from_sites(sites=sites)
  inertia = inertia_from_sites(sites=sites, pivot=center_of_mass)
  if (alignment_T is not None):
    center_of_mass = alignment_T * center_of_mass
    inertia = alignment_T.r * inertia * alignment_T.r.transpose()
  return featherstone.mcI(m=mass, c=center_of_mass, I=inertia)

def kinetic_energy(I_spatial, v_spatial):
  "RBDA Eq. 2.67"
  return 0.5 * v_spatial.dot(I_spatial * v_spatial)

def T_as_X(Tps):
  return featherstone.Xrot(Tps.r) \
       * featherstone.Xtrans(-Tps.r.transpose() * Tps.t)

class featherstone_system_model(object):

  def __init__(model, bodies):
    model.NB = len(bodies)
    model.pitch = []
    model.parent =[]
    model.Ttree = []
    model.Xtree = []
    model.I = []
    for B in bodies:
      model.pitch.append(B.J)
      model.parent.append(B.parent)
      if (B.parent == -1):
        Ttree = B.A.T0b
      else:
        Ttree = B.A.T0b * bodies[B.parent].A.Tb0
      model.Ttree.append(Ttree)
      model.Xtree.append(T_as_X(Ttree))
      model.I.append(B.I)
