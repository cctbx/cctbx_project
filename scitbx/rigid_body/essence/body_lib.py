import joint_lib
import spatial_lib
matrix = spatial_lib.matrix

class mass_points(object):

  def __init__(O, sites, masses):
    assert len(masses) == len(sites)
    O.sites = sites
    O.masses = masses
    O._sum_of_masses = None
    O._center_of_mass = None

  def sum_of_masses(O):
    if (O._sum_of_masses is None):
      sm = 0
      for mass in O.masses:
        sm += mass
      O._sum_of_masses = sm
    return O._sum_of_masses

  def center_of_mass(O):
    if (O._center_of_mass is None):
      assert len(O.masses) != 0
      assert O.sum_of_masses() != 0
      sms = matrix.col((0,0,0))
      for mass,site in zip(O.masses, O.sites):
        sms += mass * site
      O._center_of_mass = sms / O.sum_of_masses()
    return O._center_of_mass

  def inertia(O, pivot):
    m = [0] * 9
    for mass,site in zip(O.masses, O.sites):
      x,y,z = site - pivot
      m[0] += mass * (y*y+z*z)
      m[4] += mass * (x*x+z*z)
      m[8] += mass * (x*x+y*y)
      m[1] -= mass * x*y
      m[2] -= mass * x*z
      m[5] -= mass * y*z
    m[3] = m[1]
    m[6] = m[2]
    m[7] = m[5]
    return matrix.sqr(m)

  def spatial_inertia(O, alignment_cb_0b=None):
    center_of_mass = O.center_of_mass()
    inertia = O.inertia(pivot=center_of_mass)
    if (alignment_cb_0b is not None):
      center_of_mass = alignment_cb_0b * center_of_mass
      inertia = alignment_cb_0b.r * inertia * alignment_cb_0b.r.transpose()
    return spatial_lib.mci(m=O._sum_of_masses, c=center_of_mass, i=inertia)

def set_cb_tree(bodies):
  "Computes Xtree (RBDA Fig. 4.7, p. 74) for all bodies."
  for body in bodies:
    if (body.parent == -1):
      body.cb_tree = body.alignment.cb_0b
    else:
      body.cb_tree = body.alignment.cb_0b * bodies[body.parent].alignment.cb_b0

class zero_dof(object):

  def __init__(O, sites, masses):
    O.number_of_sites = len(sites)
    O.sum_of_masses = sum(masses)
    O.alignment = joint_lib.zero_dof_alignment()
    O.i_spatial = matrix.sqr([0]*36)
    O.joint = joint_lib.zero_dof()
    O.qd = O.joint.qd_zero

class six_dof(object):

  def __init__(O, sites, masses):
    O.number_of_sites = len(sites)
    mp = mass_points(sites=sites, masses=masses)
    O.sum_of_masses = mp.sum_of_masses()
    O.alignment = joint_lib.six_dof_alignment(
      center_of_mass=mp.center_of_mass())
    O.i_spatial = mp.spatial_inertia(alignment_cb_0b=O.alignment.cb_0b)
    #
    qe = matrix.col((1,0,0,0))
    qr = matrix.col((0,0,0))
    O.joint = joint_lib.six_dof(qe=qe, qr=qr)
    O.qd = O.joint.qd_zero

class spherical(object):

  def __init__(O, sites, masses, pivot):
    O.number_of_sites = len(sites)
    mp = mass_points(sites=sites, masses=masses)
    O.sum_of_masses = mp.sum_of_masses()
    O.alignment = joint_lib.spherical_alignment(pivot=pivot)
    O.i_spatial = mp.spatial_inertia(alignment_cb_0b=O.alignment.cb_0b)
    #
    qe = matrix.col((1,0,0,0))
    O.joint = joint_lib.spherical(qe=qe)
    O.qd = O.joint.qd_zero

class revolute(object):

  def __init__(O, sites, masses, pivot, normal):
    O.number_of_sites = len(sites)
    mp = mass_points(sites=sites, masses=masses)
    O.sum_of_masses = mp.sum_of_masses()
    O.alignment = joint_lib.revolute_alignment(pivot=pivot, normal=normal)
    O.i_spatial = mp.spatial_inertia(alignment_cb_0b=O.alignment.cb_0b)
    #
    O.joint = joint_lib.revolute(qe=matrix.col([0]))
    O.qd = O.joint.qd_zero

class translational(object):

  def __init__(O, sites, masses):
    O.number_of_sites = len(sites)
    mp = mass_points(sites=sites, masses=masses)
    O.sum_of_masses = mp.sum_of_masses()
    O.alignment = joint_lib.translational_alignment(
      center_of_mass=mp.center_of_mass())
    O.i_spatial = mp.spatial_inertia(alignment_cb_0b=O.alignment.cb_0b)
    #
    qr = matrix.col((0,0,0))
    O.joint = joint_lib.translational(qr=qr)
    O.qd = O.joint.qd_zero
