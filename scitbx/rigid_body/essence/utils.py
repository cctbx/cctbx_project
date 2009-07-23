import featherstone
matrix = featherstone.matrix

class mass_points(object):

  def __init__(O, masses, sites):
    assert len(masses) == len(sites)
    O.masses = masses
    O.sites = sites
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
      sm = 0
      sms = matrix.col((0,0,0))
      for mass,site in zip(O.masses, O.sites):
        sm += mass
        sms += mass * site
      O._sum_of_masses = sm
      assert O._sum_of_masses != 0
      O._center_of_mass = sms / sm
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
    return featherstone.mcI(m=O._sum_of_masses, c=center_of_mass, I=inertia)
