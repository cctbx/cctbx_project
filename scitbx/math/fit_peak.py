from __future__ import absolute_import, division, print_function
import math
from scitbx import matrix,lbfgs
from scitbx.array_family import flex
from six.moves import range

class fit_peak(object):
  """
  =============================================================================
  fit_peak Class

  This class fits a list of points (x,y,z) and their heights to a certain shape
  using the LBFGS minimizer and the function sum{(p_i - p_i0)^2}.

  -----------------------------------------------------------------------------
  Arguments:
    height_list - the list of heights (list of floats)
    xyz_list - the list of positions (list of tuples)
    shape - flag for shape to be fitted (parabola, quadratic, gaussian)
    max_iterations - the maximum number of minimization steps (int)

  Useful accessible attributes:
    self.x - the final parameters from the fit (list)
    self.vertex - the peak position (tuple)

  Notes:
    The first point should be the point closest to the peak.  If the peak
    is found outside the 26 nearest neighbors, the center point is returned.
    The order of the points is very important because the initial guess for
    minimization is obtained by solving a linear system.  See
    pick_map_neighbors for the required order.
    To use "gaussian," height_list should actually be the ln(height) since the
    least squares function is sum{[ln(p_i) - ln(p_i0)]^2} for that case.
  -----------------------------------------------------------------------------
  """
  def __init__(self,height_list=None,xyz_list=None,shape="parabola",
               max_iterations=25):

    self.height_list = height_list
    self.xyz_list = xyz_list

    # pick shape
    if (shape == "parabola"):
      self.n = parabola().n
      self.fit = parabola
      self.check_lengths(height_list=height_list,xyz_list=xyz_list)

      # construct vector and matrix for initial guess
      b_elements = []
      A_elements = []
      for i in range(self.n):
        b_elements.append(height_list[i])
        A_elements.append(xyz_list[i][0]*xyz_list[i][0])
        A_elements.append(xyz_list[i][1]*xyz_list[i][1])
        A_elements.append(xyz_list[i][2]*xyz_list[i][2])
        A_elements.append(xyz_list[i][0])
        A_elements.append(xyz_list[i][1])
        A_elements.append(xyz_list[i][2])
        A_elements.append(1.0)
      A = matrix.sqr(A_elements)
      b = matrix.col(b_elements)

    elif ((shape == "quadratic") or (shape == "gaussian")):
      self.n = quadratic().n
      if (shape == "quadratic"):
        self.fit = quadratic
      else:
        self.fit = gaussian
      self.check_lengths(height_list=height_list,xyz_list=xyz_list)

      # construct vector and matrix for initial guess
      b_elements = []
      A_elements = []
      for i in range(self.n):
        b_elements.append(height_list[i])
        A_elements.append(xyz_list[i][0]*xyz_list[i][0])
        A_elements.append(xyz_list[i][1]*xyz_list[i][1])
        A_elements.append(xyz_list[i][2]*xyz_list[i][2])
        A_elements.append(xyz_list[i][0])
        A_elements.append(xyz_list[i][1])
        A_elements.append(xyz_list[i][2])
        A_elements.append(xyz_list[i][0]*xyz_list[i][1])
        A_elements.append(xyz_list[i][0]*xyz_list[i][2])
        A_elements.append(xyz_list[i][1]*xyz_list[i][2])
        A_elements.append(1.0)
      A = matrix.sqr(A_elements)
      b = matrix.col(b_elements)

    else:
      print("fit_peak error: shape not valid")
      exit()

    # get initial guess (solve Ax = b)
    self.x = flex.double(A.inverse()*b)

    # if there are more points, minimize using the LBFGS minimizer
    if (len(height_list) > self.n):
      self.minimizer = lbfgs.run(\
        target_evaluator=self,termination_params=\
        lbfgs.termination_parameters(max_iterations=max_iterations))

    # finalize fit
    answer = self.fit(parameters=self.x)
    self.vertex = answer.vertex
    x_max = xyz_list[1][0]  # these min and max values are based on the order
    x_min = xyz_list[2][0]
    y_max = xyz_list[3][1]
    y_min = xyz_list[4][1]
    z_max = xyz_list[5][2]
    z_min = xyz_list[6][2]
    assert(x_max > x_min)
    assert(y_max > y_min)
    assert(z_max > z_min)
    if ((self.vertex[0] < x_min) or (self.vertex[0] > x_max)):
      self.vertex[0] = xyz_list[0][0]
    if ((self.vertex[1] < y_min) or (self.vertex[1] > y_max)):
      self.vertex[1] = xyz_list[0][1]
    if ((self.vertex[2] < z_min) or (self.vertex[2] > z_max)):
      self.vertex[2] = xyz_list[0][2]

  def compute_functional_and_gradients(self):
    answer = self.fit(parameters=self.x)
    f = 0.0
    g = [0.0 for i in range(self.n)]
    for i in range(len(self.height_list)):
      f_i = answer.get_height(r=self.xyz_list[i]) - self.height_list[i]
      g_i = answer.get_gradient(r=self.xyz_list[i])
      f = f + f_i*f_i
      for j in range(self.n):
        g[j] = g[j] + 2.0*f_i*g_i[j]
    return f,flex.double(g)

  def check_lengths(self,height_list=None,xyz_list=None):
    assert(len(height_list) >= self.n)
    assert(len(height_list) == len(xyz_list))

class parabola(object):
  """
  =============================================================================
  parabola Class

  This class models a parabola according to the equation,

    p = a x^2 + b y^2 + c z^2 + d x + e y + f z + g

  -----------------------------------------------------------------------------
  Arguments:
    parameters - the parameters a-g (list)

  Accesible methods:
    get_height(r) - returns p for a given r (x,y,z) (float)
    get_gradient(r) - returns the gradient of p with respect to a-g for a given
                      r (x,y,z) (tuple)

  Notes:
    This class is used in conjunction with fit_peak.  Since there are 7
    parameters, at least 7 points are needed to fit the parabola.
  -----------------------------------------------------------------------------
  """
  def __init__(self,parameters=None):
    self.n = 7
    self.p = None
    if (parameters is not None):
      assert(len(parameters) == self.n)
      self.p = parameters

      # find vertex
      if (round(self.p[0],20) == 0.0):
        rx_v = 1.0e100
      else:
        rx_v = -0.5*(self.p[3]/self.p[0])
      if (round(self.p[1],20) == 0.0):
        ry_v = 1.0e100
      else:
        ry_v = -0.5*(self.p[4]/self.p[1])
      if (round(self.p[2],20) == 0.0):
        rz_v = 1.0e100
      else:
        rz_v = -0.5*(self.p[5]/self.p[2])
      self.vertex = [rx_v,ry_v,rz_v]

  def get_height(self,r=None):
    if (self.p is None):
      print("parabaola error: parameters not initialized")
    else:
      return self.p[0]*r[0]*r[0] + self.p[1]*r[1]*r[1] + self.p[2]*r[2]*r[2] +\
             self.p[3]*r[0] + self.p[4]*r[1] + self.p[5]*r[2] + self.p[6]

  def get_gradient(self,r=None):
    return ( r[0]*r[0], r[1]*r[1], r[2]*r[2], r[0], r[1], r[2], 1.0 )

class quadratic(object):
  """
  =============================================================================
  quadratic Class

  This class models a parabola according to the equation,

    p = a x^2 + b y^2 + c z^2 + d x + e y + f z + g x y + h x z + i y z + j

  -----------------------------------------------------------------------------
  Arguments:
    parameters - the parameters a-j (list)

  Accesible methods:
    get_height(r) - returns p for a given r (x,y,z) (float)
    get_gradient(r) - returns the gradient of p with respect to a-j for a given
                      r (x,y,z) (tuple)

  Notes:
    This class is used in conjunction with fit_peak.  Since there are 10
    parameters, at least 10 points are needed to fit the quadratic.
  -----------------------------------------------------------------------------
  """
  def __init__(self,parameters=None):
    self.n = 10
    self.p = None
    if (parameters is not None):
      assert(len(parameters) == self.n)
      self.p = parameters

      # find vertex
      A = matrix.sqr([2.0*self.p[0],     self.p[6],     self.p[7],
                          self.p[6], 2.0*self.p[1],     self.p[8],
                          self.p[7],     self.p[8], 2.0*self.p[2]])
      b = matrix.col([-(self.p[3]),
                      -(self.p[4]),
                      -(self.p[5])])
      self.vertex = list(A.inverse()*b)

  def get_height(self,r=None):
    if (self.p is None):
      print("quadratic error: parameters not initialized")
    else:
      return self.p[0]*r[0]*r[0] + self.p[1]*r[1]*r[1] + self.p[2]*r[2]*r[2] +\
             self.p[3]*r[0]      + self.p[4]*r[1]      + self.p[5]*r[2]      +\
             self.p[6]*r[0]*r[1] + self.p[7]*r[0]*r[2] + self.p[8]*r[1]*r[2] +\
             self.p[9]

  def get_gradient(self,r=None):
    return ( r[0]*r[0], r[1]*r[1], r[2]*r[2], r[0], r[1], r[2],
             r[0]*r[1], r[0]*r[2], r[1]*r[2], 1.0 )

class gaussian(object):
  """
  =============================================================================
  gaussian Class

  This class models a parabola according to the equation,

    p = exp(a x^2 + b y^2 + c z^2 + d x + e y + f z + g x y + h x z + i y z + j)

  -----------------------------------------------------------------------------
  Arguments:
    parameters - the parameters a-j (list)

  Accesible methods:
    get_height(r) - returns ln(p) for a given r (x,y,z) (float)
    get_gradient(r) - returns the gradient of ln(p) with respect to a-j for a
                      given r (x,y,z) (tuple)

  Notes:
    This class is used in conjunction with fit_peak.  Since there are 10
    parameters, at least 10 points are needed to fit the quadratic.
    The positions (x,y,z) should be fractional or the exponent may return
    an overflow exception.
  -----------------------------------------------------------------------------
  """
  def __init__(self,parameters=None):
    self.n = 10
    self.p = None
    if (parameters is not None):
      assert(len(parameters) == self.n)
      self.p = parameters

      # find vertex
      A = matrix.sqr([2.0*self.p[0],     self.p[6],     self.p[7],
                          self.p[6], 2.0*self.p[1],     self.p[8],
                          self.p[7],     self.p[8], 2.0*self.p[2]])
      b = matrix.col([-(self.p[3]),
                      -(self.p[4]),
                      -(self.p[5])])
      self.vertex = list(A.inverse()*b)

  def get_height(self,r=None):
    if (self.p is None):
      print("gaussian error: parameters not initialized")
    else:
      return self.p[0]*r[0]*r[0] + self.p[1]*r[1]*r[1] + self.p[2]*r[2]*r[2] +\
             self.p[3]*r[0]      + self.p[4]*r[1]      + self.p[5]*r[2]      +\
             self.p[6]*r[0]*r[1] + self.p[7]*r[0]*r[2] + self.p[8]*r[1]*r[2] +\
             self.p[9]

  def get_gradient(self,r=None):
    h = 1.0/math.exp(self.get_height(r=r))
    return ( r[0]*r[0]*h, r[1]*r[1]*h, r[2]*r[2]*h, r[0]*h, r[1]*h, r[2]*h,
             r[0]*r[1]*h, r[0]*r[2]*h, r[1]*r[2]*h, 1.0*h )

class pick_map_neighbors(object):
  """
  =============================================================================
  pick_map_neighbors Class

  This class orders the nearest neighbors in a 3x3x3 box of gridpoints from
  nearest to farthest.

  -----------------------------------------------------------------------------
  Arguments:
    site - the center of the 3x3x3 box (tuple)
    map_in - the map of gridpoints (flex.double)

  Accessible methods:
    get_6_nearest_neighbors - returns site and 6 nearest neighbors (list of
                              floats, list of tuples)
    get_18_nearest_neighbors - returns site and 18 nearest neighbors (list of
                               floats, list of tuples)
    get_26_nearest_neighbors - returns site and 26 nearest neighbors (list of
                               floats, list of tuples)

  Notes:
    This class works in conjunction with fit_peak and supplies points in the
    proper order for that class to function.
  -----------------------------------------------------------------------------
  """
  def __init__(self,site=None,map_in=None):

    assert(site is not None)
    assert(map_in is not None)

    n_neighbors = 27
    extents = map_in.all()
    self.height_list = [0.0 for i in range(n_neighbors)]
    self.xyz_list = [(0,0,0) for i in range(n_neighbors)]

    u0 = site[0]
    up = site[0] + 1
    um = site[0] - 1
    v0 = site[1]
    vp = site[1] + 1
    vm = site[1] - 1
    w0 = site[2]
    wp = site[2] + 1
    wm = site[2] - 1

    # add initial point
    self.height_list[0] = map_in[site]
    self.xyz_list[0] = (site[0],site[1],site[2])

    # add 6 nearest neighbors (1 gridpoint away)
    self.height_list[1] = map_in[(up%extents[0],v0           ,w0           )]
    self.height_list[2] = map_in[(um%extents[0],v0           ,w0           )]
    self.height_list[3] = map_in[(u0           ,vp%extents[1],w0           )]
    self.height_list[4] = map_in[(u0           ,vm%extents[1],w0           )]
    self.height_list[5] = map_in[(u0           ,v0           ,wp%extents[2])]
    self.height_list[6] = map_in[(u0           ,v0           ,wm%extents[2])]

    self.xyz_list[1] = (up,v0,w0)
    self.xyz_list[2] = (um,v0,w0)
    self.xyz_list[3] = (u0,vp,w0)
    self.xyz_list[4] = (u0,vm,w0)
    self.xyz_list[5] = (u0,v0,wp)
    self.xyz_list[6] = (u0,v0,wm)

    # add 12 next nearest neighbors (sqrt(2) gridpoint away)
    self.height_list[7]  = map_in[(u0           ,vp%extents[1],wp%extents[2])]
    self.height_list[8]  = map_in[(um%extents[0],v0           ,wm%extents[2])]
    self.height_list[9]  = map_in[(up%extents[0],vm%extents[1],w0           )]
    self.height_list[10] = map_in[(u0           ,vm%extents[1],wp%extents[2])]
    self.height_list[11] = map_in[(up%extents[0],v0           ,wm%extents[2])]
    self.height_list[12] = map_in[(um%extents[0],vm%extents[1],w0           )]
    self.height_list[13] = map_in[(u0           ,vp%extents[1],wm%extents[2])]
    self.height_list[14] = map_in[(um%extents[0],v0           ,wp%extents[2])]
    self.height_list[15] = map_in[(up%extents[0],vp%extents[1],w0           )]
    self.height_list[16] = map_in[(u0           ,vm%extents[1],wm%extents[2])]
    self.height_list[17] = map_in[(up%extents[0],v0           ,wp%extents[2])]
    self.height_list[18] = map_in[(um%extents[0],vp%extents[1],w0           )]

    self.xyz_list[7]  = (u0,vp,wp)
    self.xyz_list[8]  = (um,v0,wm)
    self.xyz_list[9]  = (up,vm,w0)
    self.xyz_list[10] = (u0,vm,wp)
    self.xyz_list[11] = (up,v0,wm)
    self.xyz_list[12] = (um,vm,w0)
    self.xyz_list[13] = (u0,vp,wm)
    self.xyz_list[14] = (um,v0,wp)
    self.xyz_list[15] = (up,vp,w0)
    self.xyz_list[16] = (u0,vm,wm)
    self.xyz_list[17] = (up,v0,wp)
    self.xyz_list[18] = (um,vp,w0)

    # add final 8 corners (sqrt(3) gridpoint away)
    self.height_list[19] = map_in[(up%extents[0],vp%extents[1],wp%extents[2])]
    self.height_list[20] = map_in[(um%extents[0],vp%extents[1],wp%extents[2])]
    self.height_list[21] = map_in[(up%extents[0],vm%extents[1],wp%extents[2])]
    self.height_list[22] = map_in[(up%extents[0],vp%extents[1],wm%extents[2])]
    self.height_list[23] = map_in[(um%extents[0],vm%extents[1],wp%extents[2])]
    self.height_list[24] = map_in[(um%extents[0],vp%extents[1],wm%extents[2])]
    self.height_list[25] = map_in[(up%extents[0],vm%extents[1],wm%extents[2])]
    self.height_list[26] = map_in[(um%extents[0],vm%extents[1],wm%extents[2])]

    self.xyz_list[19] = (up,vp,wp)
    self.xyz_list[20] = (um,vp,wp)
    self.xyz_list[21] = (up,vm,wp)
    self.xyz_list[22] = (up,vp,wm)
    self.xyz_list[23] = (um,vm,wp)
    self.xyz_list[24] = (um,vp,wm)
    self.xyz_list[25] = (up,vm,wm)
    self.xyz_list[26] = (um,vm,wm)

    # fractionalize coordinates
    for i in range(n_neighbors):
      self.xyz_list[i] = (float(self.xyz_list[i][0])/extents[0],
                          float(self.xyz_list[i][1])/extents[1],
                          float(self.xyz_list[i][2])/extents[2])

  def get_6_nearest_neighbors(self):
    return self.height_list[0:7], self.xyz_list[0:7]

  def get_18_nearest_neighbors(self):
    return self.height_list[0:19], self.xyz_list[0:19]

  def get_26_nearest_neighbors(self):
    return self.height_list, self.xyz_list

# =============================================================================
if(__name__ == "__main__"):
  p = parabola(parameters=(-1,-1,-1,2,2,2,24))
  print(p.vertex)
