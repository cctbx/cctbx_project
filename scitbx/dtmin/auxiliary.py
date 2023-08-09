from __future__ import division
from math import log, exp, tanh
import sys

class Auxiliary:
  def get_reparameterized_parameters(self):
    "Return the vector of parameters with necessary reparameterisations performed"
    repar = self.reparameterize()
    reparameterized_x = self.get_macrocycle_parameters()[:] #the slicing is to ensure a deep copy

    if len(repar) != 0:
      for i in range(self.nmp):
        if repar[i].reparamed:
          reparameterized_x[i] = log(reparameterized_x[i] + repar[i].offset)

    return reparameterized_x

  def set_reparameterized_parameters(self, reparameterized_x):
    """Given reparameterized_x, find the original x and call set_macrocycle_parameters(x) set x in refineXYZ"""
    x = reparameterized_x[:] #the slicing is to ensure a deep copy
    repar = self.reparameterize()

    if len(repar) != 0:
      for i in range(len(repar)):
        if repar[i].reparamed:
          x[i] = exp(x[i]) - repar[i].offset

    self.set_macrocycle_parameters(x)

  def reparameterized_large_shifts(self):
    """
    derivation of the formula used:
    1d unreparameterized:
                            ls
                          |------>|
        0                 x       |
      --|-----------------|-------|---->

    reparameterized:
                rls
              |---->|
        0    rx     |
      --|-----|-----|------------------>

      rx = log(x + offset)

      rx + rls = log(x + ls + offset)

      => rls = log(x + ls + offset) - log(x + offset)

            = log(1 + ls/(x + offset) )

            ~= ls/(x + offset) for small ls/(x + offset)
    """
    x = self.get_macrocycle_parameters()
    reparameterized_ls = self.macrocycle_large_shifts()
    repar = self.reparameterize()

    if len(repar) != 0:
      for i in range(self.nmp):
        if repar[i].reparamed:
          reparameterized_ls[i] = reparameterized_ls[i]/(x[i]+repar[i].offset)

    return reparameterized_ls

  def reparameterized_bounds(self):
    """Returns the reparameterized bounds"""
    unrepar_bounds = self.bounds()
    reparameterized_bounds = self.bounds()
    repar = self.reparameterize()

    if len(reparameterized_bounds) == 0:  #no bounds case
      return reparameterized_bounds
    elif len(repar) == 0:  #bounds but no repar case
      return reparameterized_bounds
    else:  #bounds and repar case
      for i in range(self.nmp):
        if repar[i].reparamed:
          if reparameterized_bounds[i].lower_bounded():
            reparameterized_bounds[i].lower_on(log(unrepar_bounds[i].lower_limit() + repar[i].offset))
          if reparameterized_bounds[i].upper_bounded():
            reparameterized_bounds[i].upper_on(log(unrepar_bounds[i].upper_limit() + repar[i].offset))
      return reparameterized_bounds

  def reparameterized_target_gradient(self):
    r"""
    Returns the function and the reparameterized gradient

    When we are finding the gradient of f with respect to the reparameterized
    parameters rx, we can pretend that we started from rx, calculated x
    and then calculated f. In reality we started from x then calculated rx.

    x -> rx  (reparameterized x)
     \-> f

    i.e.  ' rx -> x -> f '

    rx(x) = log(x+c)  => x(rx) = e^rx - c
                        dx/drx = e^rx = x + c

    we want:  df(x(rx))   df   dx    df
              --------- = -- * --- = -- * (x+c)
                  drx      dx   drx   dx
    """
    x = self.get_macrocycle_parameters()
    reparameterized_f_g = self.target_gradient()
    repar = self.reparameterize()

    if len(repar) != 0:
      for i in range(self.nmp):
        if repar[i].reparamed:
          reparameterized_f_g[1][i] *= (x[i]+repar[i].offset)

    return reparameterized_f_g # dont need to do anything to f

  def reparameterized_target_gradient_hessian(self):
    """
    Returns the function, the reparameterized gradient,
    the reparameterized hessian and a flag to tell if
    the hessian is diagonal

     Probably easiest to understand if you write out the equation,
       all small d's are partial differentials
       rx means reparameterized x
       c_i is the reparameterisation offset used in rx_i(x_i) = log(x_i + c_i)

                        d^2 f         d      df     dx_j
    reparedHessian_ij = ----------- = ----- (---- * -----)
                        drx_i drx_j   drx_i  dx_j   drx_j

      d^2 f      dx_j     df   d^2x_j
    = ----------*-----  + ----*-----------
      drx_i dx_j drx_j    dx_j drx_i drx_j

      d^2 f     dx_i  dx_j    df   d^2 x_j
    = ---------*-----*----- + ----*-----------
      dx_i dx_j drx_i drx_j   dx_j drx_i drx_j

      d^2f                                         df
    = --------- * (x_i + c_i) * (x_j + c_j) + diag(----*(x_j + c_j) )
      dx_i dx_j                                    dx_j
    """

    x = self.get_macrocycle_parameters()
    repar = self.reparameterize()

    f = None
    reparameterized_g = None # will be populated with the unrepared gradient and ones requiring reparing will be modified accordingly
    reparameterized_h = None # will be populated with the unrepared hessian, and ones requiring reparing will be modified accordingly
    is_hessian_diagonal = None
    (f, reparameterized_g, reparameterized_h, is_hessian_diagonal) = self.target_gradient_hessian() # Note reparameterized_g and reparameterized_h are not repared yet, that is what this function does

    if len(repar) != 0:
      # repararameterize the hessian first because we need the unreparameterized gradient for its calculation
      for i in range(self.nmp):
        for j in range(self.nmp):
          if repar[i].reparamed or repar[j].reparamed:
            dxidyi,dxjdyj = 1., 1.
            if repar[i].reparamed: dxidyi = x[i]+repar[i].offset
            if repar[j].reparamed: dxjdyj = x[j]+repar[j].offset
            reparameterized_h[i,j] *= dxidyi*dxjdyj
            if i==j: reparameterized_h[i,i] += reparameterized_g[i]*dxidyi #dxidyi=d2xidyi2  # note at this point 'reparameterized_g' is not actually repared
      # repar the gradient as per reparameterized_target_gradient()
      for i in range(self.nmp):
        if repar[i].reparamed:
          reparameterized_g[i] *= (x[i]+repar[i].offset)

    return (f, reparameterized_g, reparameterized_h, is_hessian_diagonal)

  def reparameterized_adjusted_target_gradient_hessian(self):
    """
    Returns the function, the reparameterized gradient,
    the reparameterized and adjusted* hessian, a flag to tell is the hessian is diagonal
    and a flag to say whether the hessian was adjusted or not.

    *By adjusted we mean that the hessian has been through the adjust_hessian function
    """
    (f, reparameterized_g, reparameterized_h, is_hessian_diagonal) = self.reparameterized_target_gradient_hessian() # note reparameterized_g and reparameterized_h are not repared yet, that is what this function is for!

    # adjust the hessian
    (reparameterized_h, queue_new_hessian) = self.adjust_hessian(reparameterized_h)

    return (f, reparameterized_g, reparameterized_h, is_hessian_diagonal, queue_new_hessian)

  def damped_shift(self, a, old_x, p, dist, ls):
    """
    The overall idea is to shift x by a multiple of p, but avoid crossing boundaries
    and dampen any shifts of parameters that would be greater than their large_shift value.
    """

    x = [0.] * self.nmp
    bounds = self.reparameterized_bounds()
    for i in range(self.nmp):
      thisdist = a if (dist[i]<0.) else min(a,dist[i]) # neg dist (-1.) means outside of bound (should never happen) or unbounded.
      rawshift = thisdist*p[i] # Undamped shift
      relshift = abs(rawshift)/ls[i]
      dampfac = 1. if (relshift <= 0.005) else 2.*tanh(relshift/2.)/relshift
      shift = dampfac*rawshift
      x[i] = old_x[i] + shift # Shift parameter only up to its bound, damped to max of 2*largeShift

    if len(bounds) != 0:
      # correct for the possibility of minutely stepping over a boundary
      for i in range(self.nmp):
        if bounds[i].lower_bounded() and x[i] < bounds[i].lower_limit():
          x[i] = bounds[i].lower_limit()
        if bounds[i].upper_bounded() and x[i] > bounds[i].upper_limit():
          x[i] = bounds[i].upper_limit()
    return x

  #debugging
  def study_parameters(self):
    # Vary refined parameters and print out function, gradient and curvature
    parameter_names = self.macrocycle_parameter_names()
    filename = "studyParams"
    x            = self.get_reparameterized_parameters()
    old_x        = x[:] #slice to force deep_copy
    g            = [0.] * self.nmp
    unrepar_x    = self.get_macrocycle_parameters()
    unrepar_oldx = [0.] * self.nmp
    unrepar_g    = [0.] * self.nmp
    bounds = self.reparameterized_bounds()
    ls     = self.reparameterized_large_shifts()
    bounds_present = False if (len(bounds) == 0) else True

    print("")
    print("Study behaviour of parameters near current values")

    fmin = self.target() # so called because it should be at the minimum

    # First check that function values from F, FG and FGH agree
    reparameterized_f_g = self.reparameterized_target_gradient()
    gradLogLike = reparameterized_f_g[0]
    g           = reparameterized_f_g[1]
    if abs(fmin-gradLogLike) >= 0.1:
      print("Likelihoods from target() and reparameterized_f_g() disagree")
      print("Likelihood from              target(): " + str(-fmin))
      print("Likelihood from reparameterized_f_g(): " + str(-gradLogLike))

    f_g_h = self.reparameterized_target_gradient_hessian()
    hessLogLike = f_g_h[0]
    is_diagonal = f_g_h[3]
    if abs(fmin-hessLogLike) >= 0.1:
      print("Likelihoods from target() and reparameterized_f_g_h() disagree")
      print("Likelihood from                target(): " + str(-fmin))
      print("Likelihood from reparameterized_f_g_h(): " + str(-hessLogLike))

    assert(abs(fmin-gradLogLike) < 0.1)
    assert(abs(fmin-hessLogLike) < 0.1)
    fmin = gradLogLike # Use function value from gradient for consistency below

    outstream = open("studyWhatAmI", "w")
    for i in range(self.nmp):
      outstream.write(" \"" + parameter_names[i] + " \"\n")
    outstream.close()

    LOGFILE_WIDTH = 90
    ndiv = 6
    dxgh = [0.] * self.nmp # Finite diff shift for grad & hess
    for i in range(self.nmp):
      if i+1 < 10:
        outstream = open(filename + ".0" + str(i+1), "w")
      else:
        outstream = open(filename + "."  + str(i+1), "w")
      print("\nRefined Parameter #:" + str(i+1) + " " + self.macrocycle_parameter_names()[i])
      print("Centered on " + self.to_sci_format(old_x[i]))
      xmin = old_x[i] - 2.*ls[i]
      dxgh[i] = ls[i]/50 # Small compared to exploration of function value
      if bounds_present and bounds[i].lower_bounded():
        xmin = max(xmin,bounds[i].lower_limit()+dxgh[i]) # Allow room for FD calcs
        print("Lower limit: " + self.to_sci_format(bounds[i].lower_limit()))
      else: # "Unbounded" parameters may have limits depending on correlations
        g_Fake = [0.] * self.nmp
        g_Fake[i] = ls[i]/10. # Put in range of plausible shift
        max_dist = self.maximum_distance(old_x,g_Fake)[0]
        fmax = sys.float_info.max
        ftol = sys.float_info.epsilon
        specialLimit = old_x[i]-max_dist*g_Fake[i] if (max_dist < fmax-ftol) else -fmax
        if specialLimit > xmin:
          xmin = specialLimit+dxgh[i]
          print("Special lower limit: "+ self.to_sci_format(specialLimit))

      xmax = old_x[i] + 2.*ls[i]
      if bounds_present and bounds[i].upper_bounded():
        xmax = min(xmax,bounds[i].upper_limit()-dxgh[i])
        print("Upper limit: "+ self.to_sci_format(bounds[i].upper_limit()))
      else:
        g_Fake = [0.] * self.nmp
        g_Fake[i] = -ls[i]/10.
        max_dist = self.maximum_distance(old_x,g_Fake)[0]
        fmax = sys.float_info.max
        ftol = sys.float_info.epsilon
        specialLimit = old_x[i]-max_dist*g_Fake[i] if (max_dist < fmax-ftol) else fmax
        if specialLimit < xmax:
          xmax = specialLimit-dxgh[i]
          print("Special upper limit: " + self.to_sci_format(specialLimit))
      print("Large shift: "+ self.to_sci_format(ls[i]))
      print("                                                                      FD curv    FD curv")
      print(" parameter   func-min    gradient   curvature  crv*lrg^2   FD grad   from func  from grad")
      dx = (xmax-xmin)/ndiv
      for j in range(ndiv+1):
        x[i] = xmin + j*dx
        self.set_reparameterized_parameters(x)
        f_g = self.reparameterized_target_gradient()
        gradLogLike = f_g[0]
        g           = f_g[1]
        df = gradLogLike - fmin
        f_g_h = self.reparameterized_target_gradient_hessian()
        h           = f_g_h[2].deep_copy()
        thisx = x[i] # Save before finite diffs
        x[i] = x[i] + dxgh[i]
        self.set_reparameterized_parameters(x)

        f_g = self.reparameterized_target_gradient()
        fplus = f_g[0]
        gplus = f_g[1]
        x[i] = x[i] - 2*dxgh[i] # Other side of original x
        self.set_reparameterized_parameters(x)

        f_g = self.reparameterized_target_gradient()
        fminus = f_g[0]
        gminus = f_g[1]

        fdgrad = (fplus - fminus)/(2.*dxgh[i])
        fdh_from_F = (fplus+fminus-2*gradLogLike)/(dxgh[i]**2)
        fdh_from_g = (gplus[i] - gminus[i])/(2.*dxgh[i])

        x[i] = thisx # Restore

        c_str = "%+.4e %+.3e %+.4e %+.4e %+.2e %+.3e %+.3e %+.3e" % \
                (x[i],df,g[i],h[i,i],h[i,i]*(ls[i]**2),fdgrad,fdh_from_F,fdh_from_g)
        c_str = c_str[0:LOGFILE_WIDTH]
        print(c_str)
        outstream.write(str(x[i]) + " " + str(df) + " " + str(g[i]) + " " + str(h[i,i]) + " " + str(fdgrad) + " " + str(fdh_from_g) + "\n")

      x[i] = old_x[i]
      outstream.close()
    # Finite difference tests for off-diagonal Hessian elements
    if is_diagonal == False:
      first = True
      self.set_reparameterized_parameters(x)
      f_g = self.reparameterized_target_gradient()
      f_g_h = self.reparameterized_target_gradient_hessian()
      h = f_g_h[2].deep_copy()
      for i in range(self.nmp-1):
        for j in range(i+1, self.nmp):
          if h[i,j] != 0.: # Hessian may be sparse
            if first == True:
              print("\n")
              print("Off-diagonal Hessian elements at current values")
              print("Parameter #s   Hessian     FD from grad")
              first = False
            x[j] = x[j] + dxgh[j]
            self.set_reparameterized_parameters(x)

            f_g = self.reparameterized_target_gradient()
            gplus = f_g[1]
            x[j] = x[j] - 2*dxgh[j] # Other side of original x
            self.set_reparameterized_parameters(x)

            f_g = self.reparameterized_target_gradient()
            gminus = f_g[1]
            fdhess = (gplus[i] - gminus[i])/(2.*dxgh[j])

            c_str = "%4i %4i    %+.4e   %+.4e" % \
                    (i+1,j+1,h[i,j],fdhess)
            c_str = c_str[0:LOGFILE_WIDTH]
            print(c_str)
            x[j] = old_x[j]

    print("\n")
    self.set_reparameterized_parameters(old_x)

  def finite_difference_gradient(self, frac_large):
    # Fraction of large shift to use for each type of function should be
    # investigated by numerical tests, varying by, say, factors of two.
    # Gradient is with respect to original parameters, so no reparameterisation is done
    # the scheme is g_i(x) ~= (f(x + sz_i) - f(x))/sz_i, unless we are at an upper bound,
    # in which case the scheme is g_i(x) ~= (f(x) - f(x - sz_i))/sz_i
    Gradient = [0.] * self.nmp

    bounds = self.bounds()
    bounds_present = False if (len(bounds) == 0) else True
    f = fplus = fminus = sz = 0.
    x = self.get_macrocycle_parameters()
    old_x = x[:] #slice to force deep_copy
    ls = self.macrocycle_large_shifts()
    self.set_macrocycle_parameters(x) # paranoia
    f = self.target()
    for i in range(self.nmp):
      sz = frac_large*ls[i]

      x[i] = old_x[i] + sz
      if bounds_present and bounds[i].upper_bounded() and x[i] > bounds.upper_limit():
        x[i] = old_x[i] - sz
        if bounds[i].lower_bounded() and x[i] < bounds[i].lower_limit():
          # if this part of the code is reached it means that x+sz is greater than the
          # upper bound, and x-sz is lower than the lower bound, so you should reduce
          # either the large shift for the relevant parameter or frac_large. (Or change
          # your bounds)
          assert(False)
        self.set_macrocycle_parameters(x)
        fminus = self.target()
        Gradient[i] = (f - fminus)/sz
      else:
        self.set_macrocycle_parameters(x)
        fplus = self.target()
        Gradient[i] = (fplus - f)/sz

      x[i] = old_x[i]
    self.set_macrocycle_parameters(old_x)
    return (f, Gradient)

  def finite_difference_hessian_by_gradient(self, frac_large, do_repar): pass

  def finite_difference_curvatures_by_gradient(self, frac_large, do_repar): pass

  def finite_difference_hessian_by_function(self, frac_large, do_repar): pass

  def finite_difference_curvature_by_function(self, frac_large, do_repar): pass

  def to_sci_format(self, f):
    """
    Format a float as a string in scientific notation.
    Done this way as a separate function to match the c++ code layout.
    """
    return str(f)

  def maximum_distance(self, x, p):
    """
    Return:
     The multiple of the search vector p that x can be shifted by before
     hitting the furthest bound. Used to define the limit for line search.
     The vector of allowed shifts for each parameter, which can be used to
     decide which parameters can't be moved along the search direction.
    """

    max_dist, ZERO = 0., 0.
    dist = [-1.] * self.nmp # multiple of p allowed before hitting box bound, or . -1. used to indicate can't move, i.e. at a bound
    bounds = self.reparameterized_bounds()
    bounded = [False] * self.nmp # given to maximum_distance_special

    if len(bounds) != 0:
      for i in range(self.nmp):
        if (p[i] < ZERO and bounds[i].lower_bounded()) or (p[i] > ZERO and bounds[i].upper_bounded()):
          bounded[i] = True
          limit = bounds[i].lower_limit() if (p[i] < ZERO) else bounds[i].upper_limit()
          dist[i] = (limit-x[i])/p[i] # if we are in bounds will positive, on bounds, is 0.
        max_dist = max(max_dist,dist[i])

    if len(bounds) == 1 and bounded[0]: #there can't be any correlations
      return (max_dist, dist)
    else:
      # The rationale behind maximum_distance_special() is to provide a way to impose non-linear
      # equality constraints.
      # Correlated, i.e. equality constrained parameters are checked here, updating
      # bounded and dist if relevant.
      # If correlated parameters are also limited by overall upper and lower bounds,
      # those can be set as well for the normal limit checks above.
      max_dist = max(max_dist,self.maximum_distance_special(x,p,bounded,dist,max_dist))

    # If any parameters are unbounded, search can carry on essentially indefinitely
    for i in range(self.nmp):
      if bounded[i] == False:
        max_dist = 1000000.
        break
    return (max_dist, dist)

  def adjust_hessian(self, h):
    """
    Return an "adjusted" if necessary, and a flag to say whether an
    adjustment was necessary

    Eliminate non-positive curvatures (diagonal elements of the hessian).
    Return true if any curvatures were non-positive.

    Determine mean factor by which curvatures differ from 1/largeShift^2
    from parameters with positive curvatures.  This gives the average
    curvature that would be obtained if all the parameters were scaled
    so that their largeShift values were one.
    Expected curvatures (consistent with relative largeShift values) are
    computed from this and used to replace non-positive curvatures.
    Corresponding off-diagonal terms are set to zero.
    RJR Note 31/5/06: tried setting a minimum fraction of the expected
    curvature, but this degraded convergence and depended too much on
    accurate estimation of largeShift.
    DHS Note 30/4/19: if this looks hacky, that's because it is. This
    looks to be a fix implemented when only diagonal hessians were
    considered as it is still possible to have negative eigen values
    after this step for non-diagonal hessians.
    To see this consider the case [1 3] which has eigen values (4,-2)
                                  [3 1]
    """

    meanfac, ZERO = 0., 0.
    queue_new_hessian = False
    npos = 0  #number of positive diagonal hessian entries
    ls = self.reparameterized_large_shifts()
    for i in range(self.nmp):
      if h[i,i] > ZERO:
        npos += 1
        meanfac += h[i,i] * ls[i]**2
    if npos != 0:
      meanfac /= npos
      for i in range(self.nmp):
        if h[i,i] <= ZERO:
          h[i,i] = meanfac/(ls[i]**2)
          for j in range(i+1, self.nmp):
            h[i,j] = h[j,i] = ZERO
    else: # Fall back on diagonal matrix based on ls shift values
      for i in range(self.nmp):
        h[i,i] = 100./(ls[i] ** 2)
        for j in range(self.nmp):
          if i != j: h[i,j] = h[j,i] = ZERO
    if npos != self.nmp: queue_new_hessian = True
    return (h, queue_new_hessian)
