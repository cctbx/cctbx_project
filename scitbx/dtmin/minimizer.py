from __future__ import print_function
from __future__ import division
from scitbx.array_family import flex
from scitbx.dtmin.realsymmetricpseudoinverse import RealSymmetricPseudoInverse
from math import sqrt

SILENCE = 0 #
LOGFILE = 1 #
VERBOSE = 2 #
FURTHER = 3 #
TESTING = 4 #

DEF_TAB = 3

class Minimizer(object):
  def __init__(self, output_level):
    self.output_level = output_level
    self.func_count_ = 0 #for counting the number of function evaluations in the line search
    self.parameter_names_ = None

  def run(self, refine_object, protocol, ncyc, minimizer_type, study_params=False):
    """Run the minimization for refine_object leaving it the minimized state
    refine_object:  object instantiated from a class that inherits from RefineBase
                    and that implements, at a minimum, the functions in Compulsory.py
    protocol:       list of list of strings. e.g. [["first"],["second"],["all"]]
    ncyc:           maximum number of microcycles to perform each macrocycle before
                    stopping, 50 is fairly typical
    minimizer_type: string setting the method of finding the direction to perform
                    the line-search, must be either "bfgs" or "newton"
    study_params:   boolean to trigger the running of the study_parameters
                    routine at the end of each macrocycle, to allow debugging of
                    gradient and curvature calculations by comparison with values
                    calculated by finite differences
    """

    refine_object.output_level = self.output_level

    self.check_input(protocol, ncyc, minimizer_type)

    refine_object.initial_statistics()

    for macro_cycle in range(len(protocol)):
      self.log_underline(TESTING,"MACROCYCLE #" + str(macro_cycle+1) + " OF " + str(len(protocol)))

      refine_object.set_macrocycle_protocol(protocol[macro_cycle])
      self.log_protocol(LOGFILE, protocol[macro_cycle])
      if refine_object.nmp == 0:
        continue # skip this macrocycle

      refine_object.setup()

      self.parameter_names_ = refine_object.macrocycle_parameter_names()
      self.log_parameters(VERBOSE, self.parameter_names_)

      refine_object.reject_outliers()

      self.log_ellipsis_start(LOGFILE,"Performing Optimisation")
      # not implementing create_report() code in the python version

      calc_new_hessian = True                                       # for bfgs
      queue_new_hessian = False # to be set during bfgs call        # for bfgs
      ncyc_since_hessian_calc = 0                                   # for bfgs
      previous_x = flex.double(refine_object.nmp)                          # for bfgs
      previous_g = flex.double(refine_object.nmp)                          # for bfgs
      previous_h_inverse_approximation = flex.double(refine_object.nmp**2) # for bfgs
      previous_h_inverse_approximation.resize(flex.grid(refine_object.nmp,refine_object.nmp))
      termination_reason = ""

      f, initial_f, previous_f, required_decrease = 0., 0., 0., 0.
      EPS = 1.E-10
      FTOL = 1.E-6
      TWO = 2.
      too_small_shift, too_small_decrease, zero_gradient = False, False, False

      f = initial_f = previous_f = refine_object.target()

      self.log_blank(VERBOSE)
      self.log_underline(VERBOSE,"Stats before starting refinement")
      refine_object.current_statistics()
      self.log_tab(1,VERBOSE,"Optimisation statistics, macrocycle #" + str(macro_cycle+1))
      self.log_tab_printf(1,VERBOSE,"Cycle     %-15s %-18s %-18s\n",("end-this-cycle","change-from-start","change-from-last"))
      self.log_tab_printf(1,VERBOSE,"Cycle#%-2d %15.3f %18.3s %18.3s\n",(0,initial_f,"",""))

      ########################### MINIMIZATION CYCLES ######################################
      cycle = 0
      while True:
        self.log_underline(TESTING,"MACROCYCLE #" + str(macro_cycle+1) + ": CYCLE #" + str(cycle+1))

        assert(self.check_in_bounds(refine_object))

        #create_report() code goes here

        previous_f = f
        too_small_shift = too_small_decrease = zero_gradient = queue_new_hessian = False # reset
        required_decrease = FTOL*(abs(f)+abs(previous_f)+EPS)/TWO # to not converge this cycle

        if minimizer_type.lower() == "bfgs":
          (f,ncyc_since_hessian_calc, queue_new_hessian, too_small_shift, zero_gradient, previous_x, previous_g, previous_h_inverse_approximation) = \
            self.bfgs(refine_object,
                    calc_new_hessian,                 # tells bfgs to calculate a new hessian
                    ncyc_since_hessian_calc,          # reset or incremented by bfgs
                    required_decrease,                # for the line search
                    previous_x,                       # for bfgs to update h_inverse_approximation
                    previous_g,                       # for bfgs to update h_inverse_approximation
                    previous_h_inverse_approximation) # for bfgs to update h_inverse_approximation

          self.log_hessian(TESTING,previous_h_inverse_approximation, self.parameter_names_)
        elif minimizer_type.lower() == "newton":
          (f, too_small_shift) = \
            self.newton(refine_object,
                        required_decrease)    # for the line search

        refine_object.current_statistics()
        self.log_tab_printf(1,VERBOSE,"Cycle#%-2d %15.3f %18.3f %18.3f\n",(cycle+1,f,(f-initial_f),(f-previous_f)))

        too_small_decrease = (previous_f - f < required_decrease)

        (termination_reason, calc_new_hessian) = self.check_termination_criteria(
                            too_small_shift,
                            too_small_decrease,
                            zero_gradient,
                            queue_new_hessian,
                            ncyc_since_hessian_calc,
                            cycle,
                            ncyc)

        #create_report()
        #create_report()
        #create_report()
        #create_report()
        #create_report()
        #create_report()
        #create_report()
        #create_report()

        #string is empty for continuation
        if len(termination_reason) != 0:
          self.log_tab(2,FURTHER, "termination reason: " + str(termination_reason))
          break

        cycle += 1
      ################################## END OF MINIMIZATION CYCLES ####################

      self.log_ellipsis_end(LOGFILE)
      if cycle == ncyc-1:
        self.log_tab(1,LOGFILE,"--- Iteration limit reached at cycle " + str(ncyc) + " ---")
      else:
        self.log_tab(1,LOGFILE,"--- Convergence before iteration limit (" + str(ncyc) +") at cycle " + str(cycle+1) + " ---")
      self.log_tab(1,LOGFILE,"Start-this-macrocycle End-this-macrocycle Change-this-macrocycle")
      self.log_tab_printf(1,LOGFILE,"%13.3f %21.3f %21.3f\n",(initial_f,f,(f-initial_f)))
      self.log_blank(LOGFILE)
      refine_object.current_statistics()

      # study behaviour of function, gradient and curvature as parameters varied around minimum
      # using mathematica file studyParams_short.nb
      if study_params:
        refine_object.study_parameters()

      refine_object.cleanup() #clean up ready for next macrocycle
    #end of macrocycle loop

    refine_object.finalize()
    refine_object.final_statistics()

  def check_input(self, protocol, ncyc, minimizer_type):
    """Check that the protocol, ncyc and MINIMSER arguments are valid"""

    if len(protocol) == 0:
      raise RuntimeError("No protocol for refinement")

    some_refinement = False
    for macro_cycle in range(len(protocol)):
      if len(protocol[macro_cycle]) != 0:
        some_refinement = True
        break

    if some_refinement == False:
      self.log_tab(1,LOGFILE,"No refinement of parameters")
      self.log_blank(LOGFILE)
      return

    allowed_minimizers = set(["bfgs", "newton"])
    if minimizer_type.lower() not in allowed_minimizers:
      raise RuntimeError("Invalid minimizer_type value: " + minimizer_type)

    if ncyc < 1:
      raise RuntimeError("NYC must be 1 or larger")

  def bfgs(self, refine_object, calc_new_hessian, ncyc_since_hessian_calc, required_decrease, previous_x, previous_g, previous_h_inverse_approximation):
    """Move the refinable parameters by one bfgs step"""

    queue_new_hessian = False
    too_small_shift = False
    zero_gradient = False
    x = refine_object.get_reparameterized_parameters()
    g = flex.double(refine_object.nmp)
    h_inverse_approximation = flex.double(refine_object.nmp**2)
    h_inverse_approximation.resize(flex.grid(refine_object.nmp,refine_object.nmp))
    min_dx_over_esd = 0.1
    f = 0.

    self.log_tab(1, TESTING, "== BFGS ==")

    if calc_new_hessian:
      self.log_tab(1,TESTING,"== Newton Step ==")

      calc_new_hessian = False # reset
      ncyc_since_hessian_calc = 0

      self.log_tab(1,TESTING,"BFGS: Gradient and New Hessian")
      self.log_ellipsis_start(TESTING,"BFGS: Calculating Gradient and New Hessian")

      f_g_h = refine_object.reparameterized_adjusted_target_gradient_hessian()
      starting_f               = f_g_h[0]
      g                        = f_g_h[1]
      h                        = f_g_h[2].deep_copy() # does this need to be a deep copy ???
      is_hessian_diagonal      = f_g_h[3]
      queue_new_hessian = f_g_h[0]

      self.log_ellipsis_end(TESTING)
      self.log_blank(TESTING)
      self.log_vector(TESTING, "BFGS::Newton: Gradient", g, self.parameter_names_)
      self.log_blank(TESTING)

      self.log_hessian(TESTING, h, self.parameter_names_)
      self.log_blank(TESTING)

      if flex.double(g).all_eq(0):  #convert g to a flex.double to get the all_eq method
        queue_new_hessian = False
        zero_gradient = True
        return (starting_f, ncyc_since_hessian_calc, queue_new_hessian, too_small_shift, zero_gradient, x, g, previous_h_inverse_approximation)

      #Set up to repeat perturbed Newton step if insufficient progress
      hii_rms = self.hessian_diag_rms(h)
      filtered_last = 0
      f = starting_f
      previous_x = x
      ntry = 0
      while ntry < 3 and (starting_f-f < required_decrease) and (filtered_last < refine_object.nmp):
        if is_hessian_diagonal:
          for i in range(refine_object.nmp):
            h_inverse_approximation[i,i] = 1/(h[i,i]+ntry*hii_rms)
        else:
          min_to_filter = 0
          if ntry>0:
            min_to_filter = filtered_last + max(1,refine_object.nmp/10)
            if min_to_filter >= refine_object.nmp:
              min_to_filter = min(filtered_last+1, refine_object.nmp-1)
          filtered = 0
          (h_inverse_approximation, filtered) = RealSymmetricPseudoInverse(h, filtered, min_to_filter)

          if filtered:
            self.log_tab_printf(1,TESTING,"Filtered %i small eigenvectors\n",filtered)
          self.log_tab(1,TESTING,"Pseudoinverse Hessian")
          self.log_hessian(TESTING,h_inverse_approximation,self.parameter_names_)

          filtered_last = filtered

        g = flex.double(g)
        p = - h_inverse_approximation.matrix_multiply(g)
        self.log_vector(TESTING,"BFGS::Newton: Search Vector, p",p,self.parameter_names_)
        self.log_blank(TESTING)

        # Figure out how far linesearch has to go to shift at least one parameter
        # by (minimum shift/esd) * esd, using inverse Hessian to estimate covariance matrix
        required_shift = self.required_shift(min_dx_over_esd, p, h_inverse_approximation)

        (f, too_small_shift) = self.line_search(refine_object,x,f,g,p,1.,required_decrease,required_shift,too_small_shift)
        ntry += 1
      previous_h_inverse_approximation = h_inverse_approximation
    else:
      ncyc_since_hessian_calc += 1

      self.log_tab(1,TESTING,"BFGS: Gradient")
      self.log_ellipsis_start(TESTING,"BFGS: Calculating Gradient")

      (f, g) = refine_object.reparameterized_target_gradient()

      if flex.double(g).all_eq(0):  #convert g to a flex.double to get the all_eq method
        calc_new_hessian      = False
        queue_new_hessian = False
        return (f, ncyc_since_hessian_calc, queue_new_hessian, too_small_shift, zero_gradient, x, g, previous_h_inverse_approximation)

      dx = flex.double(x) - flex.double(previous_x)
      dg = flex.double(g) - flex.double(previous_g)
      Hdg = previous_h_inverse_approximation.matrix_multiply(dg)
      u = refine_object.nmp * [0.]
      dx_dot_dg  = dx.dot(dg)
      dg_dot_Hdg = dg.dot(Hdg)

      if dx_dot_dg <= 0:
        self.log_tab(1,TESTING,"BFGS: dx_dot_dg <= 0: restarting with new Hessian")
        queue_new_hessian = True
        return (f, ncyc_since_hessian_calc, queue_new_hessian, too_small_shift, zero_gradient, x, g, previous_h_inverse_approximation)
      assert(dg_dot_Hdg)
      # BFGS update of the approximation to the inverse hessian
      # see Numerical Recipes in Fortran, Second Ed. p. 420
      # h_inverse_approximation contains the approximation to the inverse hessian
      for i in range(refine_object.nmp):
        u[i] = dx[i]/dx_dot_dg - Hdg[i]/dg_dot_Hdg
      for i in range(refine_object.nmp):
        for j in range(refine_object.nmp):
          h_inverse_approximation[i,j] = previous_h_inverse_approximation[i,j] + dx[i]*dx[j]/dx_dot_dg - Hdg[i]*Hdg[j]/dg_dot_Hdg + dg_dot_Hdg*u[i]*u[j]

      self.log_blank(TESTING)
      self.log_hessian(TESTING,h_inverse_approximation, self.parameter_names_)
      self.log_blank(TESTING)

      g = flex.double(g)
      p = - h_inverse_approximation.matrix_multiply(g)

      self.log_vector(TESTING,"BFGS: Search Vector, p",p,self.parameter_names_)
      self.log_blank(TESTING)

      # Figure out how far linesearch has to go to shift at least one parameter
      # by minimum shift/esd, using inverse Hessian to estimate covariance matrix
      required_shift = self.required_shift(min_dx_over_esd, p, h_inverse_approximation)

      (f, too_small_shift) = self.line_search(refine_object,x,f,g,p,1,required_decrease,required_shift,too_small_shift)

      # Check for positive grad
      if p.dot(g) >= 0.:
        self.log_tab(1,TESTING, "BFGS: p.g >= 0: restart with new Hessian")
        queue_new_hessian = True

      previous_h_inverse_approximation = h_inverse_approximation
    return (f, ncyc_since_hessian_calc, queue_new_hessian, too_small_shift, zero_gradient, x, g, previous_h_inverse_approximation)

  def newton(self, refine_object, required_decrease):
    """Move the refinable parameters by one newton step"""

    too_small_shift = False
    f = 0.
    x = refine_object.get_reparameterized_parameters()
    min_dx_over_esd = 0.1

    # ESD is the estimated standard deviation. Why we are considering standard deviations at this
    # point in the code may not be immediately obvious.
    # The naming only really makes snse when we are minimising a log likelihood function
    # If the likelihood can locally be approxmated by a
    # multivariate normal distribution, then the log likelihood may be approximated by a multidimensional
    # quadratic with its matrix being the inverse of the covariance matrix of the distribution.

    self.log_tab(1,TESTING,"==NEWTON==")
    self.log_tab(1,TESTING,"Newton: Gradient and Hessian")
    self.log_ellipsis_start(TESTING,"Newton: Calculating the gradient and hessian")

    f_g_h = refine_object.reparameterized_adjusted_target_gradient_hessian()
    f                     = f_g_h[0]
    g                     = f_g_h[1]
    h                     = f_g_h[2]
    is_hessian_diagonal   = f_g_h[3]
    #queue_new_hessian = f_g_h[4] #unused

    h_inv = h # (shallow copy)

    if flex.double(g).all_eq(0): # convert to flex double to get all_eq method
      self.log_tab(1,TESTING,"Newton: Zero gradient: exiting")
      return (f, too_small_shift)

    self.log_ellipsis_end(TESTING)
    self.log_blank(TESTING)
    self.log_vector(TESTING,"Newton: Gradient",g,self.parameter_names_)
    self.log_blank(TESTING)
    self.log_tab(1,TESTING,"Newton: After adjusting negative diagonals, but before Inverse")
    self.log_hessian(TESTING,h,self.parameter_names_)

    # invert hessian
    if is_hessian_diagonal:
      for i in range(refine_object.nmp):
        h[i,i] = 1/h[i,i]
    else:
      filtered = 0
      (h_inv, filtered) = RealSymmetricPseudoInverse(h, filtered)
      if filtered:
        self.log_tab_printf(1,TESTING,"Newton: Filtered %i negative or small eigenvectors\n",filtered)
      self.log_tab(1,TESTING,"Newton: After Computing Pseudo-Inverse")
      self.log_hessian(TESTING,h_inv,self.parameter_names_)

    # search vector p found from the Newton equation: p = - h_inv * g
    g = flex.double(g)
    p = -h_inv.matrix_multiply(g)

    self.log_vector(TESTING,"Newton: Parameters",x,self.parameter_names_)
    self.log_vector(TESTING,"Newton: Search Vector",p,self.parameter_names_)
    self.log_blank(TESTING)

    # Figure out how far linesearch has to go to shift at least one parameter
    # by minimum shift/esd, using inverse Hessian to estimate covariance matrix
    required_shift = self.required_shift(min_dx_over_esd, p, h_inv)

    (f, too_small_shift) = self.line_search(refine_object,x,f,g,p,1.,required_decrease,required_shift,too_small_shift)
    self.log_tab_printf(1,TESTING,"Newton Function %10.6f\n", f)
    return (f, too_small_shift)

  def line_search(self, refine_object, oldx, f, g, p, starting_distance, required_decrease, required_shift, too_small_shift):
    """Step along search direction p, trying to find a good enough local minimum"""

    too_small_shift = False
    self.func_count_ = 0

    x = flex.double(refine_object.nmp)

    ZERO, HALF, ONE, TWO, FIVE = 0., .5, 1., 2., 5.
    GOLDEN = (FIVE**.5 + ONE) / 2
    GOLDFRAC = (GOLDEN-ONE)/GOLDEN
    WOLFEC1 = 1e-4
    DTOL = 1e-5

    # Check grad of function in direction of line search
    grad = p.dot(g)
    if grad >= ZERO:
      self.log_tab(1,TESTING,"Grad of target in search direction >= 0: grad = " +str(grad))
      too_small_shift = True
      #report
      return (f, too_small_shift)

    self.log_tab(1,TESTING,"Determining Stepsize")
    self.log_tab_printf(1,TESTING,"f(  Start  ) = %12.6f\n", f)
    self.log_tab_printf(1,TESTING,"|",())

    # Put distance each parameter can move before it hits a bound in dist array. If all parameters are
    # bounded in search direction, need to know maxDist we can move before hitting ultimate bound.
    (max_dist, dist) = refine_object.maximum_distance(tuple(oldx),tuple(p))
    if starting_distance > max_dist:
      if max_dist > ZERO:
        self.log_tab(1,TESTING,"To avoid exceeding bounds for all parameters, starting distance reduced from "
              + str(starting_distance) + " to " + str(max_dist))
        starting_distance = max_dist
      else:
        self.log_tab(1,TESTING,"All parameters have hit a bound -- end search")
        too_small_shift = True
        #report
        return (f, too_small_shift)

    # Get largeShifts for damping in shift_parameters
    macrocycle_large_shifts = refine_object.reparameterized_large_shifts()

    fk, flo, fhi, dk, dlo, dhi = 0., 0., 0., 0., 0., 0.

    #Sample first test point
    self.log_tab(1,TESTING,"Unit distance = " + str(starting_distance))

    dk = starting_distance
    fk = self.shift_score(refine_object, dk, x, oldx, p, dist, macrocycle_large_shifts)

    WOLFEFRAC = HALF # 0.5 slightly better in tests of 0,0.25,0.5,1
    if (dk >= WOLFEFRAC) and (fk < f+WOLFEC1*dk*grad): # Significant fraction of (quasi-)Newton shift # Wolfe condition for first step
      self.log_tab(1,TESTING,"Satisfied Wolfe condition for first step")
      self.log_tab_printf(1,TESTING,"Final Target Value: %12.6f distance: %12.6f \n",(fk, dk))
      self.log_tab(1,TESTING,"Bracketing (A) took " + str(self.func_count_) + " function evaluations")
      if dk < required_shift:
        too_small_shift = True
      #report
      return (fk, too_small_shift)

    dlo = ZERO
    flo = f

    if f <= fk: # First step is too big
      while f <= fk:
        # Shrink shift until we find a value lower than starting point
        # Fit quadratic to grad and f at 0 and dk, choose its minimum (but not too small)
        if dk < starting_distance/1.E10: # Give up if shift too small
          refine_object.set_reparameterized_parameters(oldx) # Reset x before returning
          too_small_shift = True
          self.log_tab(1,TESTING,"line_search: Shift too small")
          f = fk
          #report
          return (f, too_small_shift)

        fhi = fk
        dhi = dk
        denom = TWO * (grad*dk+(f-fk))
        if denom < 0:
          dk *= min(0.99, max(0.1,grad*dk/denom))
        else: # unlikely but possible for vanishingly small gradient
          dk *= 0.1
        fk = self.shift_score(refine_object, dk, x, oldx, p, dist, macrocycle_large_shifts)

      if fk < f+WOLFEC1*dk*grad:
        self.log_tab(1,TESTING,"Satisfied Wolfe condition after backtracking")
        self.log_tab_printf(1,TESTING,"Final Target Value: %12.6f distance: %12.6f \n", (fk, dk))
        self.log_tab(1,TESTING,"Bracketing (B) took " + str(self.func_count_) + " function evaluations")
        if dk < required_shift:
          too_small_shift = True
        #report
        return (fk, too_small_shift)

    else: # First step goes down.
      if dk >= max_dist: # First step already hit ultimate bound
        dhi = max_dist
        fhi = fk
      else:
        dhi = min(max_dist, (ONE+GOLDEN)*dk)
        fhi = self.shift_score(refine_object, dhi, x, oldx, p, dist, macrocycle_large_shifts)
        if (fhi <= fk and
            dhi >= WOLFEFRAC and       # Significant fraction of (quasi-)Newton shift
            fhi < f+WOLFEC1*dhi*grad): # Wolfe condition for first step
          self.log_tab(1,TESTING,"Satisfied Wolfe condition for extended step")
          self.log_tab_printf(1,TESTING,"Final Target Value: %12.6f distance: %12.6f \n", (fhi, dhi))
          self.log_tab(1,TESTING,"Bracketing (C) took " + str(self.func_count_) + " function evaluations")
          if dhi < required_shift:
            too_small_shift = True
          #report
          return (fhi, too_small_shift)

      if (flo > fk) and (fk >= fhi):  # Still not coming up at second step
        while (fk >= fhi) and (dhi < max_dist):
          flo = fk
          dlo = dk

          fk = fhi
          dk = dhi

          dhi = min(max_dist, dk + GOLDEN*(dk-dlo))
          fhi = self.shift_score(refine_object, dhi, x, oldx, p, dist, macrocycle_large_shifts)
          if (fhi <= fk and
              dhi >= WOLFEFRAC and # Significant fraction of (quasi-)Newton shift
              fhi < f+WOLFEC1*dhi*grad): # Wolfe condition for first step
            self.log_tab(1,TESTING,"Satisfied Wolfe condition for further extended step")
            self.log_tab_printf(1,TESTING,"Final Target Value: %12.6f distance: %12.6f \n", (fhi, dhi))
            self.log_tab(1,TESTING,"Bracketing (D) took " + str(self.func_count_) + " function evaluations")
            if dhi < required_shift:
              too_small_shift = True
            #report
            return (fhi, too_small_shift)

        if (dhi >= max_dist) and (fk >= fhi):
          # Reached boundary without defining bracket.
          # Use finite differences to test if still going down.
          # If so, stop this line search.  Otherwise, we've verified bracket.
          dk = (1.-DTOL)*dhi
          fk = self.shift_score(refine_object, dk, x, oldx, p, dist, macrocycle_large_shifts)
          if fk >= fhi:
            self.shift_parameters(refine_object,dhi,x,oldx,p,dist,macrocycle_large_shifts)
            self.log_tab(1,TESTING,"Reached boundary without defining bracket")
            if dhi < required_shift:
              too_small_shift = True
            #report
            return (fhi, too_small_shift)
      if dk >= WOLFEFRAC and fk < f+WOLFEC1*dk*grad:
        self.log_tab(1,TESTING,"Satisfied Wolfe condition for bigger than initial step")
        self.log_tab_printf(1,TESTING,"Final Target Value: %12.6f distance: %12.6f \n", (fk, dk))
        self.log_tab(1,TESTING,"Bracketing (E) took " + str(self.func_count_) + " function evaluations")

        self.shift_parameters(refine_object,dk,x,oldx,p,dist,macrocycle_large_shifts)  # Set to current best before returning
        if dk < required_shift:
          too_small_shift = True
        #report
        return (fk, too_small_shift)

    ################### End of Bracketing ###########################


    self.log_tab(1,TESTING,"Bracketing (target) took " + str(self.func_count_) + " function evaluations")
    self.log_tab(1,TESTING,"xlow= " + str(dlo) + ", xhigh= " + str(dhi))
    dmid = dk
    fmid = fk

    # BISECTING/INTERPOLATION STARTS HERE (we now know that the value is between dlo and dhi)

    # let tolerance be small but not too close to zero as to ensure our algorithm tries a new value
    # distinctly different from existing one
    XTOL = 0.01

    a, b, c = 0., 0., 0.,
    stepsize = min(dhi-dmid, dmid-dlo)
    lastlaststepsize = 0.
    laststepsize = stepsize*TWO # permit first step to actually happen
    d1, f1, d2, f2 = 0., 0., 0., 0.

    self.log_tab_printf(1,TESTING,"Line Search |",())
    for i in range(20): # fall-back upper limit on evaluations in interpolation stage
      (quadratic_coefficients, a, b, c) = self.quadratic_coefficients(flo, fmid, fhi, dlo, dmid, dhi, a, b, c)
      if quadratic_coefficients: # Make sure there will be quadratic with a minimum
        self.log_tab(1, TESTING, "Quadratic interpolation used" )
        dk = -b/(TWO*c)
        dk = max(dk,dlo)
        dk = min(dk,dhi)

      else:
        self.log_tab(1, TESTING, "No quadratic interpolation possible, taking interval midpoint..." )
        c = ZERO
        dk = dmid

      lastlaststepsize = laststepsize
      laststepsize = stepsize
      stepsize = abs(dk-dmid)

      self.log_tab_printf(0,TESTING," stepsize= %12.6f, laststepsize= %12.6f,  ",(stepsize, laststepsize))
      if (c > ZERO and                        # i.e. the interval has a local minimum, not a local maximum so proceed.
          dk-dlo > XTOL*(dhi-dlo) and         # New minimum is sufficiently far from previous points
          dhi-dk > XTOL*(dhi-dlo) and         # to make quadratic convergence probable.
          abs(dk-dmid) > XTOL*(dhi-dlo) and
          stepsize <= HALF*lastlaststepsize): # Steps are decreasing sufficiently in size
        self.log_tab(1,TESTING,"Use quadratic step")
        fk = self.shift_score(refine_object, dk, x, oldx, p, dist, macrocycle_large_shifts)

      else: # Golden section instead
        dk = dmid + GOLDFRAC*(dhi-dmid) if (dhi-dmid >= dmid-dlo) else dmid - GOLDFRAC*(dmid-dlo)
        self.log_tab(1,TESTING,"Use golden search step")
        fk = self.shift_score(refine_object, dk, x, oldx, p, dist, macrocycle_large_shifts)

      if (dk < dmid): # Get data sorted by size of shift
        d1 = dk
        f1 = fk
        d2 = dmid
        f2 = fmid

      else:
        d1 = dmid
        f1 = fmid
        d2 = dk
        f2 = fk

      if fk < fmid: # Keep current best at dmid
        dmid = dk
        fmid = fk

      self.log_tab_printf(0,TESTING,"=",())
      self.log_tab_printf(0,TESTING,"[%12.6f,%12.6f,%12.6f,%12.6f]\n",(dlo,d1,d2,dhi))
      self.log_tab_printf(0,TESTING,"[    ----    ,%12.6f,%12.6f,    ----    ]\n",(f1,f2))
      if f1 < f2:
        self.log_tab_printf(0,TESTING,"          <        ^^^^        >                     \n",())
      else:
        self.log_tab_printf(0,TESTING,"                       <       ^^^^       >          \n",())

      #convergence tests
      if fmid < f - max(required_decrease, - WOLFEC1*dmid*grad):
        self.log_tab(1,TESTING,"Stop linesearch: good enough for next cycle")
        break

      elif (d2-d1)/d2 < 0.01:
        self.log_tab(1,TESTING,"Stop linesearch: fractional change in stepsize < 0.01")
        break

      # Not converged, so prepare for next loop
      if f1 < f2:
        dhi = d2
        fhi = f2

      else:
        dlo = d1
        flo = f1

    # converged of limit in steps

    if f < fmid: # shouldn't happen, but catch possibility that function didn't improve
      refine_object.set_reparameterized_parameters(oldx)
    else:
      self.shift_parameters(refine_object,dmid,x,oldx,p,dist,macrocycle_large_shifts) # Make sure best shift has been applied.
      f = fmid

    self.log_blank(TESTING)
    self.log_tab(1,TESTING,"This line search took " + str(self.func_count_) + " function evaluations")
    self.log_blank(TESTING)
    self.log_tab_printf(1,TESTING,"Final Target Value: %12.6f distance: %12.6f \n", (f, dmid))
    self.log_tab(1,TESTING,"End distance = " +str(dmid))
    self.log_tab(1,TESTING,"Prediction ratio = " +str(dmid/starting_distance))
    self.log_blank(TESTING)
    if dmid < required_shift:
      too_small_shift = True
    #report
    return (f, too_small_shift)

  def shift_score(self, refine_object, a, x, oldx, p, dist, macrocycle_large_shifts, bcount = True):
    self.shift_parameters(refine_object,a,x,oldx,p,dist,macrocycle_large_shifts)
    f = refine_object.target()
    if bcount:
      self.func_count_ += 1
    self.log_tab_printf(1,TESTING,"f(%9.6f) = %10.6f\n", (a, f))
    return f

  def shift_parameters(self, refine_object, a, x, oldx, p, dist, macrocycle_large_shifts):
    x = refine_object.damped_shift(a,tuple(oldx),tuple(p),tuple(dist),tuple(macrocycle_large_shifts))
    refine_object.set_reparameterized_parameters(x)

  def check_termination_criteria(
    self,
    too_small_shift,
    too_small_decrease,
    zero_gradient,
    queue_new_hessian,
    ncyc_since_hessian_calc,
    cycle,
    ncyc):

    """
    Test whether we should terminate the current macrocycle

     Return a string containing the reason to terminate (empty if there is none)
     Also return a bool, calc_new_hessian if it is determined that
     the next microcycle in a bfgs minimization should start with a new hessian.
     (naturally this only makes sense if we are not terminating the macrocycle
     this time.)
    """

    calc_new_hessian = False

    if cycle == ncyc-1:
      return ("cycle limit", calc_new_hessian)

    if zero_gradient:
      self.log_tab(1,TESTING,"Zero gradient: exiting")
      return ("zero_gradient", calc_new_hessian)
    if too_small_shift or too_small_decrease:
      if (queue_new_hessian or not too_small_decrease) and (ncyc_since_hessian_calc != 0): # queue_new_hessian is only touched by bfgs
        self.log_tab(1,TESTING,"Recalculate Hessian and carry on")
        calc_new_hessian = True # tell bfgs to calculate a new hessian next time rather than doing the bfgs update to the inverse hessian approximation
        return ("", calc_new_hessian)
      else:
        self.log_blank(VERBOSE)
        self.log_tab(1,VERBOSE,"---CONVERGENCE OF MACROCYCLE---")
        self.log_blank(VERBOSE)
        return ("converged", calc_new_hessian)
    return ("",calc_new_hessian)

  def hessian_diag_rms(self, h):
    """Calculate root-mean-square of the hessian diagonal entries"""

    hii_rms = 0.
    N = int(round(sqrt(len(h))))
    assert(N)
    for i in range(N):
      hii_rms += h[i,i]**2
    hii_rms = sqrt(hii_rms/N)
    return hii_rms

  def required_shift(self, min_dx_over_esd, p, h_inv):
    """
    Figure out how far linesearch has to go to shift at least one parameter
    by minimum shift/esd, using inverse Hessian to estimate covariance matrix
    DHS: this was rather tricky to get my head around. The following is an attempt
    at a graphical explanation:

    Outer box is the esd estimated using sqrt(h(i,i))
    Inner box is min_dx_over_esd * esd  ( = min_dx_over_esd * sqrt(h(i,i)))

    We want to find the minimum multiple of p (the arrow/vector thing in the diagram)
    needed to just hit the inner box (in this case it will be less than 1 and hit the 2nd parameter's limit first)

    -------------------------------------
    |                    ^              |
    |                   /               |
    |      -----------------------      |
    |      |          /          |      |
    |      -----------------------      |
    |                                   |
    |                                   |
    -------------------------------------
    """
    required_shift = 0.
    for i in range(len(p)):
      if h_inv[i,i] > 0:
        required_shift = max(required_shift,abs(p[i])/sqrt(h_inv[i,i]))
    if required_shift > 0:
      required_shift = min_dx_over_esd/required_shift
    return required_shift

  def quadratic_coefficients(self, y1, y2, y3, x1, x2, x3, a, b, c):
    """Find the coefficients to a quadratic specified by three points on the quadratic

    Given three points in the x-y plane get the coefficients of the corresponding
    to the parabolic equation y(x) = a + b*x + c*x^2 that passes through them all
    """

    a, b, c = 0., 0., 0.

    det = x2*x3*x3 - x3*x2*x2 - x1*(x3*x3 - x2*x2) + x1*x1*(x3 - x2)
    if (abs(det) < 1e-7): # not sure how to replicate c++ version, just using a very small number for now
      return (False, a, b, c)

    invdet = 1./det

    #  a = (y1*(x2*x3*x3 - x2*x2*x3) + y2*(x1*x3*x3 - x3*x1*x1) + y3*(x1*x2*x2 - x2*x1*x1))*invdet // correct but unstable
    b = (y1*(x2*x2 - x3*x3) + y2*(x3*x3 - x1*x1) + y3*(x1*x1 - x2*x2))*invdet
    c = (y1*(x3 -x2) + y2*(x1 - x3) + y3*(x2 -x1))*invdet
    a = y1 - (b*x1 +c*x1*x1)

    if c == 0:  #this is a direct translation of the c++ code, containing a floating point equality operation
      return (False, a, b, c)
    else:
      return (True , a, b, c)

  def check_in_bounds(self, refine_object):
    bounds = refine_object.reparameterized_bounds()
    if len(bounds) == 0:
      return True
    else:
      x = refine_object.get_reparameterized_parameters()
      macrocycle_large_shifts = refine_object.reparameterized_large_shifts()
      TOL = 1.e-03 # Tolerance for deviation in bound as fraction of LargeShift
      assert(len(bounds) == len(macrocycle_large_shifts))
      for i in range(len(x)):
        thistol = macrocycle_large_shifts[i]*TOL
        if    (bounds[i].lower_bounded() and (x[i] < (bounds[i].lower_limit() - thistol) ) ) \
          or (bounds[i].upper_bounded() and (x[i] > (bounds[i].upper_limit() + thistol) ) ) :
          raise RuntimeError("Parameter #"  + str(i+1) + \
            " (Lower: " + str(bounds[i].lower_bounded()) + \
            " " + str(bounds[i].lower_limit()) + \
            ") (Value: " + str(x[i]) + \
            ") (Upper: " + str(bounds[i].upper_bounded()) + \
            " " + str(bounds[i].upper_limit()) + ")")
          return False
      return True

  ### printing methods

  def log_output(self, depth, output_level, string, add_return):
    if output_level <= self.output_level:
      print(depth*DEF_TAB*' ' + string, end='\n' if add_return else '')

  def log_tab_printf(self, t, output_level, text, fformat): #format is already taken
    if output_level <= self.output_level:
      print((t*DEF_TAB*' ' + text) % fformat, end = "")

  def log_blank(self, output_level):
    self.log_output(0,output_level,"",True)

  def log_underline(self, output_level, string):
    self.log_blank(output_level)
    self.log_tab(1,output_level, string)
    self.log_line(1,output_level, len(string),'-')

  def log_line(self, t, output_level, len, c):
    self.log_tab(t, output_level,len*c)

  def log_protocol(self, output_level, macrocycle_protocol):
    self.log_underline(output_level, "Protocol:")
    if (macrocycle_protocol[0] == "off"):
      self.log_tab(1,output_level,"No parameters to refine this macrocycle")
      self.log_blank(output_level)
    else:
      for i in range(len(macrocycle_protocol)):
        self.log_tab(2,output_level,macrocycle_protocol[i].upper() + " ON")
      self.log_blank(output_level)

  def log_tab(self, t, output_level, text, add_return=True):
    self.log_output(t, output_level, text, add_return)

  def log_ellipsis_start(self, output_level, text):
    self.log_blank(output_level)
    self.log_output(1,output_level,text+"...",True)

  def log_ellipsis_end(self, output_level):
    self.log_output(2,output_level,"...Done",True)
    self.log_blank(output_level)

  def log_parameters(self, output_level, macrocycle_parameter_names):
    self.log_tab(1,output_level, "Parameters:")
    for i in range(len(macrocycle_parameter_names)):
      self.log_tab(1,output_level, "Refined Parameter #:" + str(i+1) + " " + macrocycle_parameter_names[i])

  def log_hessian(self, output_level, hessian, macrocycle_parameter_names):
    max_dim = 100 #very generous but not infinite
    (n_rows, n_cols) = hessian.all()
    assert(n_rows == n_cols)
    self.log_tab(1,output_level,"Matrix")
    if n_rows > max_dim:
      self.log_tab(1,output_level,"Matrix is too large to write to log: size = "+str(n_rows))
    for i in range(min(n_rows, max_dim-1)):
      hess = ""
      line = ""
      for j in range(min(n_cols, max_dim-1)):
        hess = ("%+3.1e " % hessian[i,j])
        if hessian[i,j] == 0:
          line += " ---0--- "
        else:
          line += hess
      line = line[:-1] # cut the last space off for neatness
      if n_rows >= max_dim:
        self.log_tab_printf(1,output_level,"%3d [%s etc...]\n",(i+1,line))
      else:
        self.log_tab_printf(1,output_level,"%3d [%s]\n",(i+1,line))
    if n_rows >= max_dim:
      self.log_tab_printf(0,output_level," etc...\n",())
    self.log_blank(output_level)
    self.log_tab(1,output_level,"Matrix Diagonals")
    for i in range(n_rows):
      self.log_tab_printf(1,output_level,"%3d [% 9.4g] %s\n",
          (i+1,hessian[i,i],
          macrocycle_parameter_names[i] if (i < len(macrocycle_parameter_names)) else ""))
    self.log_blank(output_level)

  def log_vector(self, output_level, what, vec, macrocycle_parameter_names):
    self.log_tab(1,output_level,what)
    for i in range(len(vec)):
      self.log_tab_printf(1,output_level,"%3d [% 9.4g] %s\n",
          (i+1,vec[i],
          macrocycle_parameter_names[i] if (i < len(macrocycle_parameter_names)) else ""))
    self.log_blank(output_level)
