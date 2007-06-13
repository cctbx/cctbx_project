import cctbx.xray.structure_factors
import cctbx.xray.ext
import smtbx.refinement.ext
from cctbx import crystal
from cctbx import xray
from cctbx.xray.minimization import add_gradients as cctbx_add_gradients
from cctbx.array_family import flex
import scitbx.lbfgs
import scitbx.math
from libtbx import adopt_init_args
from stdlib import math
from cctbx import adptbx


class lbfgs(object):
  """
  A minimisation of a function of the free parameters of a X-ray structure,
  relying on the LBFGS algorithm.

  That function is represented by a functor object f which should abide to the
  following interface
    - B{f.f_obs()}
    - B{f(f_calc, compute_derivatives)}
  in the manner of L{cctbx.xray.target_functors.least_squares_residual}
  or L{cctbx.xray.target_functors.intensity_correlation}
  """

  def __init__(self, target_functor, xray_structure,
                     occupancy_penalty=None,
                     lbfgs_termination_params=None,
                     lbfgs_core_params=None,
                     correct_special_position_tolerance=1.e-2,
                     use_special_position_constraints=False,
                     cos_sin_table=True,
                     structure_factor_algorithm=None,
                     verbose=0,
                     log=None):
    """
    @type target_functor: any object like f
    @param target_functor: the numerical function of the structure to minimise
    @type xray_structure: cctbx.xray.structure
    @param xray_structure: the structure the target_functor is calculated for
    @type occupancy_penalty: ??
    @param occupancy_penalty: ??
    @type lbfgs_termination_params: like
    L{scitbx.lbfgs.termination_parameters}
    @param lbfgs_termination_params: the bunch of parameters used to decide
    where to stop the minimisation
    @type lbfgs_core_params: like L{scitbx.lbfgs.core_parameters}
    @param lbfgs_core_params: the bunch of parameters for the numerical core
    of the lbfgs algorithm
    @type correct_special_position_tolerance: number
    @param correct_special_position_tolerance: ??
    @type cos_sin_table: bool
    @param cos_sin_table: whether to use tabulated cosines and sines in
    structure factor computations
    @type structure_factor_algorithm: string
    @param structure_factor_algorithm: the name of the method to be used to
    compute structure factors; it must be one of those provided by the factory
    L{cctbx.xray.structure_factors.from_scatterers}
    @param verbose: a flag specifying the verbosity of the log printed out
    during the minimisation process; 0 means silent whereas positive numbers
    print an increasing amount of information
    """
    adopt_init_args(self, locals())
    self.scatterer_grad_flags_counts = \
        cctbx.xray.ext.scatterer_grad_flags_counts(
          self.xray_structure.scatterers())
    self.grad_flags_counts = \
        cctbx.xray.ext.scatterer_grad_flags_counts(
          self.xray_structure.scatterers())
    self.structure_factors_from_scatterers = \
      cctbx.xray.structure_factors.from_scatterers(
        miller_set=self.target_functor.f_obs(),
        cos_sin_table=cos_sin_table)
    self.structure_factor_gradients = \
      cctbx.xray.structure_factors.gradients(
        miller_set=self.target_functor.f_obs(),
        cos_sin_table=cos_sin_table)
    self.unconstrained_parameters = xray_structure.n_parameters_XXX()
    self.x = flex.double(self.n_parameters(), 0)
    xray_structure.tidy_us(u_min=1.e-6)
    self.snap_to_special_position_sites()
    self._scatterers_start = xray_structure.scatterers()
    self._d_min = self.target_functor.f_obs().d_min()
    self.first_target_value = None
    if self.verbose:
      print "begin{minimisation}"
    if log:
      try:
        print >> log
      except:
        import sys
        self.log = sys.stdout
      else:
        self.log = log
    else:
      self.log = None
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=lbfgs_termination_params,
      core_params=lbfgs_core_params,
      log=self.log
    )
    if self.verbose: print "end{minimisation}"
    self.apply_shifts()
    del self._scatterers_start
    del self._d_min
    self.compute_target(compute_gradients=False)
    self.final_target_value = self.target_result.target()

  def n_parameters(self):
    n = self.xray_structure.n_parameters_XXX()
    site_symmetry_table = self.xray_structure.site_symmetry_table()
    scatterers = self.xray_structure.scatterers()
    for i_seq in site_symmetry_table.special_position_indices():
      op = site_symmetry_table.get(i_seq)
      sc = scatterers[i_seq]
      if(sc.flags.grad_site()):
        n -= op.site_constraints().n_dependent_params()
      if(sc.flags.grad_u_aniso() and sc.flags.use_u_aniso()):
        n -= op.adp_constraints().n_dependent_params()
    return n

  def snap_to_special_position_sites(self):
    site_symmetry_table = self.xray_structure.site_symmetry_table()
    scatterers = self.xray_structure.scatterers()
    for i_seq in site_symmetry_table.special_position_indices():
      scatterers[i_seq].site = crystal.correct_special_position(
        crystal_symmetry=self.xray_structure,
        special_op=site_symmetry_table.get(i_seq).special_op(),
        site_frac=scatterers[i_seq].site,
        tolerance=self.correct_special_position_tolerance)

  def apply_shifts(self):
    apply_shifts_result = smtbx.refinement.ext.\
                        apply_special_position_constrained_shifts(
      unit_cell=self.xray_structure.unit_cell(),
      site_symmetry_table=self.xray_structure.site_symmetry_table(),
      scatterers=self._scatterers_start,
      shifts=self.x)
    shifted_scatterers = apply_shifts_result.shifted_scatterers
    if self.log is not None:
      print >> self.log, "\n******* Before applying shifts *******"
      self.xray_structure.show_scatterers(self.log)
    self.xray_structure.replace_scatterers(scatterers=shifted_scatterers)
    if self.log is not None:
      print >> self.log, "\n ********** Shifts *******************"
      self.show_shifts(self.log)
      print >> self.log, "\n******* After applying shifts ********"
      self.xray_structure.show_scatterers(self.log)
      print >> self.log, "\n**************************************\n"
    return apply_shifts_result.u_iso_refinable_params

  def show_log(self, f=None):
    import sys
    if self.log is sys.stdout: return
    if f is None: f = sys.stdout
    print >> f, self.log.getvalue()

  def show_shifts(self, log=None):
    import sys
    if log is None: log = sys.stdout
    site_symmetry_table = self.xray_structure.site_symmetry_table()
    i=0
    i_sc = 0
    for sc in self.xray_structure.scatterers():
      op = site_symmetry_table.get(i_sc)
      i_sc += 1
      print >> log, "%-4s" % sc.label
      if sc.flags.grad_site():
        n = op.site_constraints().n_independent_params()
        if n != 0:
          print >> log, ("site:" + "%7.4f, "*(n-1) + "%7.4f")\
                % tuple(self.x[i:i+n])
        i += n
      if sc.flags.grad_u_iso() and sc.flags.use_u_iso():
        if not(sc.flags.tan_u_iso() and sc.flags.param > 0):
          print >> log, "u_iso: %6.4f" % self.x[i]
          i += 1
      if sc.flags.grad_u_aniso() and sc.flags.use_u_aniso():
        n = op.adp_constraints().n_independent_params()
        print >> log, ("u_aniso:" + "%6.3f, "*(n-1) + "%6.3f") % tuple(self.x[i:i+n])
        i += n
      if sc.flags.grad_occupancy():
        print >> log, "occ: %4.2f" % self.x[i]
        i += 1
      if sc.flags.grad_fp():
        print >> log, "f': %6.4f" % self.x[i]
        i += 1
      if sc.flags.grad_fdp():
        print >> log, "f'': %6.4f" % self.x[i]
        i += 1
      print >> log

  def compute_target(self, compute_gradients):
    self.f_calc = self.structure_factors_from_scatterers(
      xray_structure=self.xray_structure,
      miller_set=self.target_functor.f_obs(),
      algorithm=self.structure_factor_algorithm).f_calc()
    self.target_result = self.target_functor(
      self.f_calc,
      compute_gradients)

  def compute_functional_and_gradients(self):
    u_iso_refinable_params = self.apply_shifts()
    self.compute_target(compute_gradients=True)
    self.f = self.target_result.target()
    if self.log is not None:
      print >> self.log, "Scale factor: %10.4f"\
            % self.target_result.scale_factor()
    if (self.first_target_value is None):
      self.first_target_value = self.f
      if self.verbose: print "initial f: %f" % self.f
    if (self.occupancy_penalty is not None
        and self.grad_flags_counts != 0):
      occupancies = self.xray_structure.scatterers().extract_occupancies()
      for occupancy in occupancies:
        self.f += self.occupancy_penalty.functional(occupancy)
    self.g = self.structure_factor_gradients(
      xray_structure=self.xray_structure,
      u_iso_refinable_params=u_iso_refinable_params,
      miller_set=self.target_functor.f_obs(),
      d_target_d_f_calc=self.target_result.derivatives(),
      n_parameters=self.unconstrained_parameters,
      algorithm=self.structure_factor_algorithm).packed()
    if (self.occupancy_penalty is not None
        and self.grad_flags_counts != 0):
      g = flex.double()
      for occupancy in occupancies:
        g.append(self.occupancy_penalty.gradient(occupancy))
      del occupancies
      cctbx_add_gradients(
        scatterers=self.xray_structure.scatterers(),
        xray_gradients=self.g,
        occupancy_gradients=g)
      del g
    reduction = smtbx.refinement.ext.special_position_constrained_gradients(
                unit_cell=self.xray_structure.unit_cell(),
                site_symmetry_table=self.xray_structure.site_symmetry_table(),
                scatterers=self.xray_structure.scatterers(),
                xray_gradients=self.g
              )
    self.g = reduction.reduced_gradients
    if (self.verbose > 1):
      print "xray.minimization line search: f,rms(g):",
      print self.f, math.sqrt(flex.mean_sq(self.g))
    return self.f, self.g

  def callback_after_step(self, minimizer):
    if (self.verbose > 0):
      print "xray.minimization step: f,iter,nfun:",
      print self.f,minimizer.iter(),minimizer.nfun()
