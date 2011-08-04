from mmtbx.scaling import absolute_scaling, relative_scaling
import iotbx.data_plots
from cctbx import adptbx
from cctbx import maptbx, miller
from cctbx.array_family import flex
from scitbx.math import phase_error
import libtbx.callbacks # import dependency
from libtbx import adopt_init_args, Auto, group_args
import libtbx
import math, sys
from cStringIO import StringIO

pi_180 = math.atan(1)/ 45

master_params_str = """
  solvent_fraction = None
    .type = float
    .style = bold
  initial_steps = 10
    .type = int
  shrink_steps = 20
    .type = int
  final_steps = 10
    .type = int
  grid_resolution_factor = 1/4
    .type = float
  d_min = None
    .type = float
    .short_caption = High resolution
  initial_d_min = None
    .type = float
    .short_caption = Initial high resolution
    .help = If this is specified, phase extension will be performed \
            from initial_d_min to d_min
  phase_extension = False
    .type = bool
  verbose = True
    .type = bool
  change_basis_to_niggli_cell = True
    .type = bool
    .short_caption = Change basis to Niggli cell
  ncs_averaging = False
    .type = bool
    .help = This functionality is not yet implemented!
    .expert_level = 3
    .style = hidden
  protein_solvent_ratio = 1.31
    .type = float
    .short_caption = Protein/solvent ratio
  density_truncation
    .style = box
  {
    fraction_min = 0.35
      .type = float
    fraction_max = None
      .type = float
  }
  solvent_modification {
    method = *flipping flattening
      .type = choice
      .short_caption = Solvent modification method
    scale_flip = True
      .type = bool
  }
  solvent_adjust = True
    .type = bool
    .short_caption = Adjust solvent
  solvent_mask {
    averaging_radius {
      initial = None
        .type = float
        .short_caption = Initial averaging radius
      final = None
        .type = float
        .short_caption = Final averaging radius
    }
  }
  anisotropic_correction = True
    .type = bool
    .help = Correct the observations for anisotropy
  asu_contents
    .help = "Defines the ASU contents"
    .short_caption = ASU contents
    .style = menu_item auto_align
  {
    n_residues=None
      .type=float
      .help="Number of residues in structural unit"
      .short_caption = Number of residues in asymmetric unit
    n_bases=None
      .type=float
      .help="Number of nucleotides in structural unit"
      .short_caption = Number of nucleotides in asymmetric unit
    n_copies_per_asu=None
      .type=float
      .help="Number of copies per ASU. If not specified, Matthews analyses is performed"
      .short_caption = Number of copies in asymmetric unit
  }

"""


def rms(flex_double):
  return math.sqrt(flex.mean(flex.pow2(flex_double)))


class local_standard_deviation_map(object):

  def __init__(self, map_coeffs, radius,
               mean_solvent_density=0,
               method=0,
               resolution_factor=1/3):
    assert map_coeffs.is_complex_array()
    self.map = map_coeffs.local_standard_deviation_map(
      radius, mean_solvent_density=mean_solvent_density,
      resolution_factor=resolution_factor)
    self.map = self.map.real_map_unpadded()

  def histogram(self, n_slots=10000):
    return flex.histogram(data=self.map.as_1d(), n_slots=n_slots)

  def mask(self, solvent_fraction):
    hist = self.histogram()
    cutoff = hist.get_cutoff(int(self.map.size()*(1-solvent_fraction)))
    mask = flex.size_t()
    mask.resize(self.map.accessor(), 1)
    mask.set_selected(self.map > cutoff, 0)
    return mask


class density_modification(object):

  def __init__(self,
               params,
               f_obs,
               hl_coeffs_start,
               map_coeffs=None,
               log=None,
               as_gui_program=False):
    if log is None: log = sys.stdout
    adopt_init_args(self, locals())
    if self.params.solvent_mask.averaging_radius.final is None:
      if self.params.initial_d_min is not None:
        self.params.solvent_mask.averaging_radius.final = self.params.initial_d_min
      elif self.params.d_min is not None:
        self.params.solvent_mask.averaging_radius.final = self.params.d_min
      else:
        self.params.solvent_mask.averaging_radius.final = self.f_obs.d_min()
    if self.params.solvent_mask.averaging_radius.initial is None:
      self.params.solvent_mask.averaging_radius.initial = \
         self.params.solvent_mask.averaging_radius.final + 1
    self.matthews_analysis()
    self.anisotropic_correction()
    self.change_of_basis_op = None
    if self.params.change_basis_to_niggli_cell:
      self.change_of_basis_op = self.f_obs.change_of_basis_op_to_niggli_cell()
      if self.change_of_basis_op.is_identity_op():
        self.change_of_basis_op = None
    if self.change_of_basis_op is not None:
      self.f_obs = self.f_obs.change_basis(self.change_of_basis_op).map_to_asu()
      self.hl_coeffs_start = self.hl_coeffs_start.change_basis(
        self.change_of_basis_op).map_to_asu()
      if self.map_coeffs is not None:
        self.map_coeffs = self.map_coeffs.change_basis(
          self.change_of_basis_op).map_to_asu()
    self.mean_solvent_density = 0
    self.phase_source_initial = None
    self.phase_source = None
    if self.params.d_min is None:
      if self.params.phase_extension:
        self.params.d_min = self.f_obs.d_min()
      else:
        self.params.d_min = self.hl_coeffs_start.d_min()
    if self.params.initial_d_min is None:
      self.params.initial_d_min = self.params.d_min
    assert self.params.initial_d_min >= self.params.d_min
    self.max_iterations = sum((self.params.initial_steps,
                               self.params.shrink_steps,
                               self.params.final_steps))
    self.i_cycle = 0
    if self.params.shrink_steps is not None and self.params.shrink_steps > 0:
      self.radius_delta = (self.params.solvent_mask.averaging_radius.initial
                           - self.params.solvent_mask.averaging_radius.final) \
          / self.params.shrink_steps
      if self.params.phase_extension:
        self.phase_extend_step = (
          self.params.initial_d_min - self.params.d_min)/self.params.shrink_steps
      else:
        self.phase_extend_step = 0
    self.complete_set = self.f_obs.complete_set()

    ref_active = (self.f_obs.sigmas() > 0) \
               & (self.f_obs.d_spacings().data() >= self.params.d_min)

    sigma_cutoff = 0
    obs_rms = 1e4
    obs_high = rms(self.f_obs.select(ref_active).data()) * obs_rms
    obs_low = flex.min(self.f_obs.select(ref_active).data())
    self.ref_flags_array = self.f_obs.array(data=(
      (self.f_obs.data() > sigma_cutoff*self.f_obs.sigmas())
      & (self.f_obs.data() >= obs_low)
      & (self.f_obs.data() <= obs_high)
      & (self.f_obs.d_spacings().data() > self.params.d_min)))
    # now setup for complete arrays
    self.ref_flags_array = self.ref_flags_array.complete_array(
      new_data_value=False, d_min=self.params.d_min)
    self.ref_flags = self.ref_flags_array.data()
    self.f_obs_complete = self.f_obs.complete_array(
      new_data_value=0, new_sigmas_value=0, d_min=self.params.d_min)
    self.hl_coeffs_start = self.hl_coeffs_start.complete_array(
      new_data_value=(0,0,0,0), d_min=self.params.d_min)

    self.hl_coeffs = self.hl_coeffs_start.select_indices(self.active_indices)
    hl_coeffs = self.hl_coeffs
    self.compute_phase_source(hl_coeffs)
    fom = flex.abs(self.phase_source.data())
    fom.set_selected(hl_coeffs.data() == (0,0,0,0), 0)
    self.fom = fom

    if self.map_coeffs is None:
      self.map_coeffs = self.f_obs_active.customized_copy(
        data=self.f_obs_active.data()*fom,
        sigmas=None).phase_transfer(phase_source=self.hl_coeffs)
      self.map_coeffs.data().set_selected(fom <= 0, 0)
      self.map = self.map_coeffs.select(fom > 0).fft_map(
        resolution_factor=self.params.grid_resolution_factor
        ).apply_volume_scaling().real_map_unpadded()
      self.map_coeffs = self.map_coeffs.select(fom > 0)
    else:
      assert self.map_coeffs.is_complex_array()
      self.map = self.map_coeffs.fft_map(
        resolution_factor=self.params.grid_resolution_factor
        ).apply_volume_scaling().real_map_unpadded()
    self.map_coeffs_start = self.map_coeffs
    self.calculate_solvent_mask()

    n_phased = (fom > 0).count(True)
    if params.verbose:
      summary = "n phased: %d\n" % n_phased
      summary += "Mean solvent density: %.4f\n" %self.mean_solvent_density
      summary += "Mean protein density: %.4f\n" %self.mean_protein_density
      summary += "RMS solvent density: %.4f\n" %self.rms_solvent_density
      summary += "RMS protein density: %.4f\n" %self.rms_protein_density
      summary += "RMS solvent/protein density ratio: %.4f\n" %(
        self.rms_solvent_density/self.rms_protein_density)
      summary += "F000/V: %.4f\n" %self.f000_over_v
      summary += "Mean FOM: %.4f\n" %flex.mean(fom.select(fom>0))
      print >> self.log, summary
      libtbx.call_back(message="summary", data=summary)
    # XXX initialize printable statistics
    self.truncate_min = None
    self.truncate_min_percent = None
    self.truncate_max = None
    self.truncate_max_percent = None
    self.k_flip = None
    self.solvent_add = None
    self.truncate_density = \
      (self.params.density_truncation.fraction_max is not None or
       self.params.density_truncation.fraction_min is not None)
    self._stats = dm_stats()
    self._stats.add_cycle(
      cycle=0,
      mean_solvent_density=self.mean_solvent_density,
      mean_protein_density=self.mean_protein_density,
      f000_over_v=self.f000_over_v,
      rms_solvent_density=self.rms_solvent_density,
      rms_protein_density=self.rms_protein_density,
      fom=flex.mean(fom.select(fom>0)))

    libtbx.call_back("start_progress_bar",
        data=group_args(label="Running %d cycles..." % self.max_iterations,
                        size=self.max_iterations))
    for self.i_cycle in range(self.max_iterations):
      self.next_cycle()
      libtbx.call_back(message="increment_progress_bar",
        data=group_args(chunk=1),
        cached=False)
    libtbx.call_back("end_progress_bar", data=None)

  def get_stats (self) :
    return self._stats

  def next_cycle(self):
    self.ncs_averaging()
    self.calculate_solvent_mask()
    self.density_truncation()
    self.solvent_flipping()
    self.solvent_flattening()
    self.solvent_adjust()
    self.compute_map_coefficients()
    self.compute_map()
    self.show_cycle_summary()

  def matthews_analysis(self):
    from mmtbx.scaling import matthews
    self.matthews_result = matthews.matthews_rupp(
      miller_array=self.f_obs,
      n_residues=self.params.asu_contents.n_residues,
      n_bases=self.params.asu_contents.n_bases,
      out=self.log, verbose=1)
    self.params.asu_contents.n_residues = self.matthews_result[0]
    self.params.asu_contents.n_bases = self.matthews_result[1]
    if self.params.asu_contents.n_copies_per_asu is None:
      self.params.asu_contents.n_copies_per_asu = self.matthews_result[2]
    if self.params.solvent_fraction is None:
      self.params.solvent_fraction = self.matthews_result[3]

  def anisotropic_correction(self):
    if not self.params.anisotropic_correction: return
    self.aniso_scale_and_b = None
    n_copies_solc = self.params.asu_contents.n_copies_per_asu
    n_residues = self.params.asu_contents.n_residues
    n_bases = self.params.asu_contents.n_bases
    self.aniso_scale_and_b = absolute_scaling.ml_aniso_absolute_scaling(
      miller_array=self.f_obs,
      n_residues=n_residues*self.f_obs.space_group().order_z()*n_copies_solc,
      n_bases=n_bases*self.f_obs.space_group().order_z()*n_copies_solc)
    self.aniso_scale_and_b.show(out=self.log,verbose=1)
    b_cart = self.aniso_scale_and_b.b_cart
    trace = sum(b_cart[:3])/3
    b_cart = [b_cart[0]-trace, b_cart[1]-trace, b_cart[2]-trace,
              b_cart[3], b_cart[4], b_cart[5]]
    u_star = adptbx.u_cart_as_u_star(self.f_obs.unit_cell(), adptbx.b_as_u(b_cart))
    self.f_obs = absolute_scaling.anisotropic_correction(
      self.f_obs, 0.0, u_star).set_observation_type(self.f_obs)

  def compute_phase_source(self, hl_coeffs, n_steps=72):
    integrator = miller.phase_integrator(n_steps=n_steps)
    self.phase_source_previous = self.phase_source
    self.phase_source = miller.array(
      miller_set=hl_coeffs,
      data=integrator(
        space_group=hl_coeffs.space_group(),
        miller_indices=hl_coeffs.indices(),
        hendrickson_lattman_coefficients=hl_coeffs.data()))
    if self.phase_source_initial is None:
      self.phase_source_initial = self.phase_source
    return self.phase_source

  def compute_map(self):
    self.map = self.map_coeffs.fft_map(
      resolution_factor=self.params.grid_resolution_factor
      ).apply_volume_scaling().real_map_unpadded()

  def calculate_solvent_mask(self):
    # calculate mask
    lsd = local_standard_deviation_map(
      self.map_coeffs,
      self.radius,
      mean_solvent_density=self.mean_solvent_density,
      resolution_factor=self.params.grid_resolution_factor,
      method=2)
    self.rms_map = lsd.map
    self.mask = lsd.mask(self.params.solvent_fraction)
    # setup solvent/protein selections
    self.solvent_selection = (self.mask == 1)
    self.protein_selection = (self.mask == 0)
    self.solvent_iselection = self.solvent_selection.iselection()
    self.protein_iselection = self.protein_selection.iselection()
    self.n_solvent_grid_points = self.mask.count(1)
    self.n_protein_grid_points = self.mask.count(0)
    # map statistics
    self.mean_protein_density = self.mean_protein_density_start = flex.mean(
      self.map.select(self.protein_iselection))
    self.mean_solvent_density = self.mean_solvent_density_start = flex.mean(
      self.map.select(self.solvent_iselection))
    self.mask_percent = self.n_solvent_grid_points/(self.mask.size()) * 100
    self.f000_over_v = ((
      (1/self.params.protein_solvent_ratio) * self.mean_protein_density)
                        - self.mean_solvent_density) \
        * (self.params.protein_solvent_ratio/(self.params.protein_solvent_ratio-1))
    self.rms_protein_density = rms(self.map.select(self.protein_iselection))
    self.rms_solvent_density = rms(self.map.select(self.solvent_iselection))
    self.standard_deviation_local_rms = flex.mean_and_variance(
      lsd.map.as_1d()).unweighted_sample_standard_deviation()

  def density_truncation(self):
    min_fraction = self.params.density_truncation.fraction_min
    max_fraction = self.params.density_truncation.fraction_max
    if min_fraction is None and max_fraction is None: return
    if min_fraction is Auto:
      min_fraction = self.mean_protein_density-self.f000_over_v
    hist = flex.histogram(
      self.map.select(self.protein_iselection), n_slots=10000)
    if max_fraction is not None:
      self.truncate_max = hist.get_cutoff(
        int(self.n_protein_grid_points * (1-max_fraction)))
      truncate_max_sel = (self.map > self.truncate_max) & self.protein_selection
      self.map.set_selected(truncate_max_sel, self.truncate_max)
      self.truncate_max_percent = (
        truncate_max_sel.count(True) / self.n_protein_grid_points) * 100
    if min_fraction is not None:
      self.truncate_min = hist.get_cutoff(
        int(self.n_protein_grid_points * (1-min_fraction)))
      truncate_min_sel = (self.map < self.truncate_min) & self.protein_selection
      self.map.set_selected(truncate_min_sel, self.truncate_min)
      self.truncate_min_percent = (
        truncate_min_sel.count(True) / self.n_protein_grid_points) * 100
    self.mean_protein_density = flex.mean(
      self.map.select(self.protein_iselection))

  def solvent_flipping(self):
    if not self.params.solvent_modification.method == "flipping": return
    if (self.i_cycle + 1) == self.max_iterations:
      self.k_flip = 0
    else:
      self.k_flip = -(1-self.params.solvent_fraction)/self.params.solvent_fraction
      if self.params.solvent_modification.scale_flip:
        rms_protein_density_new = math.sqrt(
          flex.mean(flex.pow2(self.map.select(self.protein_iselection))))
        self.k_flip *= math.pow(
          rms_protein_density_new/self.rms_protein_density, 2)
    self.map.as_1d().copy_selected(
      self.solvent_iselection,
      (self.mean_solvent_density
       + self.k_flip * (self.map - self.mean_solvent_density)).as_1d())
    self.mean_solvent_density = flex.mean(
      self.map.select(self.solvent_iselection))

  def solvent_flattening(self):
    if not self.params.solvent_modification.method == "flattening": return
    self.map.set_selected(self.solvent_selection, self.mean_solvent_density)

  def solvent_adjust(self):
    if not self.params.solvent_adjust: return
    min_solvent_density = flex.min(self.map.select(self.solvent_iselection))
    min_protein_density = flex.min(self.map.select(self.protein_iselection))
    self.solvent_add = ((self.mean_protein_density-min_protein_density)
                   /self.params.protein_solvent_ratio) \
                + min_protein_density - self.mean_solvent_density
    self.map.as_1d().copy_selected(
      self.solvent_iselection, (self.map + self.solvent_add).as_1d())
    #self.mean_solvent_density = flex.mean(self.map.select(self.solvent_iselection))
    self.mean_solvent_density = (1-self.params.solvent_fraction) \
        * (self.mean_solvent_density+self.solvent_add-self.mean_protein_density)

  def compute_map_coefficients(self):
    f_obs = self.f_obs_complete.select(self.f_obs_complete.d_spacings().data() >= self.d_min)
    f_calc = f_obs.structure_factors_from_map(self.map, use_sg=True)
    f_obs_active = f_obs.select_indices(self.active_indices)
    minimized = relative_scaling.ls_rel_scale_driver(
      f_obs_active,
      f_calc.as_amplitude_array().select_indices(self.active_indices),
      use_intensities=False,
      use_weights=False)
    #minimized.show()
    f_calc = f_calc.customized_copy(data=f_calc.data()\
                                    * math.exp(-minimized.p_scale)\
                                    * adptbx.debye_waller_factor_u_star(
                                      f_calc.indices(), minimized.u_star))
    f_calc_active = f_calc.common_set(f_obs_active)
    matched_indices = f_obs.match_indices(self.f_obs_active)
    lone_indices_selection = matched_indices.single_selection(0)
    from mmtbx.max_lik import maxlik
    alpha_beta_est = maxlik.alpha_beta_est_manager(
      f_obs=f_obs_active,
      f_calc=f_calc_active,
      free_reflections_per_bin=140,
      flags=flex.bool(f_obs_active.size()),
      interpolation=True)
    alpha, beta = alpha_beta_est.alpha_beta_for_each_reflection(
      f_obs=self.f_obs_complete.select(self.f_obs_complete.d_spacings().data() >= self.d_min))
    f_obs.data().copy_selected(
      lone_indices_selection.iselection(), flex.abs(f_calc.data()))
    t = maxlik.fo_fc_alpha_over_eps_beta(
      f_obs=f_obs,
      f_model=f_calc,
      alpha=alpha,
      beta=beta)
    hl_coeff = flex.hendrickson_lattman(
      t * flex.cos(f_calc.phases().data()),
      t * flex.sin(f_calc.phases().data()))
    dd = alpha.data()
    #
    hl_array = f_calc.array(
      data=self.hl_coeffs_start.common_set(f_calc).data()+hl_coeff)
    self.compute_phase_source(hl_array)
    fom = flex.abs(self.phase_source.data())
    mFo = hl_array.array(data=f_obs.data()*self.phase_source.data())
    DFc = hl_array.array(data=dd*f_calc.as_amplitude_array().phase_transfer(
        self.phase_source).data())
    centric_flags = f_obs.centric_flags().data()
    acentric_flags = ~centric_flags
    fo_scale = flex.double(centric_flags.size())
    fc_scale = flex.double(centric_flags.size())
    fo_scale.set_selected(acentric_flags, 2)
    fo_scale.set_selected(centric_flags, 1)
    fc_scale.set_selected(acentric_flags, 1)
    fc_scale.set_selected(centric_flags, 0)
    fo_scale.set_selected(lone_indices_selection, 0)
    fc_scale.set_selected(lone_indices_selection, -1)
    self.map_coeffs = hl_array.array(
      data=mFo.data()*fo_scale - DFc.data()*fc_scale)
    self.fom = hl_array.array(data=fom)
    self.hl_coeffs = hl_array
    # statistics
    self.r1_factor = f_obs_active.r1_factor(f_calc_active)
    fom = fom.select(matched_indices.pair_selection(0))
    self.r1_factor_fom = flex.sum(
      fom * flex.abs(f_obs_active.data() - f_calc_active.as_amplitude_array().data())) \
        / flex.sum(fom * f_obs_active.data())
    phase_source, phase_source_previous = self.phase_source.common_sets(
      self.phase_source_previous)
    self.mean_delta_phi = phase_error(
      flex.arg(phase_source.data()), flex.arg(phase_source_previous.data()))
    phase_source, phase_source_initial = self.phase_source.common_sets(
      self.phase_source_initial)
    self.mean_delta_phi_initial = phase_error(
      flex.arg(phase_source.data()), flex.arg(phase_source_initial.data()))
    self.mean_fom = flex.mean(fom)
    fom = f_obs_active.array(data=fom)
    fom.setup_binner(reflections_per_bin=1000)
    self.mean_fom_binned = fom.mean(use_binning=True)

  def show_cycle_summary(self, out=None):
    if not self.params.verbose: return
    if out is None: out = sys.stdout
    self.more_statistics = maptbx.more_statistics(self.map)
    self._stats.add_cycle(
      cycle=self.i_cycle+1,
      radius=self.radius,
      d_min=self.d_min,
      mask_percent=self.mask_percent,
      mean_solvent_density=self.mean_solvent_density,
      mean_protein_density=self.mean_protein_density,
      f000_over_v=self.f000_over_v,
      truncate_density=self.truncate_density,
      truncate_min=self.truncate_min,
      truncate_min_percent=self.truncate_min_percent,
      truncate_max=self.truncate_max,
      truncate_max_percent=self.truncate_max_percent,
      k_flip=self.k_flip,
      solvent_add=self.solvent_add,
      rms_solvent_density=self.rms_solvent_density,
      rms_protein_density=self.rms_protein_density,
      standard_deviation_local_rms=self.standard_deviation_local_rms,
      mean_delta_phi=flex.mean(self.mean_delta_phi)/pi_180,
      mean_delta_phi_initial=flex.mean(self.mean_delta_phi_initial)/pi_180,
      r1_factor=self.r1_factor,
      r1_factor_fom=self.r1_factor_fom,
      fom=self.mean_fom,
      fom_binned=self.mean_fom_binned,
      skewness=self.more_statistics.skewness())
    summary = self._stats.format_summary()
    print >> self.log, summary
    self.log.flush()
    if (not self.as_gui_program) :
      libtbx.call_back(message="summary",
        data=summary,
        accumulate=True)
    else :
      libtbx.call_back(message="plot_current_stats",
        data=self._stats.get_fom_for_plot())

  def ncs_averaging(self):
    if not self.params.ncs_averaging: return
    else: raise NotImplementedError

  class f_obs_active(libtbx.property):
    def fget(self):
      return self.f_obs_complete.select_indices(self.active_indices)

  class active_indices(libtbx.property):
    def fget(self):
      if self.i_cycle == 0 and not hasattr(self, "_active_indices"):
        self._active_indices = None
      if (self.params.d_min < self.params.initial_d_min and
          self.i_cycle > self.params.initial_steps and
          self.i_cycle < (self.params.initial_steps + self.params.shrink_steps)):
        self._active_indices = None
      if self._active_indices is None:
        sel = (self.ref_flags_array.d_spacings().data() >= self.d_min) & self.ref_flags
        self._active_indices = self.ref_flags_array.select(sel).indices()
      return self._active_indices

  def miller_array_in_original_setting(self, miller_array):
    if self.change_of_basis_op is not None:
      return miller_array.change_basis(self.change_of_basis_op.inverse())
    return miller_array

  class map_coeffs_in_original_setting(libtbx.property):
    def fget(self):
      return self.miller_array_in_original_setting(self.map_coeffs)

  class hl_coeffs_in_original_setting(libtbx.property):
    def fget(self):
      return self.miller_array_in_original_setting(self.hl_coeffs)

  class radius(libtbx.property):
    def fget(self):
      if self.i_cycle == 0 or self.i_cycle < self.params.initial_steps:
        return self.params.solvent_mask.averaging_radius.initial
      elif self.i_cycle < (self.params.initial_steps + self.params.shrink_steps):
        return (self.params.solvent_mask.averaging_radius.initial -
                (self.radius_delta * (self.i_cycle - self.params.initial_steps + 1)))
      else:
        return self.params.solvent_mask.averaging_radius.final

  class d_min(libtbx.property):
    def fget(self):
      if self.i_cycle == 0 or self.i_cycle < self.params.initial_steps:
        return self.params.initial_d_min
      elif self.i_cycle < (self.params.initial_steps + self.params.shrink_steps):
        return (self.params.initial_d_min -
                (self.phase_extend_step * (self.i_cycle - self.params.initial_steps + 1)))
      else:
        return self.params.d_min

class dm_stats (object) :
  def __init__ (self) :
    self._stats = []

  def add_cycle (self, **kwds) :
    cycle_stats = group_args(**kwds)
    self._stats.append(cycle_stats)

  def get_cycle_stats (self, i_cycle=-1) :
    return self._stats[i_cycle]

  def extract_loggraph (self) :
    table = iotbx.data_plots.table_data(
      title="Density modification statistics by cycle",
      column_names=["cycle","fom","mean_protein_density","mean_solvent_density",
                    "rms_protein_density","rms_solvent_density"],
      column_labels=["Cycle", "Figure of Merit", "Mean protein density",
        "Mean solvent density", "RMS protein density", "RMS solvent density"],
      graph_names=["FOM vs. cycle", "Mean density vs. cycle",
        "RMS density vs. cycle"],
      graph_columns=[[0,1],[0,2,3],[0,4,5]])
    for stats in self._stats :
      table.add_row([
        stats.cycle,
        stats.fom,
        stats.mean_protein_density,
        stats.mean_solvent_density,
        stats.rms_protein_density,
        stats.rms_solvent_density])
    return table

  def format_summary (self, i_cycle=-1) :
    stats = self._stats[i_cycle]
    summary = "#"*80 + "\n"
    summary += "Cycle %i\n" %(stats.cycle)
    summary += "d min: %.2f\n" % stats.d_min
    summary += "Mask averaging radius: %.2f\n" % stats.radius
    summary += "Solvent mask volume (%%): %.4f\n" % stats.mask_percent
    summary += "Mean solvent density: %.4f\n" % stats.mean_solvent_density
    summary += "Mean protein density: %.4f\n" % stats.mean_protein_density
    summary += "F000/V: %.4f\n" % stats.f000_over_v
    if (stats.truncate_density) :
      summary += "Protein density truncation:\n"
      if (stats.truncate_min is not None) :
        summary += "  min = %7.4f (%.2f%%)\n" %(
          stats.truncate_min, stats.truncate_min_percent)
      if (stats.truncate_max is not None) :
        summary += "  max = %7.4f (%.2f%%)\n" %(
          stats.truncate_max, stats.truncate_max_percent)
    if (stats.k_flip is not None) :
      summary += "Solvent flipping factor: %.4f\n" %stats.k_flip
    if (stats.solvent_add is not None) :
      summary += "Solvent level raised by: %.4f\n" % stats.solvent_add
    summary += "RMS solvent density: %.4f\n" % stats.rms_solvent_density
    summary += "RMS protein density: %.4f\n" % stats.rms_protein_density
    summary += "RMS solvent/protein density ratio: %.4f\n" %(
      stats.rms_solvent_density/stats.rms_protein_density)
    summary += "Standard deviation (local RMS): %.4f\n" %(
      stats.standard_deviation_local_rms)
    summary += "Mean delta phi: %.4f\n" % stats.mean_delta_phi
    summary += "Mean delta phi (initial): %.4f\n" %stats.mean_delta_phi_initial
    summary += "R1-factor:       %.2f\n" % stats.r1_factor
    summary += "R1-factor (fom): %.2f\n" % stats.r1_factor_fom
    summary += "Skewness: %.4f\n" % stats.skewness
    summary += "Mean figure of merit (FOM): %.3f\n" % stats.fom
    s = StringIO()
    summary += "Mean figure of merit (FOM) in resolution bins:\n"
    stats.fom_binned.show(f=s, data_fmt="%.4f")
    summary += s.getvalue()
    summary += "#"*80 + "\n"
    summary += "\n"
    return summary

  def get_fom_for_plot (self) :
    return [ stats.fom for stats in self._stats ]
