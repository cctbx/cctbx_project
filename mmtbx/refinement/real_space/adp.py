from __future__ import absolute_import, division, print_function
import mmtbx.refinement.real_space.utils
import mmtbx.refinement.utils
from scitbx.array_family import flex
from cctbx import adptbx
from mmtbx import bulk_solvent
from libtbx.test_utils import approx_equal
from six.moves import range
from cctbx import crystal
from cctbx import xray
from cctbx import maptbx
from mmtbx.refinement import adp_refinement
from cctbx import adp_restraints
from libtbx import group_args
from libtbx.str_utils import format_value
import scitbx.minimizers

import boost_adaptbx.boost.python as bp
cctbx_maptbx_ext = bp.import_ext("cctbx_maptbx_ext")

import random
random.seed(0)
flex.set_random_seed(0)

# Termination of L-BFGS-B minimization of B-factors: stop when the relative
# decrease of the target in one iteration is below lbfgsb_factr*(machine
# precision), about 2.e-6. Target and gradients are computed via FFT, with
# limited accuracy; with the L-BFGS-B default (1.e+7) the minimizer went on
# for several times more evaluations in numerical noise, changing neither
# R-factor (differences below 0.0002) nor B-factors (within a fraction of A**2).
lbfgsb_factr = 1.e+10
# No termination on the size of the gradient: it is an absolute number, and
# target and gradient are scaled (see tg).
lbfgsb_pgtol = 0

class tg_vrm(object):
  def __init__(self,
               vrm,
               map_data,
               x,
               restraints_weight,
               bound_flags,
               lower_bound,
               upper_bound):
    self.bound_flags = bound_flags
    self.lower_bound = lower_bound
    self.upper_bound = upper_bound
    self.restraints_weight = restraints_weight
    self.vrm = vrm
    self.map_data = map_data
    self.plain_pair_sym_table = None
    self.adp_iso_params = None
    if(self.restraints_weight is not None):
      self.plain_pair_sym_table = \
        mmtbx.refinement.real_space.adp.get_plain_pair_sym_table(
          crystal_symmetry = self.vrm.crystal_symmetry,
          sites_frac       = self.vrm.xray_structure.sites_frac())
      self.adp_iso_params = \
        adp_refinement.adp_restraints_master_params.extract().iso
    # required fields
    self.x = x
    self.t = None
    self.g = None
    self.d = None
    self.use_curvatures=False
    #
    self.weight = None
    if(self.restraints_weight is not None):
      self.weight = self._weight()
      if abs(self.weight)<1.e-6: self.weight = 1  #XXX
    self.n_evaluations = 0
    self.update_target_and_grads(x=x)

  def _weight(self):
    R = self._restraints()
    num = R.gradients.norm()
    den = self._data().gradients(map_data = self.map_data).grad_uiso.norm()
    if(den==0): return 1
    return num/den

  def _restraints(self):
    return adp_restraints.energies_iso(
      plain_pair_sym_table = self.plain_pair_sym_table,
      xray_structure       = self.vrm.xray_structure,
      parameters           = self.adp_iso_params,
      compute_gradients    = True,
      use_u_local_only     = self.adp_iso_params.use_u_local_only,
      use_hd               = False)

  def _data(self):
    return self.vrm

  def _compute(self, x):
    self.n_evaluations += 1
    for bs, xi in zip(self.vrm.bcr_scatterers, x):
      bs.scatterer.u_iso = xi
    self.vrm.update()
    D = self._data()
    DG = D.gradients(map_data = self.map_data).grad_uiso
    DT = D.target()
    if(self.restraints_weight is not None):
      R = self._restraints()
      self.tgo = group_args(
        target   = DT*self.weight + R.target   *self.restraints_weight,
        gradient = DG*self.weight + R.gradients*self.restraints_weight)
    else:
      self.tgo = group_args(
        target   = DT,
        gradient = DG)
    return self.tgo

  def update(self, x):
    self.update_target_and_grads(x = x)

  def update_target_and_grads(self, x):
    self.x = x
    self.tgo = self._compute(x=self.x)
    self.t = self.tgo.target
    self.g = self.tgo.gradient

  def target(self): return self.t

  def gradients(self): return self.g

  def gradient(self): return self.gradients()

# Amplitudes ("Fobs") for refinement of B-factors and occupancies are made
# from the map in a box around the model, with the map outside atom_radius of
# the atoms masked out (density of other molecules, of what is not modelled,
# noise). The mask also cuts off what the atoms of the model contribute
# outside it, while structure factors calculated from the model have all of
# it, and refined parameters make up for the difference (B-factors off by
# tens of A**2 at 4-6 A, occupancies off by 0.1-0.2). With
# fill_outside_mask the map outside the mask is replaced with the map
# calculated from the present model, scaled to the observed one inside the
# mask: the two sides are then compared like with like, and the data are
# updated when the model has changed (update_f_obs).
fill_outside_mask = True
# Group B refinement without individual one after it: stop when a round
# lowers the R-factor by less than this
group_min_gain = 0.002
box_cushion_extra = 0.

class data_from_map(object):
  """
  Map in a box around the model, mask around atoms, and amplitudes made from
  them: see fill_outside_mask.
  """
  def __init__(self, map_data, xray_structure, atom_radius, d_min):
    self.d_min = d_min
    self.box = mmtbx.utils.extract_box_around_model_and_map(
      xray_structure = xray_structure,
      map_data       = map_data,
      box_cushion    = atom_radius + box_cushion_extra)
    xb = self.box.xray_structure_box
    # (the structure in the box has the scattering table of the input one)
    self.xray_structure = xb.deep_copy_scatterers()
    self.map_box = self.box.map_box
    assert self.map_box.origin() == (0,0,0)
    self.mask = cctbx_maptbx_ext.mask(
      sites_frac                  = xb.sites_frac(),
      unit_cell                   = xb.unit_cell(),
      n_real                      = self.map_box.all(),
      mask_value_inside_molecule  = 1,
      mask_value_outside_molecule = 0,
      radii                       = flex.double(xb.scatterers().size(),
                                                atom_radius))
    self.inside = self.mask.as_1d() > 0.5
    self.crystal_gridding = maptbx.crystal_gridding(
      unit_cell             = xb.unit_cell(),
      space_group_info      = xb.space_group_info(),
      pre_determined_n_real = self.map_box.all())
    self.scale_and_offset = None
    # Observed map limited to d_min, as the calculated one is: masking mixes
    # Fourier coefficients, so what the map has beyond d_min would end up in
    # the amplitudes within d_min, with no counterpart on the model side.
    if(fill_outside_mask):
      self.map_box = self._coefficients(map_box = self.map_box).fft_map(
        crystal_gridding = self.crystal_gridding).real_map_unpadded()

  def _coefficients(self, map_box):
    self.box.map_box = map_box
    return self.box.box_map_coefficients(d_min = self.d_min)

  def f_obs_complex(self, f_calc = None):
    """
    Map coefficients of the masked map; with f_calc (of the model in the
    box, any scale) the map outside the mask is the one calculated from
    them, scaled to the observed map inside the mask.
    """
    if(f_calc is None or not fill_outside_mask):
      return self._coefficients(map_box = self.map_box * self.mask)
    rho_c = f_calc.fft_map(
      crystal_gridding = self.crystal_gridding).real_map_unpadded()
    x = rho_c.as_1d().select(self.inside)
    y = self.map_box.as_1d().select(self.inside)
    if(x.size() < 2):
      return self._coefficients(map_box = self.map_box * self.mask)
    xm, ym = flex.mean(x), flex.mean(y)
    den = flex.sum((x-xm)*(x-xm))
    k = 0
    if(den > 0): k = flex.sum((x-xm)*(y-ym))/den
    if(k <= 0): # no positive correlation: nothing sensible to fill with
      return self._coefficients(map_box = self.map_box * self.mask)
    c = ym - k*xm
    self.scale_and_offset = (k, c)
    filled = (self.map_box - c)*self.mask + rho_c*(1 - self.mask)*k
    return self._coefficients(map_box = filled)

def update_f_obs(fmodel):
  """
  Renew amplitudes in fmodel (made by map_and_model_to_fmodel) after its
  model changed: see fill_outside_mask.
  """
  data = getattr(fmodel, "data_from_map", None)
  if(data is None or not fill_outside_mask): return
  f_obs = abs(data.f_obs_complex(f_calc = fmodel.f_calc()))
  assert f_obs.indices().all_eq(fmodel.f_obs().indices())
  fmodel.update(f_obs = f_obs)
  fmodel.update_all_scales(update_f_part1=False, apply_back_trace=True,
    remove_outliers=False)

def map_and_model_to_fmodel(map_data, xray_structure, atom_radius, d_min,
                            reset_adp=True):
  data = data_from_map(
    map_data       = map_data,
    xray_structure = xray_structure,
    atom_radius    = atom_radius,
    d_min          = d_min)
  f_obs_complex = data.f_obs_complex()
  f_obs = abs(f_obs_complex)
  if(flex.mean(f_obs.data())<1.e-6): return None
  xrs = data.xray_structure
  if(reset_adp):
    vals_init = xrs.extract_u_iso_or_u_equiv()
    xrs = xrs.set_b_iso(value=0)
    assert approx_equal(flex.mean(xrs.extract_u_iso_or_u_equiv()),0.)
    f_calc = f_obs.structure_factors_from_scatterers(
      xray_structure = xrs).f_calc()
    sc = flex.sum(abs(f_obs).data()*abs(f_calc).data())/ \
         flex.sum(abs(f_calc).data()*abs(f_calc).data())
    f_calc = f_calc.array(data = f_calc.data()*sc)
    o = bulk_solvent.complex_f_kb_scaled(
      f1      = f_obs_complex.data(),
      f2      = f_calc.data(),
      b_range = flex.double(range(5,505,5)),
      ss      = 1./flex.pow2(f_calc.d_spacings().data()) / 4.)
    xrs = xrs.set_b_iso(value=o.b())
    k_isotropic = flex.double(f_calc.data().size(), o.k())
    if(o.k()<1.e-6):
      k_isotropic = flex.double(f_calc.data().size(), 1)
      xrs.set_u_iso(values = vals_init)
  fmodel = mmtbx.f_model.manager(f_obs = f_obs, xray_structure = xrs)
  if(reset_adp):
    fmodel.update_core(k_isotropic = k_isotropic)
  fmodel.update(target_name="ls_wunit_k1")
  fmodel.update_all_scales(update_f_part1=False, apply_back_trace=True,
    remove_outliers=False)
  fmodel.data_from_map = data
  update_f_obs(fmodel = fmodel)
  return fmodel

def get_plain_pair_sym_table(crystal_symmetry, sites_frac, plain_pairs_radius=5):
  asu_mappings = crystal.symmetry.asu_mappings(crystal_symmetry,
    buffer_thickness = plain_pairs_radius)
  special_position_settings = crystal.special_position_settings(
    crystal_symmetry = crystal_symmetry)
  sites_cart = crystal_symmetry.unit_cell().orthogonalize(sites_frac)
  site_symmetry_table = special_position_settings.site_symmetry_table(
    sites_cart = sites_cart)
  asu_mappings.process_sites_frac(
    original_sites      = sites_frac,
    site_symmetry_table = site_symmetry_table)
  pair_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
  pair_asu_table.add_all_pairs(distance_cutoff = plain_pairs_radius)
  return pair_asu_table.extract_pair_sym_table()

# Hydrogens are not seen in the map and are not restrained: refined as
# atoms of their own their B is noise (many end at zero). They ride instead,
# as in reciprocal-space refinement: B(H) = hydrogen_b_scale * B(parent).
hydrogen_b_scale = 1.2

def hydrogen_riding(model):
  """
  Hydrogens of model and the atoms they are bonded to: from restraints if
  the model has them, otherwise the nearest non-hydrogen atom of the same
  residue. None if there is nothing to ride.
  """
  if(not model.has_hd()): return None
  table = model.xh_connectivity_table()
  if(table is None):
    table = []
    for rg in model.get_hierarchy().residue_groups():
      atoms = rg.atoms()
      hd = flex.bool([a.element_is_hydrogen() for a in atoms])
      if(hd.count(True) == 0 or hd.count(False) == 0): continue
      heavy = atoms.select((~hd).iselection())
      sites = heavy.extract_xyz()
      for h in atoms.select(hd.iselection()):
        d = (sites - h.xyz).dot()
        table.append([heavy[flex.min_index(d)].i_seq, h.i_seq])
  if(len(table) == 0): return None
  return group_args(
    i_x   = flex.size_t([t[0] for t in table]),
    i_h   = flex.size_t([t[1] for t in table]),
    scale = hydrogen_b_scale)

class tg(object):
  def __init__(self, fmodel, x, restraints_weight, bound_flags,
               lower_bound, upper_bound, riding = None):
    self.riding = riding
    self.bound_flags = bound_flags
    self.lower_bound = lower_bound
    self.upper_bound = upper_bound
    self.restraints_weight = restraints_weight
    self.fmodel = fmodel
    self.plain_pair_sym_table = None
    self.adp_iso_params = None
    if(self.restraints_weight is not None):
      self.plain_pair_sym_table = get_plain_pair_sym_table(
        crystal_symmetry = self.fmodel.xray_structure.crystal_symmetry(),
        sites_frac       = self.fmodel.xray_structure.sites_frac())
      self.adp_iso_params = \
        adp_refinement.adp_restraints_master_params.extract().iso
    self.fmodel.xray_structure.scatterers().flags_set_grads(state=False)
    self.fmodel.xray_structure.scatterers().flags_set_grad_u_iso(
      iselection = self.fmodel.xray_structure.all_selection().iselection())
    # required fields
    self.x = x
    self.t = None
    self.g = None
    self.d = None
    self.use_curvatures=False
    #
    self.weight = None
    if(self.restraints_weight is not None):
      self.weight = self._weight()
    self.n_evaluations = 0
    # Scale of target and gradients. All variables are bounded, and then
    # the first step of L-BFGS-B is minus the gradient itself: its size is
    # set by the units of the target. Too small, and the decrease of the
    # target in the first iteration is below the termination threshold (the
    # minimizer stops having done nothing); so the target is scaled to make
    # the largest first shift of a B-factor 1 A**2.
    self.scale = 1
    self.update_target_and_grads(x=x)
    g_max = flex.max(flex.abs(self.g)) if self.g.size() > 0 else 0
    if(g_max > 0):
      self.scale = 1./g_max
      self.t = self.t*self.scale
      self.g = self.g*self.scale

  def _weight(self):
    """
    Weight of the data term that makes its gradient as large as that of the
    restraints at the start. With all B equal restraints have zero gradient:
    no such weight (it would be zero, and nothing would be refined), 1 is
    used.
    """
    num = self._restraints().gradients.norm()
    den = self._data().gradient_xray.norm()
    if(den==0): return 1
    result = num/den
    if(result < 1.e-6): result = 1
    return result

  def _restraints(self):
    return adp_restraints.energies_iso(
      plain_pair_sym_table = self.plain_pair_sym_table,
      xray_structure       = self.fmodel.xray_structure,
      parameters           = self.adp_iso_params,
      compute_gradients    = True,
      use_u_local_only     = self.adp_iso_params.use_u_local_only,
      use_hd               = False)

  def _data(self):
    fmodels = mmtbx.fmodels(fmodel_xray = self.fmodel)
    return fmodels.target_and_gradients(compute_gradients=True)

  def _compute(self, x):
    self.n_evaluations += 1
    if(self.riding is not None):
      x = x.deep_copy().set_selected(self.riding.i_h,
        x.select(self.riding.i_x)*self.riding.scale)
    self.fmodel.xray_structure.set_b_iso(values = x)
    self.fmodel.update_xray_structure(update_f_calc = True)
    D = self._data()
    # Variables are B, gradients of data and restraints terms are with
    # respect to U: dT/dB = dT/dU / (8*pi**2). (Without this factor the
    # gradient given to the minimizer was 79 times too large for the target,
    # and its line search failed at times: trial minimizations that did not
    # move at all.)
    u_per_b = adptbx.b_as_u(1.)
    if(self.restraints_weight is not None):
      R = self._restraints()
      self.tgo = group_args(
        target   = D.target()*self.weight + R.target*self.restraints_weight,
        gradient = (D.gradient_xray*self.weight +
                    R.gradients*self.restraints_weight)*u_per_b)
    else:
      self.tgo = group_args(
        target   = D.target(),
        gradient = D.gradient_xray*u_per_b)
    if(self.riding is not None):
      # B(H) is a function of B(parent): its gradient goes to the parent
      g = self.tgo.gradient
      g.add_selected(self.riding.i_x,
        g.select(self.riding.i_h)*self.riding.scale)
      g.set_selected(self.riding.i_h, 0)
    return self.tgo

  def update(self, x):
    self.update_target_and_grads(x = x)

  def update_target_and_grads(self, x):
    self.x = x
    self.tgo = self._compute(x=self.x)
    self.t = self.tgo.target*self.scale
    self.g = self.tgo.gradient*self.scale

  def target(self): return self.t

  def gradients(self): return self.g

  def gradient(self): return self.gradients()

class ncs_aware_refinement(object):
  def __init__(self, map_model_manager, d_min, atom_radius, nproc=1,
               log = None, individual = True, restraints_weight = 1,
               group_mode = None):
    self.mmm               = map_model_manager
    self.nproc             = nproc
    self.d_min             = d_min
    self.atom_radius       = atom_radius
    self.log               = log
    self.individual        = individual
    self.restraints_weight = restraints_weight
    self.group_mode        = group_mode
    #
    ncs_groups = self.mmm.model().get_ncs_groups()
    if(ncs_groups is None or len(ncs_groups)==0):
      values = self.run_one()
      self.mmm.model().set_b_iso(values = values)
    else:
      values = self.mmm.model().get_b_iso()
      for i, g in enumerate(ncs_groups):
        values_g = self.run_one(selection = g.master_iselection)
        values = values.set_selected(g.master_iselection, values_g)
        for j, c in enumerate(g.copies):
          values = values.set_selected(c.iselection, values_g)
      self.mmm.model().set_b_iso(values = values)

  def run_one(self, selection=None):
    """
    Refine ADP of the whole model (or of the selected part, an NCS master)
    against the map in one calculation. nproc is not used: refining chain by
    chain in parallel gave different, and for some models clearly worse,
    results (each chain was fitted to density that also holds its neighbours),
    so there is one code path and the result does not depend on nproc.
    """
    model = self.mmm.model()
    if(selection is not None): model = model.select(selection)
    model.get_hierarchy().atoms().reset_i_seq()
    return self.run_one_one(args = [model,])

  def run_one_one(self, args):
    model = args[0]
    log = self.log
    #
    fmodel = map_and_model_to_fmodel(
      map_data       = self.mmm.map_data().deep_copy(),
      xray_structure = model.get_xray_structure(),
      atom_radius    = self.atom_radius,
      d_min          = self.d_min)
    if(fmodel is None):
      return model.get_xray_structure().extract_u_iso_or_u_equiv()*adptbx.u_as_b(1.)
    # selections for group ADP
    ph_box = model.get_hierarchy()
    ph_box.atoms().reset_i_seq()
    group_adp_sel = []
    if(self.group_mode=="one_per_residue"):
      for rg in ph_box.residue_groups():
        group_adp_sel.append(rg.atoms().extract_i_seq())
    elif(self.group_mode=="one_per_chain"):
      for chain in ph_box.chains():
        group_adp_sel.append(chain.atoms().extract_i_seq())
    elif(self.group_mode=="two_per_residue"):
      group_adp_sel = \
        mmtbx.refinement.utils.two_group_selections_per_residue(model)
    else:
      assert self.group_mode is None
    b_isos = fmodel.xray_structure.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1.)
    if self.group_mode is not None:
      #
      if(flex.max(b_isos)<1.e-2):
        b_isos = flex.random_double(model.size())*10
        model.set_b_iso(values = b_isos)
        fmodel.xray_structure.set_b_iso(values = b_isos)
        fmodel.update_xray_structure(xray_structure = fmodel.xray_structure,
          update_f_calc = True)
      #
      # One round if individual refinement follows. Otherwise up to five,
      # with the data renewed before each (update_f_obs); rounds after the
      # first seldom help (R-factor changes in the fourth digit, or goes
      # up) and each costs as much as the first, so they stop as soon as a
      # round gains less than group_min_gain in R-factor.
      number_of_rounds = 5
      if(self.individual): number_of_rounds = 1
      for i_round in range(number_of_rounds):
        if(i_round > 0): update_f_obs(fmodel = fmodel)
        r_before = fmodel.r_work()
        group_b_manager = mmtbx.refinement.group.manager(
          fmodel                   = fmodel,
          selections               = group_adp_sel,
          convergence_test         = False,
          max_number_of_iterations = 50,
          number_of_macro_cycles   = 1,
          run_finite_differences_test = False,
          use_restraints           = True,
          refine_adp               = True,
          log                      = log)
        if(r_before - fmodel.r_work() < group_min_gain): break
      fmodel.update_all_scales(update_f_part1=False, apply_back_trace=True,
        remove_outliers=False)
      b_isos = fmodel.xray_structure.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1.)
      model.set_b_iso(values = b_isos)
    #
    if(self.individual):
      riding = hydrogen_riding(model = model)
      if(log is not None):
        rmsb = format_value("%4.2f", model.rms_b_iso_or_b_equiv()).strip()
        print("r_work (start): %6.4f rms_B_bonded: %s"%(fmodel.r_work(),
          rmsb), file=log)
      rw = self.restraints_weight
      if(rw is not None and rw > 0): # None, 0 or negative: no restraints
        b_start = b_isos # start of every tightening pass, never the result
        flipped = False
        for it in range(1,20):
          # after a tightening pass restart from b_start (tg normalizes there)
          if(flipped):
            fmodel.xray_structure.set_b_iso(values = b_start)
            fmodel.update_xray_structure(update_f_calc = True)
          update_f_obs(fmodel = fmodel)
          x = fmodel.xray_structure.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1.)
          lower = flex.double(x.size(), 0)
          upper = flex.double(x.size(), flex.max(x)*2)
          calculator = tg(
            fmodel            = fmodel,
            x                 = x,
            restraints_weight = rw,
            bound_flags       = flex.int(x.size(), 2),
            lower_bound       = lower,
            upper_bound       = upper,
            riding            = riding)
          b_isos_prev = b_isos
          rms_b_prev = model.rms_b_iso_or_b_equiv()
          scitbx.minimizers.lbfgs(
            mode='lbfgsb', max_iterations=100, calculator=calculator,
            factr=lbfgsb_factr, pgtol=lbfgsb_pgtol)
          b_isos = fmodel.xray_structure.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1.)
          model.set_b_iso(values = b_isos)
          if(rms_b_prev is not None):
            rms_b = model.rms_b_iso_or_b_equiv()
            if(log is not None):
              print(
                "r_work: %6.4f rms_B_bonded: %4.2f restraints_weight: %6.4f"%(
                fmodel.r_work(), rms_b, rw),
                "evaluations: %d"%calculator.n_evaluations, file=log)
            if(rms_b<5):
              rw = rw/2
              if(flipped): # bracketed: keep the result closer to rms_B=5
                if(abs(rms_b_prev-5) < abs(rms_b-5)):
                  b_isos = b_isos_prev
                  model.set_b_iso(values = b_isos)
                break
              # Loosening the restraints did not change anything: they do not
              # matter any more, further trials would give the same result.
              # (Data are renewed before each trial, so results are never
              # identical: 0.1 A**2 on average is no change.)
              if(it > 1 and
                 flex.mean(flex.abs(b_isos - b_isos_prev)) < 0.1): break
            else:
              # it=1: b_isos_prev is the group result, do not fall back to it
              if(rms_b > rms_b_prev and it > 1):
                b_isos = b_isos_prev
                model.set_b_iso(values = b_isos)
                break
              rw = rw*2
              flipped = True
          if(rms_b_prev is None): break
      else:
        for it in range(1,2):
          update_f_obs(fmodel = fmodel)
          x = fmodel.xray_structure.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1.)
          ma = min(999,flex.max(x)*2)
          if ma<1: ma = 100
          lower = flex.double(x.size(), 0)
          upper = flex.double(x.size(), ma )
          calculator = tg(
            fmodel            = fmodel,
            x                 = x,
            restraints_weight = None,
            bound_flags       = flex.int(x.size(), 2),
            lower_bound       = lower,
            upper_bound       = upper,
            riding            = riding)
          m = scitbx.minimizers.lbfgs(
            mode='lbfgsb', max_iterations=100, calculator=calculator,
            factr=lbfgsb_factr, pgtol=lbfgsb_pgtol)
          b_isos = fmodel.xray_structure.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1.)
          model.set_b_iso(values = b_isos)
          if(log is not None):
            print("r_work: %6.4f"%fmodel.r_work(), file=log)
    #
    fmodel.xray_structure.set_b_iso(values = b_isos)
    fmodel.update_xray_structure(xray_structure = fmodel.xray_structure,
      update_f_calc = True)
    if(log is not None):
      print("r_work (final): %6.4f"%fmodel.r_work(), file=log)
    #
    return b_isos
