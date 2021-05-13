from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from scitbx import matrix
import math
from libtbx import adopt_init_args
import scitbx.lbfgs
from mmtbx.bulk_solvent import kbu_refinery
from cctbx import maptbx
import mmtbx.masks
import boost_adaptbx.boost.python as bp
asu_map_ext = bp.import_ext("cctbx_asymmetric_map_ext")
from libtbx import group_args
from mmtbx import bulk_solvent
from mmtbx.ncs import tncs
from collections import OrderedDict
import mmtbx.f_model
import sys
from libtbx.test_utils import approx_equal

from mmtbx import masks
from cctbx.masks import vdw_radii_from_xray_structure
ext = bp.import_ext("mmtbx_masks_ext")
mosaic_ext = bp.import_ext("mmtbx_mosaic_ext")

APPLY_SCALE_K1_TO_FOBS = False

def moving_average(x, n):
  r = []
  for i, xi in enumerate(x):
    s = 0
    cntr = 0
    for j in range(max(0,i-n), min(i+n+1, len(x))):
      s+=x[j]
      cntr+=1
    s = s/cntr
    r.append(s)
  return r

# Utilities used by algorithm 2 ------------------------------------------------

class minimizer(object):
  def __init__(self, max_iterations, calculator):
    adopt_init_args(self, locals())
    self.x = self.calculator.x
    self.cntr=0
    exception_handling_params = scitbx.lbfgs.exception_handling_parameters(
      ignore_line_search_failed_step_at_lower_bound=True,
      )
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      exception_handling_params=exception_handling_params,
      termination_params=scitbx.lbfgs.termination_parameters(
        max_iterations=max_iterations))

  def compute_functional_and_gradients(self):
    self.cntr+=1
    self.calculator.update_target_and_grads(x=self.x)
    t = self.calculator.target()
    g = self.calculator.gradients()
    #print "step: %4d"%self.cntr, "target:", t, "params:", \
    #  " ".join(["%10.6f"%i for i in self.x]), math.log(t)
    return t,g

class minimizer2(object):

  def __init__(self, calculator, min_iterations=0, max_iterations=2000):
    adopt_init_args(self, locals())
    self.x = self.calculator.x
    self.n = self.x.size()
    self.cntr=0

  def run(self, use_curvatures=0):
    self.minimizer = kbu_refinery.lbfgs_run(
      target_evaluator=self,
      min_iterations=self.min_iterations,
      max_iterations=self.max_iterations,
      use_curvatures=use_curvatures)
    self(requests_f_and_g=True, requests_diag=False)
    return self

  def __call__(self, requests_f_and_g, requests_diag):
    self.cntr+=1
    self.calculator.update_target_and_grads(x=self.x)
    if (not requests_f_and_g and not requests_diag):
      requests_f_and_g = True
      requests_diag = True
    if (requests_f_and_g):
      self.f = self.calculator.target()
      self.g = self.calculator.gradients()
      self.d = None
    if (requests_diag):
      self.d = self.calculator.curvatures()
      #assert self.d.all_ne(0)
      if(self.d.all_eq(0)): self.d=None
      else:
        self.d = 1 / self.d
    #print "step: %4d"%self.cntr, "target:", self.f, "params:", \
    #  " ".join(["%10.6f"%i for i in self.x]) #, math.log(self.f)
    return self.x, self.f, self.g, self.d

class tg(object):
  def __init__(self, x, i_obs, F, use_curvatures):
    self.x = x
    self.i_obs = i_obs
    self.F = F
    self.t = None
    self.g = None
    self.d = None
    # Needed to do sums from small to large to prefent loss
    s = flex.sort_permutation(self.i_obs.data())
    self.i_obs = self.i_obs.select(s)
    self.F = [f.select(s) for f in self.F]
    #
    self.sum_i_obs = flex.sum(self.i_obs.data()) # needed for Python version
    self.use_curvatures=use_curvatures
    self.tgo = mosaic_ext.alg2_tg(
      F     = [f.data() for f in self.F],
      i_obs = self.i_obs.data())
    self.update_target_and_grads(x=x)

  def update(self, x):
    self.update_target_and_grads(x = x)

  def update_target_and_grads(self, x):
    self.x = x
    self.tgo.update(self.x)
    self.t = self.tgo.target()
    self.g = self.tgo.gradient()
#
# Reference implementation in Python
#    s = 1 #180/math.pi
#    i_model = flex.double(self.i_obs.data().size(),0)
#    for n, kn in enumerate(self.x):
#      for m, km in enumerate(self.x):
#        tmp = self.F[n].data()*flex.conj(self.F[m].data())
#        i_model += kn*km*flex.real(tmp)
#        #pn = self.F[n].phases().data()*s
#        #pm = self.F[m].phases().data()*s
#        #Fn = flex.abs(self.F[n].data())
#        #Fm = flex.abs(self.F[m].data())
#        #i_model += kn*km*Fn*Fm*flex.cos(pn-pm)
#    diff = i_model - self.i_obs.data()
#    #print (flex.min(diff), flex.max(diff))
#    t = flex.sum(diff*diff)/4
#    #
#    g = flex.double()
#    for j in range(len(self.F)):
#      tmp = flex.double(self.i_obs.data().size(),0)
#      for m, km in enumerate(self.x):
#        tmp += km * flex.real( self.F[j].data()*flex.conj(self.F[m].data()) )
#        #pj = self.F[j].phases().data()*s
#        #pm = self.F[m].phases().data()*s
#        #Fj = flex.abs(self.F[j].data())
#        #Fm = flex.abs(self.F[m].data())
#        #tmp += km * Fj*Fm*flex.cos(pj-pm)
#      g.append(flex.sum(diff*tmp))
#    self.t = t/self.sum_i_obs
#    self.g = g/self.sum_i_obs
#    #print (self.t,t1)
#    #print (list(self.g))
#    #print (list(g1))
#    #print ()
#    #assert approx_equal(self.t, t1, 5)
#    #assert approx_equal(self.g, g1, 1.e-6)
#
    if self.use_curvatures:
      d = flex.double()
      for j in range(len(self.F)):
        tmp1 = flex.double(self.i_obs.data().size(),0)
        tmp2 = flex.double(self.i_obs.data().size(),0)
        for m, km in enumerate(self.x):
          zz = flex.real( self.F[j].data()*flex.conj(self.F[m].data()) )
          tmp1 += km * zz
          tmp2 += zz
          #pj = self.F[j].phases().data()*s
          #pm = self.F[m].phases().data()*s
          #Fj = flex.abs(self.F[j].data())
          #Fm = flex.abs(self.F[m].data())
          #tmp += km * Fj*Fm*flex.cos(pj-pm)
        d.append(flex.sum(tmp1*tmp1 + tmp2))
      self.d=d

  def target(self): return self.t

  def gradients(self): return self.g

  def gradient(self): return self.gradients()

  def curvatures(self): return self.d/self.sum_i_obs
#-------------------------------------------------------------------------------

def write_map_file(crystal_symmetry, map_data, file_name):
  from iotbx import mrcfile
  mrcfile.write_ccp4_map(
    file_name   = file_name,
    unit_cell   = crystal_symmetry.unit_cell(),
    space_group = crystal_symmetry.space_group(),
    map_data    = map_data,
    labels      = flex.std_string([""]))

class refinery(object):
  def __init__(self, fmodel, fv, alg, anomaly=True, log = sys.stdout):
    assert alg in ["alg0", "alg2", "alg4", None]
    self.log            = log
    self.f_obs          = fmodel.f_obs()
    self.r_free_flags   = fmodel.r_free_flags()
    d_spacings          = self.f_obs.d_spacings().data()
    dsel                = d_spacings > 3
    k_mask_overall      = fmodel.k_masks()[0]
    self.bin_selections = fmodel.bin_selections
    #
    k_total = fmodel.k_total()
    self.f_calc         = fmodel.f_model()
    self.F              = [self.f_calc.deep_copy()] + fv.keys()
    #
    for it in range(3):
      #
      if alg is not None: ALG = alg
      else:
        if it ==0: ALG = "alg4"
        else:      ALG = "alg2"
      #
      if it>0:
        self.F = [self.fmodel.f_model().deep_copy()] + self.F[1:]
      self._print("cycle: %2d"%it)
      self._print("  volumes: "+" ".join([str(fv[f]) for f in self.F[1:]]))
      f_obs   = self.f_obs.deep_copy()
      if it==0: k_total = fmodel.k_total()
      else:     k_total = self.fmodel.k_total()
      i_obs   = f_obs.customized_copy(data = f_obs.data()*f_obs.data())
      K_MASKS = OrderedDict()

      self.bin_selections = self.f_obs.log_binning(
        n_reflections_in_lowest_resolution_bin = 100*len(self.F))

      for i_bin, sel in enumerate(self.bin_selections):
        d_max, d_min = f_obs.select(sel).d_max_min()
        if d_max<3: continue
        bin = "  bin %2d: %5.2f-%-5.2f: "%(i_bin, d_max, d_min)
        F = [f.select(sel) for f in self.F]
        k_total_sel = k_total.select(sel)
        F_scaled = [F[0].deep_copy()]+[f.customized_copy(data=f.data()*k_total_sel) for f in F[1:]]

        #r00=bulk_solvent.r_factor(f_obs.select(sel).data()*k_total_sel, F[0].data()*k_total_sel)

        # algorithm_0
        if(ALG=="alg0"):
          k_masks = algorithm_0(
            f_obs = f_obs.select(sel),
            F     = F_scaled,
            kt=k_total_sel)

        #fd = flex.complex_double(F[0].data().size())
        #for i,f in enumerate(F):
        #  fd = fd + f.data()*k_masks[i]
        #r0=bulk_solvent.r_factor(f_obs.select(sel).data()*k_total_sel, fd*k_total_sel)

        # algorithm_4
        if(ALG=="alg4"):
          if it==0: phase_source = fmodel.f_model().select(sel)
          else:     phase_source = self.fmodel.f_model().select(sel)
          k_masks = algorithm_4(
            f_obs             = self.f_obs.select(sel),
            F                 = F_scaled,
            auto_converge_eps = 0.0001,
            phase_source = phase_source)

        #fd = flex.complex_double(F[0].data().size())
        #for i,f in enumerate(F):
        #  fd = fd + f.data()*k_masks[i]
        #r4=bulk_solvent.r_factor(f_obs.select(sel).data()*k_total_sel, fd*k_total_sel)

        # algorithm_2
        if(ALG=="alg2"):
          k_masks = algorithm_2(
            i_obs          = i_obs.select(sel),
            F              = F_scaled,
            x              = self._get_x_init(i_bin),
            use_curvatures = False)

        #fd = flex.complex_double(F[0].data().size())
        #for i,f in enumerate(F):
        #  fd = fd + f.data()*k_masks[i]
        #r2=bulk_solvent.r_factor(f_obs.select(sel).data()*k_total_sel, fd*k_total_sel)

        #self._print(bin+" ".join(["%6.2f"%k for k in k_masks])+" %6.4f %6.4f %6.4f %6.4f"%(r00,r0,r4, r2))
        k_mean = flex.mean(k_mask_overall.select(sel))
        k_masks_plus = [k_masks[0]]+[k_mean + k for k in k_masks[1:]]
        self._print(bin+" ".join(["%6.2f"%k for k in k_masks_plus]) )
        K_MASKS[sel] = [k_masks, k_masks_plus]
      #
      if(len(self.F)==2): break # stop and fall back onto using largest mask
      #
      #
      #print()
      #self.update_k_masks(K_MASKS)
      #for k_masks in K_MASKS.values():
      #  self._print(bin+" ".join(["%6.2f"%k for k in k_masks]))
      #
      f_calc_data = self.f_calc.data().deep_copy()
      f_bulk_data = flex.complex_double(fmodel.f_calc().data().size(), 0)
      for sel, k_masks in zip(K_MASKS.keys(), K_MASKS.values()):
        k_masks = k_masks[0] # 1 is shifted!
        f_bulk_data_ = flex.complex_double(sel.count(True), 0)
        for i_mask, k_mask in enumerate(k_masks):
          if i_mask==0:
            f_calc_data = f_calc_data.set_selected(sel,
              f_calc_data.select(sel)*k_mask)
            continue
          f_bulk_data_ += self.F[i_mask].data().select(sel)*k_mask
        f_bulk_data = f_bulk_data.set_selected(sel,f_bulk_data_)
      #
      self.update_F(K_MASKS)
      f_bulk = fmodel.f_calc().customized_copy(data = f_bulk_data)

      if(len(self.F)==2):
        self.fmodel = mmtbx.f_model.manager(
          f_obs          = self.f_obs,
          r_free_flags   = self.r_free_flags,
          f_calc         = fmodel.f_calc(),
          f_mask         = self.F[1],
          k_mask         = flex.double(f_obs.data().size(),1)
          )
        self.fmodel.update_all_scales(remove_outliers=False,
          apply_scale_k1_to_f_obs = APPLY_SCALE_K1_TO_FOBS)
      else:
        self.fmodel = mmtbx.f_model.manager(
          f_obs          = self.f_obs,
          r_free_flags   = self.r_free_flags,
          f_calc         = self.f_obs.customized_copy(data = f_calc_data),
          bin_selections = self.bin_selections,
          f_mask         = f_bulk,
          k_mask         = flex.double(f_obs.data().size(),1)
          )
        self.fmodel.update_all_scales(remove_outliers=False,
          apply_scale_k1_to_f_obs = APPLY_SCALE_K1_TO_FOBS)
        #
        self.fmodel = mmtbx.f_model.manager(
          f_obs          = self.f_obs,
          r_free_flags   = self.r_free_flags,
          #f_calc         = self.f_obs.customized_copy(data = f_calc_data),
          f_calc         = self.fmodel.f_calc(),
          f_mask         = self.fmodel.f_bulk(),
          k_mask         = flex.double(f_obs.data().size(),1)
          )
        self.fmodel.update_all_scales(remove_outliers=False,
          apply_scale_k1_to_f_obs = APPLY_SCALE_K1_TO_FOBS)
        self._print(self.fmodel.r_factors(prefix="  "))

      #self._print(self.fmodel.r_factors(prefix="  "))
      self.mc = self.fmodel.electron_density_map().map_coefficients(
        map_type   = "mFobs-DFmodel",
        isotropize = True,
        exclude_free_r_reflections = False)

  #def update_k_masks(self, K_MASKS):
  #  tmp = []
  #  for i_mask, F in enumerate(self.F):
  #    k_masks = [k_masks_bin[i_mask] for k_masks_bin in K_MASKS.values()]
  #    found = False
  #    for i_bin, k_masks_bin in enumerate(K_MASKS.values()):
  #      if(not found and k_masks_bin[i_mask]<=0.009):
  #        found = True
  #        K_MASKS.values()[i_bin][i_mask]=0
  #      elif found:
  #        K_MASKS.values()[i_bin][i_mask]=0

  def _print(self, m):
    if(self.log is not None):
      print(m, file=self.log)

  def update_F(self, K_MASKS):
    tmp = []
    for i_mask, F in enumerate(self.F):
      k_masks = [k_masks_bin[1][i_mask] for k_masks_bin in K_MASKS.values()]
      if(i_mask == 0):                        tmp.append(self.F[0])
      elif moving_average(k_masks,2)[0]>=0.03: tmp.append(F)
      self.F = tmp[:]

  def _get_x_init(self, i_bin):
    return flex.double([1] + [1]*len(self.F[1:]))
    #k_maks1_init = 0.35 - i_bin*0.35/len(self.bin_selections)
    #x = flex.double([1,k_maks1_init])
    #x.extend( flex.double(len(self.F)-2, 0.1))
    #return x

def get_f_mask(xrs, ma, step, option = 2, r_shrink = None, r_sol = None):
  crystal_gridding = maptbx.crystal_gridding(
    unit_cell        = xrs.unit_cell(),
    space_group_info = xrs.space_group_info(),
    symmetry_flags   = maptbx.use_space_group_symmetry,
    step             = step)
  n_real = crystal_gridding.n_real()
  atom_radii = vdw_radii_from_xray_structure(xray_structure = xrs)
  mask_params = masks.mask_master_params.extract()
  grid_step_factor = ma.d_min()/step
  if(r_shrink is not None): mask_params.shrink_truncation_radius = r_shrink
  if(r_sol is not None):    mask_params.solvent_radius           = r_sol
  mask_params.grid_step_factor = grid_step_factor
  # 1
  if(option==1):
    asu_mask = ext.atom_mask(
      unit_cell                = xrs.unit_cell(),
      group                    = xrs.space_group(),
      resolution               = ma.d_min(),
      grid_step_factor         = grid_step_factor,
      solvent_radius           = mask_params.solvent_radius,
      shrink_truncation_radius = mask_params.shrink_truncation_radius)
    asu_mask.compute(xrs.sites_frac(), atom_radii)
    fm_asu = asu_mask.structure_factors(ma.indices())
    f_mask = ma.set().array(data = fm_asu)
  # 2
  elif(option==2):
    asu_mask = ext.atom_mask(
      unit_cell                = xrs.unit_cell(),
      space_group              = xrs.space_group(),
      gridding_n_real          = n_real,
      solvent_radius           = mask_params.solvent_radius,
      shrink_truncation_radius = mask_params.shrink_truncation_radius)
    asu_mask.compute(xrs.sites_frac(), atom_radii)
    fm_asu = asu_mask.structure_factors(ma.indices())
    f_mask = ma.set().array(data = fm_asu)
  # 3
  elif(option==3):
    mask_p1 = mmtbx.masks.mask_from_xray_structure(
      xray_structure           = xrs,
      p1                       = True,
      for_structure_factors    = True,
      solvent_radius           = mask_params.solvent_radius,
      shrink_truncation_radius = mask_params.shrink_truncation_radius,
      n_real                   = n_real,
      in_asu                   = False).mask_data
    maptbx.unpad_in_place(map=mask_p1)
    mask = asu_map_ext.asymmetric_map(
      xrs.crystal_symmetry().space_group().type(), mask_p1).data()
    f_mask = ma.structure_factors_from_asu_map(
      asu_map_data = mask, n_real = n_real)
  # 4
  elif(option==4):
    f_mask = masks.bulk_solvent(
      xray_structure              = xrs,
      ignore_zero_occupancy_atoms = False,
      solvent_radius              = mask_params.solvent_radius,
      shrink_truncation_radius    = mask_params.shrink_truncation_radius,
      ignore_hydrogen_atoms       = False,
      grid_step                   = step,
      atom_radii                  = atom_radii).structure_factors(
        miller_set = ma)
  elif(option==5):
    o = mmtbx.masks.bulk_solvent(
      xray_structure              = xrs,
      ignore_zero_occupancy_atoms = False,
      solvent_radius              = mask_params.solvent_radius,
      shrink_truncation_radius    = mask_params.shrink_truncation_radius,
      ignore_hydrogen_atoms       = False,
      gridding_n_real             = n_real,
      atom_radii                  = atom_radii)
    assert approx_equal(n_real, o.data.accessor().all())
    f_mask = o.structure_factors(ma)
  elif(option==6):
    # XXX No control over n_real, so results with others don't match
    mask_manager = masks.manager(
      miller_array      = ma,
      miller_array_twin = None,
      mask_params       = mask_params)
    f_mask = mask_manager.shell_f_masks(xray_structure=xrs, force_update=True)[0]
  else: assert 0
  #
  return f_mask

def filter_mask(mask_p1, volume_cutoff, crystal_symmetry,
                for_structure_factors = False):
  co = maptbx.connectivity(
    map_data                   = mask_p1,
    threshold                  = 0.01,
    preprocess_against_shallow = True,
    wrapping                   = True)
  mi, ma = flex.min(mask_p1), flex.max(mask_p1)
  print (mask_p1.size(), (mask_p1<0).count(True))
  assert mi == 0, mi
  assert ma == 1, ma
  a,b,c    = crystal_symmetry.unit_cell().parameters()[:3]
  na,nb,nc = mask_p1.accessor().all()
  step = flex.mean(flex.double([a/na, b/nb, c/nc]))
  if(crystal_symmetry.space_group_number() != 1):
    co.merge_symmetry_related_regions(space_group=crystal_symmetry.space_group())
  conn = co.result().as_double()
  z = zip(co.regions(),range(0,co.regions().size()))
  sorted_by_volume = sorted(z, key=lambda x: x[0], reverse=True)
  for i_seq, p in enumerate(sorted_by_volume):
    v, i = p
    if(i==0): continue # skip macromolecule
    # skip small volume
    volume = v*step**3
    if volume < volume_cutoff:
      conn = conn.set_selected(conn==i, 0)
  conn = conn.set_selected(conn>0, 1)
  if for_structure_factors:
    conn = conn / crystal_symmetry.space_group().order_z()
  return conn

class mosaic_f_mask(object):
  def __init__(self,
               xray_structure,
               step,
               volume_cutoff=None,
               mean_diff_map_threshold=None,
               compute_whole=False,
               preprocess_against_shallow=True,
               largest_only=False,
               wrapping=True,
               f_obs=None,
               r_sol=1.1,
               r_shrink=0.9,
               f_calc=None,
               log = None,
               write_masks=False):
    adopt_init_args(self, locals())
    #
    self.dsel = f_obs.d_spacings().data()>=0 # XXX WHY????????????
    self.miller_array = f_obs.select(self.dsel)
    #
    self.crystal_symmetry = self.xray_structure.crystal_symmetry()
    # compute mask in p1 (via ASU)
    self.crystal_gridding = maptbx.crystal_gridding(
      unit_cell        = xray_structure.unit_cell(),
      space_group_info = xray_structure.space_group_info(),
      symmetry_flags   = maptbx.use_space_group_symmetry,
      step             = step)
    self.n_real = self.crystal_gridding.n_real()
    # XXX Where do we want to deal with H and occ==0?
    mask_p1 = mmtbx.masks.mask_from_xray_structure(
      xray_structure           = xray_structure,
      p1                       = True,
      for_structure_factors    = True,
      solvent_radius           = r_sol,
      shrink_truncation_radius = r_shrink,
      n_real                   = self.n_real,
      in_asu                   = False).mask_data
    maptbx.unpad_in_place(map=mask_p1)
    self.f_mask_whole = None
    if(compute_whole):
      mask = asu_map_ext.asymmetric_map(
        xray_structure.crystal_symmetry().space_group().type(), mask_p1).data()
      self.f_mask_whole = self.miller_array.structure_factors_from_asu_map(
        asu_map_data = mask, n_real = self.n_real)
    self.solvent_content = 100.*mask_p1.count(1)/mask_p1.size()
    if(write_masks):
      write_map_file(crystal_symmetry=xray_structure.crystal_symmetry(),
        map_data=mask_p1, file_name="mask_whole.mrc")
    # conn analysis
    co = maptbx.connectivity(
      map_data                   = mask_p1,
      threshold                  = 0.01,
      preprocess_against_shallow = preprocess_against_shallow,
      wrapping                   = wrapping)
    co.merge_symmetry_related_regions(space_group=xray_structure.space_group())
    del mask_p1
    self.conn = co.result().as_double()
    z = zip(co.regions(),range(0,co.regions().size()))
    sorted_by_volume = sorted(z, key=lambda x: x[0], reverse=True)
    #
    f_mask_data_0 = flex.complex_double(f_obs.data().size(), 0)
    f_mask_data   = flex.complex_double(f_obs.data().size(), 0)
    self.FV = OrderedDict()
    self.mc = None
    diff_map = None
    mean_diff_map = None
    self.regions = OrderedDict()
    self.f_mask_0 = None
    self.f_mask = None
    #
    if(log is not None):
      print("  #    volume_p1    uc(%) mFo-DFc: min,max,mean,sd", file=log)
    #
    for i_seq, p in enumerate(sorted_by_volume):
      v, i = p
      # skip macromolecule
      if(i==0): continue
      # skip small volume
      volume = v*step**3
      uc_fraction = v*100./self.conn.size()
      if(volume_cutoff is not None):
        if volume < volume_cutoff: continue

      selection = self.conn==i
      mask_i_asu = self.compute_i_mask_asu(selection = selection, volume = volume)
      volume_asu = (mask_i_asu>0).count(True)*step**3

      if(uc_fraction >= 1):
        f_mask_i = self.compute_f_mask_i(mask_i_asu)
        f_mask_data_0 += f_mask_i.data()
      elif(largest_only): break

      if(uc_fraction < 1 and diff_map is None):
        diff_map = self.compute_diff_map(f_mask_data = f_mask_data_0)

      mi,ma,me,sd = None,None,None,None
      if(diff_map is not None):
        blob = diff_map.select(selection.iselection())
        mean_diff_map = flex.mean(diff_map.select(selection.iselection()))
        mi,ma,me = flex.min(blob), flex.max(blob), flex.mean(blob)
        sd = blob.sample_standard_deviation()

      if(log is not None):
        print("%3d"%i_seq,"%12.3f"%volume, "%8.4f"%round(uc_fraction,4),
            "%7s"%str(None) if diff_map is None else "%7.3f %7.3f %7.3f %7.3f"%(
              mi,ma,me,sd), file=log)

      if(mean_diff_map_threshold is not None and
         mean_diff_map is not None and mean_diff_map<=mean_diff_map_threshold):
        continue

      self.regions[i_seq] = group_args(
        id          = i,
        i_seq       = i_seq,
        volume      = volume,
        uc_fraction = uc_fraction,
        diff_map    = group_args(mi=mi, ma=ma, me=me, sd=sd))

      f_mask_i = self.compute_f_mask_i(mask_i_asu)
      f_mask_data += f_mask_i.data()

      self.FV[f_mask_i] = [round(volume, 3), round(uc_fraction,1)]
    #
    self.f_mask_0 = f_obs.customized_copy(data = f_mask_data_0)
    self.f_mask   = f_obs.customized_copy(data = f_mask_data)
    self.do_mosaic = False
    self.n_regions = len(self.FV.keys())
    if(self.n_regions>1):
      self.do_mosaic = True

  def compute_f_mask_i(self, mask_i_asu):
    f_mask_i = self.miller_array.structure_factors_from_asu_map(
      asu_map_data = mask_i_asu, n_real = self.n_real)
    data = flex.complex_double(self.dsel.size(), 0)
    data = data.set_selected(self.dsel, f_mask_i.data())
    return self.f_obs.set().array(data = data)

  def compute_diff_map(self, f_mask_data):
    if(self.f_calc is None): return None
    f_mask = self.f_obs.customized_copy(data = f_mask_data)
    fmodel = mmtbx.f_model.manager(
      f_obs  = self.f_obs,
      f_calc = self.f_calc,
      f_mask = f_mask)
    fmodel = fmodel.select(self.dsel)
    fmodel.update_all_scales(remove_outliers=True,
      apply_scale_k1_to_f_obs = APPLY_SCALE_K1_TO_FOBS)
    self.mc = fmodel.electron_density_map().map_coefficients(
      map_type   = "mFobs-DFmodel",
      isotropize = True,
      exclude_free_r_reflections = False)
    fft_map = self.mc.fft_map(crystal_gridding = self.crystal_gridding)
    fft_map.apply_sigma_scaling()
    return fft_map.real_map_unpadded()

  def compute_i_mask_asu(self, selection, volume):
    mask_i = flex.double(flex.grid(self.n_real), 0)
    mask_i = mask_i.set_selected(selection, 1)
    if(self.write_masks):
      write_map_file(
        crystal_symmetry = self.crystal_symmetry,
        map_data         = mask_i,
        file_name        = "mask_%s.mrc"%str(round(volume,3)))
    tmp = asu_map_ext.asymmetric_map(
      self.crystal_symmetry.space_group().type(), mask_i).data()
    return tmp

def algorithm_0(f_obs, F, kt):
  """
  Grid search
  """
  fc, f_masks = F[0], F[1:]
  k_mask_trial_range=[]
  s = -1
  while s<1:
    k_mask_trial_range.append(s)
    s+=0.0001
  r = []
  fc_data = fc.data()
  for i, f_mask in enumerate(f_masks):
    #print("mask ",i)
    assert f_obs.data().size() == fc.data().size()
    assert f_mask.data().size() == fc.data().size()
    #print (bulk_solvent.r_factor(f_obs.data(),fc_data))
    kmask_, k_ = \
      bulk_solvent.k_mask_and_k_overall_grid_search(
        f_obs.data()*kt,
        fc_data*kt,
        f_mask.data()*kt,
        flex.double(k_mask_trial_range),
        flex.bool(fc.data().size(),True))
    r.append(kmask_)
    fc_data += fc_data*k_ + kmask_*f_mask.data()
    #print (bulk_solvent.r_factor(f_obs.data(),fc_data + kmask_*f_mask.data(),k_))
  r = [1,]+r
  return r

def algorithm_2(i_obs, F, x, use_curvatures=True, macro_cycles=10):
  """
  Unphased one-step search
  """
  calculator = tg(i_obs = i_obs, F=F, x = x, use_curvatures=use_curvatures)
  for it in range(macro_cycles):
    if(use_curvatures):
      m = minimizer(max_iterations=100, calculator=calculator)
    else:
      upper = flex.double([1.1] + [1]*(x.size()-1))
      lower = flex.double([0.9] + [-1]*(x.size()-1))
      #upper = flex.double([10] + [5]*(x.size()-1))
      #lower = flex.double([0.1] + [-5]*(x.size()-1))
      #upper = flex.double([10] + [0.65]*(x.size()-1))
      #lower = flex.double([0.1] + [0]*(x.size()-1))

      #upper = flex.double([1] + [0.65]*(x.size()-1))
      #lower = flex.double([1] + [0]*(x.size()-1))
      #upper = flex.double([1] + [5.65]*(x.size()-1))
      #lower = flex.double([1] + [-5]*(x.size()-1))
      m = tncs.minimizer(
        potential       = calculator,
        use_bounds      = 2,
        lower_bound     = lower,
        upper_bound     = upper,
        initial_values  = x).run()
    calculator = tg(i_obs = i_obs, F=F, x = m.x, use_curvatures=use_curvatures)
  if(use_curvatures):
    for it in range(10):
      m = minimizer(max_iterations=100, calculator=calculator)
      calculator = tg(i_obs = i_obs, F=F, x = m.x, use_curvatures=use_curvatures)
      m = minimizer2(max_iterations=100, calculator=calculator).run(use_curvatures=True)
      calculator = tg(i_obs = i_obs, F=F, x = m.x, use_curvatures=use_curvatures)
  return m.x

def algorithm_3(i_obs, fc, f_masks):
  """
  Unphased two-step search
  """
  F = [fc]+f_masks
  Gnm = []
  cs = {}
  cntr=0
  nm=[]
  # Compute and store Gnm
  for n, Fn in enumerate(F):
    for m, Fm in enumerate(F):
      if m < n:
        continue
      Gnm.append( flex.real( Fn.data()*flex.conj(Fm.data()) ) )
      cs[(n,m)] = cntr
      cntr+=1
      nm.append((n,m))
  # Keep track of indices for "upper triangular matrix vs full"
  for k,v in zip(list(cs.keys()), list(cs.values())):
    i,j=k
    if i==j: continue
    else: cs[(j,i)]=v
  # Generate and solve system Ax=b, x = A_1*b
  A = []
  b = []
  for u, Gnm_u in enumerate(Gnm):
    for v, Gnm_v in enumerate(Gnm):
      scale = 2
      n,m=nm[v]
      if n==m: scale=1
      A.append( flex.sum(Gnm_u*Gnm_v)*scale )
    b.append( flex.sum(Gnm_u * i_obs.data()) )
  A = matrix.sqr(A)
  A_1 = A.inverse()
  b = matrix.col(b)
  x = A_1 * b
  # Expand Xmn from solution x
  Xmn = []
  for n, Fn in enumerate(F):
    rows = []
    for m, Fm in enumerate(F):
      x_ = x[cs[(n,m)]]
      rows.append(x_)
    Xmn.append(rows)
  # Do formula (19)
  lnK = []
  for j, Fj in enumerate(F):
    t1 = flex.sum( flex.log( flex.double(Xmn[j]) ) )
    t2 = 0
    for n, Fn in enumerate(F):
      for m, Fm in enumerate(F):
        t2 += math.log(Xmn[n][m])
    t2 = t2 / (2*len(F))
    lnK.append( 1/len(F)*(t1-t2) )
  return [math.exp(x) for x in lnK]

def algorithm_4(f_obs, F, phase_source, max_cycles=100, auto_converge_eps=1.e-7,
                use_cpp=True):
  """
  Phased simultaneous search (alg4)
  """
  fc, f_masks = F[0], F[1:]
  fc = fc.deep_copy()
  F = [fc]+F[1:]
  # C++ version
  if(use_cpp):
    return mosaic_ext.alg4(
      [f.data() for f in F],
      f_obs.data(),
      phase_source.data(),
      max_cycles,
      auto_converge_eps)
  # Python version (1.2-3 times slower, but much more readable!)
  cntr = 0
  x_prev = None
  while True:
    f_obs_cmpl = f_obs.phase_transfer(phase_source = phase_source)
    A = []
    b = []
    for j, Fj in enumerate(F):
      A_rows = []
      for n, Fn in enumerate(F):
        Gjn = flex.real( Fj.data()*flex.conj(Fn.data()) )
        A_rows.append( flex.sum(Gjn) )
      Hj = flex.real( Fj.data()*flex.conj(f_obs_cmpl.data()) )
      b.append(flex.sum(Hj))
      A.extend(A_rows)
    A = matrix.sqr(A)
    A_1 = A.inverse()
    b = matrix.col(b)
    x = A_1 * b
    #
    fc_d = flex.complex_double(phase_source.indices().size(), 0)
    for i, f in enumerate(F):
      fc_d += f.data()*x[i]
    phase_source = phase_source.customized_copy(data = fc_d)
    x_ = x[:]
    #
    cntr+=1
    if(cntr>max_cycles): break
    if(x_prev is None): x_prev = x_[:]
    else:
      max_diff = flex.max(flex.abs(flex.double(x_prev)-flex.double(x_)))
      if(max_diff<=auto_converge_eps): break
      x_prev = x_[:]
  return x_
