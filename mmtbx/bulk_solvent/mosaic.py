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
from collections import OrderedDict
import mmtbx.f_model
import sys
from libtbx.test_utils import approx_equal
from libtbx.utils import Sorry
from cctbx import miller

from mmtbx import masks
from cctbx.masks import vdw_radii_from_xray_structure
ext = bp.import_ext("mmtbx_masks_ext")
mosaic_ext = bp.import_ext("mmtbx_mosaic_ext")

from scitbx import minimizers

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

class minimizer2(object):

  def __init__(self, calculator, min_iterations=0, max_iterations=2000):
    adopt_init_args(self, locals())
    self.x = self.calculator.x
    self.n = self.x.size()
    self.cntr=0

  def run(self, use_curvatures):
    self.minimizer = kbu_refinery.lbfgs_run(
      target_evaluator=self,
      min_iterations=self.min_iterations,
      max_iterations=self.max_iterations,
      use_curvatures=use_curvatures)
    self(requests_f_and_g=True, requests_diag=use_curvatures)
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
  def __init__(self, x, i_obs, F, use_curvatures,
               bound_flags=None, lower_bound=None, upper_bound=None):
    self.x = x
    self.i_obs = i_obs
    self.F = F
    self.t = None
    self.g = None
    self.d = None
    self.bound_flags=bound_flags
    self.lower_bound=lower_bound
    self.upper_bound=upper_bound
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

def thiken_bins(bins, n, ds):
  result = []
  group = []
  cntr = 0
  for bin in bins:
    bs = bin.count(True)
    if bs>=n and len(group)==0:
      result.append(bin)
      cntr=0
    else:
      group.append(bin)
      cntr += bs
    if cntr>=n:
      result.append(group)
      group = []
      cntr = 0
  #
  for i, r in enumerate(result):
    if(not isinstance(r, flex.bool)):
      r0 = r[0]
      for r_ in r:
        r0 = r0 | r_
      result[i] = r0
  #
  if(len(result)==0):
    result = [flex.bool(bins[0].size(), True)]
  #
  tmp=[]
  for i, bin in enumerate(result):
    ds_ = ds.select(bin)
    mi, ma = flex.min(ds_), flex.max(ds_)
    if(mi<3 and ma>=3 and len(tmp)>0):
      tmp[i-1] = tmp[i-1] | bin
    else:
      tmp.append(bin)
  result = tmp
  #
  return result

class refinery(object):
  def __init__(self, fmodel, fv, alg, log = sys.stdout):
    assert alg in ["alg0", "alg2", "alg4", "alg4a"]
    self.log             = log
    self.fv              = fv
    self.f_obs           = fmodel.f_obs()
    self.r_free_flags    = fmodel.r_free_flags()
    k_mask_overall       = fmodel.k_masks()[0]
    bin_selections_input = fmodel.bin_selections
    self.f_calc          = fmodel.f_calc()
    phase_source_init    = fmodel.f_model()
    k_total              = fmodel.k_total()
    self.F               = [self.f_calc.deep_copy()] + list(self.fv.keys())
    n_zones_start        = len(self.F)
    r4_start             = fmodel.r_work4()
    r4_best              = r4_start
    del fmodel
    self.fmodel          = None
    r4_start             = None
    r4_best              = None
    self.fmodel_best     = None
    self.n_regions       = None
    #

    #alg_save = alg

    for it in range(5):

      #if it in [0,1]:
      #  alg = "alg4a"
      #else:
      #  alg = alg_save

      #
      if(it==1):
        r4_start         = self.fmodel.r_work4()
        r4_best          = r4_start
        self.fmodel_best = self.fmodel.deep_copy()
      if(it>1):
        r4 = self.fmodel.r_work4()
        if(abs(round(r4-r4_start,4))<1.e-4):
          break
        r4_start = r4

      self._print("cycle: %2d (regions: %d)"%(it, len(self.F[1:])))
      self._print("  Region #:                  "+"".join(["%7d"%fv[f] for f in self.F[1:]]))

      f_obs   = self.f_obs.deep_copy()
      if it>0: k_total = self.fmodel.k_total()
      i_obs   = f_obs.customized_copy(data = f_obs.data()*f_obs.data())
      K_MASKS = OrderedDict()

      if(alg.endswith('a')):
        self.bin_selections = [self.f_calc.d_spacings().data()>4]
      else:
        self.bin_selections = thiken_bins(
          bins=bin_selections_input, n=50*len(self.F), ds=self.f_calc.d_spacings().data())
      for i_bin, sel in enumerate(self.bin_selections):
        d_max, d_min = f_obs.select(sel).d_max_min()
        #if d_max<3 or d_min<3: continue
        if d_max<3: continue
        bin = "  bin %2d: %5.2f-%-5.2f: "%(i_bin, d_max, d_min)
        F = [f.select(sel) for f in self.F]
        k_total_sel = k_total.select(sel)
        #
        F_scaled = [f.customized_copy(data=f.data()*k_total_sel) for f in F]
        # algorithm_0
        if(alg.startswith("alg0")):
          k_masks = algorithm_0(
            f_obs = f_obs.select(sel),
            F     = F_scaled,
            kt=k_total_sel)
        # algorithm_4
        if(alg.startswith("alg4")):
          if it==0: phase_source = phase_source_init.select(sel)
          else:     phase_source = self.fmodel.f_model().select(sel)
          k_masks = algorithm_4(
            f_obs             = self.f_obs.select(sel),
            F                 = F_scaled,
            auto_converge_eps = 0.0001,
            phase_source = phase_source)
        # algorithm_2
        if(alg.startswith("alg2")):
          k_masks = algorithm_2(
            i_obs          = i_obs.select(sel),
            F              = F_scaled,
            x              = self._get_x_init(i_bin),
            use_curvatures = False)

        self._print(bin+" ".join(["%6.2f"%k for k in k_masks]) )
        K_MASKS[sel] = [k_masks, k_masks]
      #
      if(len(self.F)==2): break # stop and fall back onto using largest mask
      #
      f_calc_data = self.f_calc.data().deep_copy()
      f_bulk_data = flex.complex_double(self.f_calc.data().size(), 0)
      k_mask_0 = None
      for sel, k_masks in zip(K_MASKS.keys(), K_MASKS.values()):
        k_masks = k_masks[0] # 1 is shifted!
        f_bulk_data_ = flex.complex_double(sel.count(True), 0)
        for i_mask, k_mask in enumerate(k_masks):
          if i_mask==0:
            f_calc_data = f_calc_data.set_selected(sel,
              f_calc_data.select(sel)*k_mask)
            k_mask_0 = k_mask
            continue
          if k_mask_0 is not None: k_mask = k_mask/k_mask_0
          f_bulk_data_ += self.F[i_mask].data().select(sel)*k_mask
        f_bulk_data = f_bulk_data.set_selected(sel,f_bulk_data_ )
      #

      #if it in [0,1]:
      #  self.update_F(K_MASKS, cutoff=0)
      #else:
      #  self.update_F(K_MASKS, cutoff=0.03)
      self.update_F(K_MASKS, cutoff=0.03)


      f_bulk = self.f_calc.customized_copy(data = f_bulk_data)

      if(len(self.F)==2):
        self.fmodel = mmtbx.f_model.manager(
          f_obs          = self.f_obs,
          r_free_flags   = self.r_free_flags,
          f_calc         = self.f_calc,
          f_mask         = self.F[1],
          k_mask         = flex.double(f_obs.data().size(),1)
          )
        self.fmodel.update_all_scales(remove_outliers=False,
          apply_scale_k1_to_f_obs = APPLY_SCALE_K1_TO_FOBS)
        self._print(self.fmodel.r_factors(prefix="  "))
      else:
        self.fmodel = mmtbx.f_model.manager(
          f_obs          = self.f_obs,
          r_free_flags   = self.r_free_flags,
          #f_calc         = self.f_obs.customized_copy(data = f_calc_data),
          f_calc         = self.f_calc,
          bin_selections = bin_selections_input,
          #bin_selections = self.bin_selections,
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
          bin_selections = bin_selections_input,
          f_calc         = self.fmodel.f_calc(),
          f_mask         = self.fmodel.f_bulk(),
          k_mask         = flex.double(f_obs.data().size(),1)
          )
        self.fmodel.update_all_scales(remove_outliers=False,
          apply_scale_k1_to_f_obs = APPLY_SCALE_K1_TO_FOBS)
        self._print(self.fmodel.r_factors(prefix="  "))
      #
      r4 = self.fmodel.r_work4()
      if(r4_best is not None and r4<=r4_best):
        r4_best = r4
        self.fmodel_best = self.fmodel.deep_copy()
    #
    self.fmodel = self.fmodel_best
    self.n_regions = len(self.F)-1


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

  def update_F(self, K_MASKS, cutoff=0.03):
    tmp = []
    for i_mask, F in enumerate(self.F):
      k_masks = [k_masks_bin[1][i_mask] for k_masks_bin in K_MASKS.values()]
      #print(i_mask, flex.mean(flex.double(k_masks)), [round(i,2) for i in k_masks])
      if(i_mask == 0):                        tmp.append(self.F[0])
      #elif moving_average(k_masks,2)[0]>=cutoff: tmp.append(F)
      elif flex.mean(flex.double(k_masks))>0.01: tmp.append(F)
      self.F = tmp[:]

  def _get_x_init(self, i_bin):
    return flex.double([1] + [1]*len(self.F[1:]))
    #k_maks1_init = 0.35 - i_bin*0.35/len(self.bin_selections)
    #x = flex.double([1,k_maks1_init])
    #x.extend( flex.double(len(self.F)-2, 0.1))
    #return x

def get_f_mask(xrs, ma, step, option = 2, r_shrink = None, r_sol = None):
  #
  result = ma.deep_copy()
  sel = ma.d_spacings().data()>=3
  ma = ma.select(sel)
  #
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
  data = flex.complex_double(result.indices().size(), 0)
  data = data.set_selected(sel, f_mask.data())
  result = result.array(data = data)
  return result

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


def modify(fofc, R, cutoff, crystal_gridding, mask=None):
  G = fofc.g_function(R=R)

  fft_map = fofc.fft_map(crystal_gridding = crystal_gridding)
  fft_map.apply_volume_scaling()
  map_data = fft_map.real_map_unpadded()

  if mask is not None: map_data = map_data*mask

  sel = map_data<cutoff
  map_data = map_data.set_selected(sel, cutoff)

  B = miller.structure_factor_box_from_map(
    crystal_symmetry=fofc.crystal_symmetry(),
    map=map_data, n_real=None,
     anomalous_flag=False, include_000=True)
  GB = B.g_function(R=R)
  B = B.customized_copy(data = B.data()*GB)
  fft_map = B.fft_map(crystal_gridding = crystal_gridding)
  fft_map.apply_volume_scaling()
  map_data = fft_map.real_map_unpadded()

  print(flex.min(map_data), flex.max(map_data), flex.mean(map_data))

  map_data = map_data.set_selected(map_data<abs(flex.min(map_data)), 0)

  print(flex.min(map_data), flex.max(map_data), flex.mean(map_data))

  F=B


  #F = fofc.structure_factors_from_map(map = map_data, in_place_fft=False,
  #  use_scale=True, use_sg=False)
  #
  #F = F.customized_copy(data = F.data()*G)
  #
  #fft_map = F.fft_map(crystal_gridding = crystal_gridding)
  #fft_map.apply_volume_scaling()
  #map_data = fft_map.real_map_unpadded()

  mtz_dataset = F.as_mtz_dataset(column_root_label = "F")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "map.mtz")

  write_map_file(crystal_symmetry=fofc.crystal_symmetry(),
    map_data=map_data, file_name="map.ccp4")

  return map_data

class mask_and_regions(object):
  def __init__(self,
      xray_structure,
      crystal_gridding,
      step=0.6,
      r_sol=1.1,
      r_shrink=0.9,
      volume_cutoff=50,
      wrapping=True,
      force_symmetry=True,
      log=None,
      modifier=None
      ):
    """
    Split 0/1 traditional bulk-solvent mask into a series of isolated masks
    covering the whole unit cell in P1.
    """
    self.crystal_gridding = crystal_gridding
    self.crystal_symmetry = xray_structure.crystal_symmetry()
    self.step = step
    # Compute mask in P1
    self.mask_p1 = mmtbx.masks.mask_from_xray_structure(
      xray_structure           = xray_structure,
      p1                       = True,
      for_structure_factors    = True,
      solvent_radius           = r_sol,
      shrink_truncation_radius = r_shrink,
      n_real                   = self.crystal_gridding.n_real(),
      in_asu                   = False).mask_data
    maptbx.unpad_in_place(map=self.mask_p1)
    #
    # TO-DO: instead of excluding, add them as mosaic regions, may be?
    #
    tmp  = flex.double(flex.grid(self.mask_p1.all()), 1)
    if(modifier is not None):
      modifier = modifier * self.mask_p1
      sel = modifier <  -2.0 # XXX Depending on this, results can vary a lot
      modifier = modifier.set_selected(sel, 1)
      modifier = modifier.set_selected(~sel, 0)
      co = maptbx.connectivity(
        map_data                   = modifier,
        threshold                  = 0.01,
        preprocess_against_shallow = False,
        wrapping                   = wrapping)
      if(force_symmetry and xray_structure.space_group().type().number() != 1):
        co.merge_symmetry_related_regions(
          space_group = xray_structure.space_group())
      conn = co.result().as_double()
      z = zip(co.regions(),range(0,co.regions().size()))
      sorted_by_volume = sorted(z, key=lambda x: x[0], reverse=True)
      for i_seq, p in enumerate(sorted_by_volume):
        v, i = p
        if(i==0): continue
        # Skip small volume regions
        volume = v*self.step**3
        tmp = tmp.set_selected(conn==i, 0)
        if(volume_cutoff is not None and volume < volume_cutoff):
          break
    modifier = tmp
    #
    #
    # Solvent fraction
    self.solvent_content=100.*(self.mask_p1!=0).count(True)/self.mask_p1.size()
    if(log is not None):
      print("Solvent content: %7.3f"%self.solvent_content, file=log)
    # Connectivity
    co = maptbx.connectivity(
      map_data                   = self.mask_p1 * modifier, #!!!!!!!!!!!!!!!!!!!
      threshold                  = 0.01,
      preprocess_against_shallow = False, # XXX WHY False?
      wrapping                   = wrapping)
    if(force_symmetry and xray_structure.space_group().type().number() != 1):
      co.merge_symmetry_related_regions(
        space_group = xray_structure.space_group())
    # Regions
    self.conn = co.result().as_double() # 0 = protein, 1 = solvent, >1 = other
    z = zip(co.regions(),range(0,co.regions().size()))
    self.sorted_by_volume = sorted(z, key=lambda x: x[0], reverse=True)
    self.n_regions_total = co.regions().size()
    print("Total number of regions: %d"%self.n_regions_total, file=log)
    #
    if(log is not None):
      print("Regions (V > %d):"%volume_cutoff, file=log)
      print("  #   #   Vol. (A^3)    uc(%)", file=log)
    self.small_selection = flex.size_t()
    self.regions = []
    cntr = 0
    for i_seq, p in enumerate(self.sorted_by_volume):
      v, i = p
      # skip macromolecule
      if(i==0): continue
      # Region (i)selection
      iselection = (self.conn==i).iselection()
      # Skip small volume regions
      volume = v*self.step**3
      uc_fraction = v*100./self.conn.size()
      if(volume_cutoff is not None and volume < volume_cutoff):
        self.small_selection.extend(iselection)
        continue
      # Region
      region = group_args(
        i_seq       = cntr,
        i           = i,
        volume      = volume,
        uc_fraction = uc_fraction,
        iselection  = iselection)
      cntr+=1
      self.regions.append(region)
      if(log is not None):
        print("%3d %3d %12.3f %8.4f"%(region.i_seq, region.i, region.volume,
          region.uc_fraction), file=log)
    #
    self.n_regions = len(self.regions)
    print("Number of regions (V > %d): %d"%(
      volume_cutoff, self.n_regions), file=log)

  def get_selection(self, region):
    result = flex.bool(self.conn.size(), region.iselection)
    result.resize(self.conn.accessor())
    return result

  def get_region_asu_mask(self, region, map_data=None, write=False):
    # This is ASU mask, NOT P1 !
    selection = self.get_selection(region = region)
    mask_i = flex.double(flex.grid(self.crystal_gridding.n_real()), 0)
    mask_i = mask_i.set_selected(selection, 1)
    if(map_data is not None):
      mask_i = mask_i * map_data # may beed to divide by max value of map_data
    if(write): # This write P1 mask for region
      write_map_file(
        crystal_symmetry = self.crystal_symmetry,
        map_data         = mask_i,
        file_name        = "mask_%s.mrc"%str(round(region.volume,3)))
    return asu_map_ext.asymmetric_map(
      self.crystal_symmetry.space_group().type(), mask_i).data()

  def _mask_to_sf(self, mask_data, miller_array, d_min_cutoff=3):
    # Cutoff incoming Miller array
    result = miller_array.deep_copy()
    d_spacings    = miller_array.d_spacings().data()
    sel           = d_spacings >= 3
    miller_array_ = miller_array.select(sel)
    # Compute SF for truncated Miller array
    mask_asu = asu_map_ext.asymmetric_map(
      self.crystal_symmetry.space_group().type(), mask_data).data()
    sf = miller_array_.structure_factors_from_asu_map(
        asu_map_data = mask_asu, n_real = self.crystal_gridding.n_real())
    # Inflate Miller array to original hkl set before returning
    data = flex.complex_double(d_spacings.size(), 0)
    data = data.set_selected(sel, sf.data())
    return result.set().array(data=data)

  def f_mask_whole(self, miller_array):
    return self._mask_to_sf(mask_data=self.mask_p1, miller_array=miller_array)

  def f_mask_whole_filtered(self, miller_array):
    mask_p1_corrected = self.mask_p1.deep_copy()
    sel = flex.bool(self.conn.size(), self.small_selection)
    sel.resize(self.conn.accessor())
    mask_p1_corrected = mask_p1_corrected.set_selected(sel, 0)
    return self._mask_to_sf(
      mask_data=mask_p1_corrected, miller_array=miller_array)

class f_masks(object):
  def __init__(self,
               mask_and_regions,
               crystal_gridding,
               xray_structure,
               f_obs,
               f_calc=None,
               r_free_flags=None,
               mean_diff_map_threshold=0.5,
               log = None,
               write_masks=False):
    """
    Compute Fmask(whole), Fmask(main) and list of mosaic contributions.
    Fmask match input miller array set but zeroed out in 3A-better resolution.
    """
    adopt_init_args(self, locals())
    #
    self.d_spacings   = f_obs.d_spacings().data()
    self.sel_gte3     = self.d_spacings >= 3
    self.miller_array = f_obs.select(self.sel_gte3)
    #
    self.crystal_symmetry = self.xray_structure.crystal_symmetry()
    self.n_real = self.crystal_gridding.n_real()
    #
    f_mask_data_main = flex.complex_double(f_obs.data().size(), 0)
    self.f_mask_main = None
    self.FV          = OrderedDict()
    self.mFoDFc_main = None
    self.diff_map    = None # mFo-DFc map computed using F_mask_0 (main mask)
    #
    if(log is not None):
      print("Regions satisfying volume and map criteria:", file=log)
      print("                                   mFo-DFc (sigma)", file=log)
      print("  #   Vol. (A^3)    uc(%)     min     max    mean      sd", file=log)
    #
    for region in self.mask_and_regions.regions:
      f_mask_i = None # must be here inside the loop!
      # Compute i-th region mask
      mask_i_asu = self.mask_and_regions.get_region_asu_mask(region = region)
      # Compute F_mask_0 (F_mask for main mask)
      if(region.uc_fraction >= 1):
        f_mask_i = self.compute_f_mask_i(mask_i_asu)
        f_mask_data_main += f_mask_i.data()
      # Compute mFo-DFc map using main mask (once done computing main mask!)
      if(region.uc_fraction < 1 and self.diff_map is None and
         self.f_calc is not None and mean_diff_map_threshold is not None):
        self.diff_map = self.compute_diff_map(f_mask_data = f_mask_data_main)
      if(self.diff_map is not None):
        # Analyze mFo-DFc map in the i-th region
        map_info = self._get_map_info(iselection = region.iselection)
        # Skip regions with weak mean density
        if(map_info.mean is not None and map_info.mean < mean_diff_map_threshold):
          continue
        if(log is not None):
          print("%3d"%region.i_seq,"%12.3f"%region.volume,
                "%8.4f"%round(region.uc_fraction,4), map_info.string, file=log)
      else:
        print("%3d"%region.i_seq,"%12.3f"%region.volume,
              "%8.4f"%round(region.uc_fraction,4), file=log)
      # Compute F_mask for i-th region
      if(f_mask_i is None):
        f_mask_i = self.compute_f_mask_i(mask_i_asu)
      # Compose result object
      self.FV[f_mask_i] = region.i_seq
    # END OF LOOP OVER REGIONS



    """

    #self.FV          = OrderedDict()
    #self.FV[f_obs.customized_copy(data = f_mask_data_main)] = 0

    self.diff_map = self.compute_diff_map(f_mask_data = f_mask_data_main)

    mask_ = mmtbx.masks.mask_from_xray_structure(
      xray_structure           = self.xray_structure,
      p1                       = True,
      for_structure_factors    = True,
      #solvent_radius           = 0.5,
      #shrink_truncation_radius = 0.5,
      n_real                   = self.n_real,
      in_asu                   = False).mask_data
    maptbx.unpad_in_place(map=mask_)

    mc = self.mFoDFc_main
    fft_map = mc.fft_map(
      symmetry_flags   = maptbx.use_space_group_symmetry,
      crystal_gridding = self.crystal_gridding)
    fft_map.apply_volume_scaling()
    map_data = fft_map.real_map_unpadded() #* mask_

    map_data = map_data.set_selected(map_data<0,0)

    #sgt = self.f_obs.space_group().type()
    #asu_map = asu_map_ext.asymmetric_map(sgt, map_data)
    #map_data_asu = asu_map.data()
    #map_data_asu = map_data_asu.shift_origin()
    #
    #asu_map = asu_map_ext.asymmetric_map(sgt, map_data_asu,
    #  self.crystal_gridding.n_real())
    #f_mask_i = self.f_calc.customized_copy(
    #  indices = self.f_calc.indices(),
    #  data    = asu_map.structure_factors(self.f_calc.indices() ))

    def _sf_from_map(map_data):
      return self.f_calc.structure_factors_from_map(
        map            = map_data,
        use_scale      = True,
        anomalous_flag = False,
        use_sg         = False)

    #f_mask_i = _sf_from_map(map_data = map_data)
    #
    #
    #self.FV[f_mask_i] = -99

    OFFSET = region.i_seq+1
    print("<><><><><><><><><><>")
    dmd = modify(fofc=self.mFoDFc_main, R=2, cutoff=0.3, #0.1/10,
       crystal_gridding=self.crystal_gridding,
       mask = self.mask_and_regions.mask_p1)
    #dmd = map_data
    #STOP()

    mask_p1 = self.mask_and_regions.mask_p1.deep_copy()

    #mask_ = mmtbx.masks.mask_from_xray_structure(
    #  xray_structure           = self.xray_structure,
    #  p1                       = True,
    #  for_structure_factors    = True,
    #  solvent_radius           = 0.5,
    #  shrink_truncation_radius = 0.5,
    #  n_real                   = self.n_real,
    #  in_asu                   = False).mask_data
    #maptbx.unpad_in_place(map=mask_)
    #dmd = dmd * mask_

    co = maptbx.connectivity(
      map_data                   = dmd,
      threshold                  = 0.04, #0.001/10,
      preprocess_against_shallow = False,
      wrapping                   = True)
    if(xray_structure.space_group().type().number() != 1): # not P1
      co.merge_symmetry_related_regions(
        space_group = xray_structure.space_group())
    conn = co.result().as_double()
    z = zip(co.regions(),range(0,co.regions().size()))
    sorted_by_volume = sorted(z, key=lambda x: x[0], reverse=True)
    for i_seq, p in enumerate(sorted_by_volume):
      i_seq = OFFSET + i_seq
      v, i = p

      #if(i==0 or i==1): continue
      if i==0: continue

      step=0.6
      volume = v*step**3
      volume_cutoff = 2500
      if(volume_cutoff is not None and volume < volume_cutoff): continue
      uc_fraction = v*100./conn.size()
      print(i, "%3d"%i_seq,"%12.3f"%volume, "%8.4f"%round(uc_fraction,4),
            "%7s"%str(None), file=log)
      iselection = (conn==i).iselection()
      selection = flex.bool(conn.size(), iselection)
      selection.resize(conn.accessor())

      mask_p1 = mask_p1.set_selected(selection, 0)

      region = group_args(
        i           = i,
        i_seq       = i_seq,
        volume      = volume,
        uc_fraction = uc_fraction,
        iselection  = iselection)
      self.mask_and_regions.regions.append( region )

      dmd_i = dmd.deep_copy()
      mask_i_asu = dmd_i.set_selected(~selection, 0)
      #mask_i_asu = mask_i_asu.set_selected(mask_i_asu>0, 1)
      f_mask_i = _sf_from_map(map_data=mask_i_asu)

      #mask_i_asu = mask_i_asu/flex.max(mask_i_asu)
      #mask_i_asu = self.mask_and_regions.get_region_asu_mask(region = region, map_data=dmd_i)
      #mask_i_asu = mask_i_asu/flex.max(mask_i_asu)
      #
      #f_mask_i = self.compute_f_mask_i(mask_i_asu)
      self.FV[f_mask_i] = -1*i_seq
    print("<><><><><><><><><><>")


    fmodel = mmtbx.f_model.manager(
      f_obs        = self.f_obs,
      f_calc       = self.f_calc,
      r_free_flags = self.r_free_flags,
      f_mask       = _sf_from_map(map_data=mask_p1))
    fmodel.update_all_scales(remove_outliers=True,
      apply_scale_k1_to_f_obs = APPLY_SCALE_K1_TO_FOBS)
    mc = fmodel.electron_density_map().map_coefficients(
      map_type   = "mFobs-DFmodel",
      isotropize = False,
      exclude_free_r_reflections = True)
    mtz_dataset = mc.as_mtz_dataset(column_root_label = "F")
    mtz_object = mtz_dataset.mtz_object()
    mtz_object.write(file_name = "map2.mtz")
    """


    #
    # Determine number of secondary regions. Must happen here!
    # Preliminarily if need to do mosaic.

    self.n_regions = len(self.FV.values())
    print("Number of regions (mean(mFo-DFc) >= %s): %d"%(
      mean_diff_map_threshold, self.n_regions), file=log)
    if(self.n_regions==0):
      raise Sorry("Region selection criteria lead to no regions")
    self.do_mosaic = False
    if(self.n_regions>1 and flex.max(self.d_spacings)>6):
      self.do_mosaic = True
    # Finalize main Fmask
    self.f_mask_main = f_obs.customized_copy(data = f_mask_data_main)
    # Delete large objects from memory
    del self.diff_map

  def _get_sum_f_masks(self):
    f_masks = list(self.FV.keys())
    data = f_masks[0].data().deep_copy()
    for f in f_masks[1:]:
      data += f.data()
    return self.f_obs.set().array(data = data)

  def _get_map_info(self, iselection):
    if(self.diff_map is None):
      return group_args(min=None, max=None, mean=None, sd=None,
                        string="   None    None    None    None")
    blob = self.diff_map.select(iselection)
    mi,ma,me = flex.min(blob), flex.max(blob), flex.mean(blob)
    sd = blob.sample_standard_deviation()
    return group_args(min=mi, max=ma, mean=me, sd=sd,
                      string="%7.3f %7.3f %7.3f %7.3f"%(mi,ma,me,sd))

  def _inflate(self, f):
    data = flex.complex_double(self.d_spacings.size(), 0)
    data = data.set_selected(self.sel_gte3, f.data())
    return self.f_obs.set().array(data = data)

  def compute_f_mask_i(self, mask_i_asu):
    return self._inflate(self.miller_array.structure_factors_from_asu_map(
      asu_map_data = mask_i_asu, n_real = self.n_real))

  def compute_diff_map(self, f_mask_data):
    if(self.f_calc is None): return None
    f_mask = self.f_obs.customized_copy(data = f_mask_data)
    fmodel = mmtbx.f_model.manager(
      f_obs        = self.f_obs,
      f_calc       = self.f_calc,
      r_free_flags = self.r_free_flags,
      f_mask       = f_mask)
    fmodel.update_all_scales(remove_outliers=True,
      apply_scale_k1_to_f_obs = APPLY_SCALE_K1_TO_FOBS)
    self.mFoDFc_main = fmodel.electron_density_map().map_coefficients(
      map_type   = "mFobs-DFmodel",
      isotropize = False,
      exclude_free_r_reflections = True)
    fft_map = self.mFoDFc_main.fft_map(crystal_gridding = self.crystal_gridding)
    fft_map.apply_sigma_scaling()
    return fft_map.real_map_unpadded()

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

def algorithm_2(i_obs, F, x, use_curvatures, use_lbfgsb, macro_cycles=10,
                max_iterations=100):
  """
  Unphased one-step search
  """
  assert use_curvatures in [True, False, None]
  assert use_lbfgsb in [True, False, None]
  upper = flex.double([1.1] + [12]*(x.size()-1))
  lower = flex.double([0.9] + [0]*(x.size()-1))
  for it in range(macro_cycles):
    x = x.set_selected(x<0,1.e-6)
    if(use_curvatures is True): # Curvatures only
      calculator = tg(i_obs = i_obs, F=F, x = x, use_curvatures=True)
      m = minimizer2(
        max_iterations=max_iterations, calculator=calculator).run(use_curvatures=True)
      x = m.x
      x = x.set_selected(x<0,1.e-6)
    elif(use_curvatures is False): # No curvatures at all
      calculator = tg(i_obs = i_obs, F=F, x = x, use_curvatures=False)
      if(use_lbfgsb is True):
        m = scitbx.minimizers.lbfgs(
           mode='lbfgsb', max_iterations=max_iterations, calculator=calculator)
      else:
        m = scitbx.minimizers.lbfgs(
           mode='lbfgs', max_iterations=max_iterations, calculator=calculator)
      x = m.x
      x = x.set_selected(x<0,1.e-6)
    elif(use_curvatures is None):
      calculator = tg(i_obs = i_obs, F=F, x = x, use_curvatures=True)
      m = minimizer2(
        max_iterations=max_iterations, calculator=calculator).run(use_curvatures=True)
      x = m.x
      x = x.set_selected(x<0,1.e-6)
      calculator = tg(i_obs = i_obs, F=F, x = x, use_curvatures=False,
        bound_flags=flex.int(x.size(),2), lower_bound=lower, upper_bound=upper)
      if(use_lbfgsb is True):
        m = scitbx.minimizers.lbfgs(
           mode='lbfgsb', max_iterations=max_iterations, calculator=calculator)
      else:
        m = scitbx.minimizers.lbfgs(
           mode='lbfgs', max_iterations=max_iterations, calculator=calculator)
      x = m.x
      x = x.set_selected(x<0,1.e-6)
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
                use_cpp=True, x_init=None):
  """
  Phased simultaneous search (alg4)
  """
  fc, f_masks = F[0], F[1:]
  if(x_init is not None):
    assert x_init.size() == len(f_masks)
    # Apply scale
    tmp = []
    for fm, s in zip(f_masks, x_init):
      fm = fm.array(data = fm.data()*s)
      tmp.append(fm)
    f_masks = tmp
  #
  fc = fc.deep_copy()
  F = [fc]+f_masks
  # C++ version
  if(use_cpp):
    result = mosaic_ext.alg4(
      [f.data() for f in F],
      f_obs.data(),
      phase_source.data(),
      max_cycles,
      auto_converge_eps)
  else:
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
    result = x_
  if(x_init is not None):
    for i,tmp in enumerate(x_init):
      result[i+1] = result[i+1] * x_init[i]
  return result
