
"""
Methods and utilities for generating composite omit maps.  The actual
end-user application, phenix.composite_omit_map, is part of the phenix
sources.
"""

from __future__ import absolute_import, division, print_function
import mmtbx.f_model
from cctbx import miller
from cctbx import maptbx
import cctbx.miller
from scitbx.array_family import flex
import boost_adaptbx.boost.python as bp
from six.moves import zip
from six.moves import range
asu_map_ext = bp.import_ext("cctbx_asymmetric_map_ext")
from libtbx import slots_getstate_setstate, \
  slots_getstate_setstate_default_initializer
from libtbx.utils import null_out
import libtbx.phil
import operator
import sys
from libtbx import adopt_init_args
import mmtbx.real_space

omit_map_phil = libtbx.phil.parse("""
n_debias_cycles = 2
  .type = int
neutral_volume_box_cushion_width = 1
  .type = float
full_resolution_map = True
  .type = bool
""")

class run(object):
  """
  Composite residual OMIT map. Details:
    Afonine et al. (2014). FEM: Feature Enhanced Map
  """
  def __init__(
        self,
        fmodel,
        crystal_gridding,
        box_size_as_fraction=0.1,
        max_boxes=2000,
        neutral_volume_box_cushion_width=2,
        full_resolution_map=True,
        log=sys.stdout):
    adopt_init_args(self, locals())
    self.map_type="mFo-DFc" # using other map types is much worse
    self.sd = self.get_p1_map_unscaled(fmodel = self.fmodel,
      map_type=self.map_type).sample_standard_deviation()
    xrs_p1 = fmodel.xray_structure.expand_to_p1(sites_mod_positive=True)
    self.sgt = fmodel.f_obs().space_group().type()
    self.r = flex.double() # for tests only
    self.acc_asu = None
    # bulk-solvent
    mp = mmtbx.masks.mask_master_params.extract()
    mp.n_real = crystal_gridding.n_real()
    mp.step=None
    mmtbx_masks_asu_mask_obj = mmtbx.masks.asu_mask(
      xray_structure = xrs_p1,
      mask_params = mp)
    self.bulk_solvent_mask = mmtbx_masks_asu_mask_obj.mask_data_whole_uc()
    self.bulk_solvent_mask_asu=self.to_asu_map(map_data=self.bulk_solvent_mask)
    # atom map
    self.atom_map = mmtbx.real_space.sampled_model_density(
      xray_structure = xrs_p1,
      n_real         = crystal_gridding.n_real()).data()
    self.n_real = self.atom_map.focus()
    self.atom_map_asu=self.to_asu_map(map_data=self.atom_map)
    # extras
    self.zero_cmpl_ma = fmodel.f_calc().customized_copy(
      data = flex.complex_double(fmodel.f_calc().size(), 0))
    self.r_factor_omit = flex.double() # for regression test only
    # result map
    self.map_result_asu = flex.double(flex.grid(self.atom_map_asu.focus()))
    # iterate over boxes
    self.box_iterator()

  def map_coefficients(self, filter_noise=True):
    map_coefficients = self.result_as_sf()
    if(filter_noise):
      return self.noise_filtered_map_coefficients(
        map_coefficients=map_coefficients)
    else:
      return map_coefficients

  def noise_filtered_map_coefficients(self, map_coefficients):
    from mmtbx.maps import fem
    f_map = self.fmodel.electron_density_map().map_coefficients(
      map_type                   = "2mFo-DFc",
      isotropize                 = True,
      exclude_free_r_reflections = False,
      fill_missing               = True)
    selection = fem.good_atoms_selection(
      crystal_gridding = self.crystal_gridding,
      map_coeffs       = f_map,#map_coefficients,
      xray_structure   = self.fmodel.xray_structure)
    fft_map = cctbx.miller.fft_map(
      crystal_gridding     = self.crystal_gridding,
      fourier_coefficients = map_coefficients)
    fft_map.apply_sigma_scaling()
    m = fft_map.real_map_unpadded()
    m = fem.low_volume_density_elimination(m=m, fmodel=self.fmodel,
      selection=selection, end=11)
    ### Filter by 2mFo-DFc filled map
    fft_map = cctbx.miller.fft_map(
      crystal_gridding     = self.crystal_gridding,
      fourier_coefficients = f_map)
    fft_map.apply_sigma_scaling()
    filter_mask = fft_map.real_map_unpadded()
    filter_mask = fem.low_volume_density_elimination(m=filter_mask,
      fmodel=self.fmodel, selection=selection, end=6)#11)
    sel = filter_mask<0.25
    filter_mask = filter_mask.set_selected(sel, 0)
    sel = filter_mask>=0.25
    filter_mask = filter_mask.set_selected(sel, 1)
    ###
    return map_coefficients.structure_factors_from_map(
      map            = m*filter_mask,
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = False)

  def box_iterator(self):
    b = maptbx.boxes(
      n_real   = self.atom_map_asu.focus(),
      fraction = self.box_size_as_fraction,
      max_boxes= self.max_boxes,
      log      = self.log)
    def get_wide_box(s,e): # define wide box: neutral + phased volumes
      if(self.neutral_volume_box_cushion_width>0):
        sh = self.neutral_volume_box_cushion_width
        ss = [max(s[i]-sh,0) for i in [0,1,2]]
        ee = [min(e[i]+sh,n_real_asu[i]) for i in [0,1,2]]
      else: ss,ee = s,e
      return ss,ee
    n_real_asu = b.n_real
    n_boxes = len(b.starts)
    i_box = 0
    for s,e in zip(b.starts, b.ends):
      i_box+=1
      sw,ew = get_wide_box(s=s,e=e)
      fmodel_omit = self.omit_box(start=sw, end=ew)
      r = fmodel_omit.r_work()
      self.r.append(r) # for tests only
      if(self.log):
        print("r(curr,min,max,mean)=%6.4f %6.4f %6.4f %6.4f"%(r,
          flex.min(self.r), flex.max(self.r), flex.mean(self.r)), i_box, n_boxes, file=self.log)
      omit_map_data = self.asu_map_from_fmodel(
        fmodel=fmodel_omit, map_type=self.map_type)
      maptbx.copy_box(
        map_data_from = omit_map_data,
        map_data_to   = self.map_result_asu,
        start         = s,
        end           = e)
    self.map_result_asu.reshape(self.acc_asu)

  def result_as_sf(self):
    asu_map_omit = asu_map_ext.asymmetric_map(
      self.sgt, self.map_result_asu, self.n_real)
    map_coefficients = self.zero_cmpl_ma.customized_copy(
      indices = self.zero_cmpl_ma.indices(),
      data    = asu_map_omit.structure_factors(self.zero_cmpl_ma.indices()))
    # full resolution map (reflections in sphere, not in box!)
    if(self.full_resolution_map):
      full_set = self.zero_cmpl_ma.complete_set(d_min=self.zero_cmpl_ma.d_min())
      fill = self.zero_cmpl_ma.customized_copy(
        indices = full_set.indices(),
        data    = asu_map_omit.structure_factors(full_set.indices()))
      map_coefficients = map_coefficients.complete_with(
        other=fill, scale=True)
    return map_coefficients

  def to_asu_map(self, map_data):
    r = asu_map_ext.asymmetric_map(self.sgt, map_data).data()
    if(self.acc_asu is None): self.acc_asu = r.accessor()
    # Why I'm doing this? This may be a bug for some SG!
    # See: exercise_structure_factors_from_map_and_asu_map in miller tests.
    return r.shift_origin()

  def omit_box(self, start, end):
    def asu_omit_map_as_sf(m_asu, s, e):
      m_asu_omit = maptbx.set_box_copy(value = 0, map_data_to = m_asu,
        start = s, end = e)
      m_asu_omit.reshape(self.acc_asu)
      asu_m_omit = asu_map_ext.asymmetric_map(self.sgt, m_asu_omit, self.n_real)
      return self.zero_cmpl_ma.customized_copy(
        indices = self.zero_cmpl_ma.indices(),
        data    = asu_m_omit.structure_factors(self.zero_cmpl_ma.indices()))
    f_calc_omit = asu_omit_map_as_sf(m_asu=self.atom_map_asu, s=start, e=end)
    f_mask_omit = asu_omit_map_as_sf(m_asu=self.bulk_solvent_mask_asu,
      s=start, e=end)
    fmodel_omit = mmtbx.f_model.manager(
      f_obs        = self.fmodel.f_obs(),
      r_free_flags = self.fmodel.r_free_flags(),
      k_isotropic  = self.fmodel.k_isotropic(),
      k_anisotropic= self.fmodel.k_anisotropic(),
      k_mask       = self.fmodel.k_masks(),
      f_calc       = f_calc_omit,
      f_mask       = f_mask_omit)
    return fmodel_omit

  def as_p1_map(self, map_asu=None):
    if(map_asu is None): map_asu = self.map_result_asu
    asu_map_omit = asu_map_ext.asymmetric_map(self.sgt, map_asu, self.n_real)
    p1_map = asu_map_omit.symmetry_expanded_map()
    maptbx.unpad_in_place(map=p1_map)
    sd = p1_map.sample_standard_deviation()
    if(sd != 0):
      p1_map = p1_map/sd
    return p1_map

  def get_p1_map_unscaled(self, fmodel, map_type):
    f_map = fmodel.electron_density_map().map_coefficients(
      map_type                   = map_type,
      isotropize                 = True,
      exclude_free_r_reflections = True,
      fill_missing               = False)
    fft_map = cctbx.miller.fft_map(
      crystal_gridding     = self.crystal_gridding,
      fourier_coefficients = f_map)
    return fft_map.real_map_unpadded()

  def asu_map_from_fmodel(self, fmodel, map_type):
    map_data = self.get_p1_map_unscaled(fmodel=fmodel, map_type=map_type)
    return asu_map_ext.asymmetric_map(self.sgt, map_data).data().shift_origin()

################################################################################

########################################################################
# XXX LEGACY IMPLEMENTATION
#
# The code below is similar to the method used by CNS, which operates on
# xray scatterers rather than the map itself.  As a result it requires
# refinement (w/ annealing) of the partial model to thoroughly de-bias
# the phases, and is therefore much slower than the protocol above.  It
# is preserved for maximum flexibility and because the same underlying
# method is also useful for analyzing static disorder at high resolution.
#
# Note that most of the program logic lives in the phenix tree due to
# its use of phenix.refine.

class omit_box(slots_getstate_setstate_default_initializer):
  """
  Defines a region in fractional coordinates containing a selection of atoms
  to be omitted.
  """
  __slots__ = [ "frac_min", "frac_max", "selection", "serial" ]

  @property
  def n_atoms(self):
    return len(self.selection)

  def show(self, out=sys.stdout, prefix=""):
    print(prefix + \
      "box %d: atoms: %6d  extents: (%.4f, %.4f, %.4f) to (%.4f, %.4f, %.4f)" \
      % tuple([ self.serial, len(self.selection) ] +
        list(self.frac_min) + list(self.frac_max)), file=out)

class omit_regions(slots_getstate_setstate):
  """
  Groups together multiple omit_box objects (without any implicit spatial
  relationship); this is used to reduce the total number of bins of omitted
  atoms and to make the fraction omitted per bin approximately similar.
  """
  __slots__ = [ "boxes", "selection", "serial" ]
  def __init__(self, serial, selection=None):
    self.serial = serial
    self.boxes = []
    self.selection = selection
    if (selection is None):
      self.selection = flex.size_t()

  @property
  def n_boxes(self):
    return len(self.boxes)

  @property
  def n_atoms(self):
    return len(self.selection)

  def add_box(self, box):
    self.boxes.append(box)
    self.selection = join_selections(box.selection, self.selection)
    return self

  def combine_with(self, other):
    for box in other.boxes :
      self.add_box(box)
    return self

  def show(self, out=sys.stdout, prefix=""):
    if (len(self.boxes) == 0):
      print(prefix + "Region %d: empty" % self.serial, file=out)
    else :
      print(prefix + "Region %d:" % self.serial, file=out)
      for box in self.boxes :
        box.show(out=out, prefix=prefix+"  ")

class omit_region_results(slots_getstate_setstate_default_initializer):
  __slots__ = [ "serial", "map_coeffs_list", "r_work", "r_free" ]

def create_omit_regions(xray_structure,
    selection=None,
    fraction_omit=0.05,
    optimize_binning=True,
    box_cushion_radius=2.5,
    even_boxing=False,
    asu_buffer_thickness=1.e-5,
    log=None):
  """
  Divide the asymmetric unit into boxes of atoms to omit.  This will include
  a cushion around the actual omit region.  Although the step size is uniform,
  by default boxes will be grouped as needed to make the distribution roughly
  even and keep the total number of omit regions to a minimum.
  If a region does not contain any atoms, it will be skipped and filled in
  with the background map.

  Returns a list of omit_regions objects, which specify the selection of atoms
  to omit along with the boxes of interest (in fractional coordinates).
  """
  if (log is None) : log = null_out()
  if (selection is None):
    selection = flex.bool(xray_structure.scatterers().size(), True)
  iselection = selection.iselection()
  boxes = []
  xrs = xray_structure.sites_mod_positive()
  unit_cell = xrs.unit_cell()
  asu_mappings = xrs.asu_mappings(buffer_thickness=asu_buffer_thickness)
  asu = asu_mappings.asu()
  #print asu.box_min(), asu.box_max()
  sites_asu = flex.vec3_double()
  for i_seq, mappings in enumerate(asu_mappings.mappings()):
    if (not selection[i_seq]) : continue
    for mapping in mappings :
      site_cart = mapping.mapped_site()
      site_frac = unit_cell.fractionalize(site_cart=site_cart)
      sites_asu.append(site_frac)
      break
  x_cushion, y_cushion, z_cushion = unit_cell.fractionalize(
    site_cart=[ box_cushion_radius ] * 3)
  x_min, y_min, z_min = sites_asu.min()
  x_max, y_max, z_max = sites_asu.max()
  #print "min", x_min, y_min, z_min
  #print "max", x_max, y_max, z_max
  non_hd_sel = ~(xray_structure.hd_selection())
  n_omit_heavy_atoms = non_hd_sel.select(iselection).count(True)
  n_omit_per_box = n_omit_heavy_atoms * fraction_omit
  assert (n_omit_heavy_atoms > 0)
  if (even_boxing) : # XXX remove?
    x_step = asu_buffer_thickness + (x_max - x_min) / 4
    y_step = asu_buffer_thickness + (y_max - y_min) / 4
    z_step = asu_buffer_thickness + (z_max - z_min) / 2
  else :
    volume = n_omit_heavy_atoms * 9.0
    cube_length = (volume * fraction_omit) ** (1/3.)
    x_step, y_step, z_step = unit_cell.fractionalize(site_cart=[cube_length]*3)
  # to facilitate using the C++ loops to determine which sites fall inside the
  # box, we first split the fractional coordinates into separate x, y, z 1D
  # arrays.  not sure if there's a better way.
  sites_1d = sites_asu.as_double()
  sel_base = flex.size_t(range(len(sites_asu)))
  sites_x = sites_1d.select(sel_base*3)
  sites_y = sites_1d.select(sel_base*3+1)
  sites_z = sites_1d.select(sel_base*3+2)
  x_start = x_min
  serial = 1
  while (x_start < x_max):
    y_start = y_min
    while (y_start < y_max):
      z_start = z_min
      while (z_start < z_max):
        x_end = x_start + x_step
        y_end = y_start + y_step
        z_end = z_start + z_step
        x_min_box = x_start - x_cushion
        x_max_box = x_end + x_cushion
        y_min_box = y_start - y_cushion
        y_max_box = y_end + y_cushion
        z_min_box = z_start - z_cushion
        z_max_box = z_end + z_cushion
        box_selection = ((sites_x >= x_min_box) & (sites_x < x_max_box) &
                         (sites_y >= y_min_box) & (sites_y < y_max_box) &
                         (sites_z >= z_min_box) & (sites_z < z_max_box))
        box_iselection = box_selection.iselection()
        n_atoms_current_box = len(box_iselection)
        if (n_atoms_current_box > 0):
          frac_max = [min(x_end, 1.0), min(y_end, 1.0), min(z_end, 1.0)]
          new_box = omit_box(
            frac_min=[x_start, y_start, z_start],
            frac_max=frac_max,
            selection=iselection.select(box_iselection),
            serial=serial)
          boxes.append(new_box)
          serial += 1
        z_start += z_step
      y_start += y_step
    x_start += x_step
  groups = []
  for box in boxes :
    group = omit_regions(serial=box.serial).add_box(box)
    groups.append(group)
  if (optimize_binning):
    # http://en.wikipedia.org/wiki/Bin_packing_problem
    groups = sorted(groups, key=operator.attrgetter("n_atoms"), reverse=True)
    while True :
      n_combined = 0
      i_group = 0
      while i_group < len(groups):
        groups[i_group].serial = i_group + 1
        j_group = 0
        while j_group < i_group :
          n_atoms_combined = groups[i_group].n_atoms + groups[j_group].n_atoms
          if (n_atoms_combined <= n_omit_per_box):
            groups[j_group].combine_with(groups[i_group])
            del groups[i_group]
            n_combined += 1
            break
          j_group += 1
        i_group += 1
      if (n_combined == 0):
        break
  assert (len(set([ g.serial for g in groups ])) == len(groups))
  return groups

def join_selections(sel1, sel2):
  intersections = sel1.intersection_i_seqs(sel2)
  unique_sel = flex.bool(len(sel1), True)
  unique_sel.set_selected(intersections[0], False)
  sel1 = sel1.select(unique_sel)
  sel1.extend(sel2)
  return flex.sorted(sel1)

def combine_maps(
    map_arrays,
    omit_groups,
    background_map_coeffs,
    resolution_factor,
    flatten_background=False,
    sigma_scaling=False,
    control_map=False):
  """
  For each box, FFT the corresponding omit map coefficients, extract the
  omit regions, and copy them to the combined map, using Marat's asymmetric
  map class to handle symmetry expansion internally.
  """
  assert len(map_arrays) == len(omit_groups)
  space_group = background_map_coeffs.space_group()
  fft_map = background_map_coeffs.fft_map(
    symmetry_flags=maptbx.use_space_group_symmetry,
    resolution_factor=resolution_factor).apply_volume_scaling()
  if (sigma_scaling):
    fft_map.apply_sigma_scaling()
  background_map = fft_map.real_map_unpadded()
  #print "full map:", background_map.focus()
  if (flatten_background):
    sel_all = flex.bool(background_map.as_1d().size(), True)
    background_map.as_1d().set_selected(sel_all, 0)
  asym_map = asu_map_ext.asymmetric_map(space_group.type(), background_map)
  origin = asym_map.data().origin()
  n_real_all = background_map.focus()
  n_real_asu = asym_map.data().focus()
  m_real_asu = asym_map.data().all()
  #print "ORIGIN:", origin
  #print "N_REAL:", n_real_asu
  #print "M_REAL:", m_real_asu
  if (not control_map):
    f2g = maptbx.frac2grid(n_real_all)
    n = 0
    for group, map_coeffs in zip(omit_groups, map_arrays):
      space_group = map_coeffs.space_group()
      omit_fft_map = map_coeffs.fft_map(
        resolution_factor=resolution_factor,
        symmetry_flags=maptbx.use_space_group_symmetry).apply_volume_scaling()
      if (sigma_scaling):
        omit_fft_map.apply_sigma_scaling()
      omit_map = omit_fft_map.real_map_unpadded()
      omit_asym_map = asu_map_ext.asymmetric_map(space_group.type(), omit_map)
      assert omit_asym_map.data().focus() == n_real_asu
      for box in group.boxes :
        if (control_map) : continue
        grid_start = tuple(f2g(box.frac_min))
        grid_end = grid_max(f2g(box.frac_max), n_real_asu)
        #print grid_start, grid_end
        for u in range(grid_start[0], grid_end[0]+1):
          for v in range(grid_start[1], grid_end[1]+1):
            for w in range(grid_start[2], grid_end[2]+1):
              try :
                asym_map.data()[(u,v,w)] = omit_asym_map.data()[(u,v,w)]
              except IndexError :
                raise IndexError("n_real: %s  origin: %s  index: %s" %
                  (str(n_real_asu), str(origin), str((u,v,w))))
  return asym_map.map_for_fft()

def grid_max(grid_coords, n_real):
  return (min(grid_coords[0], n_real[0]-1), min(grid_coords[1], n_real[1]-1),
          min(grid_coords[2], n_real[2]-1))
