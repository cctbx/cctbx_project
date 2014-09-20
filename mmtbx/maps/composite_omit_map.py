
"""
Methods and utilities for generating composite omit maps.  The actual
end-user application, phenix.composite_omit_map, is part of the phenix
sources.
"""

from __future__ import division
import mmtbx.f_model
from cctbx import miller
from cctbx import maptbx
import cctbx.miller
from scitbx.array_family import flex
import boost.python
asu_map_ext = boost.python.import_ext("cctbx_asymmetric_map_ext")
from libtbx import slots_getstate_setstate, \
  slots_getstate_setstate_default_initializer
from libtbx.utils import null_out
import libtbx.phil
import sys

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
  Composite full-omit map: omit entire box in real-space corresponding to
  Fmodel, which includs atomic model and non-atomic model (bulk-solvent and
  scales).  This is much faster than the refinement-based method, and is used
  by default in phenix.composite_omit_map.
  """
  def __init__(
        self,
        crystal_gridding,
        fmodel,
        map_type,
        box_size_as_fraction=0.03,
        n_debias_cycles=2,
        neutral_volume_box_cushion_width=1,
        full_resolution_map=True,
        log=sys.stdout):
    self.crystal_gridding = crystal_gridding
    # assert compatibility of symops with griding
    assert self.crystal_gridding._symmetry_flags is not None
    self.sgt = fmodel.f_obs().space_group().type()
    self.zero_cmpl_ma = fmodel.f_calc().customized_copy(
      data = flex.complex_double(fmodel.f_calc().size(), 0))
    # embedded utility functions
    def get_map(fmodel, map_type, crystal_gridding, asu=True):
      f_map = fmodel.electron_density_map().map_coefficients(
        map_type                   = map_type,
        isotropize                 = True,
        exclude_free_r_reflections = True,
        fill_missing               = False)
      fft_map = cctbx.miller.fft_map(
        crystal_gridding     = crystal_gridding,
        fourier_coefficients = f_map)
      if(asu): return asu_map_ext.asymmetric_map(self.sgt,
        fft_map.real_map_unpadded()).data()
      else:
        return fft_map.real_map_unpadded()
    # f_model map
    f_model_map_data = fmodel.f_model_scaled_with_k1().fft_map(
      symmetry_flags   = maptbx.use_space_group_symmetry,
      crystal_gridding = self.crystal_gridding).real_map_unpadded()
    self.n_real = f_model_map_data.focus()
    # extract asu map from full P1
    f_model_map_data_asu=asu_map_ext.asymmetric_map(
      self.sgt, f_model_map_data).data()
    self.acc = f_model_map_data_asu.accessor()
    f_model_map_data_asu = f_model_map_data_asu.shift_origin()
    # set up boxes
    b = maptbx.boxes(
      n_real   = f_model_map_data_asu.focus(),
      fraction = box_size_as_fraction,
      log      = log)
    self.map_result_asu = flex.double(flex.grid(b.n_real))
    assert f_model_map_data_asu.focus()==b.n_real
    assert b.n_real==self.map_result_asu.focus()
    n_real_asu = b.n_real
    self.r = flex.double() # for regression test only
    n_boxes = len(b.starts)
    i_box = 0
    for s,e in zip(b.starts, b.ends):
      i_box+=1
      # define wide box: neutral + phased volumes
      if(neutral_volume_box_cushion_width>0):
        sh = neutral_volume_box_cushion_width
        ss = [max(s[i]-sh,0) for i in [0,1,2]]
        ee = [min(e[i]+sh,n_real_asu[i]) for i in [0,1,2]]
      else: ss,ee = s,e
      # omit wide box from f_model map, repeat n_debias_cycles times
      f_model_map_data_asu_ = f_model_map_data_asu.deep_copy()
      for i in xrange(n_debias_cycles):
        f_model_omit, f_model_map_data_asu_ = self.omit_box(s=ss, e=ee,
          md_asu=f_model_map_data_asu_)
      # get fmodel for omit map calculation
      fmodel_ = mmtbx.f_model.manager(
        f_obs        = fmodel.f_obs(),
        r_free_flags = fmodel.r_free_flags(),
        f_calc       = f_model_omit,
        f_mask       = self.zero_cmpl_ma)
      rw = fmodel_.r_work()
      self.r.append(rw) # for regression test only
      f_map_data_asu = get_map(fmodel=fmodel_, map_type=map_type,
        crystal_gridding=self.crystal_gridding)
      f_map_data_asu = f_map_data_asu.shift_origin()
      if(log):
        print >> log, "box %2d of %2d:"%(i_box, n_boxes), s, e, "%6.4f"%rw
      assert f_map_data_asu.focus() == self.map_result_asu.focus()
      maptbx.copy_box(
        map_data_from = f_map_data_asu,
        map_data_to   = self.map_result_asu,
        start         = s,
        end           = e)
    # result
    self.map_result_asu.reshape(self.acc)
    self.asu_map_omit = asu_map_ext.asymmetric_map(
      self.sgt, self.map_result_asu, self.n_real)
    self.map_coefficients = self.zero_cmpl_ma.customized_copy(
      indices = self.zero_cmpl_ma.indices(),
      data    = self.asu_map_omit.structure_factors(self.zero_cmpl_ma.indices()))
    # full resolution map (reflections in sphere, not in box!)
    if(full_resolution_map):
      cs = self.zero_cmpl_ma.complete_set(d_min=self.zero_cmpl_ma.d_min())
      asu_map_omit = asu_map_ext.asymmetric_map(
        self.sgt,self.map_result_asu,self.n_real)
      fill = self.zero_cmpl_ma.customized_copy(
        indices = cs.indices(),
        data    = asu_map_omit.structure_factors(cs.indices()))
      self.map_coefficients = self.map_coefficients.complete_with(
        other=fill, scale=True)

  def omit_box(self, s, e, md_asu):
    md_asu_omit = maptbx.set_box_copy(value = 0, map_data_to = md_asu,
      start = s, end = e)
    md_asu_omit.reshape(self.acc)
    asu_map_omit = asu_map_ext.asymmetric_map(self.sgt, md_asu_omit, self.n_real)
    ma_omit = self.zero_cmpl_ma.customized_copy(
      indices = self.zero_cmpl_ma.indices(),
      data    = asu_map_omit.structure_factors(self.zero_cmpl_ma.indices()))
    fft_map = cctbx.miller.fft_map(
      crystal_gridding = self.crystal_gridding, fourier_coefficients = ma_omit)
    md = fft_map.real_map_unpadded()
    asu_map = asu_map_ext.asymmetric_map(self.sgt, md)
    md_asu_omit = asu_map.data()
    md_asu_omit = md_asu_omit.shift_origin()
    return ma_omit, md_asu_omit

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

class omit_box (slots_getstate_setstate_default_initializer) :
  """
  Defines a region in fractional coordinates containing a selection of atoms
  to be omitted.
  """
  __slots__ = [ "frac_min", "frac_max", "selection", "serial" ]

  @property
  def n_atoms (self) :
    return len(self.selection)

  def show (self, out=sys.stdout, prefix="") :
    print >> out, prefix + \
      "box %d: atoms: %6d  extents: (%.4f, %.4f, %.4f) to (%.4f, %.4f, %.4f)" \
      % tuple([ self.serial, len(self.selection) ] +
        list(self.frac_min) + list(self.frac_max))

class omit_regions (slots_getstate_setstate) :
  """
  Groups together multiple omit_box objects (without any implicit spatial
  relationship); this is used to reduce the total number of bins of omitted
  atoms and to make the fraction omitted per bin approximately similar.
  """
  __slots__ = [ "boxes", "selection", "serial" ]
  def __init__ (self, serial, selection=None) :
    self.serial = serial
    self.boxes = []
    self.selection = selection
    if (selection is None) :
      self.selection = flex.size_t()

  @property
  def n_boxes (self) :
    return len(self.boxes)

  @property
  def n_atoms (self) :
    return len(self.selection)

  def add_box (self, box) :
    self.boxes.append(box)
    self.selection = join_selections(box.selection, self.selection)
    return self

  def combine_with (self, other) :
    for box in other.boxes :
      self.add_box(box)
    return self

  def show (self, out=sys.stdout, prefix="") :
    if (len(self.boxes) == 0) :
      print >> out, prefix + "Region %d: empty" % self.serial
    else :
      print >> out, prefix + "Region %d:" % self.serial
      for box in self.boxes :
        box.show(out=out, prefix=prefix+"  ")

class omit_region_results (slots_getstate_setstate_default_initializer) :
  __slots__ = [ "serial", "map_coeffs_list", "r_work", "r_free" ]

def create_omit_regions (xray_structure,
    selection=None,
    fraction_omit=0.05,
    optimize_binning=True,
    box_cushion_radius=2.5,
    even_boxing=False,
    asu_buffer_thickness=1.e-5,
    log=None) :
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
  if (selection is None) :
    selection = flex.bool(xray_structure.scatterers().size(), True)
  iselection = selection.iselection()
  boxes = []
  xrs = xray_structure.sites_mod_positive()
  unit_cell = xrs.unit_cell()
  asu_mappings = xrs.asu_mappings(buffer_thickness=asu_buffer_thickness)
  asu = asu_mappings.asu()
  #print asu.box_min(), asu.box_max()
  sites_asu = flex.vec3_double()
  for i_seq, mappings in enumerate(asu_mappings.mappings()) :
    if (not selection[i_seq]) : continue
    for mapping in mappings :
      site_cart = mapping.mapped_site()
      print site_cart
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
  while (x_start < x_max) :
    y_start = y_min
    while (y_start < y_max) :
      z_start = z_min
      while (z_start < z_max) :
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
        if (n_atoms_current_box > 0) :
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
  if (optimize_binning) :
    # http://en.wikipedia.org/wiki/Bin_packing_problem
    groups = sorted(groups, lambda a,b: cmp(b.n_atoms, a.n_atoms))
    while True :
      n_combined = 0
      i_group = 0
      while i_group < len(groups) :
        groups[i_group].serial = i_group + 1
        j_group = 0
        while j_group < i_group :
          n_atoms_combined = groups[i_group].n_atoms + groups[j_group].n_atoms
          if (n_atoms_combined <= n_omit_per_box) :
            groups[j_group].combine_with(groups[i_group])
            del groups[i_group]
            n_combined += 1
            break
          j_group += 1
        i_group += 1
      if (n_combined == 0) :
        break
  assert (len(set([ g.serial for g in groups ])) == len(groups))
  return groups

def join_selections (sel1, sel2) :
  intersections = sel1.intersection_i_seqs(sel2)
  unique_sel = flex.bool(len(sel1), True)
  unique_sel.set_selected(intersections[0], False)
  sel1 = sel1.select(unique_sel)
  sel1.extend(sel2)
  return flex.sorted(sel1)

def combine_maps (
    map_arrays,
    omit_groups,
    background_map_coeffs,
    resolution_factor,
    flatten_background=False,
    sigma_scaling=False,
    control_map=False) :
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
  if (sigma_scaling) :
    fft_map.apply_sigma_scaling()
  background_map = fft_map.real_map_unpadded()
  #print "full map:", background_map.focus()
  if (flatten_background) :
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
  if (not control_map) :
    f2g = maptbx.frac2grid(n_real_all)
    n = 0
    for group, map_coeffs in zip(omit_groups, map_arrays) :
      space_group = map_coeffs.space_group()
      omit_fft_map = map_coeffs.fft_map(
        resolution_factor=resolution_factor,
        symmetry_flags=maptbx.use_space_group_symmetry).apply_volume_scaling()
      if (sigma_scaling) :
        omit_fft_map.apply_sigma_scaling()
      omit_map = omit_fft_map.real_map_unpadded()
      omit_asym_map = asu_map_ext.asymmetric_map(space_group.type(), omit_map)
      assert omit_asym_map.data().focus() == n_real_asu
      for box in group.boxes :
        if (control_map) : continue
        grid_start = tuple(f2g(box.frac_min))
        grid_end = grid_max(f2g(box.frac_max), n_real_asu)
        #print grid_start, grid_end
        for u in range(grid_start[0], grid_end[0]+1) :
          for v in range(grid_start[1], grid_end[1]+1) :
            for w in range(grid_start[2], grid_end[2]+1) :
              try :
                asym_map.data()[(u,v,w)] = omit_asym_map.data()[(u,v,w)]
              except IndexError :
                raise IndexError("n_real: %s  origin: %s  index: %s" %
                  (str(n_real_asu), str(origin), str((u,v,w))))
  return asym_map.map_for_fft()

def grid_max (grid_coords, n_real) :
  return (min(grid_coords[0], n_real[0]-1), min(grid_coords[1], n_real[1]-1),
          min(grid_coords[2], n_real[2]-1))
