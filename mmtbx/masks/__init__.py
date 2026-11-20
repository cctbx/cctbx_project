from __future__ import absolute_import, division, print_function
import boost_adaptbx.boost.python as bp
from cctbx.array_family import flex
from six.moves import zip
from six.moves import range
ext = bp.import_ext("mmtbx_masks_ext")
from mmtbx_masks_ext import *

from cctbx.masks import around_atoms, vdw_radii_from_xray_structure
from cctbx import maptbx
import sys
import iotbx.xplor.map
import iotbx.phil
from libtbx import introspection
from libtbx import adopt_init_args
from copy import deepcopy
import mmtbx.masks
asu_map_ext = bp.import_ext("cctbx_asymmetric_map_ext")

number_of_mask_calculations = 0

# Default mask radii solvent_radius=0.95 / shrink_truncation_radius=0.7 A
# (historical values 1.1 / 0.9) were derived 2026-11 on 55 low-resolution and
# 10 high-resolution refinements with and without riding H. With explicit H
# electrons the historical mask mis-fits the d > 6 A data (Rfree +0.002 on the
# 55 models); the optimum moves by about -0.15 A in solvent_radius along the
# ridge solvent_radius - shrink_truncation_radius = 0.2-0.3 A. Without H the
# new values are neutral (Rfree -0.0005, Rwork -0.0015 on the same sets).
# Radii of the model mask used by real-space (cryo-EM) tools: CC_mask and
# related map-model correlations, the mtriage model-based soft mask. They
# define published metrics and are not tied to the calibration of the
# bulk-solvent mask above, so they keep the historical values.
real_space_solvent_radius = 1.1
real_space_shrink_truncation_radius = 0.9

mask_master_params = iotbx.phil.parse("""\
  use_asu_masks = True
    .type = bool
    .help = Compute the mask in the asymmetric unit (faster for high \
            symmetry). If False the mask is computed in the whole unit \
            cell, with the same result within the accuracy of the grid. \
            More than one radial shell requires True.
  solvent_radius = 0.95
    .type = float
    .help = Probe radius (A) added to the atomic van der Waals radii when \
            building the bulk-solvent mask.
  shrink_truncation_radius = 0.7
    .type = float
    .help = Shrink radius (A) applied to the solvent region.
  use_resolution_based_gridding = False
    .type = bool
    .help = If enabled, mask grid step = d_min / grid_step_factor (kept \
            between 0.15 and 0.8 A), and step and n_real are not used.
  step = 0.6
    .type = float
    .help = Grid step (in A) for solvent mask calculation. Set to None to \
            use n_real.
  n_real = None
    .type = ints(size=3)
    .help = Gridding (three integers) for mask calculation. Requires \
            step=None.
  grid_step_factor = 4.0
    .type = float
    .help = Grid step is resolution divided by grid_step_factor. Used only \
            if use_resolution_based_gridding=True.
  mean_shift_for_mask_update = 0.1
    .type = float
    .help = The mask is re-calculated if the mean shift of atoms (in A) \
            since the previous mask calculation is larger than this value.
    .expert_level=2
  ignore_zero_occupancy_atoms = True
    .type = bool
    .short_caption = Ignore zero-occupancy atoms when calculating bulk solvent mask (recommended)
    .help = Do not include atoms with zero occupancy into mask calculation
  ignore_hydrogens = True
    .type = bool
    .help = Ignore H or D atoms in mask calculation. Always False for \
            neutron data.
  n_radial_shells = 1
    .type = int
    .help = Number of shells in a radial shell bulk solvent model. More \
            than one shell requires bulk_solvent_and_scale.mode=slow.
  radial_shell_width = 1.3
    .type = float
    .help = Radial shell width (A). Used only if n_radial_shells > 1.
  Fmask_res_high = 3.0
    .type = float
    .help = Assume Fmask=0 for reflections outside the range inf>d>=Fmask_res_high
""")

def _sorry(problem, options):
  """
  Raise Sorry that states the problem and the ways out.
  """
  from libtbx.utils import Sorry
  lines = ["mask: %s"%problem]
  for i, o in enumerate(options):
    if(i == 0): lines.append("  Either %s"%o)
    else:       lines.append("  or %s"%o)
  raise Sorry("\n".join(lines))

def validate_mask_params(mask_params, d_min_given=False):
  """
  Stop with a clear message if mask parameters contradict each other or
  have meaningless values.

  Parameters
  ----------
  mask_params : libtbx.phil.scope_extract
      Extract of mask_master_params.
  d_min_given : bool
      True if the caller defines the gridding by resolution itself
      (asu_mask(d_min=...)); step and n_real are not checked then.

  Raises
  ------
  libtbx.utils.Sorry
  """
  p = mask_params
  for name in ["solvent_radius", "shrink_truncation_radius"]:
    value = getattr(p, name)
    if(value is None or value < 0):
      _sorry("mask.%s=%s is not possible."%(name, str(value)),
        ["set mask.%s to a non-negative number."%name])
  if(p.grid_step_factor is not None and p.grid_step_factor <= 0):
    _sorry("mask.grid_step_factor=%s is not possible."%str(
      p.grid_step_factor), ["set mask.grid_step_factor to a positive number."])
  if(p.use_resolution_based_gridding or d_min_given):
    if(p.grid_step_factor is None):
      _sorry("mask.use_resolution_based_gridding=True requires "
        "mask.grid_step_factor.", ["set mask.grid_step_factor (default 4),",
        "set mask.use_resolution_based_gridding=False."])
  else:
    if(p.step is not None and p.n_real is not None):
      _sorry("mask.step=%s and mask.n_real=%s are both set."%(str(p.step),
        ",".join([str(n) for n in p.n_real])),
        ["set mask.step=None to use mask.n_real,",
         "remove mask.n_real to use mask.step."])
    if(p.step is None and p.n_real is None):
      _sorry("mask.step=None and mask.n_real=None: the grid is not defined.",
        ["set mask.step (default 0.6),", "set mask.n_real,",
         "set mask.use_resolution_based_gridding=True."])
  if(p.step is not None and p.step <= 0):
    _sorry("mask.step=%s is not possible."%str(p.step),
      ["set mask.step to a positive number (default 0.6)."])
  if(p.n_real is not None and min(p.n_real) < 1):
    _sorry("mask.n_real=%s is not possible."%",".join(
      [str(n) for n in p.n_real]), ["set three positive integers."])
  if(p.n_radial_shells is None or p.n_radial_shells < 1):
    _sorry("mask.n_radial_shells=%s is not possible."%str(p.n_radial_shells),
      ["set mask.n_radial_shells=1 or more."])
  if(p.n_radial_shells > 1):
    if(not p.use_asu_masks):
      _sorry("mask.n_radial_shells=%d requires mask.use_asu_masks=True."%
        p.n_radial_shells, ["set mask.use_asu_masks=True,",
        "set mask.n_radial_shells=1."])
    if(p.radial_shell_width is None or p.radial_shell_width <= 0):
      _sorry("mask.radial_shell_width=%s is not possible."%str(
        p.radial_shell_width),
        ["set mask.radial_shell_width to a positive number (default 1.3),",
         "set mask.n_radial_shells=1."])
  for name in ["mean_shift_for_mask_update", "Fmask_res_high"]:
    value = getattr(p, name)
    if(value is None or value < 0):
      _sorry("mask.%s=%s is not possible."%(name, str(value)),
        ["set mask.%s to a non-negative number."%name])

def gridding_source(mask_params, d_min=None):
  """
  One-line description of what defines the mask grid, for logs.
  """
  p = mask_params
  if(p.use_resolution_based_gridding):
    if(d_min is None): return "d_min/grid_step_factor"
    return "d_min/grid_step_factor=%.2f/%.2f"%(d_min, p.grid_step_factor)
  if(p.step is not None): return "step=%.2f"%p.step
  return "n_real=%s"%",".join([str(n) for n in p.n_real])

class bulk_solvent(around_atoms):

  def __init__(self,
        xray_structure,
        ignore_zero_occupancy_atoms,
        solvent_radius,
        shrink_truncation_radius,
        ignore_hydrogen_atoms=True,
        gridding_n_real=None,
        grid_step=None,
        atom_radii=None):
     global number_of_mask_calculations
     number_of_mask_calculations += 1
     assert [gridding_n_real, grid_step].count(None) == 1
     self.xray_structure = xray_structure
     if (gridding_n_real is None):
       gridding_n_real = maptbx.crystal_gridding(
         unit_cell=xray_structure.unit_cell(),
         step=grid_step).n_real()
     if(atom_radii is None):
       atom_radii = vdw_radii_from_xray_structure(xray_structure =
         self.xray_structure)
     sites_frac = xray_structure.sites_frac()
     self.n_atoms_excluded = 0
     selection = flex.bool(xray_structure.scatterers().size(), True)
     if(ignore_zero_occupancy_atoms):
       selection &= xray_structure.scatterers().extract_occupancies() > 0
     if(ignore_hydrogen_atoms):
       selection &= (~xray_structure.hd_selection())
     sites_frac = sites_frac.select(selection)
     atom_radii = atom_radii.select(selection)
     self.n_atoms_excluded = selection.count(False)
     around_atoms.__init__(self,
       unit_cell           = xray_structure.unit_cell(),
       space_group_order_z = xray_structure.space_group().order_z(),
       sites_frac          = sites_frac,
       atom_radii          = atom_radii,
       gridding_n_real     = gridding_n_real,
       solvent_radius      = solvent_radius,
       shrink_truncation_radius = shrink_truncation_radius)
     introspection.virtual_memory_info().update_max()

  def show_summary(self, out=None):
    if (out is None): out = sys.stdout
    print("|- mask summary -----------------------------------------"\
                  "---------------------|", file=out)
    sr = "solvent radius:            %4.2f A" % self.solvent_radius
    st = "shrink truncation radius:  %4.2f A" % self.shrink_truncation_radius
    na = "number of atoms: %d" % self.xray_structure.scatterers().size()
    n_real = self.data.accessor().focus()
    gr = "gridding:   %s" % str(n_real)
    gs = "grid steps: (%s) A" % ", ".join(["%.2f" % (l/n) for l,n in zip(
      self.xray_structure.unit_cell().parameters()[:3],
      n_real)])
    sc = "estimated solvent content: %.1f %%" % (
      self.contact_surface_fraction*100)
    l = max(len(sr), len(st), len(na))
    sr += " "*(l-len(sr))
    st += " "*(l-len(st))
    na += " "*(l-len(na))
    def show_line(s):
      print("| " + s + " " * max(0, 75-len(s)) + " |", file=out)
    gap = 6
    show_line(s=sr+" "*gap+gr)
    show_line(s=st+" "*gap+gs)
    show_line(s=na+" "*gap+sc)
    print("|"+"-"*77+"|", file=out)

  def mask_as_xplor_map(self, file_name):
    gridding = iotbx.xplor.map.gridding(n     = self.data.focus(),
                                        first = (0,0,0),
                                        last  = self.data.focus())
    iotbx.xplor.map.writer(
                          file_name          = file_name,
                          is_p1_cell         = True,
                          title_lines        = [' REMARKS Mask map""',],
                          unit_cell          = self.xray_structure.unit_cell(),
                          gridding           = gridding,
                          data               = self.data.as_double(),
                          average            = -1,
                          standard_deviation = -1)

  def structure_factors(self, miller_set, zero_high=True):
    if(zero_high):
      sel = miller_set.d_spacings().data()>=3
      miller_set_sel = miller_set.select(sel)
      result = miller_set_sel.structure_factors_from_map(
        map = self.data, use_scale = True, anomalous_flag = False, use_sg = True)
      data = flex.complex_double(miller_set.indices().size(), 0)
      data = data.set_selected(sel, result.data())
      return miller_set.array(data=data)
    else:
      result = miller_set.structure_factors_from_map(
        map = self.data, use_scale = True, anomalous_flag = False, use_sg = True)
      return miller_set.array(data=result.data())

  def subtract_non_uniform_solvent_region_in_place(self, non_uniform_mask):
    assert non_uniform_mask.accessor() == self.data.accessor()
    self.data.set_selected(non_uniform_mask > 0, 0)

def gridding_n_real(mask_params, xray_structure, d_min=None):
  """
  Mask gridding from mask parameters: d_min/grid_step_factor if d_min is
  given, else step if given, else n_real.

  Parameters
  ----------
  mask_params : libtbx.phil.scope_extract
  xray_structure : cctbx.xray.structure
  d_min : float, optional
      Resolution; give it only for resolution-based gridding.

  Returns
  -------
  tuple of int
  """
  if(d_min is not None):
    # Same as atom_mask::determine_gridding (atom_mask.cpp): step limits and
    # symmetry flags.
    from cctbx import sgtbx
    step = min(0.8, max(0.15, d_min/mask_params.grid_step_factor))
    return maptbx.crystal_gridding(
      unit_cell         = xray_structure.unit_cell(),
      space_group_info  = xray_structure.space_group_info(),
      symmetry_flags    = sgtbx.search_symmetry_flags(
        use_space_group_symmetry = True,
        use_space_group_ltr      = 0,
        use_seminvariants        = True,
        use_normalizer_k2l       = True,
        use_normalizer_l2n       = True),
      d_min             = 2*step,
      resolution_factor = 0.5).n_real()
  if(mask_params.step is None): return tuple(mask_params.n_real)
  return maptbx.crystal_gridding(
    unit_cell        = xray_structure.unit_cell(),
    space_group_info = xray_structure.space_group_info(),
    symmetry_flags   = maptbx.use_space_group_symmetry,
    step             = mask_params.step).n_real()

class asu_mask(object):
  def __init__(self,
               xray_structure,
               mask_params=None,
               d_min = None,
               atom_radius = None):
    """
    If d_min is not None then it defines the gridding. Otherwise mask_params is
    used for that.
    """
    adopt_init_args(self, locals())
    if(mask_params is None):
      mask_params = mmtbx.masks.mask_master_params.extract()
      self.mask_params = mask_params
    validate_mask_params(mask_params, d_min_given = d_min is not None)
    if self.mask_params.ignore_hydrogens:
      hd_sel = self.xray_structure.hd_selection()
      if hd_sel.count(True) > 0:
        self.xray_structure = self.xray_structure.select(~hd_sel)
    if(atom_radius is None):
      self.atom_radii = vdw_radii_from_xray_structure(xray_structure =
        self.xray_structure)
    else:
      self.atom_radii = flex.double(self.xray_structure.scatterers().size(),
        atom_radius)
    self.n_real = None
    if(d_min is None):
      self.n_real = gridding_n_real(
        mask_params    = mask_params,
        xray_structure = self.xray_structure)
    if(d_min is not None):
      self.asu_mask = atom_mask(
        unit_cell                = self.xray_structure.unit_cell(),
        group                    = self.xray_structure.space_group(),
        resolution               = self.d_min,
        grid_step_factor         = self.mask_params.grid_step_factor,
        solvent_radius           = self.mask_params.solvent_radius,
        shrink_truncation_radius = self.mask_params.shrink_truncation_radius)
    else:
      self.asu_mask = atom_mask(
        unit_cell                = self.xray_structure.unit_cell(),
        space_group              = self.xray_structure.space_group(),
        gridding_n_real          = self.n_real,
        solvent_radius           = self.mask_params.solvent_radius,
        shrink_truncation_radius = self.mask_params.shrink_truncation_radius)
    sites_frac = self.xray_structure.sites_frac()
    if(self.mask_params.n_radial_shells > 1):
      # number of shell radii is one less than number of shells
      # last shell is of unknown radius
      shell_rads = [self.mask_params.radial_shell_width] * \
        (self.mask_params.n_radial_shells-1)
      # TODO: Should first shell width be:
      shell_rads[0] -= self.mask_params.solvent_radius
      if( shell_rads[0]<0. ):
        shell_rads[0] = 0.
      self.asu_mask.compute(sites_frac, self.atom_radii, shell_rads)
    else:
      self.asu_mask.compute(sites_frac, self.atom_radii)

  def mask_data_whole_uc(self):
    mask_data = self.asu_mask.mask_data_whole_uc() / \
      self.xray_structure.space_group().order_z()
    return maptbx.copy(mask_data, flex.grid(mask_data.focus()))

class manager(object):
  def __init__(self, miller_array,
                     xray_structure = None,
                     miller_array_twin = None,
                     mask_params = None,
                     compute_mask = True):
    adopt_init_args(self, locals())
    if(self.mask_params is not None): self.mask_params = mask_params
    else: self.mask_params = mask_master_params.extract()
    validate_mask_params(self.mask_params)
    Fmask_res_high = self.mask_params.Fmask_res_high
    self.sel_Fmask_res = miller_array.d_spacings().data() >= Fmask_res_high
    self.sel_Fmask_res_twin = None
    if(miller_array_twin is not None):
       self.sel_Fmask_res_twin = miller_array_twin.d_spacings().data() >= Fmask_res_high
    self.xray_structure = self._load_xrs(xray_structure)
    self.atom_radii = None
    self._f_mask = None
    self._f_mask_twin = None
    self.solvent_content_via_mask = None
    self.layer_volume_fractions = None
    self.sites_cart = None
    if(self.xray_structure is not None):
      self.atom_radii = vdw_radii_from_xray_structure(xray_structure =
       self.xray_structure)
      self.sites_cart = self.xray_structure.sites_cart()
      twin=False
      if(self.miller_array_twin is not None): twin=True
      if(compute_mask): self.compute_f_mask()
    if( not (self._f_mask is None) ):
      assert self._f_mask[0].data().size() == self.miller_array.indices().size()

  def _load_xrs(self, xray_structure):
    if(xray_structure is None): return None
    selection = flex.bool(xray_structure.scatterers().size(), True)
    if(self.mask_params.ignore_zero_occupancy_atoms):
      selection &= xray_structure.scatterers().extract_occupancies() > 0
    if(self.mask_params.ignore_hydrogens):
      selection &= (~xray_structure.hd_selection())
    return xray_structure.select(selection)

  def deep_copy(self):
    return self.select(flex.bool(self.miller_array.indices().size(),True))

  def select(self, selection):
    miller_array_twin = None
    if(self.miller_array_twin is not None):
      miller_array_twin = self.miller_array_twin.select(selection)
    new_manager = manager(
      miller_array      = self.miller_array.select(selection),
      miller_array_twin = miller_array_twin,
      xray_structure    = self.xray_structure,
      mask_params       = deepcopy(self.mask_params),
      compute_mask      = False)
    if(self._f_mask is not None):
      assert self._f_mask[0].data().size() == self.miller_array.indices().size()
      new_manager._f_mask = []
      for fm in self._f_mask:
        assert fm.data().size() == selection.size(), (fm.data().size(), "!=", selection.size())
        new_manager._f_mask.append(fm.select(selection = selection))
    if(self._f_mask_twin is not None):
      new_manager._f_mask_twin = []
      for fm in self._f_mask_twin:
        new_manager._f_mask_twin.append( fm.select(selection = selection) )
    new_manager.solvent_content_via_mask = self.solvent_content_via_mask
    return new_manager

  def shell_f_masks(self, xray_structure = None, force_update=False):
    if(xray_structure is not None):
      if(force_update or self._f_mask is None):
        self.xray_structure = self._load_xrs(xray_structure)
        self.sites_cart = self.xray_structure.sites_cart()
        self.compute_f_mask()
      else:
        xray_structure = self._load_xrs(xray_structure)
        flag = self._need_update_mask(sites_cart_new =
          xray_structure.sites_cart())
        if(flag):
          self.xray_structure = xray_structure
          self.sites_cart = self.xray_structure.sites_cart()
          self.compute_f_mask()
    elif(self._f_mask is None and self.xray_structure is not None):
      self.compute_f_mask()
    # TODO: should return self._f_mask[:], to avoid modification in upper levels?
    return self._f_mask

  # TODO: this seems to be unused
  def shell_f_masks_twin(self):
    return self._f_mask_twin

  def _need_update_mask(self, sites_cart_new):
    if(self.sites_cart is not None and
       self.sites_cart.size() != sites_cart_new.size()): return True
    if(self.sites_cart is not None):
      atom_atom_distances = flex.sqrt((sites_cart_new - self.sites_cart).dot())
      mean_shift = flex.mean_default(atom_atom_distances,0)
      if(mean_shift > self.mask_params.mean_shift_for_mask_update):
        return True
      else: return False
    else: return True

  def compute_f_mask(self):
    d_min = None
    if(self.mask_params.use_resolution_based_gridding):
      d_min = self.miller_array.d_min()
      assert self.mask_params.grid_step_factor is not None
    if(not self.mask_params.use_asu_masks):
      assert self.mask_params.n_radial_shells <= 1
      mask_obj = self.bulk_solvent_mask(d_min = d_min)
      self._f_mask = self._compute_f(
        mask_obj=mask_obj, ma=self.miller_array, sel=self.sel_Fmask_res)
      if(self.miller_array_twin is not None):
        self._f_mask_twin = self._compute_f(
          mask_obj=mask_obj, ma=self.miller_array_twin, sel=self.sel_Fmask_res_twin)
      self.solvent_content_via_mask = mask_obj \
        .contact_surface_fraction
    else:
      global number_of_mask_calculations
      number_of_mask_calculations += 1
      asu_mask_obj = asu_mask(
        xray_structure = self.xray_structure,
        d_min          = d_min,
        mask_params    = self.mask_params).asu_mask
      self._f_mask = self._compute_f(mask_obj=asu_mask_obj,
        ma=self.miller_array, sel=self.sel_Fmask_res)
      if(self.miller_array_twin is not None):
        assert self.miller_array.indices().size() == \
               self.miller_array_twin.indices().size()
        self._f_mask_twin = self._compute_f(mask_obj=asu_mask_obj,
          ma=self.miller_array_twin, sel=self.sel_Fmask_res_twin)
      self.solvent_content_via_mask = asu_mask_obj.contact_surface_fraction
      self.layer_volume_fractions = asu_mask_obj.layer_volume_fractions()
    assert self._f_mask[0].data().size() == self.miller_array.data().size()

  def _compute_f(self, mask_obj, ma, sel):
    result = []
    ma_sel = ma.select(sel)
    for i in range(0, self.mask_params.n_radial_shells):
      if(self.mask_params.use_asu_masks):
        fm_asu = mask_obj.structure_factors(ma_sel.indices(),i+1)
      else:
        assert self.mask_params.n_radial_shells==1
        # The mask covers the whole unit cell: no symmetry summation. sel
        # (Fmask_res_high) is the only resolution cutoff.
        ma_sel_p1 = ma_sel.customized_copy(
          crystal_symmetry = mask_obj.xray_structure.crystal_symmetry())
        fm_asu = mask_obj.structure_factors(ma_sel_p1, zero_high=False).data()
      data = flex.complex_double(ma.indices().size(), 0)
      data = data.set_selected(sel, fm_asu)
      result.append( ma.set().array(data = data) )
    return result

  def bulk_solvent_mask(self, d_min=None):
    """
    Mask in the whole unit cell, on the same grid as the asu mask. d_min is
    given for resolution-based gridding only.

    mmtbx.masks.bulk_solvent marks the grid around the given sites only, so
    it gets the structure expanded to P1.
    """
    mp = self.mask_params
    n_real = gridding_n_real(mask_params = mp,
      xray_structure = self.xray_structure, d_min = d_min)
    return bulk_solvent(
      xray_structure = self.xray_structure.expand_to_p1(
        sites_mod_positive = True),
      gridding_n_real             = n_real,
      ignore_zero_occupancy_atoms = mp.ignore_zero_occupancy_atoms,
      ignore_hydrogen_atoms       = mp.ignore_hydrogens,
      solvent_radius              = mp.solvent_radius,
      shrink_truncation_radius    = mp.shrink_truncation_radius)

class mask_from_xray_structure(object):
  def __init__(
        self,
        xray_structure,
        p1,
        for_structure_factors,
        n_real,
        atom_radii=None,
        atom_radius=None,
        solvent_radius=None,
        shrink_truncation_radius=None,
        in_asu=False,
        rad_extra=None):
    if([solvent_radius, shrink_truncation_radius].count(None)>0):
      mask_params = mask_master_params.extract()
      if(solvent_radius is None): solvent_radius = mask_params.solvent_radius
      if(shrink_truncation_radius is None):
        shrink_truncation_radius = mask_params.shrink_truncation_radius
    xrs = xray_structure
    sgt = xrs.space_group().type() # must be BEFORE going to P1, obviously!
    if(p1): xrs = xrs.expand_to_p1(sites_mod_positive=True)
    if([atom_radii,atom_radius].count(None)==2):
      atom_radii = vdw_radii_from_xray_structure(xray_structure = xrs)
    elif atom_radius is not None:
      assert atom_radii is None
      atom_radii = flex.double(xrs.scatterers().size(), atom_radius)
    if(rad_extra is not None):
      atom_radii = atom_radii+rad_extra
    self.asu_mask = mmtbx.masks.atom_mask(
      unit_cell                = xrs.unit_cell(),
      space_group              = xrs.space_group(),
      gridding_n_real          = n_real,
      solvent_radius           = solvent_radius,
      shrink_truncation_radius = shrink_truncation_radius)
    self.asu_mask.compute(xrs.sites_frac(), atom_radii)
    self.mask_data = self.asu_mask.mask_data_whole_uc()
    if(for_structure_factors):
      self.mask_data = self.mask_data / xrs.space_group().order_z()
    if(in_asu):
      assert p1
      asu_map_ = asu_map_ext.asymmetric_map(sgt, self.mask_data)
      self.mask_data = asu_map_.data()

class smooth_mask(object):
  def __init__(self, xray_structure, n_real, rad_smooth, atom_radii=None,
                     solvent_radius=None, shrink_truncation_radius=None):
    mask_binary = mask_from_xray_structure(
      xray_structure           = xray_structure,
      p1                       = True,
      for_structure_factors    = False,
      n_real                   = n_real,
      atom_radii               = atom_radii,
      solvent_radius           = solvent_radius,
      shrink_truncation_radius = shrink_truncation_radius,
      rad_extra                = rad_smooth).mask_data
    s = mask_binary>0.5
    mask_binary = mask_binary.set_selected(s, 0)
    mask_binary = mask_binary.set_selected(~s, 1)
    maptbx.unpad_in_place(map=mask_binary)
    self.mask_smooth = maptbx.smooth_map(
      map              = mask_binary,
      crystal_symmetry = xray_structure.crystal_symmetry(),
      rad_smooth       = rad_smooth)
