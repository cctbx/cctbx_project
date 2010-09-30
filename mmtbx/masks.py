import boost.python
from cctbx.array_family import flex
ext = boost.python.import_ext("mmtbx_masks_ext")
from mmtbx_masks_ext import *

from cctbx.masks import around_atoms, vdw_radii_from_xray_structure
from cctbx import maptbx
from scitbx import fftpack
from scitbx import matrix
import sys
import iotbx.xplor.map
import iotbx.phil
from libtbx.utils import Sorry
from libtbx import introspection
from libtbx import adopt_init_args
from copy import deepcopy
import masks

number_of_mask_calculations = 0

mask_master_params = iotbx.phil.parse("""\
  use_asu_masks = True
    .type = bool
  solvent_radius = 1.11
    .type = float
  shrink_truncation_radius = 0.9
    .type = float
  grid_step_factor = 4.0
    .type = float
    .help = The grid step for the mask calculation is determined as \
            highest_resolution divided by grid_step_factor. This is considered \
            as suggested value and may be adjusted internally based on the \
            resolution.
  verbose = 1
    .type = int
    .expert_level=3
  mean_shift_for_mask_update = 0.1
    .type = float
    .help = Value of overall model shift in refinement to updates the mask.
    .expert_level=2
  ignore_zero_occupancy_atoms = True
    .type = bool
    .help = Include atoms with zero occupancy into mask calculation
  ignore_hydrogens = False
    .type = bool
    .help = Ignore H or D atoms in mask calculation
""")

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
    print >> out, "|- mask summary -----------------------------------------"\
                  "---------------------|"
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
      print >> out, "| " + s + " " * max(0, 75-len(s)) + " |"
    gap = 6
    show_line(s=sr+" "*gap+gr)
    show_line(s=st+" "*gap+gs)
    show_line(s=na+" "*gap+sc)
    print >> out, "|"+"-"*77+"|"

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

  def structure_factors(self, miller_set):
    result = miller_set.structure_factors_from_map(
      map = self.data, use_scale = True, anomalous_flag = False, use_sg = True)
    return miller_set.array(data=result.data())

  def subtract_non_uniform_solvent_region_in_place(self, non_uniform_mask):
    assert non_uniform_mask.accessor() == self.data.accessor()
    self.data.set_selected(non_uniform_mask > 0, 0)

class manager(object):
  def __init__(self, miller_array,
                     xray_structure,
                     miller_array_twin = None,
                     mask_params = None):
    adopt_init_args(self, locals())
    if(self.mask_params is not None): self.mask_params = mask_params
    else: self.mask_params = mask_master_params.extract()
    self.grid_step = self._get_grid_step()
    self.atom_radii = None
    self._f_mask = None
    self._f_mask_twin = None
    self.solvent_content_via_mask = None
    self.sites_cart = None
    if(xray_structure is not None):
      self.atom_radii = vdw_radii_from_xray_structure(xray_structure =
       self.xray_structure)
      self.xray_structure = self.xray_structure.deep_copy_scatterers()
      self.sites_cart = self.xray_structure.sites_cart()
      twin=False
      if(self.miller_array_twin is not None): twin=True
      self.compute_f_mask(twin = twin)

  def deep_copy(self):
    return self.select(flex.bool(self.miller_array.indices().size(),True))

  def select(self, selection):
    miller_array_twin = None
    if(self.miller_array_twin is not None):
      miller_array_twin = self.miller_array_twin.select(selection)
    new_manager = manager(
      miller_array      = self.miller_array.select(selection),
      miller_array_twin = miller_array_twin,
      xray_structure    = None,
      mask_params       = deepcopy(self.mask_params))
    if(self._f_mask is not None):
      new_manager._f_mask = self._f_mask.select(selection = selection)
    if(self._f_mask_twin is not None):
      new_manager._f_mask_twin = self._f_mask_twin.select(selection = selection)
    new_manager.solvent_content_via_mask = self.solvent_content_via_mask
    return new_manager

  def _get_grid_step(self):
    if not (self.mask_params.grid_step_factor > 0) :
      raise Sorry("Inappropriate value for grid_step_factor: must be "+
        "positive and non-zero.")
    assert self.mask_params.grid_step_factor > 0
    step = self.miller_array.d_min()/self.mask_params.grid_step_factor
    if(step < 0.15): step = 0.15
    step = min(0.8, step)
    return step

  def f_mask(self, xray_structure_new = None, force_update = False,
             twin = False):
    if(twin): f_mask = self._f_mask_twin
    else: f_mask = self._f_mask
    if(xray_structure_new is None): return f_mask
    else:
      if(force_update or f_mask is None):
        self.xray_structure = xray_structure_new.deep_copy_scatterers()
        self.sites_cart = xray_structure_new.sites_cart()
        return self.compute_f_mask()
      else:
        flag = self._need_update_mask(sites_cart_new =
          xray_structure_new.sites_cart())
        if(flag):
          self.xray_structure = xray_structure_new.deep_copy_scatterers()
          self.sites_cart = xray_structure_new.sites_cart()
          return self.compute_f_mask(twin=twin)
        else:
          return f_mask

  def _need_update_mask(self, sites_cart_new):
    if(self.sites_cart is not None and
       self.sites_cart.size() != sites_cart_new.size()): return True
    if(self.sites_cart is not None):
      atom_atom_distances = flex.sqrt((sites_cart_new - self.sites_cart).dot())
      mean_shift = flex.mean(atom_atom_distances)
      if(mean_shift > self.mask_params.mean_shift_for_mask_update):
        return True
      else: return False
    else: return True

  def compute_f_mask(self, twin=False):
    if(not self.mask_params.use_asu_masks):
      bulk_solvent_mask_obj = self.bulk_solvent_mask()
      self._f_mask = bulk_solvent_mask_obj.structure_factors(
        miller_set = self.miller_array)
      if(self.miller_array_twin is not None):
        self._f_mask_twin = bulk_solvent_mask_obj.structure_factors(
          miller_set = self.miller_array_twin)
      self.solvent_content_via_mask = bulk_solvent_mask_obj \
        .contact_surface_fraction
    else:
      self.atom_radii = vdw_radii_from_xray_structure(xray_structure =
        self.xray_structure)
      asu_mask = masks.atom_mask(
        unit_cell                = self.xray_structure.unit_cell(),
        group                    = self.xray_structure.space_group(),
        resolution               = self.miller_array.d_min(),
        grid_step_factor         = self.mask_params.grid_step_factor,
        solvent_radius           = self.mask_params.solvent_radius,
        shrink_truncation_radius = self.mask_params.shrink_truncation_radius)
      # XXX duplication from old mask calculation code, see above
      selection = flex.bool(self.xray_structure.scatterers().size(), True)
      if(self.mask_params.ignore_zero_occupancy_atoms):
        selection &= self.xray_structure.scatterers().extract_occupancies() > 0
      if(self.mask_params.ignore_hydrogens): # it is very essential
        selection &= (~self.xray_structure.hd_selection())
      sites_frac = self.xray_structure.sites_frac()
      sites_frac = sites_frac.select(selection)
      atom_radii = self.atom_radii.select(selection)
      #
      asu_mask.compute(sites_frac, atom_radii)
      fm_asu = asu_mask.structure_factors(self.miller_array.indices())
      self._f_mask = self.miller_array.set().array(data = fm_asu)
      if(self.miller_array_twin is not None):
        assert self.miller_array.indices().size() == \
               self.miller_array_twin.indices().size()
        fm_asu = asu_mask.structure_factors(self.miller_array_twin.indices())
        self._f_mask_twin = self.miller_array_twin.set().array(data = fm_asu)
      self.solvent_content_via_mask = asu_mask.contact_surface_fraction
    if(twin): return self._f_mask_twin
    else: return self._f_mask

  def bulk_solvent_mask(self):
    mp = self.mask_params
    return bulk_solvent(
      xray_structure              = self.xray_structure,
      grid_step                   = self.grid_step,
      ignore_zero_occupancy_atoms = mp.ignore_zero_occupancy_atoms,
      ignore_hydrogen_atoms       = mp.ignore_hydrogens,
      solvent_radius              = mp.solvent_radius,
      shrink_truncation_radius    = mp.shrink_truncation_radius)
