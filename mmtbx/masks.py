import boost.python
from cctbx.array_family import flex
ext = boost.python.import_ext("mmtbx_masks_ext")
from mmtbx_masks_ext import *

from cctbx.masks import around_atoms, vdw_radii_from_xray_structure
from cctbx import maptbx
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
  n_radial_shells = 1
    .type = int
    .help = Number of shells in a radial shell bulk solvent model
  radial_shell_width = 1.3
    .type = float
    .help = Radial shell width
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

class asu_mask(object):
  def __init__(self, xray_structure, d_min=None, n_real=None, mask_params = None):
    adopt_init_args(self, locals())
    assert [d_min, n_real].count(None) == 1
    if(self.mask_params is None):
      self.mask_params = mask_master_params.extract()
    self.atom_radii = vdw_radii_from_xray_structure(xray_structure =
      self.xray_structure)
    if(d_min is not None):
      self.asu_mask = masks.atom_mask(
        unit_cell                = self.xray_structure.unit_cell(),
        group                    = self.xray_structure.space_group(),
        resolution               = self.d_min,
        grid_step_factor         = self.mask_params.grid_step_factor,
        solvent_radius           = self.mask_params.solvent_radius,
        shrink_truncation_radius = self.mask_params.shrink_truncation_radius)
    else:
      self.asu_mask = masks.atom_mask(
        unit_cell                = self.xray_structure.unit_cell(),
        space_group              = self.xray_structure.space_group(),
        gridding_n_real          = self.n_real,
        solvent_radius           = self.mask_params.solvent_radius,
        shrink_truncation_radius = self.mask_params.shrink_truncation_radius)
    selection = flex.bool(self.xray_structure.scatterers().size(), True)
    if(self.mask_params.ignore_zero_occupancy_atoms):
      selection &= self.xray_structure.scatterers().extract_occupancies() > 0
    if(self.mask_params.ignore_hydrogens):
      selection &= (~self.xray_structure.hd_selection())
    sites_frac = self.xray_structure.sites_frac()
    sites_frac = sites_frac.select(selection)
    atom_radii = self.atom_radii.select(selection)
    if(self.mask_params.n_radial_shells > 1):
      # number of shell radii is one less than number of shells
      # last shell is of unknown radius
      shell_rads = [self.mask_params.radial_shell_width] * \
        (self.mask_params.n_radial_shells-1)
      # TODO: Should first shell width be:
      shell_rads[0] -= self.mask_params.solvent_radius
      if( shell_rads[0]<0. ):
        shell_rads[0] = 0.
      self.asu_mask.compute(sites_frac, atom_radii, shell_rads)
    else:
      self.asu_mask.compute(sites_frac, atom_radii)

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
    self.grid_step = self._get_grid_step()
    self.atom_radii = None
    self._f_mask = None
    self._f_mask_twin = None
    self.solvent_content_via_mask = None
    self.layer_volume_fractions = None
    self.sites_cart = None
    if(xray_structure is not None):
      self.atom_radii = vdw_radii_from_xray_structure(xray_structure =
       self.xray_structure)
      self.xray_structure = self.xray_structure.deep_copy_scatterers()
      self.sites_cart = self.xray_structure.sites_cart()
      twin=False
      if(self.miller_array_twin is not None): twin=True
      if(compute_mask): self.compute_f_mask()
    if( not (self._f_mask is None) ):
      assert self._f_mask[0].data().size() == self.miller_array.indices().size()

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

  def _get_grid_step(self):
    if not (self.mask_params.grid_step_factor > 0) :
      raise Sorry("Inappropriate value for grid_step_factor: must be "+
        "positive and non-zero.")
    assert self.mask_params.grid_step_factor > 0
    step = self.miller_array.d_min()/self.mask_params.grid_step_factor
    if(step < 0.15): step = 0.15 # XXX arbitrary, see also masks/atom_mask.cpp
    step = min(0.8, step) # XXX arbitrary, can lead to very large maps
    return step

  def shell_f_masks(self, xray_structure = None, force_update=False):
    if(xray_structure is not None):
      if(force_update or self._f_mask is None):
        self.xray_structure = xray_structure.deep_copy_scatterers()
        self.sites_cart = xray_structure.sites_cart()
        self.compute_f_mask()
      else:
        flag = self._need_update_mask(sites_cart_new =
          xray_structure.sites_cart())
        if(flag):
          self.xray_structure = xray_structure.deep_copy_scatterers()
          self.sites_cart = xray_structure.sites_cart()
          self.compute_f_mask()
    elif(self._f_mask is None and self.xray_structure is not None):
      self.compute_f_mask()
    # TODO: should return self._f_mask[:], to avoid modification in upper levels?
    return self._f_mask

  def f_mask(self, xray_structure = None, force_update=False):
    assert self.mask_params.n_radial_shells <= 1
    fmsks = self.shell_f_masks(xray_structure, force_update)
    if( fmsks is None ):
      return None
    assert len(fmsks) == 1
    return fmsks[0]

  def shell_f_masks_twin(self):
    return self._f_mask_twin

  def f_mask_twin(self):
    assert self.mask_params.n_radial_shells <= 1
    if( self._f_mask_twin is None ):
      return None
    assert len(self._f_mask_twin) == 1
    return self._f_mask_twin[0]

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
    if(not self.mask_params.use_asu_masks):
      assert self.mask_params.n_radial_shells <= 1
      bulk_solvent_mask_obj = self.bulk_solvent_mask()
      self._f_mask = [bulk_solvent_mask_obj.structure_factors(
        miller_set = self.miller_array)]
      if(self.miller_array_twin is not None):
        self._f_mask_twin = [bulk_solvent_mask_obj.structure_factors(
          miller_set = self.miller_array_twin)]
      self.solvent_content_via_mask = bulk_solvent_mask_obj \
        .contact_surface_fraction
    else:
      asu_mask_obj = asu_mask(
        xray_structure = self.xray_structure,
        d_min          = self.miller_array.d_min(),
        mask_params    = self.mask_params).asu_mask
      self._f_mask = []
      for i in range(0, self.mask_params.n_radial_shells):
        fm_asu = asu_mask_obj.structure_factors(self.miller_array.indices(),i+1)
        self._f_mask.append( self.miller_array.set().array(data = fm_asu) )
      if(self.miller_array_twin is not None):
        assert self.miller_array.indices().size() == \
               self.miller_array_twin.indices().size()
        self._f_mask_twin = []
        for i in range(0,self.mask_params.n_radial_shells):
          fm_asu = asu_mask_obj.structure_factors(self.miller_array_twin.indices())
          self._f_mask_twin.append( self.miller_array_twin.set().array(data = fm_asu) )
      self.solvent_content_via_mask = asu_mask_obj.contact_surface_fraction
      self.layer_volume_fractions = asu_mask_obj.layer_volume_fractions()
    assert self._f_mask[0].data().size() == self.miller_array.data().size()

  def bulk_solvent_mask(self):
    mp = self.mask_params
    return bulk_solvent(
      xray_structure              = self.xray_structure,
      grid_step                   = self.grid_step,
      ignore_zero_occupancy_atoms = mp.ignore_zero_occupancy_atoms,
      ignore_hydrogen_atoms       = mp.ignore_hydrogens,
      solvent_radius              = mp.solvent_radius,
      shrink_truncation_radius    = mp.shrink_truncation_radius)
