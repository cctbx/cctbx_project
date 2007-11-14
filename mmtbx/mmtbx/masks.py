import boost.python
ext = boost.python.import_ext("mmtbx_masks_ext")
from mmtbx_masks_ext import *

from cctbx import maptbx
from cctbx.eltbx import van_der_waals_radii
from cctbx.array_family import flex
from scitbx import fftpack
from scitbx import matrix
import sys
import iotbx.xplor.map
import iotbx.phil
from libtbx import introspection
from libtbx import adopt_init_args
from copy import deepcopy

number_of_mask_calculations = 0

mask_master_params = iotbx.phil.parse("""\
  solvent_radius = 1.11
    .type = float
  shrink_truncation_radius = 0.9
    .type = float
  grid_step_factor = 4.0
    .type = float
    .help = The grid step for the mask calculation is determined as \
            highest_resolution devided by grid_step_factor. This is considered \
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
""")

class bulk_solvent(around_atoms):

  def __init__(self,
        xray_structure,
        ignore_zero_occupancy_atoms,
        solvent_radius,
        shrink_truncation_radius,
        gridding_n_real=None,
        grid_step=None):
     global number_of_mask_calculations
     number_of_mask_calculations += 1
     assert [gridding_n_real, grid_step].count(None) == 1
     self.xray_structure = xray_structure
     if (gridding_n_real is None):
       gridding_n_real = maptbx.crystal_gridding(
         unit_cell=xray_structure.unit_cell(),
         step=grid_step).n_real()
     atom_radii = flex.double()
     # XXX use scattering dictionary and set_selected
     # XXX use monomer library definitions for radii
     unknown = []
     for i_seq, scatterer in enumerate(xray_structure.scatterers()):
       try:
         atom_radii.append(
           van_der_waals_radii.vdw.table[scatterer.element_symbol()])
       except:
         unknown.append(scatterer.element_symbol())
     sites_frac = xray_structure.sites_frac()
     if(len(unknown) > 0):
        raise RuntimeError("Atoms with unknown van der Waals radius: ",unknown)
     self.n_atoms_excluded = 0
     if(ignore_zero_occupancy_atoms):
       selection = xray_structure.scatterers().extract_occupancies() > 0
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
    fft_manager = fftpack.real_to_complex_3d(self.data.focus())
    padded_data = maptbx.copy(
      self.data.as_double(),
      flex.grid(fft_manager.m_real()).set_focus(fft_manager.n_real()))
    map_of_coeff = fft_manager.forward(padded_data)
    scale = miller_set.unit_cell().volume() \
          / matrix.col(fft_manager.n_real()).product()
    map_of_coeff *= scale # XXX scale from_map.data() instead
    from_map = maptbx.structure_factors.from_map(
      space_group=miller_set.space_group(),
      anomalous_flag=False,
      miller_indices=miller_set.indices(),
      complex_map=map_of_coeff,
      conjugate_flag=True)
    return miller_set.array(data=from_map.data())

  def subtract_non_uniform_solvent_region_in_place(self, non_uniform_mask):
    assert non_uniform_mask.accessor() == self.data.accessor()
    self.data.set_selected(non_uniform_mask > 0, 0)

class manager(object):
  def __init__(self, miller_array,
                     xray_structure,
                     mask_params = None):
    adopt_init_args(self, locals())
    if(self.mask_params is not None): self.mask_params = mask_params
    else: self.mask_params = mask_master_params.extract()
    self.grid_step = self._get_grid_step()
    if(xray_structure is not None):
      self.xray_structure = self.xray_structure.deep_copy_scatterers()
      self.sites_cart = self.xray_structure.sites_cart()
      self._f_mask = self.compute_f_mask()
    else: self._f_mask = None

  def deep_copy(self):
    new_manager = manager(miller_array   = self.miller_array.deep_copy(),
                          xray_structure = None,
                          mask_params    = deepcopy(self.mask_params))
    if(self.xray_structure is not None):
      new_manager.xray_structure = self.xray_structure.deep_copy_scatterers()
      new_manager.sites_cart     = new_manager.xray_structure.sites_cart()
    if(self._f_mask is not None):
      new_manager._f_mask = self._f_mask.deep_copy()
    return new_manager

  def select(self, selection):
    new_manager = self.deep_copy()
    new_manager.miller_array = self.miller_array.select(selection = selection)
    if(self._f_mask is not None):
      new_manager._f_mask = self._f_mask.select(selection = selection)
    return new_manager

  def _get_grid_step(self):
    assert self.mask_params.grid_step_factor > 0
    step = self.miller_array.d_min()/self.mask_params.grid_step_factor
    if(step < 0.15): step = 0.15
    step = min(0.8, step)
    return step

  def f_mask(self, xray_structure_new = None, force_update = False):
    if(xray_structure_new is None): return self._f_mask
    else:
      if(force_update or self._f_mask is None):
        self.xray_structure = xray_structure_new.deep_copy_scatterers()
        self.sites_cart = xray_structure_new.sites_cart()
        return self.compute_f_mask()
      else:
        flag = self._need_update_mask(sites_cart_new =
          xray_structure_new.sites_cart())
        if(flag):
          self.xray_structure = xray_structure_new.deep_copy_scatterers()
          self.sites_cart = xray_structure_new.sites_cart()
          return self.compute_f_mask()
        else:
          return self._f_mask

  def _need_update_mask(self, sites_cart_new):
    if(self.sites_cart.size() != sites_cart_new.size()): return True
    atom_atom_distances = flex.sqrt((sites_cart_new - self.sites_cart).dot())
    mean_shift = flex.mean(atom_atom_distances)
    if(mean_shift > self.mask_params.mean_shift_for_mask_update):
      return True
    else: return False

  def compute_f_mask(self):
    bulk_solvent_mask_obj = self.bulk_solvent_mask()
    self._f_mask = bulk_solvent_mask_obj.structure_factors(
      miller_set = self.miller_array)
    return self._f_mask

  def bulk_solvent_mask(self):
    mp = self.mask_params
    return bulk_solvent(
      xray_structure              = self.xray_structure,
      grid_step                   = self.grid_step,
      ignore_zero_occupancy_atoms = mp.ignore_zero_occupancy_atoms,
      solvent_radius              = mp.solvent_radius,
      shrink_truncation_radius    = mp.shrink_truncation_radius)
