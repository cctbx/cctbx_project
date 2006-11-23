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

mask_master_params = iotbx.phil.parse("""\
  solvent_radius = 1.11
    .type = float
  shrink_truncation_radius = 0.9
    .type = float
  grid_step_factor = 4.0
    .type = float
  verbose = 1
    .type = int
  mean_shift_for_mask_update = 0.1
    .type = float
""")

class bulk_solvent(around_atoms):

  def __init__(self,
        xray_structure,
        gridding_n_real=None,
        grid_step =None,
        solvent_radius=1.0,
        shrink_truncation_radius=1.0):
     assert [gridding_n_real, grid_step].count(None) == 1
     self.xray_structure = xray_structure
     self.grid_step = grid_step
     if(self.grid_step is not None): self.grid_step = min(0.8, self.grid_step)
     if(gridding_n_real is None):
       gridding_n_real = maptbx.crystal_gridding(
         unit_cell=xray_structure.unit_cell(),
         step= self.grid_step).n_real()
     atom_radii = flex.double()
     # XXX use scattering dictionary and set_selected
     for scatterer in xray_structure.scatterers():
       try:
         atom_radii.append(
           van_der_waals_radii.vdw.table[scatterer.element_symbol()])
       except Exception, e:
         raise RuntimeError,\
               "scatterer.element_symbol()= %s"%str(scatterer.element_symbol())
     around_atoms.__init__(self,
       unit_cell=xray_structure.unit_cell(),
       space_group_order_z=xray_structure.space_group().order_z(),
       sites_frac=xray_structure.sites_frac(),
       atom_radii=atom_radii,
       gridding_n_real=gridding_n_real,
       solvent_radius=solvent_radius,
       shrink_truncation_radius=shrink_truncation_radius)
     introspection.virtual_memory_info().update_max()

  def show_summary(self, out=None):
    if (out is None): out = sys.stdout
    print >> out, "|-mask parameters----------------------------------------"\
                  "---------------------|"
    print >> out, "| solvent radius        = %4.2f A           shrink truncat"\
                  "ion radius = %4.2f A  |"%(self.solvent_radius,\
                  self.shrink_truncation_radius)
    print >> out, "| solvent content       = %.1f %%           grid_step     "\
                  "           = %5.3f A | "%(self.contact_surface_fraction*100,
                  self.grid_step)
    part_1 = "| number of grid points = %s"%(str(self.data.accessor().focus()))
    n = 78 - len(part_1)
    print >> out, part_1 + " "*n +"|"
    print >> out, "|"+"-"*77+"|"
    print >> out

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
    anomalous_flag = False
    conjugate_flag = True
    from_map = maptbx.structure_factors.from_map(
      miller_set.space_group(),
      anomalous_flag,
      miller_set.indices(),
      map_of_coeff,
      conjugate_flag)
    return miller_set.array(data=from_map.data())

  def subtract_non_uniform_solvent_region_in_place(self, non_uniform_mask):
    assert non_uniform_mask.accessor() == self.data.accessor()
    self.data.set_selected(non_uniform_mask > 0, 0)
