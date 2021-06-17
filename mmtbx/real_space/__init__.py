from __future__ import absolute_import, division, print_function
import cctbx.array_family.flex
from scitbx.array_family import flex
from cctbx.xray import ext
from libtbx.test_utils import approx_equal

class sampled_model_density(object):

  def __init__(self, xray_structure,
                     n_real    = None,
                     grid_step = None):
    assert [n_real, grid_step].count(None) == 1
    self.xray_structure = xray_structure
    # XXX Use crystal gridding
    if(n_real is None):
       a_cell, b_cell, c_cell = self.xray_structure.unit_cell().parameters()[:3]
       nx,ny,nz = \
              int(a_cell/grid_step),int(b_cell/grid_step),int(c_cell/grid_step)
    else:
       nx,ny,nz = n_real[0],n_real[1],n_real[2]
    #
    occs = self.xray_structure.scatterers().extract_occupancies()
    sel = occs>0
    self.xray_structure = self.xray_structure.select(sel)
    #
    self.xray_structure.tidy_us()
    us = flex.min_default(self.xray_structure.extract_u_iso_or_u_equiv(),0)
    u_base = self.xray_structure.min_u_cart_eigenvalue()
    u_base = min(us, u_base)
    self.sampled_model_density = ext.sampled_model_density(
       unit_cell                        = self.xray_structure.unit_cell(),
       scatterers                       = self.xray_structure.scatterers(),
       scattering_type_registry =self.xray_structure.scattering_type_registry(),
       fft_n_real                       = (nx,ny,nz),
       fft_m_real                       = (nx,ny,nz),
       u_base                           = u_base,
       wing_cutoff                      = 1.e-6,
       exp_table_one_over_step_size     = -1000,
       force_complex                    = False,
       sampled_density_must_be_positive = False,
       tolerance_positive_definite      = 1.e-5)
    assert approx_equal(self.sampled_model_density.u_extra(), 0.0)

  def data(self):
    return self.sampled_model_density.real_map_unpadded()

  def write_as_xplor_map(self, file_name):
    import iotbx.xplor.map
    model_map = self.data()
    gridding = iotbx.xplor.map.gridding(n     = model_map.focus(),
                                        first = (0,0,0),
                                        last  = model_map.focus())
    iotbx.xplor.map.writer(file_name   = file_name,
                           is_p1_cell  = True,
                           title_lines = [' REMARKS FILENAME=""',
                                          ' REMARKS '],
                           unit_cell   = self.xray_structure.unit_cell(),
                           gridding    = gridding,
                           data        = model_map,
                           average     = -1,
                           standard_deviation=-1)

class around_atom(object):

  def __init__(self, unit_cell,
                     map_data,
                     radius,
                     shell,
                     site_frac):
    from cctbx import maptbx
    obj = maptbx.grid_points_in_sphere_around_atom_and_distances(
                                                         unit_cell = unit_cell,
                                                         data      = map_data,
                                                         radius    = radius,
                                                         shell     = shell,
                                                         site_frac = site_frac)
    data = obj.data_at_grid_points()
    dist = obj.distances()
    p = flex.sort_permutation(dist)
    self.data_ = data.select(p)
    self.dist_ = dist.select(p)

  def data(self):
    return self.data_

  def distances(self):
    return self.dist_
