from __future__ import absolute_import, division, print_function
from libtbx.math_utils import ifloor, iceil
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
from cctbx import crystal
from cctbx import maptbx
from cctbx import uctbx
from cctbx import xray

class around_model(object):
  """
  Extract map box around model. Model can be either xray_structure or
  (pdb_hierarchy and crystal_symmetry) or both (xray_structure and pdb_hierarchy).
  """
  def __init__(self,
               map_data,
               cushion,
               xray_structure=None,
               pdb_hierarchy=None,
               crystal_symmetry=None):
    # safeguards
    if(xray_structure is None):
      assert [pdb_hierarchy, crystal_symmetry].count(None) == 0
      cs = crystal_symmetry
      uc = cs.unit_cell()
      sites_frac = uc.fractionalize(pdb_hierarchy.atoms().extract_xyz())
    else:
      if(crystal_symmetry is not None):
        assert crystal_symmetry.is_similar_symmetry(
          xray_structure.crystal_symmetry())
      if(pdb_hierarchy is not None):
        assert approx_equal(xray_structure.sites_cart(),
          pdb_hierarchy.atoms().extract_xyz())
      cs = xray_structure.crystal_symmetry()
      uc = cs.unit_cell()
      sites_frac = xray_structure.sites_frac()
    assert cushion >= 0
    # convert cushion into fractional vector
    cushion_frac = flex.double(uc.fractionalize((cushion,)*3))
    # find fractional corners
    frac_min = sites_frac.min()
    frac_max = sites_frac.max()
    frac_max = list(flex.double(frac_max)+cushion_frac)
    frac_min = list(flex.double(frac_min)-cushion_frac)
    # find corner grid nodes
    na = map_data.all()
    self.gridding_first = [ifloor(f*n) for f,n in zip(frac_min, na)]
    self.gridding_last  = [ iceil(f*n) for f,n in zip(frac_max, na)]
    # extract map box
    self.map_data = maptbx.copy(map_data,self.gridding_first,self.gridding_last)
    # get shift vector as result of boxing
    self.shift_frac = [-self.map_data.origin()[i]/na[i] for i in range(3)]
    self.shift_cart = cs.unit_cell().orthogonalize(self.shift_frac)
    # this makes the box 'forget' about its parent -- the original box
    self.map_data.reshape(flex.grid(self.map_data.all()))
    # get crystal symmetry of the box
    p = uc.parameters()
    abc = [p[i] * self.map_data.all()[i]/na[i] for i in range(3)]
    box_uc = uctbx.unit_cell(parameters=(abc[0],abc[1],abc[2],p[3],p[4],p[5]))
    self.crystal_symmetry = crystal.symmetry(
      unit_cell = box_uc, space_group="P1")
    # new xray_structure or pdb_hierarchy
    sites_frac_new = box_uc.fractionalize(
      uc.orthogonalize(sites_frac+self.shift_frac))
    sites_cart_new = box_uc.orthogonalize(sites_frac_new)
    self.xray_structure, self.pdb_hierarchy = None,None
    if(xray_structure is not None):
      scatterers = xray_structure.scatterers().deep_copy()
      scatterers.set_sites(sites_frac_new)
      sp = crystal.special_position_settings(self.crystal_symmetry)
      self.xray_structure = xray.structure(sp, scatterers)
    if(pdb_hierarchy is not None):
      self.pdb_hierarchy = pdb_hierarchy.deep_copy()
      self.pdb_hierarchy.atoms().set_xyz(sites_cart_new)
