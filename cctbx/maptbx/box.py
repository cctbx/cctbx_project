from __future__ import absolute_import, division, print_function
from libtbx.math_utils import ifloor, iceil
from cctbx import crystal, maptbx, uctbx
from scitbx.array_family import flex
from libtbx.utils import null_out
import iotbx.map_manager
import mmtbx.model

class around_model(object):
  """
  Extract map box around atomic model. Box is in P1 and has origin at (0,0,0).
  """
  def __init__(self, map_manager, model, cushion):
    self.map_manager = map_manager
    self.model = model
    # safeguards
    assert isinstance(map_manager, iotbx.map_manager.map_manager)
    assert isinstance(model, mmtbx.model.manager)
    assert self.map_manager.map_data().accessor().origin() == (0,0,0)
    assert map_manager.crystal_symmetry().is_similar_symmetry(
      model.crystal_symmetry())
    assert cushion >= 0
    # get items needed to do the shift
    cs = model.crystal_symmetry()
    uc = cs.unit_cell()
    sites_frac = model.get_sites_frac()
    map_data = map_manager.map_data()
    # convert cushion into fractional vector
    cushion_frac = flex.double(uc.fractionalize((cushion,)*3))
    # find fractional corners
    frac_min = sites_frac.min()
    frac_max = sites_frac.max()
    frac_max = list(flex.double(frac_max)+cushion_frac)
    frac_min = list(flex.double(frac_min)-cushion_frac)
    # find corner grid nodes
    all_orig = map_data.all()
    self.gridding_first = [ifloor(f*n) for f,n in zip(frac_min, all_orig)]
    self.gridding_last  = [ iceil(f*n) for f,n in zip(frac_max, all_orig)]
    #
    self.box_all=[j-i+1 for i,j in zip(self.gridding_first,self.gridding_last)]
    # get shift vector as result of boxing
    self.shift_frac = [-self.gridding_first[i]/all_orig[i] for i in range(3)]
    self.shift_cart = cs.unit_cell().orthogonalize(self.shift_frac)
    # get crystal symmetry of the box
    p = uc.parameters()
    abc = [p[i] * self.box_all[i]/all_orig[i] for i in range(3)]
    box_uc = uctbx.unit_cell(parameters=(abc[0],abc[1],abc[2],p[3],p[4],p[5]))
    self.crystal_symmetry = crystal.symmetry(
      unit_cell = box_uc, space_group="P1")

  def apply_to_map(self, map_manager=None):
    if(map_manager is None):
      map_data = self.map_manager.map_data()
    else:
      assert map_manager.crystal_symmetry().is_similar_symmetry(
        self.map_manager.crystal_symmetry())
      ma1 = map_manager.map_data().accessor()
      ma2 = self.map_manager.map_data().accessor()
      assert ma1.all()    == ma2.all()
      assert ma1.origin() == ma2.origin()
      assert ma1.focus()  == ma2.focus()
      map_data = map_manager.map_data()
    map_box = maptbx.copy(map_data, self.gridding_first, self.gridding_last)
    map_box.reshape(flex.grid(self.box_all))
    return iotbx.map_manager.map_manager(
      unit_cell_grid             = map_box.accessor().all(),
      unit_cell_crystal_symmetry = self.crystal_symmetry,
      origin_shift_grid_units    = (0,0,0),
      map_data                   = map_box)

  def apply_to_model(self, model=None):
    if(model is None): model = self.model
    else:
      assert model.crystal_symmetry().is_similar_symmetry(
        self.model.crystal_symmetry())
    box_uc = self.crystal_symmetry.unit_cell()
    uc     = model.crystal_symmetry().unit_cell()
    sites_frac = model.get_sites_frac()
    sites_frac_new = box_uc.fractionalize(
      uc.orthogonalize(sites_frac+self.shift_frac))
    sites_cart_new = box_uc.orthogonalize(sites_frac_new)
    pdb_hierarchy = model.get_hierarchy().deep_copy()
    pdb_hierarchy.atoms().set_xyz(sites_cart_new)
    return mmtbx.model.manager(
      model_input      = None,
      pdb_hierarchy    = pdb_hierarchy,
      crystal_symmetry = self.crystal_symmetry,
      log              = null_out())
