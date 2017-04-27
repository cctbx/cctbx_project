from __future__ import division
import sys
from libtbx import group_args
from cctbx import maptbx
import iotbx.phil
from libtbx import adopt_init_args
from cctbx.maptbx import resolution_from_map_and_model
import mmtbx.utils
from libtbx import group_args
from cctbx import miller
from mmtbx.maps import correlation

master_params_str = """
  map_file_name = None
    .type = str
    .help = Map file name
  model_file_name = None
    .type = str
    .help = Model file name
  resolution = None
    .type = float
    .help = Data (map) resolution
  scattering_table = wk1995  it1992  *n_gaussian  neutron electron
    .type = choice
    .help = Scattering table (X-ray, neutron or electron)
  atom_radius = None
    .type = float
    .help = Atom radius for masking. If undefined then calculated automatically
  compute_fsc_curve_model = True
    .type = bool
    .help = Compute model-map CC in reciprocal space: FSC(model map, data map)
  compute_d_model = True
    .type = bool
    .help = Resolution estimate using model and map
  compute_d99 = True
    .type = bool
    .help = Resolution estimate d99
  mask_maps = True
    .type = bool
    .help = Mask out region outside molecule
"""

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=False)

def get_box(map_data, pdb_hierarchy, xray_structure):
  if(pdb_hierarchy is not None):
    box = mmtbx.utils.extract_box_around_model_and_map(
      xray_structure         = xray_structure,
      map_data               = map_data,
      box_cushion            = 5.0,
      selection              = None,
      density_select         = None,
      threshold              = None)
    pdb_hierarchy.adopt_xray_structure(box.xray_structure_box)
    return group_args(
      map_data       = box.map_box,
      xray_structure = box.xray_structure_box,
      pdb_hierarchy  = pdb_hierarchy)
  else:
    return None

class mtriage(object):
  def __init__(self,
               map_data,
               crystal_symmetry,
               params=master_params(),
               half_map_data_1=None,
               half_map_data_2=None,
               pdb_hierarchy=None,
               nproc=1):
    adopt_init_args(self, locals())
    assert [half_map_data_1, half_map_data_2].count(None) in [0,2]
    # Results
    self.d9              = None
    self.d99             = None
    self.d999            = None
    self.d99_1           = None
    self.d99_2           = None
    self.d_model         = None
    self.b_iso_overall   = None
    self.d_fsc           = None
    self.d_fsc_model     = None
    self.fsc_curve       = None
    self.fsc_curve_model = None
    # Internal work objects
    self.f   = None
    self.f1  = None
    self.f2  = None
    self.box = None
    self.xray_structure = None

  def validate(self):
    if(not [self.half_map_data_1, self.half_map_data_2].count(None) in [0,2]):
      raise Sorry("None or two half-maps are required.")
    if(self.half_map_data_1 is not None):
      correlation.assert_same_gridding(
        map_1 = self.half_map_data_1,
        map_2 = self.half_map_data_2,
        Sorry_message="Half-maps have different gridding.")
      correlation.assert_same_gridding(
        map_1 = self.map_data,
        map_2 = self.half_map_data_2,
        Sorry_message="Half-maps and full map have different gridding.")
    if(self.crystal_symmetry.space_group().type().number()!=1):
      raise Sorry("Symmetry must be P1")

  def run(self):
    # Extract xrs from pdb_hierarchy
    self._get_xray_structure()
    # Compute mask
#    self._compute_mask()
    # Extract box around model with map
    self._get_box()
    # Compute d99
    self._compute_d99()
    # Compute d_model
    self._compute_d_model()
    # Compute half-map FSC
    self._compute_half_map_fsc()
    # Map-model FSC and d_fsc_model
    self._compute_model_map_fsc()

  def _get_xray_structure(self):
    if(self.pdb_hierarchy is not None):
      self.pdb_hierarchy.atoms().reset_i_seq()
      self.xray_structure = self.pdb_hierarchy.extract_xray_structure(
        crystal_symmetry = self.crystal_symmetry)

  def _get_box(self):
    if(self.pdb_hierarchy is not None):
      self.box = get_box(
        map_data       = self.map_data,
        pdb_hierarchy  = self.pdb_hierarchy,
        xray_structure = self.xray_structure)

  def _compute_d99(self):
    if(not self.params.compute_d99): return
    d99_obj = maptbx.d99(
      map              = self.map_data,
      crystal_symmetry = self.crystal_symmetry)
    self.d9   = d99_obj.result.d9
    self.d99  = d99_obj.result.d99
    self.d999 = d99_obj.result.d999
    self.f = d99_obj.f
    d99_obj_1, d99_obj_2 = None,None
    if(self.half_map_data_1 is not None):
      d99_obj_1 = maptbx.d99(
        map              = self.half_map_data_1,
        crystal_symmetry = self.crystal_symmetry)
      d99_obj_2 = maptbx.d99(
        map              = self.half_map_data_2,
        crystal_symmetry = self.crystal_symmetry)
      self.d99_1 = d99_obj_1.result.d99
      self.d99_2 = d99_obj_2.result.d99
      self.f1 = d99_obj_1.f
      self.f2 = d99_obj_2.f

  def _compute_d_model(self):
    if(not self.params.compute_d_model): return
    if(self.pdb_hierarchy is not None):
      o = resolution_from_map_and_model.run(
        map_data         = self.box.map_data,
        xray_structure   = self.box.xray_structure,
        pdb_hierarchy    = self.box.pdb_hierarchy,
        d_min_min        = 1.7,
        nproc            = self.nproc)
      self.d_model       = o.d_min
      self.b_iso_overall = o.b_iso

  def _compute_half_map_fsc(self):
    if(self.half_map_data_1 is not None):
      self.fsc_curve = self.f1.d_min_from_fsc(
        other = self.f2, bin_width=100, fsc_cutoff=0.143)
      self.d_fsc = self.fsc_curve.d_min

  def _compute_model_map_fsc(self):
    if(not self.params.compute_fsc_curve_model): return
    if(self.pdb_hierarchy is not None):
      f_obs = miller.structure_factor_box_from_map(
        map              = self.box.map_data,
        crystal_symmetry = self.box.xray_structure.crystal_symmetry())
      f_calc = f_obs.structure_factors_from_scatterers(
        xray_structure = self.box.xray_structure).f_calc()
      self.fsc_curve_model = f_calc.d_min_from_fsc(
        other=f_obs, bin_width=100, fsc_cutoff=0.0)
      self.d_fsc_model = self.fsc_curve_model.d_min

  def show_summary(self):
    r = self.get_results()
    print "d99          : ", r.d99
    print "d99_1        : ", r.d99_1
    print "d99_2        : ", r.d99_2
    print "d_model      : ", r.d_model
    print "b_iso_overall: ", r.b_iso_overall
    print "d_fsc        : ", r.d_fsc
    print "d_fsc_model  : ", r.d_fsc_model
    #
    of = open("fsc_curve_model","w")
    for a,b in zip(r.fsc_curve_model.fsc.d_inv, r.fsc_curve_model.fsc.fsc):
      print >>of, "%15.9f %15.9f"%(a,b)
    of.close()
    #
    if(r.fsc_curve is not None):
      of = open("fsc_curve","w")
      for a,b in zip(r.fsc_curve.fsc.d_inv, r.fsc_curve.fsc.fsc):
        print >>of, "%15.9f %15.9f"%(a,b)
      of.close()

  def get_results(self):
    return group_args(
      d9              = self.d9,
      d99             = self.d99,
      d999            = self.d999,
      d99_1           = self.d99_1,
      d99_2           = self.d99_2,
      d_model         = self.d_model,
      b_iso_overall   = self.b_iso_overall,
      d_fsc           = self.d_fsc,
      d_fsc_model     = self.d_fsc_model,
      fsc_curve       = self.fsc_curve,
      fsc_curve_model = self.fsc_curve_model)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
