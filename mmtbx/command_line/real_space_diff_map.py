"""Computer real-space difference map"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.real_space_diff_map

from scitbx.array_family import flex
import sys
import iotbx.pdb
from libtbx.utils import Sorry
import mmtbx.utils
from cctbx import maptbx
from cctbx import miller
from cctbx import uctbx
from cctbx import crystal
from libtbx import adopt_init_args
from six.moves import zip
from six.moves import range

import boost_adaptbx.boost.python as bp
cctbx_maptbx_ext = bp.import_ext("cctbx_maptbx_ext")

legend = """phenix.real_space_diff_map:
  Given PDB file and a map file calculate difference map:
    ExperimentalMap-ModelMap (like Fo-Fc map, in reciprocal space).

How to run:
  phenix.real_space_diff_map model.pdb map.ccp4 resolution=3.5

Feedback:
  PAfonine@lbl.gov
  phenixbb@phenix-online.org"""

master_params_str = """
  map_type = *vector *obs_phase mask_model
    .type = choice(multi=True)
    .help = vector: (Fobs,Pobs)-(Fcalc,Pcalc), obs_phase: (Fobs-Fcalc,Pobs), \
            mask_model: mask out interpreted density and eliminate weak and \
            low-volume map values
  map_file_name = None
    .type = str
  model_file_name = None
    .type = str
  resolution = None
    .type = float
  scattering_table = wk1995  it1992  *n_gaussian  neutron electron
    .type = choice
  wrapping = False
    .type = bool
    .short_caption = Wrapping
    .help = Wrapping defines whether values outside map boundary can be mapped \
            inside with unit cell translations. Normally True for crystal \
            structures and False for cryo-EM
  ignore_symmetry_conflicts = None
    .type = bool
    .short_caption = Ignore symmetry conflicts
    .help = Ignore symmetry differences between input model and map (use values \
            from map).
"""

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=False)

def broadcast(m, log):
  print("-"*79, file=log)
  print(m, file=log)
  print("*"*len(m), file=log)

def run(args, log=sys.stdout):
  print("-"*79, file=log)
  print(legend, file=log)
  print("-"*79, file=log)
  inputs = mmtbx.utils.process_command_line_args(args = args,
    master_params = master_params(),
    suppress_symmetry_related_errors = True)
  params = inputs.params.extract()
  # model
  broadcast(m="Input PDB:", log=log)
  file_names = inputs.pdb_file_names
  if(len(file_names) != 1): raise Sorry("PDB file has to given.")
  from iotbx.data_manager import DataManager
  dm = DataManager()
  dm.set_overwrite(True)
  model = dm.get_model(file_names[0])

  # map
  broadcast(m="Input map:", log=log)
  if(inputs.ccp4_map is None): raise Sorry("Map file has to given.")

  from iotbx.map_model_manager import map_model_manager
  mam = map_model_manager(model = model, map_manager = inputs.ccp4_map,
     wrapping = params.wrapping,
     ignore_symmetry_conflicts = params.ignore_symmetry_conflicts)

  mam.model().setup_scattering_dictionaries(
     scattering_table=params.scattering_table)
  mam.model().get_xray_structure().show_summary(f=log, prefix="  ")
  inputs.ccp4_map.show_summary(prefix="  ")

  # estimate resolution
  d_min = params.resolution
  if(d_min is None):
    raise Sorry("Map resolution must be given.")
  print("  d_min: %6.4f"%d_min, file=log)
  #
  if("obs_phase" in params.map_type):
    result_obj = compdiff(
      map_data_obs = mam.map_manager().map_data(), # NOTE this will always wrap map
      xrs          = mam.model().get_xray_structure(),
      d_min        = d_min,
      vector_map   = False)
    output_map_manager=mam.map_manager().customized_copy(
        map_data=result_obj.map_result)
    dm.write_real_map_file(output_map_manager, "map_model_difference_1.ccp4")
  #
  if("vector" in params.map_type):
    result_obj = compdiff(
      map_data_obs = mam.map_manager().map_data(),
      xrs          = mam.model().get_xray_structure(),
      d_min        = d_min,
      vector_map   = True)
    output_map_manager=mam.map_manager().customized_copy(
        map_data=result_obj.map_result)
    dm.write_real_map_file(output_map_manager, "map_model_difference_2.ccp4")
  #
  if("mask_model" in params.map_type):
    map_data_result = remove_model_density(
      map_data = mam.map_manager().map_data(),
      xrs      = mam.model().get_xray_structure())
    output_map_manager=mam.map_manager().customized_copy(
        map_data=map_data_result)
    dm.write_real_map_file(output_map_manager, "map_model_difference_3.ccp4")

def remove_model_density(map_data, xrs, rad_inside=2):
  #
  map_data = map_data - flex.mean(map_data)
  map_data = map_data.set_selected(map_data < 0, 0)
  sd = map_data.sample_standard_deviation()
  assert sd != 0
  map_data = map_data / sd
  #
  map_at_atoms = flex.double()
  for site_frac in xrs.sites_frac():
    mv = map_data.tricubic_interpolation(site_frac)
    map_at_atoms.append( mv )
  print (flex.mean(map_at_atoms), flex.max(map_at_atoms))
  mmax = flex.max(map_at_atoms)
  cut = 0
  print (dir(map_data))
  while cut<mmax:
    map_data_ = map_data.deep_copy()
    map_data_ = map_data_.set_selected(map_data<cut, 0)
    map_data_ = map_data_.set_selected(map_data>=cut, 1)
    cut+=1

    zz = flex.double()
    for site_frac in xrs.sites_frac():
      mv = map_data_.value_at_closest_grid_point(site_frac)
      zz.append( mv )
    print(cut,  (zz==1).count(True)/zz.size()*100. )

  #
  #radii = flex.double(xrs.sites_frac().size(), rad_inside)
  #mask = cctbx_maptbx_ext.mask(
  #  sites_frac                  = xrs.sites_frac(),
  #  unit_cell                   = xrs.unit_cell(),
  #  n_real                      = map_data.all(),
  #  mask_value_inside_molecule  = 0,
  #  mask_value_outside_molecule = 1,
  #  radii                       = radii)

  mask = mmtbx.masks.mask_from_xray_structure(
    xray_structure           = xrs,
    p1                       = True,
    for_structure_factors    = True,
    solvent_radius           = None,
    shrink_truncation_radius = None,
    n_real                   = map_data.accessor().all(),
    in_asu                   = False).mask_data
  maptbx.unpad_in_place(map=mask)


  map_data = map_data * mask
  map_data = map_data.set_selected(map_data < flex.mean(map_at_atoms)/6, 0)
  #
  n = map_data.accessor().all()
  abc = xrs.unit_cell().parameters()[:3]
  print(abc[0]/n[0], abc[1]/n[1], abc[2]/n[2])

  step = abc[0]/n[0]

  co = maptbx.connectivity(
    map_data                   = map_data.deep_copy(),
    threshold                  = 0.0,
    preprocess_against_shallow = True,
    wrapping                   = False)
  conn = co.result().as_double()
  z = zip(co.regions(),range(0,co.regions().size()))
  sorted_by_volume = sorted(z, key=lambda x: x[0], reverse=True)
  mask_ = flex.double(flex.grid(n), 0)
  for i_seq, p in enumerate(sorted_by_volume):
    v, i = p
    if i_seq==0: continue
    volume = v*step**3
    print(v, volume)
    if 1:#(volume<3):
      sel = conn==i
      mask_ = mask_.set_selected(sel, 1)

  #
  return map_data*mask_

def scale_k1(x,y):
  x = x.as_1d()
  y = y.as_1d()
  den=flex.sum(y*y)
  if(abs(den)<1.e-9): return 0
  return flex.sum(x*y)/den

def write_ccp4_map(map_data, unit_cell, space_group, file_name):
  iotbx.mrcfile.write_ccp4_map(
    file_name      = file_name,
    unit_cell      = unit_cell,
    space_group    = space_group,
    map_data       = map_data.as_double(),
    labels=flex.std_string([" "]))

def scale_two_real_maps_in_fourier_space(m1, m2, cs, d_min, vector_map):
  f1 = maptbx.map_to_map_coefficients(m=m1, cs=cs, d_min=d_min)
  f2 = maptbx.map_to_map_coefficients(m=m2, cs=cs, d_min=d_min)
  if(vector_map):
    f2 = f2.phase_transfer(phase_source=f1)
  ss = 1./flex.pow2(f1.d_spacings().data()) / 4.
  bs = flex.double([i for i in range(0,100)])
  mc = mmtbx.bulk_solvent.complex_f_minus_f_kb_scaled(
    f1.data(),f2.data(),bs,ss)
  crystal_gridding = maptbx.crystal_gridding(
    unit_cell             = cs.unit_cell(),
    space_group_info      = cs.space_group_info(),
    pre_determined_n_real = m1.all())
  fft_map = miller.fft_map(
    crystal_gridding     = crystal_gridding,
    fourier_coefficients = f1.array(data=mc))
  return fft_map.real_map_unpadded()

class compdiff(object):
  def __init__(
        self,
        map_data_obs,
        xrs,
        d_min,
        vector_map,
        box_dimension=30):
    adopt_init_args(self, locals())
    self.crystal_gridding = maptbx.crystal_gridding(
      unit_cell             = self.xrs.unit_cell(),
      space_group_info      = self.xrs.space_group_info(),
      pre_determined_n_real = self.map_data_obs.all())
    self.n_real = self.crystal_gridding.n_real()
    crystal_gridding = maptbx.crystal_gridding(
      unit_cell             = self.xrs.unit_cell(),
      space_group_info      = self.xrs.space_group_info(),
      pre_determined_n_real = self.map_data_obs.all())
    mc = xrs.structure_factors(d_min=d_min).f_calc()
    fft_map = miller.fft_map(
      crystal_gridding     = crystal_gridding,
      fourier_coefficients = mc)
    fft_map.apply_sigma_scaling()
    self.map_data_calc = fft_map.real_map_unpadded()
    scale = scale_k1(x=self.map_data_obs, y=self.map_data_calc)
    self.map_data_calc = self.map_data_calc * scale
    #
    # result map
    self.map_result = flex.double(flex.grid(self.map_data_obs.all()))
    # iterate over boxes
    self.box_iterator()

  def box_iterator(self):
    p = self.xrs.unit_cell().parameters()
    b = maptbx.boxes_by_dimension(
      n_real = self.n_real,
      dim    = self.box_dimension,
      abc    = p[:3])
    i_box = 0
    for s,e in zip(b.starts, b.ends):
      i_box+=1
      map_box_obs  = maptbx.copy(self.map_data_obs,  s, e)
      map_box_calc = maptbx.copy(self.map_data_calc, s, e)
      map_box_obs.reshape(flex.grid(map_box_obs.all()))
      map_box_calc.reshape(flex.grid(map_box_calc.all()))
      #######
      # XXX Copy-paste from map_box
      abc = []
      for i in range(3):
        abc.append( p[i] * map_box_calc.all()[i]/self.n_real[i] )
      ucb = uctbx.unit_cell(
        parameters=(abc[0],abc[1],abc[2],p[3],p[4],p[5]))
      cs = crystal.symmetry(unit_cell=ucb, space_group="P1")
      #######
      diff_map = scale_two_real_maps_in_fourier_space(
        m1         = map_box_obs,
        m2         = map_box_calc,
        cs         = cs,
        d_min      = self.d_min,
        vector_map = self.vector_map)
      maptbx.set_box(
        map_data_from = diff_map,
        map_data_to   = self.map_result,
        start         = s,
        end           = e)
    sd = self.map_result.sample_standard_deviation()
    if(sd!=0):
      self.map_result = self.map_result/sd


if (__name__ == "__main__"):
  run(args=sys.argv[1:])

