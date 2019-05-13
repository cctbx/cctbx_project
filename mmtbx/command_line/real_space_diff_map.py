from __future__ import division
from __future__ import print_function
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

legend = """phenix.real_space_diff_map:
  Given PDB file and a map file calculate difference map:
    ExperimentalMap-ModelMap (like Fo-Fc map, in reciprocal space).

How to run:
  phenix.real_space_diff_map model.pdb map.ccp4 resolution=3.5

Feedback:
  PAfonine@lbl.gov
  phenixbb@phenix-online.org"""

master_params_str = """
  map_file_name = None
    .type = str
  model_file_name = None
    .type = str
  resolution = None
    .type = float
  scattering_table = wk1995  it1992  *n_gaussian  neutron electron
    .type = choice
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
    master_params = master_params())
  params = inputs.params.extract()
  # model
  broadcast(m="Input PDB:", log=log)
  file_names = inputs.pdb_file_names
  if(len(file_names) != 1): raise Sorry("PDB file has to given.")
  pi = iotbx.pdb.input(file_name = file_names[0])
  h = pi.construct_hierarchy()
  xrs = pi.xray_structure_simple(crystal_symmetry=inputs.crystal_symmetry)
  xrs.scattering_type_registry(table = params.scattering_table)
  xrs.show_summary(f=log, prefix="  ")
  # map
  broadcast(m="Input map:", log=log)
  if(inputs.ccp4_map is None): raise Sorry("Map file has to given.")
  inputs.ccp4_map.show_summary(prefix="  ")
  map_data = inputs.ccp4_map.map_data()
  # shift origin if needed
  soin = maptbx.shift_origin_if_needed(map_data=map_data,
    sites_cart=xrs.sites_cart(), crystal_symmetry=xrs.crystal_symmetry())
  map_data = soin.map_data
  xrs.set_sites_cart(soin.sites_cart)
  # estimate resolution
  d_min = params.resolution
  if(d_min is None):
    raise Sorry("Map resolution must be given.")
  print("  d_min: %6.4f"%d_min, file=log)
  #
  result_obj = compdiff(
    map_data_obs = map_data,
    xrs          = xrs,
    d_min        = d_min,
    vector_map   = False)
  write_ccp4_map(
    map_data    = result_obj.map_result,
    unit_cell   = xrs.unit_cell(),
    space_group = xrs.space_group(),
    file_name   = "map_model_difference_1.ccp4")
  #
  result_obj = compdiff(
    map_data_obs = map_data,
    xrs          = xrs,
    d_min        = d_min,
    vector_map   = True)
  write_ccp4_map(
    map_data    = result_obj.map_result,
    unit_cell   = xrs.unit_cell(),
    space_group = xrs.space_group(),
    file_name   = "map_model_difference_2.ccp4")

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
  bs = flex.double([i for i in xrange(0,100)])
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
