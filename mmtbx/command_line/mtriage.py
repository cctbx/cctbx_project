from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.mtriage

import sys
import iotbx.pdb
from libtbx import group_args
from libtbx.utils import Sorry
import mmtbx.utils
from cctbx import crystal
from scitbx.array_family import flex
import mmtbx.maps.mtriage

master_params_str = """\
  map_file_name = None
    .type = str
    .help = Map file name
  half_map_file_name_1 = None
    .type = str
    .help = Half map file name
  half_map_file_name_2 = None
    .type = str
    .help = Half map file name
  model_file_name = None
    .type = str
    .help = Model file name
  include scope mmtbx.maps.mtriage.master_params
"""

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=True)

def broadcast(m, log):
  print >> log, "-"*79
  print >> log, m
  print >> log, "*"*len(m)

def get_map(map_file):
  ccp4_map = iotbx.ccp4_map.map_reader(file_name=map_file)
  map_data = ccp4_map.data.as_double()
  space_group_number = ccp4_map.space_group_number
  cs = crystal.symmetry(ccp4_map.unit_cell().parameters(), space_group_number)
  return map_data, cs

def get_inputs(args, log, master_params):
  """
  Eventually, this will be centralized.
  """
  inputs = mmtbx.utils.process_command_line_args(
    args          = args,
    master_params = master_params)
  e = inputs.params.extract()
  #
  map_data=None
  crystal_symmetry=None
  half_map_data_1=None
  half_map_data_2=None
  pdb_hierarchy=None
  css = []
  # Model
  if(e.model_file_name is not None):
    pdb_hierarchy = iotbx.pdb.input(
      file_name = e.model_file_name).construct_hierarchy()
  # Map
  if(e.map_file_name is not None):
    map_data, cs = get_map(map_file=e.map_file_name)
    css.append(cs)
  # Half-maps
  assert [e.half_map_file_name_1, e.half_map_file_name_2].count(None) in [0,2]
  if(e.half_map_file_name_1 is not None):
    half_map_data_1, cs = get_map(map_file=e.half_map_file_name_1)
    css.append(cs)
  if(e.half_map_file_name_2 is not None):
    half_map_data_2, cs = get_map(map_file=e.half_map_file_name_2)
    css.append(cs)
  # Crystal symmetry
  if(len(css)>0):
    crystal_symmetry = css[0]
  if(len(css)>1):
    for cs in css[1:]:
      assert crystal_symmetry.is_similar_symmetry(cs)
  return group_args(
    params           = inputs.params.extract(),
    pdb_hierarchy    = pdb_hierarchy,
    map_data         = map_data,
    half_map_data_1  = half_map_data_1,
    half_map_data_2  = half_map_data_2,
    crystal_symmetry = crystal_symmetry)

def show_histogram(data, n_slots, log, data_1=None, data_2=None):
  hm = flex.histogram(data = data.as_1d(), n_slots = n_slots)
  if(data_1 is None):
    print >> log, "                   Values                 Map"
    lc_1 = hm.data_min()
    s_1 = enumerate(hm.slots())
    for (i_1,n_1) in s_1:
      hc_1 = hm.data_min() + hm.slot_width() * (i_1+1)
      print >> log, "%15.4f - %-15.4f : %d" % (lc_1, hc_1, n_1)
      lc_1 = hc_1
  else:
    print >> log, \
      "                   Values                 Map Half-map1 Half-map2"
    data_min = min(flex.min(data), flex.min(data_1), flex.min(data_2))
    data_max = min(flex.max(data), flex.max(data_1), flex.max(data_2))
    h0 = flex.histogram(data = data.as_1d(), data_min=data_min,
      data_max=data_max, n_slots = n_slots)
    h1 = flex.histogram(data = data_1.as_1d(), data_min=data_min,
      data_max=data_max, n_slots = n_slots)
    h2 = flex.histogram(data = data_2.as_1d(), data_min=data_min,
      data_max=data_max, n_slots = n_slots)
    lc_1 = data_min
    s_0 = enumerate(h0.slots())
    s_1 = h1.slots()
    s_2 = h2.slots()
    for (i_1,n_1) in s_0:
      hc_1 = data_min + h0.slot_width() * (i_1+1)
      print >> log, "%15.4f - %-15.4f : %9d %9d %9d" % (
        lc_1, hc_1, n_1, s_1[i_1], s_2[i_1])
      lc_1 = hc_1

def run(args, log=sys.stdout):
  """phenix.mtriage:
  Given map file and optionally model and half-map files compute map statistics.

How to run:
  phenix.mtriage model_file_name=m.pdb map_file_name=m.map half_map_file_name_1=m1.map half_map_file_name_2=m2.map

Optional: model_file_name=, half_map_file_name_1=, half_map_file_name_2=

Feedback:
  PAfonine@lbl.gov
  phenixbb@phenix-online.org
  """
  assert len(locals().keys()) == 2 # intentional
  print >> log, "-"*79
  print >> log, run.__doc__
  print >> log, "-"*79
  # Get inputs
  inputs = get_inputs(
    args          = args,
    log           = log,
    master_params = master_params())
  # Map
  map_data = inputs.map_data
  if(map_data is None):
    raise Sorry("Map is not defined.")
  # Model
  broadcast(m="Input model:", log=log)
  if(inputs.pdb_hierarchy is None):
    print >> log, "No model specified."
  else:
    inputs.pdb_hierarchy.show(level_id="chain")
  # Crystal symmetry
  broadcast(m="Box (unit cell) info:", log=log)
  if(inputs.crystal_symmetry is None):
    raise Sorry("Box (unit cell) parameters are not defined.")
  inputs.crystal_symmetry.show_summary(f=log)
  #
  # Map statistics
  #
  broadcast(m="Map statistics:", log=log)
  print >> log, "Map:"
  a = map_data.accessor()
  print >> log, "  origin:      ", a.origin()
  print >> log, "  last:        ", a.last()
  print >> log, "  focus:       ", a.focus()
  print >> log, "  all:         ", a.all()
  print >> log, "  min,max,mean:", map_data.as_1d().min_max_mean().as_tuple()
  print >> log, "Half-maps:"
  if(inputs.half_map_data_1 is None):
    print >> log, "  Half-maps are not provided."
  else:
    for i, m in enumerate([inputs.half_map_data_1,inputs.half_map_data_2]):
      print >> log, "  half-map:", i+1
      print >> log, "    origin:", m.accessor().origin()
      print >> log, "    last:  ", m.accessor().last()
      print >> log, "    focus: ", m.accessor().focus()
      print >> log, "    all:   ", m.accessor().all()
      print >> log, "    min,max,mean:", m.as_1d().min_max_mean().as_tuple()
  # Map
  print >> log, "Histogram of map values:"
  show_histogram(
    data    = map_data,
    data_1  = inputs.half_map_data_1,
    data_2  = inputs.half_map_data_2,
    n_slots = 20,
    log     = log)
  #
  task_obj = mmtbx.maps.mtriage.mtriage(
    map_data         = inputs.map_data,
    crystal_symmetry = inputs.crystal_symmetry,
    params           = inputs.params,
    half_map_data_1  = inputs.half_map_data_1,
    half_map_data_2  = inputs.half_map_data_2,
    pdb_hierarchy    = inputs.pdb_hierarchy)
  task_obj.validate()
  task_obj.run()
  results = task_obj.get_results()
  # show results
  print >> log, "Map resolution estimates:"
  print >> log, "  using map alone (d99)             :", results.d99
  print >> log, "  comparing with model (d_model)    :", results.d_model
  print >> log, "    b_iso_overall                   :", results.b_iso_overall
  print >> log, "  FSC(map,model map)=0 (d_fsc_model):", results.d_fsc_model
  print >> log, "  d99 (half map 1)                  :", results.d99_1
  print >> log, "  d99 (half map 2)                  :", results.d99_2
  print >> log, "  FSC(half map 1,2)=0.143 (d_fsc)   :", results.d_fsc
  print >> log
  #
  if(results.fsc_curve_model is not None):
    fn = "fsc_model.mtriage.log"
    of = open(fn,"w")
    for a,b in zip(results.fsc_curve_model.fsc.d_inv,
                   results.fsc_curve_model.fsc.fsc):
      print >>of, "%15.9f %15.9f"%(a,b)
    of.close()
    print >> log, "FSC(model map, map) is written to %s"%fn
  #
  if(results.fsc_curve is not None):
    fn = "fsc_half_maps.mtriage.log"
    of = open("fsc_half_maps.mtriage.log","w")
    for a,b in zip(results.fsc_curve.fsc.d_inv, results.fsc_curve.fsc.fsc):
      print >>of, "%15.9f %15.9f"%(a,b)
    of.close()
    print >> log, "FSC(half map 1, half map 1) is written to %s"%fn
  #
  return None

if (__name__ == "__main__"):
  assert run(args=sys.argv[1:]) is None # assert here is intentional
