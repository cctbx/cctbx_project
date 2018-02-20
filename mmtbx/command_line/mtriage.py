from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.mtriage

import sys
import iotbx.pdb
from libtbx import group_args
import mmtbx.utils
import mmtbx.maps.mtriage
from iotbx import ccp4_map
from scitbx.array_family import flex
from libtbx.str_utils import format_value
from libtbx import introspection
from libtbx.utils import null_out

master_params_GUI_str = """\
  include scope libtbx.phil.interface.tracking_params
  gui
    .help = "GUI-specific parameter required for output directory"
  {
    output_dir = None
    .type = path
    .style = output_dir
  }
"""

master_params_str = """\
  map_file_name = None
    .type = path
    .help = Map file name
    .short_caption = Map file
    .style = file_type:ccp4_map bold input_file
  half_map_file_name_1 = None
    .type = path
    .help = Half map file name
    .short_caption = Half map file 1
    .style = file_type:ccp4_map input_file
  half_map_file_name_2 = None
    .type = path
    .help = Half map file name
    .short_caption = Half map file 2
    .style = file_type:ccp4_map input_file
  model_file_name = None
    .type = path
    .help = Model file name
    .short_caption = Model file
    .style = file_type:pdb input_file
  include scope mmtbx.maps.mtriage.master_params
  fsc_model_plot_file_name_prefix = fsc_model
    .type = str
  fsc_half_maps_file_name_prefix = fsc_half_maps
    .type = str
  write_mask_file = True
    .type = bool
  mask_file_name = mask.ccp4
    .type = str
  %(master_params_GUI_str)s
""" % vars()

master_params = iotbx.phil.parse(master_params_str, process_includes=True)

def broadcast(m, log):
  print >> log, "-"*79
  print >> log, m
  print >> log, "*"*len(m)

def get_inputs(args, log, master_params):
  """
  Eventually, this will be centralized.
  """
  inputs = mmtbx.utils.process_command_line_args(
    args          = args,
    master_params = master_params)
  e = inputs.params.extract()
  # Model
  pdb_inp = None
  if(e.model_file_name is not None):
    pdb_inp = iotbx.pdb.input(file_name = e.model_file_name)
  # Map
  map_inp = None
  if(e.map_file_name is not None):
    map_inp = iotbx.ccp4_map.map_reader(file_name=e.map_file_name)
  # Half-maps
  map_inp_1 = None
  if(e.half_map_file_name_1 is not None):
    map_inp_1 = iotbx.ccp4_map.map_reader(file_name=e.half_map_file_name_1)
  map_inp_2 = None
  if(e.half_map_file_name_2 is not None):
    map_inp_2 = iotbx.ccp4_map.map_reader(file_name=e.half_map_file_name_2)
  return group_args(
    map_inp   = map_inp,
    map_inp_1 = map_inp_1,
    map_inp_2 = map_inp_2,
    pdb_inp   = pdb_inp,
    params    = inputs.params.extract())

def show_histogram(map_histograms, log):
  if(map_histograms.h_half_map_1 is None):
    hm = map_histograms.h_map
    print >> log, "                   Values                 Map"
    lc_1 = hm.data_min()
    s_1 = enumerate(hm.slots())
    for (i_1,n_1) in s_1:
      hc_1 = hm.data_min() + hm.slot_width() * (i_1+1)
      print >> log, "%8.4f - %-8.4f : %d" % (lc_1, hc_1, n_1)
      lc_1 = hc_1
  else:
    print >> log, \
      "                Full map                   Half-map1 Half-map2"
    h0 = map_histograms.h_map
    h1 = map_histograms.h_half_map_1
    h2 = map_histograms.h_half_map_2
    data_min = map_histograms._data_min
    lc_2 = data_min
    lc_1 = h0.data_min()
    s_0 = enumerate(h0.slots())
    s_1 = h1.slots()
    s_2 = h2.slots()
    for (i_1,n_1) in s_0:
      hc_1 = h0.data_min() + h0.slot_width() * (i_1+1)
      hc_2 = data_min + h2.slot_width() * (i_1+1)
      print >> log, "%8.4f - %-8.4f : %9d %8.4f - %-8.4f : %9d %9d" % (
        lc_1, hc_1, n_1, lc_2, hc_2, s_1[i_1], s_2[i_1])
      lc_1 = hc_1
      lc_2 = hc_2
    print >> log, "  Half-maps, correlation of histograms: ", \
      map_histograms.half_map_histogram_cc

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
  introspection.virtual_memory_info().show_if_available(out=null_out(),
    show_max=True) # just to initialize something
  # Get inputs
  inputs = get_inputs(
    args          = args,
    log           = log,
    master_params = master_params)
  #
  task_obj = mmtbx.maps.mtriage.mtriage(
    map_inp   = inputs.map_inp,
    map_inp_1 = inputs.map_inp_1,
    map_inp_2 = inputs.map_inp_2,
    pdb_inp   = inputs.pdb_inp,
    params    = inputs.params)
  results = task_obj.get_results()
  #
  # Map statistics
  #
  broadcast(m="Map statistics:", log=log)
  print >> log, "Map:"
  print >> log, "  origin:      ", results.counts.origin
  print >> log, "  last:        ", results.counts.last
  print >> log, "  focus:       ", results.counts.focus
  print >> log, "  all:         ", results.counts.all
  print >> log, "  min,max,mean:", results.counts.min_max_mean
  print >> log, "  d_min_corner:", "%7.3f"%results.counts.d_min_corner
  #
  print >> log, "Half-maps:"
  if(inputs.map_inp_1 is None):
    print >> log, "  Half-maps are not provided."
  #
  print >> log, "Histogram(s) of map values (masked):"
  show_histogram(map_histograms = results.histograms, log = log)
  # show results
  fv = format_value
  fs = "%8.2f"
  rm = results.masked
  ru = results.unmasked
  if([rm,ru].count(None)==0):
    print >> log, "Map resolution estimates:              masked unmasked"
    print >> log, "  using map alone (d99)            :", fv(fs,rm.d99)          , fv(fs,ru.d99)
    print >> log, "  using map alone (d9999)          :", fv(fs,rm.d9999)        , fv(fs,ru.d9999)
    print >> log, "  comparing with model (d_model)   :", fv(fs,rm.d_model)      , fv(fs,ru.d_model)
    print >> log, "    b_iso_overall                  :", fv(fs,rm.b_iso_overall), fv(fs,ru.b_iso_overall)
    print >> log, "  comparing with model (d_model_b0):", fv(fs,rm.d_model_b0)   , fv(fs,ru.d_model_b0)
    print >> log, "    b_iso_overall=0"
    print >> log, "  d_fsc_model:"
    print >> log, "    FSC(map,model map)=0           :", fv(fs,rm.d_fsc_model_0)   , fv(fs,ru.d_fsc_model_0)
    print >> log, "    FSC(map,model map)=0.143       :", fv(fs,rm.d_fsc_model_0143), fv(fs,ru.d_fsc_model_0143)
    print >> log, "    FSC(map,model map)=0.5         :", fv(fs,rm.d_fsc_model_05)  , fv(fs,ru.d_fsc_model_05)
    print >> log, "  d99 (half map 1)                 :", fv(fs,rm.d99_1)           , fv(fs,ru.d99_1)
    print >> log, "  d99 (half map 2)                 :", fv(fs,rm.d99_2)           , fv(fs,ru.d99_2)
    print >> log, "  FSC(half map 1,2)=0.143 (d_fsc)  :", fv(fs,rm.d_fsc)           , fv(fs,ru.d_fsc)
    print >> log
    #
    print >> log, "Radius used for mask smoothing:", format_value("%6.2f", results.masked.radius_smooth)
    print >> log
  else:
    r = rm
    if(r is None): r=ru
    print >> log, "Map resolution estimates:              masked unmasked"
    print >> log, "  using map alone (d99)            :", fv(fs,r.d99)
    print >> log, "  using map alone (d9999)          :", fv(fs,r.d9999)
    print >> log, "  comparing with model (d_model)   :", fv(fs,r.d_model)
    print >> log, "    b_iso_overall                  :", fv(fs,r.b_iso_overall)
    print >> log, "  comparing with model (d_model_b0):", fv(fs,r.d_model_b0)
    print >> log, "    b_iso_overall=0"
    print >> log, "  d_fsc_model:"
    print >> log, "    FSC(map,model map)=0           :", fv(fs,r.d_fsc_model_0)
    print >> log, "    FSC(map,model map)=0.143       :", fv(fs,r.d_fsc_model_0143)
    print >> log, "    FSC(map,model map)=0.5         :", fv(fs,r.d_fsc_model_05)
    print >> log, "  d99 (half map 1)                 :", fv(fs,r.d99_1)
    print >> log, "  d99 (half map 2)                 :", fv(fs,r.d99_2)
    print >> log, "  FSC(half map 1,2)=0.143 (d_fsc)  :", fv(fs,r.d_fsc)
    print >> log
    #
    print >> log, "Radius used for mask smoothing:", format_value("%s", str(r.radius_smooth))
    print >> log
  #
  for r in [(results.masked, "masked"),(results.unmasked, "unmasked")]:
    if(r[0] is None): continue
    # FSC_model curve
    if(r[0].fsc_curve_model is not None):
      file_name = "%s.%s.mtriage.log"%(
        inputs.params.fsc_model_plot_file_name_prefix, r[1])
      of = open(file_name,"w")
      for a,b in zip(r[0].fsc_curve_model.d_inv, r[0].fsc_curve_model.fsc):
        print >> of, "%15.9f %15.9f"%(a,b)
      of.close()
      print >> log, "FSC(model map, map) is written to %s"%file_name
    # Mask
    if(inputs.params.write_mask_file and r[0].mask is not None):
      print >> log, "Mask is written to %s"%inputs.params.mask_file_name
      ccp4_map.write_ccp4_map(
        file_name   = inputs.params.mask_file_name,
        unit_cell   = inputs.crystal_symmetry.unit_cell(),
        space_group = inputs.crystal_symmetry.space_group(),
        map_data    = r[0].mask,
        labels      = flex.std_string(["mask"]))
    # FSC (half-maps) curve
    if(r[0].fsc_curve is not None):
      file_name = "%s.%s.mtriage.log"%(
        inputs.params.fsc_half_maps_file_name_prefix, r[1])
      of = open(file_name,"w")
      for a,b in zip(r[0].fsc_curve.fsc.d_inv, r[0].fsc_curve.fsc.fsc):
        print >> of, "%15.9f %15.9f"%(a,b)
      of.close()
      print >> log, "FSC(half map 1, half map 1) is written to %s"%file_name
  #
  return results # required for GUI

# GUI-specific class for running command
# required for GUI

from libtbx import runtime_utils
from wxGUI2 import utils
import os

def validate_params(params):
  return True

class launcher (runtime_utils.target_with_save_result) :
  def run (self) :
    utils.safe_makedirs(self.output_dir)
    os.chdir(self.output_dir)
    result = run(args=self.args, log=sys.stdout)
    return result

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
