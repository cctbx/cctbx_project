from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.map_model_cc
# LIBTBX_SET_DISPATCHER_NAME phenix.model_map_cc

import sys
import time
import iotbx.pdb
from libtbx import group_args, easy_pickle
from libtbx.utils import Sorry
import mmtbx.utils
import mmtbx.maps.map_model_cc
from libtbx.str_utils import format_value
from mmtbx.maps import mtriage

master_params_str = """\
  map_file_name = None
    .type = str
    .help = Map file name
  model_file_name = None
    .type = str
    .help = Model file name
  pkl_file_name = None
    .type = path
    .help = File name for pickle file with results
    .expert_level = 3
  include scope mmtbx.maps.map_model_cc.master_params
  output_file_name_prefix=None
    .type=str
"""

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=True)


def broadcast(m, log):
  print >> log, "-"*79
  print >> log, m
  print >> log, "*"*len(m)

def get_inputs(args,
               log,
               master_params,
               need_map=True,
               need_model_hierarchy=True,
               need_crystal_symmetry=True):
  """
  Eventually, this will be centralized.
  """
  inputs = mmtbx.utils.process_command_line_args(
    args          = args,
    master_params = master_params)
  # Model
  pdb_file_name, pdb_hierarchy = None, None
  if(need_model_hierarchy):
    file_names = inputs.pdb_file_names
    if(len(file_names) != 1):
      raise Sorry("One model (PDB or mmCIF) required.")
    pdb_file_name = file_names[0]
    pdb_inp = iotbx.pdb.input(
      file_name = pdb_file_name)
    pdb_hierarchy = pdb_inp.construct_hierarchy()
  # Map
  ccp4_map_object = None
  if(need_map):
    if(inputs.ccp4_map is None):
      raise Sorry("Map file has to given.")
    ccp4_map_object = inputs.ccp4_map
  # Crystal symmetry
  crystal_symmetry = None
  if(need_crystal_symmetry):
    crystal_symmetry = inputs.crystal_symmetry
    if(crystal_symmetry is None):
      raise Sorry("No box (unit cell) info found.")
  #
  return group_args(
    params           = inputs.params.extract(),
    pdb_file_name    = pdb_file_name,
    pdb_inp          = pdb_inp,
    pdb_hierarchy    = pdb_hierarchy,
    ccp4_map_object  = ccp4_map_object,
    crystal_symmetry = crystal_symmetry)

def get_fsc(map_inp, pdb_inp, params):
  result = None
  if(params.compute.fsc):
    mtriage_params = mtriage.master_params().extract()
    mtriage_params.scattering_table = params.scattering_table
    mtriage_params.compute.map_counts = False
    mtriage_params.compute.fsc_curve_model = True
    mtriage_params.compute.d_fsc_model_05 = False
    mtriage_params.compute.d_fsc_model_0 = False
    mtriage_params.compute.d_fsc_model_0143 = False
    mtriage_params.compute.d_model = False
    mtriage_params.compute.d_model_b0 = False
    mtriage_params.compute.d99 = False
    mtriage_params.mask_maps = True
    #mtriage_params.radius_smooth = self.atom_radius
    mtriage_params.resolution = params.resolution
    result = mtriage.mtriage(
      map_inp = map_inp,
      pdb_inp = pdb_inp,
      params  = mtriage_params).get_results().masked.fsc_curve_model
  return result

def run(args, log=sys.stdout):
  """phenix.map_model_cc or phenix.model_map_cc:
  Given PDB file and a map file calculate model-map coorelation.

How to run:
  phenix.map_model_cc model.pdb map.ccp4 resolution=3
  phenix.model_map_cc model.pdb map.ccp4 resolution=3

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
  # Model
  broadcast(m="Input PDB:", log=log)
  print >> log, inputs.pdb_file_name # ideally this should not be available here
  inputs.pdb_hierarchy.show(level_id="chain")
  # Crystal symmetry
  broadcast(m="Box (unit cell) info:", log=log)
  inputs.crystal_symmetry.show_summary(f=log)
  # Map
  broadcast(m="Input map:", log=log)
  inputs.ccp4_map_object.show_summary(prefix="  ")
  # Run task in 4 separate steps
  task_obj = mmtbx.maps.map_model_cc.map_model_cc(
    map_data         = inputs.ccp4_map_object.data.as_double(),
    pdb_hierarchy    = inputs.pdb_hierarchy,
    crystal_symmetry = inputs.crystal_symmetry,
    params           = inputs.params.map_model_cc)
  task_obj.validate()
  task_obj.run()
  results = task_obj.get_results()
  #results.fsc=None
  results.fsc = get_fsc(
    map_inp = inputs.ccp4_map_object,
    pdb_inp = inputs.pdb_inp,
    params  = inputs.params.map_model_cc)
  if inputs.params.pkl_file_name is not None:
    easy_pickle.dump(file_name=inputs.params.pkl_file_name, obj=group_args(
      cc_mask        = results.cc_mask,
      cc_volume      = results.cc_volume,
      cc_peaks       = results.cc_peaks,
      cc_per_chain   = results.cc_per_chain,
      # cc_per_residue = results.cc_per_residue,
      fsc            = results.fsc,
      ))
  #
  broadcast(m="Map resolution:", log=log)
  print >> log, "  Resolution:", inputs.params.map_model_cc.resolution
  broadcast(m="Map-model CC (overall):", log=log)
  print >> log, "  CC_mask  : %s" % format_value("%6.4f", results.cc_mask)
  print >> log, "  CC_volume: %s" % format_value("%6.4f", results.cc_volume)
  print >> log, "  CC_peaks : %s" % format_value("%6.4f", results.cc_peaks)
  log.flush()
  if results.fsc is not None:
    broadcast(m="Model-map FSC:", log=log)
    if(inputs.params.output_file_name_prefix is None):
      fsc_fn = "map_model_fsc.log"
    else:
      fsc_fn = "%s_fsc.log"%inputs.params.output_file_name_prefix
    out = open(fsc_fn, "w")
    print >> log, "  saved to:", fsc_fn
    print >> out, "    1/resolution    CC"
    for a,b in zip(results.fsc.d_inv, results.fsc.fsc):
      print >> out, "%15.9f %15.9f"%(a,b)
    out.close()
  if len(results.cc_per_chain) + len(results.cc_per_residue) > 0:
    broadcast(m="Map-model CC (local):", log=log)
  # Per chain
  if len(results.cc_per_chain) > 0:
    print >> log, "Per chain:"
    print >> log, "chain ID  CC       <B>    <occ>   N atoms"
    fmt = "%s        %7.4f %8.3f %4.2f    %d"
    for r in results.cc_per_chain:
      print fmt%(r.chain_id, r.cc, r.b_iso_mean, r.occ_mean, r.n_atoms)
  # Per residue
  if len(results.cc_per_residue) > 0:
    print >> log, "Per residue:"
    if(inputs.params.output_file_name_prefix is None):
      fn = "cc_per_residue.log"
    else:
      fn = "%s_cc_per_residue.log"%inputs.params.output_file_name_prefix
    out = open(fn, "w")
    print >> log, "  saved to:", fn
    fmt = "%s %s %s %7.4f %8.3f %4.2f"
    for r in results.cc_per_residue:
      print >> out, fmt%(r.chain_id, r.resname, r.resseq, r.cc, r.b_iso_mean,
        r.occ_mean)
    out.close()
    #
  return None

if (__name__ == "__main__"):
  t0 = time.time()
  assert run(args=sys.argv[1:]) is None # assert here is intentional
  print "Total time: %8.3f"%(time.time()-t0)
