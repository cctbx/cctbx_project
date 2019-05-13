
from __future__ import division
from __future__ import print_function
import iotbx.phil
from libtbx.utils import Sorry
import time
import sys

master_phil = iotbx.phil.parse("""
include scope mmtbx.utils.cmdline_input_phil_str
selection = None
  .type = atom_selection
nproc = Auto
  .type = int
output_file_name = refined.pdb
  .type = path
include scope mmtbx.building.alternate_conformations.real_space_annealing.master_params_str
""", process_includes=True)

def run(args, out=None):
  if (out is None) : out = sys.stdout
  from mmtbx.building.alternate_conformations import real_space_annealing
  import mmtbx.maps.utils
  import mmtbx.utils
  cmdline = mmtbx.utils.cmdline_load_pdb_and_data(
    args=args,
    master_phil=master_phil,
    process_pdb_file=True,
    scattering_table="n_gaussian",
    out=out)
  params = cmdline.params
  working_phil = master_phil.format(python_object=params)
  master_phil.fetch_diff(source=working_phil).show(out=out)
  fmodel = cmdline.fmodel
  hierarchy = cmdline.pdb_hierarchy
  sele_cache = hierarchy.atom_selection_cache()
  if (params.selection is None):
    raise Sorry("Please specify the selection parameter.")
  selection = sele_cache.selection(params.selection)
  assert (selection.count(True) > 0)
  t1 = time.time()
  models = []
  refined_hierarchy = real_space_annealing.refine_into_difference_density(
    fmodel=fmodel,
    pdb_hierarchy=hierarchy,
    processed_pdb_file=cmdline.processed_pdb_file,
    selection=selection,
    params=params,
    nproc=params.nproc,
    out=out).as_pdb_ensemble(log=out)
  t2 = time.time()
  print("refinement time: %.3fs" % (t2-t1), file=out)
  if (refined_hierarchy is not None):
    refined_hierarchy.write_pdb_file(file_name=params.output_file_name,
      crystal_symmetry = fmodel.xray_structure.crystal_symmetry())
    print("wrote %s" % params.output_file_name)
  else :
    raise Sorry("No conformations surviving filtering step!")

if (__name__ == "__main__"):
  run(sys.argv[1:])
