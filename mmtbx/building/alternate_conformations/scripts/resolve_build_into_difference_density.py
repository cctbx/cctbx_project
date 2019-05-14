
from __future__ import absolute_import, division, print_function
import iotbx.phil
from libtbx.utils import Sorry
import time
import sys

master_phil = iotbx.phil.parse("""
include scope mmtbx.utils.cmdline_input_phil_str
selection = None
  .type = atom_selection
partial_occupancy = 0.6
  .type = float
nproc = Auto
  .type = int
output_file_name = resolve.pdb
  .type = path
resolve_build {
  include scope mmtbx.building.alternate_conformations.resolve_build.resolve_build_params_str
}
""", process_includes=True)

def run(args, out=None):
  if (out is None) : out = sys.stdout
  from mmtbx.building.alternate_conformations import resolve_build
  from mmtbx.building import disorder
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
  two_fofc_map, fofc_map = disorder.get_partial_omit_map(
    fmodel=fmodel,
    selection=selection.iselection(),
    selection_delete=None,#nearby_water_selection,
    negate_surrounding=False,
    partial_occupancy=params.partial_occupancy)
  build = resolve_build.resolve_builder(
    params=params.resolve_build,
    pdb_hierarchy=hierarchy,
    xray_structure=fmodel.xray_structure,
    processed_pdb_file=cmdline.processed_pdb_file,
    target_map=fofc_map,
    selection=selection,
    d_min=fmodel.f_obs().d_min(),
    out=out)

if (__name__ == '__main__'):
  run(sys.argv[1:])
