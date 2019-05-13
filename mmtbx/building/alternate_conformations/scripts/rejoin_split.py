
from __future__ import division
from __future__ import print_function
import libtbx.phil
import sys

master_phil = libtbx.phil.parse("""
file_name = None
  .type = path
model_error_ml = None
  .type = float
include scope mmtbx.building.alternate_conformations.rejoin_phil
""", process_includes=True)

def run(args, out=sys.stdout):
  from mmtbx.building import alternate_conformations
  import iotbx.phil
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=master_phil,
    float_def="model_error_ml",
    pdb_file_def="file_name")
  params = cmdline.work.extract()
  assert (params.file_name is not None)
  pdb_in = cmdline.get_file(params.file_name)
  hierarchy = pdb_in.file_object.construct_hierarchy()
  xrs = pdb_in.file_object.xray_structure_simple()
  n_modified = alternate_conformations.rejoin_split_single_conformers(
    pdb_hierarchy=hierarchy,
    crystal_symmetry=xrs,
    model_error_ml=params.model_error_ml,
    params=params,
    verbose=True,
    log=out)
  alternate_conformations.finalize_model(
    pdb_hierarchy=hierarchy,
    xray_structure=hierarchy.extract_xray_structure(
      crystal_symmetry=pdb_in.file_object.xray_structure_simple()),
    set_b_iso=True,
    convert_to_isotropic=True)
  if (n_modified == 0):
    print("No residues modified.", file=out)
  else :
    open("rejoined.pdb", "w").write(hierarchy.as_pdb_string())
    print("wrote rejoined.pdb", file=out)

if (__name__ == "__main__"):
  run(sys.argv[1:])
