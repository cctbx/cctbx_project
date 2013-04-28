
from __future__ import division
import iotbx.phil
from libtbx.str_utils import make_sub_header
from libtbx.utils import null_out
from libtbx import Auto
import os.path
import sys

master_phil = iotbx.phil.parse("""
include scope mmtbx.utils.cmdline_input_phil_str
include scope mmtbx.building.extend_sidechains.master_params
output_model = None
  .type = path
output_map_coeffs = None
  .type = path
""", process_includes=True)

def run (args, out=sys.stdout, verbose=True) :
  import mmtbx.building.extend_sidechains
  import mmtbx.utils
  input_out = out
  if (not verbose) :
    input_out = null_out()
  cmdline = mmtbx.utils.cmdline_load_pdb_and_data(
    args=args,
    master_phil=master_phil,
    process_pdb_file=False,
    out=input_out,
    usage_string="""\
mmtbx.extend_sidechains model.pdb data.mtz [restraints.cif] [options]

Rebuild sidechains with missing non-hydrogen atoms.  Includes real-space
refinement (but needs work).""")
  pdb_hierarchy = cmdline.pdb_hierarchy
  params = cmdline.params
  fmodel = cmdline.fmodel
  xrs = cmdline.xray_structure
  if (params.build_hydrogens is Auto) :
    params.build_hydrogens = xrs.hd_selection().count(True) > 0
  make_sub_header("Filling in partial sidechains", out=out)
  prefilter_callback = mmtbx.building.extend_sidechains.prefilter(
    fmodel=fmodel,
    out=out)
  n_extended = mmtbx.building.extend_sidechains.extend_protein_model(
    pdb_hierarchy=pdb_hierarchy,
    hydrogens=params.build_hydrogens,
    max_atoms_missing=params.max_atoms_missing,
    prefilter_callback=prefilter_callback,
    log=out)
  print >> out, "  %d sidechains extended." % n_extended
  if (n_extended > 0) and (not params.skip_rsr) :
    make_sub_header("Real-space refinement", out=out)
    new_hierarchy = mmtbx.building.extend_sidechains.refit_residues(
      pdb_hierarchy=pdb_hierarchy,
      cif_objects=[ co for fn, co in cmdline.cif_objects ],
      fmodel=cmdline.fmodel,
      use_rotamers=params.use_rotamers,
      anneal=params.anneal_residues,
      out=out)
  else :
    new_hierarchy = pdb_hierarchy
  base = os.path.splitext(params.input.pdb.file_name[0])[0]
  if (params.output_model is None) :
    params.output_model = base + "_extended.pdb"
  f = open(params.output_model, "w")
  f.write(new_hierarchy.as_pdb_string(fmodel.xray_structure))
  f.close()
  print >> out, "  wrote new model to %s" % params.output_model
  if (params.output_map_coeffs is None) :
    params.output_map_coeffs = base + "_maps.mtz"
  from mmtbx.maps.utils import get_maps_from_fmodel
  import iotbx.map_tools
  two_fofc_map, fofc_map = get_maps_from_fmodel(fmodel)
  iotbx.map_tools.write_map_coeffs(
    fwt_coeffs=two_fofc_map,
    delfwt_coeffs=fofc_map,
    file_name=params.output_map_coeffs)
  print >> out, "  wrote map coefficients to %s" % params.output_map_coeffs

if (__name__ == "__main__") :
  run(sys.argv[1:])
