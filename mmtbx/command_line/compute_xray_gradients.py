
# TODO: tests

import mmtbx.f_model
import mmtbx.utils
import iotbx.phil
from scitbx.array_family import flex
from libtbx.str_utils import make_header
from libtbx.utils import Usage
import sys
import os

master_phil = iotbx.phil.parse("""
include scope mmtbx.utils.cmdline_input_phil_str
output {
  prefix = gradients
    .type = str
  serial = 1
    .type = int
  output_dir = None
    .type = path
}
options {
  bulk_solvent_and_scale = True
    .type = bool
  force_bss = False
    .type = bool
  remove_outliers = True
    .type = bool
}
include scope mmtbx.command_line.fmodel.fmodel_from_xray_structure_params
bss {
  include scope mmtbx.bulk_solvent.bulk_solvent_and_scaling.master_params
}
""", process_includes=True)

def compute_gradients (params,
                       xray_structure,
                       pdb_hierarchy,
                       f_obs,
                       r_free_flags,
                       out=sys.stdout) :
  params.bss.apply_back_trace_of_b_cart=False
  if (not params.options.bulk_solvent_and_scale) :
    assert (not None in [params.fmodel.k_sol,params.fmodel.b_sol,
            params.fmodel.b_cart])
  fmodel = mmtbx.f_model.manager(
    xray_structure               = xray_structure,
    sf_and_grads_accuracy_params = params.structure_factors_accuracy,
    r_free_flags                 = r_free_flags,
    mask_params                  = params.mask,
    f_obs                        = f_obs,
    target_name                  = "ml",
    k_sol                        = params.fmodel.k_sol,
    b_sol                        = params.fmodel.b_sol,
    b_cart                       = params.fmodel.b_cart)
  #fmodel_info = fmodel.info()
  #fmodel_info.show_rfactors_targets_scales_overall(out=out)
  if params.options.bulk_solvent_and_scale :
    make_header("Bulk solvent and scaling", out=out)
    fmodel.update_solvent_and_scale(params=params.bss)
    fmodel_info = fmodel.info()
    fmodel_info.show_rfactors_targets_scales_overall(out=out)
  else :
    make_header("Update X-ray structure", out=out)
    fmodel.apply_back_b_iso()
    fmodel.update_xray_structure(update_f_mask=True)
  if params.options.remove_outliers :
    print >> out, ""
    fmodel = fmodel.remove_outliers(show=True, log=out)
  fmodel.xray_structure.scatterers().flags_set_grads(state=False)
  selection = flex.bool(xray_structure.sites_cart().size(), True)
  #print selection.iselection().size()
  fmodel.xray_structure.scatterers().flags_set_grad_site(
    iselection=selection.iselection())
  make_header("Calculating gradients", out=out)
  tf = fmodel.target_functor()
  tf.prepare_for_minimization()
  tfx_r = tf(compute_gradients=True)
  target_work = tfx_r.target_work()
  #print target_work
  sf = tfx_r.gradients_wrt_atomic_parameters().packed()
  assert (sf.size() == (pdb_hierarchy.atoms().size() * 3))
  xyz_grads = flex.vec3_double(sf)
  info = fmodel.info()
  info.show_targets(out=out, text="Target values")
  params.fmodel.b_cart = fmodel.b_cart()
  params.fmodel.k_sol = fmodel.k_sol()
  params.fmodel.b_sol = fmodel.b_sol()
  if (params.output.output_dir is None) :
    params.output.output_dir = os.getcwd()
  file_base = os.path.join(params.output.output_dir,
    "%s_%d" % (params.output.prefix, params.output.serial))
  params.output.serial += 1
  if (not params.options.force_bss) :
    params.options.bulk_solvent_and_scale = False
  f = open("%s.eff" % file_base, "w")
  final_phil = master_phil.format(python_object=params)
  final_phil.show(out=f)
  f.close()
  pdb_out = open("%s.pdb" % file_base, "w")
  pdb_out.write(pdb_hierarchy.as_pdb_string())
  pdb_out.close()
  grads_out = open("%s.xyz" % file_base, "w")
  print >> grads_out, "# TARGET = %.12g" % target_work
  for site_grads in xyz_grads :
    print >> grads_out, "%.12g %.12g %.12g" % site_grads
  grads_out.close()
  print >> out, "PDB file:"
  print >> out, "  %s.pdb" % file_base
  print >> out, "Gradients:"
  print >> out, "  %s.xyz" % file_base
  print >> out, "Parameters for next run:"
  print >> out, "  %s.eff" % file_base

def run (args, out=sys.stdout) :
  if (len(args) == 0) :
    raise Usage("""
phenix.compute_xray_gradients [model.pdb] [data.mtz] [params.eff] [options...]

This program will write out a parameter (.eff) file when it is finished,
which can be used for the next round without bulk-solvent correction.
""")
  cmdline = mmtbx.utils.cmdline_load_pdb_and_data(
    args=args,
    master_phil=master_phil,
    out=out,
    process_pdb_file=False,
    create_fmodel=False)
  params = cmdline.params
  f_obs = cmdline.f_obs
  r_free_flags = cmdline.r_free_flags
  xray_structure = cmdline.xray_structure
  compute_gradients(
    params=cmdline.params,
    xray_structure=cmdline.xray_structure,
    pdb_hierarchy=cmdline.pdb_hierarchy,
    f_obs=cmdline.f_obs,
    r_free_flags=cmdline.r_free_flags,
    out=out)

if __name__ == "__main__" :
  run(sys.argv[1:])
