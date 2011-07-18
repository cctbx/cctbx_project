from libtbx import easy_run
from iotbx.cns.space_group_symbols import cns_format
import os, shutil, time

def run_cns_density_modification(params, fo, hl_coeffs):
  cur_dir = os.path.abspath(os.path.curdir)
  if os.path.exists("tmp_cns"):
    shutil.rmtree("tmp_cns")
  os.mkdir("tmp_cns")
  os.chdir(os.path.abspath("tmp_cns"))
  try:
    result = do_dirty_work(params, fo, hl_coeffs)
  finally:
    # make sure we always end up back here no matter what happens
    os.chdir(cur_dir)
  result.show_stdout()

def do_dirty_work(params, fo, hl_coeffs):
  fo.export_as_cns_hkl(
    file_object=open("f_obs.hkl", "wb"),
    file_name="f_obs.hkl")
  hl_coeffs.export_as_cns_hkl(
    file_object=open("hl_coeffs.hkl", "wb"),
    file_name="hl_coeffs.hkl")

  cns_solve_dir = os.environ.get("CNS_SOLVE")
  if cns_solve_dir is None:
    raise RuntimeError("Environment variable CNS_SOLVE is not defined")

  shutil.copyfile(
    "%s/inputs/xtal_phase/density_modify.inp" % cns_solve_dir, "tmp.inp")

  cns_params = {}
  cns_params.setdefault("space_group", cns_format(
    space_group_info=params.input.space_group))
  uc_params = params.input.unit_cell.parameters()
  cns_params["a"] = uc_params[0]
  cns_params["b"] = uc_params[1]
  cns_params["c"] = uc_params[2]
  cns_params["alpha"] = uc_params[3]
  cns_params["beta"] = uc_params[4]
  cns_params["gamma"] = uc_params[5]

  if params.d_min is None:
    if params.phase_extension:
      params.d_min = fo.d_min()
    else:
      params.d_min = hl_coeffs_start.d_min()
  if params.initial_d_min is None:
    params.initial_d_min = params.d_min
  assert params.initial_d_min >= params.d_min
  cns_params["d_min"] = params.d_min
  cns_params["initial_d_min"] = params.initial_d_min
  cns_params["phase_extend"] = str(params.phase_extension).lower()
  cns_params["prot_to_solv"] = params.protein_solvent_ratio
  cns_params["trunc_min"] = params.density_truncation.fraction_min
  cns_params["trunc_max"] = params.density_truncation.fraction_max
  cns_params["trunc_max"] = 1
  cns_params["solcon"] = params.solvent_fraction
  if params.solvent_mask.averaging_radius.final is None:
    params.solvent_mask.averaging_radius.final = params.d_min
  if params.solvent_mask.averaging_radius.initial is None:
    params.solvent_mask.averaging_radius.initial \
          = params.solvent_mask.averaging_radius.final + 1
  cns_params["start_ave_radius"] = params.solvent_mask.averaging_radius.initial
  cns_params["finish_ave_radius"] = params.solvent_mask.averaging_radius.final
  cns_params["initial_steps"] = params.initial_steps
  cns_params["shrink_steps"] = params.shrink_steps
  cns_params["final_steps"] = params.final_steps
  cns_params["grid_resolution_factor"] = params.grid_resolution_factor

  params_file = open("params.inp", "wb")
  params_file.write(cns_density_modify_inp_template % cns_params)
  params_file.close()

  cmd = "%s/bin/cns_transfer -def params.inp -inp tmp.inp -out density_modify.inp" %(
    cns_solve_dir)
  print cmd
  result = easy_run.fully_buffered(command=cmd).raise_if_errors_or_output()

  cmd = "cns_solve < density_modify.inp > density_modify.out"
  print cmd
  t0 = time.time()
  result = easy_run.fully_buffered(command=cmd).raise_if_errors()
  print "CNS time: %.2f" %(time.time()-t0)
  return result

cns_density_modify_inp_template = """\
{- begin block parameter definition -} define(
{===>} sg="%(space_group)s";

{===>} a=%(a)f;
{===>} b=%(b)f;
{===>} c=%(c)f;
{===>} alpha=%(alpha)f;
{===>} beta=%(beta)f;
{===>} gamma=%(gamma)f;

{===>} reflection_infile_1="f_obs.hkl";
{===>} reflection_infile_2="hl_coeffs.hkl";

{===>} obs_f="FOBS";
{===>} obs_sigf="SIGMA";
{===>} obs_pa="PA";
{===>} obs_pb="PB";
{===>} obs_pc="PC";
{===>} obs_pd="PD";
{===>} low_res=500.0;
{===>} high_res=%(d_min)f;

{================== non-crystallographic symmetry ====================}

{===>} averaging=false;
{===>} ncs_mask_infile="";
{===>} ncs_infile="";

{======================= density truncation ==========================}

{===>} truncation=true;
{===>} prot_to_solv=%(prot_to_solv)f;
{===>} trunc_min=%(trunc_min)f;
{===>} trunc_max=%(trunc_max)f;

{====================== solvent modification =========================}

{===>} solvent_mod=true;
{===>} solvent_method="flip";
{===>} scale_flip=true;

{========================== solvent adjust ===========================}

{===>} solvent_adjust=true;

{=========================== solvent mask ============================}

{===>} mask_sol_type="rms";
{===>} solvent_mask_infile="";
{===>} mask_complete=true;
{===>} solcon=%(solcon)f;
{===>} start_ave_radius=%(start_ave_radius)f;
{===>} finish_ave_radius=%(finish_ave_radius)f;

{======================= modification scheme =========================}

{===>} initial_steps=%(initial_steps)i;
{===>} shrink_steps=%(shrink_steps)i;
{===>} final_steps=%(final_steps)i;

{========================= phase extension ===========================}

{===>} phase_extend=%(phase_extend)s;
{===>} initial_highres=%(initial_d_min)f;

{===================== modification parameters =======================}

{===>} map_mode="unit_cell";

{========================= output arrays =============================}

{===>} out_f="full_mod";
{===>} out_fom="fom_mod";
{===>} out_pa="pa_mod";
{===>} out_pb="pb_mod";
{===>} out_pc="pc_mod";
{===>} out_pd="pd_mod";

{=========================== output files ============================}

{===>} output_root="density_modify";
{===>} write_map=true;
{===>} write_mask=false;

 ) {- end block parameter definition -}
"""
