# LIBTBX_SET_DISPATCHER_NAME phenix.molprobity
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1

from __future__ import division
from libtbx.program_utils.result import program_result
from libtbx.utils import Sorry, multi_out
from libtbx import runtime_utils
from libtbx import easy_pickle
from libtbx import Auto
import libtbx.load_env
import os.path
import sys

def get_master_phil () :
  from mmtbx.command_line import generate_master_phil_with_inputs
  return generate_master_phil_with_inputs(
    enable_automatic_twin_detection=True,
    enable_pdb_interpretation_params=True,
    enable_stop_for_unknowns=False,
    enable_unmerged_data=True,
    enable_cdl=Auto,
    phil_string="""
molprobity {
  outliers_only = True
    .type = bool
  keep_hydrogens = Auto
    .type = bool
    .help = Keep hydrogens in input file (instead of re-generating them with \
      Reduce).  If set to Auto, the behavior will depend on whether the \
      neutron scattering table is used (regardless of whether we actually \
      have experimental data).
  nuclear = False
    .type = bool
    .short_caption = "Use nuclear hydrogen positions"
  min_cc_two_fofc = 0.8
    .type = float
  n_bins = 10
    .type = int
    .short_caption = Number of resolution bins
  use_pdb_header_resolution_cutoffs = False
    .type = bool
    .short_caption = Use resolution cutoffs in PDB header
  count_anomalous_pairs_separately = False
    .type = bool
    .expert_level = 2
  rotamer_library = *500 8000
    .type = choice
    .help = Library of rotamer probabilities (Top500 or Top8000)
    .expert_level = 2
  flags
    .expert_level = 3
  {
    include scope mmtbx.validation.molprobity.master_phil_str
  }
  ligand_selection = None
    .type = atom_selection
    .expert_level = 3
}
output {
  quiet = False
    .type = bool
  probe_dots = True
    .type = bool
    .short_caption = Save Probe dots for Coot
  kinemage = False
    .type = bool
    .short_caption = Save Kinemage file for KiNG
  percentiles = False
    .type = bool
    .help = Show percentile rankings for summary statistics
  coot = True
    .type = bool
    .help = Write Coot script
  maps = Auto
    .type = bool
    .short_caption = Save map coefficients
    .help = Write map coefficients (if experimental data supplied)
  map_options
    .short_caption = Advanced options for map coefficients
  {
    fill_missing_f_obs = True
      .type = bool
    exclude_free_r_reflections = False
      .type = bool
  }
  prefix = None
    .type = str
    .style = hidden
  pickle = False
    .type = bool
    .style = hidden
  wxplots = False
    .type = bool
    .help = Display plots in wxPython
    .style = hidden
  gui_dir = None
    .type = path
    .short_caption = Output directory
    .help = Output directory (Phenix GUI only).
    .style = output_dir
  include scope libtbx.phil.interface.tracking_params
}
""")

usage_string = """\
phenix.molprobity model.pdb [data.mtz] [options ...]

Run comprehensive MolProbity validation plus R-factor calculation (if data
supplied).
"""

def run (args,
    out=sys.stdout,
    program_name="phenix.molprobity",
    ignore_missing_modules=False,
    return_input_objects=False) : # for testing
  rotarama_dir = libtbx.env.find_in_repositories(
    relative_path="chem_data/rotarama_data",
    test=os.path.isdir)
  if (rotarama_dir is None) :
    raise ImportError("Rotamer and Ramachandran distributions not available; "+
      "you will need these to run MolProbity.")
  elif (((not libtbx.env.has_module("reduce")) or
         (not libtbx.env.has_module("probe"))) and
         (not ignore_missing_modules)) :
    raise ImportError("Reduce and/or Probe not configured.")
  import mmtbx.validation.molprobity
  import mmtbx.command_line
  cmdline = mmtbx.command_line.load_model_and_data(
    args=args,
    master_phil=get_master_phil(),
    require_data=False,
    create_fmodel=True,
    process_pdb_file=True,
    usage_string=usage_string,
    prefer_anomalous=True,
    out=out)
  params = cmdline.params
  fmodel = cmdline.fmodel
  if (params.output.maps is Auto) and (fmodel is not None) :
    params.output.maps = True
  elif (params.output.maps == True) and (fmodel is None) :
    raise Sorry("Map output requires experimental data.")
  if (params.molprobity.keep_hydrogens is Auto) :
    params.molprobity.keep_hydrogens = \
      (params.input.scattering_table == "neutron")
  header_info = mmtbx.validation.molprobity.pdb_header_info(
    pdb_file=params.input.pdb.file_name[0],
    pdb_hierarchy=cmdline.pdb_hierarchy)
  pdb_prefix = os.path.splitext(os.path.basename(
    params.input.pdb.file_name[0]))[0]
  if (params.output.prefix is None) :
    params.output.prefix = "molprobity"
  probe_file = None
  if (params.output.probe_dots) or (params.output.kinemage) :
    probe_file = params.output.prefix + "_probe.txt"
  raw_data = cmdline.raw_data
  validation = mmtbx.validation.molprobity.molprobity(
    pdb_hierarchy=cmdline.pdb_hierarchy,
    xray_structure=cmdline.xray_structure,
    fmodel=fmodel,
    flags=params.molprobity.flags,
    crystal_symmetry=cmdline.crystal_symmetry,
    geometry_restraints_manager=cmdline.geometry,
    raw_data=cmdline.raw_data,
    unmerged_data=cmdline.unmerged_i_obs,
    all_chain_proxies=cmdline.processed_pdb_file.all_chain_proxies,
    header_info=header_info,
    keep_hydrogens=params.molprobity.keep_hydrogens,
    nuclear=params.molprobity.nuclear,
    save_probe_unformatted_file=probe_file,
    min_cc_two_fofc=params.molprobity.min_cc_two_fofc,
    n_bins_data=params.molprobity.n_bins,
    outliers_only=params.molprobity.outliers_only,
    use_pdb_header_resolution_cutoffs=\
      params.molprobity.use_pdb_header_resolution_cutoffs,
    count_anomalous_pairs_separately=\
      params.molprobity.count_anomalous_pairs_separately,
    file_name=params.input.pdb.file_name[0],
    ligand_selection=params.molprobity.ligand_selection,
    rotamer_library=params.molprobity.rotamer_library)
  map_file = None
  if (not params.output.quiet) :
    out2 = multi_out()
    out2.register("stdout", out)
    f = open(params.output.prefix + ".out", "w")
    out2.register("txt_out", f)
    validation.show(out=out2,
      outliers_only=params.molprobity.outliers_only,
      show_percentiles=params.output.percentiles)
    f.close()
    print >> out, ""
    print >> out, "Results written to %s.out" % params.output.prefix
    if (params.output.kinemage) :
      if (cmdline.pdb_hierarchy.models_size() == 1) :
        assert (probe_file is not None)
        import mmtbx.kinemage.validation
        kin_file = "%s.kin" % params.output.prefix
        kin_out = \
          mmtbx.kinemage.validation.export_molprobity_result_as_kinemage(
            result=validation,
            pdb_hierarchy=cmdline.pdb_hierarchy,
            geometry=cmdline.geometry,
            probe_file=probe_file,
            keep_hydrogens=params.molprobity.keep_hydrogens,
            pdbID=pdb_prefix)
        f = open(kin_file, "w")
        f.write(kin_out)
        f.close()
        if (not params.output.quiet) :
          print >> out, "Wrote kinemage to %s" % kin_file
      else :
        print >> out, "Kinemage output not available for multiple MODELs."
    if (params.output.pickle) :
      easy_pickle.dump("%s.pkl" % params.output.prefix, validation)
      if (not params.output.quiet) :
        print >> out, "Saved result to %s.pkl" % params.output.prefix
    if (params.output.coot) :
      coot_file = "%s_coot.py" % params.output.prefix
      validation.write_coot_script(coot_file)
      if (not params.output.quiet) :
        print >> out, "Wrote script for Coot: %s" % coot_file
    if (params.output.maps == True) :
      import mmtbx.maps.utils
      import iotbx.map_tools
      map_file = "%s_maps.mtz" % params.output.prefix
      two_fofc_map, fofc_map = mmtbx.maps.utils.get_maps_from_fmodel(
        fmodel=fmodel,
        fill_missing_f_obs=params.output.map_options.fill_missing_f_obs,
        exclude_free_r_reflections=\
          params.output.map_options.exclude_free_r_reflections)
      anom_map = None
      if (fmodel.f_obs().anomalous_flag()) :
        anom_map = mmtbx.maps.utils.get_anomalous_map(fmodel)
      iotbx.map_tools.write_map_coeffs(
        file_name=map_file,
        fwt_coeffs=two_fofc_map,
        delfwt_coeffs=fofc_map,
        anom_coeffs=anom_map)
      print >> out, "Wrote map coefficients to %s" % map_file
  else :
    print >> out, ""
    validation.show_summary(out=out, show_percentiles=params.output.percentiles)
  if (params.output.wxplots) :
    try :
      import wxtbx.app
    except ImportError, e :
      raise Sorry("wxPython not available.")
    else :
      app = wxtbx.app.CCTBXApp(0)
      validation.display_wx_plots()
      app.MainLoop()
  if (return_input_objects) :
    return validation, cmdline
  return result(
    program_name="phenix.molprobity",
    job_title=params.output.job_title,
    directory=os.getcwd(),
    map_file=map_file,
    other_result=validation,
    other_files=[ params.input.pdb.file_name[0] ])

class launcher (runtime_utils.target_with_save_result) :
  def run (self) :
    os.mkdir(self.output_dir)
    os.chdir(self.output_dir)
    return run(args=self.args, out=sys.stdout)

class result (program_result) :
  """
  Wrapper object for Phenix GUI.
  """
  @property
  def validation (self) :
    return self.other_result

  @property
  def pdb_file (self) :
    return self.other_files[0]

  def get_final_stats (self) :
    return self.validation.get_statistics_for_phenix_gui()

if (__name__ == "__main__") :
  run(sys.argv[1:])
