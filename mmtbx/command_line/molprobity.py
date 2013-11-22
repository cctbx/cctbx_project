# LIBTBX_SET_DISPATCHER_NAME phenix.molprobity
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1

from __future__ import division
from libtbx.utils import Sorry, multi_out
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
    .help = '''Use nuclear hydrogen positions'''
  min_cc_two_fofc = 0.8
    .type = float
}
output {
  quiet = False
    .type = bool
  probe_dots = True
    .type = bool
  kinemage = False
    .type = bool
  coot = True
    .type = bool
    .help = Write Coot script
  maps = Auto
    .type = bool
    .help = Write maps (if experimental data supplied)
  prefix = None
    .type = str
  pickle = False
    .type = bool
    .style = hidden
  wxplots = False
    .type = bool
    .help = Display plots in wxPython
    .style = hidden
}
pdb_interpretation {
  include scope mmtbx.monomer_library.pdb_interpretation.master_params
  stop_for_unknowns = True
    .type = bool
}""")

usage_string = """\
phenix.molprobity model.pdb [data.mtz] [options ...]

Run comprehensive MolProbity validation plus R-factor calculation (if data
supplied).
"""

def run (args, out=sys.stdout, return_model_fmodel_objects=False) :
  rotarama_dir = libtbx.env.find_in_repositories(
    relative_path="chem_data/rotarama_data",
    test=os.path.isdir)
  if (rotarama_dir is None) :
    raise ImportError("Rotamer and Ramachandran distributions not available; "+
      "you will need these to run MolProbity.")
  elif ((not libtbx.env.has_module("reduce")) or
        (not libtbx.env.has_module("probe"))) :
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
    use_conformation_dependent_library=True,
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
    pdb_file=params.input.pdb.file_name[0])
  pdb_prefix = os.path.splitext(os.path.basename(
    params.input.pdb.file_name[0]))[0]
  if (params.output.prefix is None) :
    params.output.prefix = "molprobity"
  probe_file = None
  if (params.output.probe_dots) or (params.output.kinemage) :
    probe_file = params.output.prefix + "_probe.txt"
  result = mmtbx.validation.molprobity.molprobity(
    pdb_hierarchy=cmdline.pdb_hierarchy,
    xray_structure=cmdline.xray_structure,
    fmodel=fmodel,
    crystal_symmetry=cmdline.crystal_symmetry,
    geometry_restraints_manager=cmdline.geometry,
    header_info=header_info,
    keep_hydrogens=params.molprobity.keep_hydrogens,
    nuclear=params.molprobity.nuclear,
    save_probe_unformatted_file=probe_file,
    min_cc_two_fofc=params.molprobity.min_cc_two_fofc)
  if (not params.output.quiet) :
    out2 = multi_out()
    out2.register("stdout", out)
    f = open(params.output.prefix + ".out", "w")
    out2.register("txt_out", f)
    result.show(out=out2, outliers_only=params.molprobity.outliers_only)
    f.close()
    print >> out, ""
    print >> out, "Results written to %s.out" % params.output.prefix
    if (params.output.kinemage) :
      assert (probe_file is not None)
      import mmtbx.kinemage.validation
      kin_file = "%s.kin" % params.output.prefix
      kin_out = mmtbx.kinemage.validation.export_molprobity_result_as_kinemage(
        result=result,
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
    if (params.output.pickle) :
      easy_pickle.dump("%s.pkl" % params.output.prefix, result)
      if (not params.output.quiet) :
        print >> out, "Saved result to %s.pkl" % params.output.prefix
    if (params.output.coot) :
      coot_file = "%s_coot.py" % params.output.prefix
      result.write_coot_script(coot_file)
      if (not params.output.quiet) :
        print >> out, "Wrote script for Coot: %s" % coot_file
    if (params.output.maps) :
      pass
  else :
    print >> out, ""
    result.show_summary(out=out)
  if (params.output.wxplots) :
    try :
      import wxtbx.app
    except ImportError, e :
      raise Sorry("wxPython not available.")
    else :
      app = wxtbx.app.CCTBXApp(0)
      result.display_wx_plots()
      app.MainLoop()
  if (return_model_fmodel_objects) :
    return result, cmdline.pdb_hierarchy, cmdline.fmodel
  return result

if (__name__ == "__main__") :
  run(sys.argv[1:])
