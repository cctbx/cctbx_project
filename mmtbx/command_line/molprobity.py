# LIBTBX_SET_DISPATCHER_NAME phenix.molprobity
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1

from __future__ import division
import mmtbx.validation.molprobity
import iotbx.phil
from libtbx import easy_pickle
from libtbx.utils import Sorry
import sys

def get_master_phil () :
  return iotbx.phil.parse("""
include scope mmtbx.utils.cmdline_input_phil_str
molprobity {
  outliers_only = True
    .type = bool
  keep_hydrogens = True
    .type = bool
    .help = '''Keep hydrogens in input file'''
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
}
""", process_includes=True)

usage_string = """\
phenix.molprobity model.pdb [data.mtz] [options ...]

Run comprehensive MolProbity validation plus R-factor calculation (if data
supplied).
"""

def run (args, out=sys.stdout) :
  import mmtbx.utils
  cmdline = mmtbx.utils.cmdline_load_pdb_and_data(
    args=args,
    master_phil=get_master_phil(),
    require_data=False,
    create_fmodel=True,
    process_pdb_file=True,
    usage_string=usage_string,
    out=out)
  params = cmdline.params
  header_info = mmtbx.validation.molprobity.pdb_header_info(
    pdb_file=params.input.pdb.file_name[0])
  if (params.output.prefix is None) :
    params.output.prefix = "molprobity"
  probe_file = None
  if (params.output.probe_dots) :
    probe_file = params.output.prefix + "_probe.txt"
  result = mmtbx.validation.molprobity.molprobity(
    pdb_hierarchy=cmdline.pdb_hierarchy,
    xray_structure=cmdline.xray_structure,
    fmodel=cmdline.fmodel,
    geometry_restraints_manager=cmdline.geometry,
    header_info=header_info,
    keep_hydrogens=params.molprobity.keep_hydrogens,
    nuclear=params.molprobity.nuclear,
    save_probe_unformatted_file=probe_file,
    min_cc_two_fofc=params.molprobity.min_cc_two_fofc)
  if (not params.output.quiet) :
    result.show(out=out, outliers_only=params.molprobity.outliers_only)
  else :
    print >> out, ""
    result.show_summary(out=out)
  print >> out, ""
  if (params.output.kinemage) :
    kin_file = "%s.kin" % params.output.prefix
    f = open(kin_file, "w")
    f.write(result.as_kinemage())
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
  if params.output.wxplots :
    try :
      import wxtbx.app
    except ImportError, e :
      raise Sorry("wxPython not available.")
    else :
      app = wxtbx.app.CCTBXApp(0)
      result.display_wx_plots()
      app.MainLoop()
  return result

if (__name__ == "__main__") :
  run(sys.argv[1:])
