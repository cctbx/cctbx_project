# LIBTBX_SET_DISPATCHER_NAME phenix.rotalyze
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1

from __future__ import division
import mmtbx.validation.rotalyze
import iotbx.phil
from libtbx.utils import Usage
import sys

def get_master_phil():
  return iotbx.phil.parse(input_string="""
    include scope mmtbx.validation.molprobity_cmdline_phil_str
    data_version = 8000
      .type = str
      .help = '''Choose which dataset to use, Top \'8000\' or \'500\' '''
    show_errors = False
      .type = bool
      .help = '''Print out errors'''
    wxplot = False
      .type = bool
      .help = Display interactive plots (requires wxPython and Matplotlib)
""", process_includes=True)

usage_string = """\
phenix.rotalyze file.pdb [params.eff] [options ...]

Options:

  model=input_file        input PDB file
  outliers_only=False   only print outliers
  verbose=False         verbose text output
  data_version=8000     Give \'500\' to use the Top500 data and
                        \'8000\' to use the Top8000 data.

Example:

  phenix.rotalyze model=1ubq.pdb outliers_only=True
"""

def run (args, out=sys.stdout, quiet=False) :
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=get_master_phil(),
    pdb_file_def="model",
    usage_string=usage_string)
  params = cmdline.work.extract()
  if (params.model is None) :
    raise Usage(usage_string)
  pdb_in = cmdline.get_file(params.model, force_type="pdb")
  hierarchy = pdb_in.file_object.hierarchy
  result = mmtbx.validation.rotalyze.rotalyze(
    pdb_hierarchy=hierarchy,
    data_version=params.data_version,
    show_errors=params.show_errors,
    outliers_only=params.outliers_only,
    out=out,
    quiet=quiet)
  if params.verbose:
    result.show_old_output(out=out, verbose=True)
  if params.wxplot :
    try :
      import wxtbx.app
    except ImportError, e :
      raise Sorry("wxPython not available.")
    else :
      app = wxtbx.app.CCTBXApp(0)
      result.display_wx_plots()
      app.MainLoop()

if (__name__ == "__main__") :
  run(sys.argv[1:])
