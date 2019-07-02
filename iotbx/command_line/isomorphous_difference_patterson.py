# LIBTBX_SET_DISPATCHER_NAME cctbx.isomorphous_difference_patterson

from __future__ import absolute_import, division, print_function
from libtbx.utils import Sorry, Usage, show_development_warning
import libtbx.callbacks # import dependency
import libtbx.phil
from six.moves import cStringIO as StringIO
import os
import sys

master_phil = libtbx.phil.parse("""
data_file_1 = None
  .type = path
data_file_2 = None
  .type = path
labels_1 = None
  .type = str
labels_2 = None
  .type = str
include scope iotbx.command_line.patterson_map.patterson_map_phil_str
relative_diff_limit = 0.1
  .type = float
map_file_name = iso_diff.ccp4
  .type = path
""", process_includes=True)

def run(args, out=sys.stdout):
  show_development_warning(out=out)
  if (len(args) == 0) or ("--help" in args):
    raise Usage("""\
cctbx.isomorphous_difference_patterson xtal1.mtz xtal2.mtz [options]

Calculate an isomorphous difference Patterson map from two datasets.  You may
specify both in a single MTZ file if you want, but you must manually specify
the data_file_1 and data_file_2 parameters.  Output is a CCP4-format map.

Full parameters:
%s
""" % master_phil.as_str(prefix="  "))
  from iotbx.command_line import patterson_map
  from iotbx import file_reader
  import iotbx.symmetry
  sources = []
  data_file_1 = data_file_2 = None
  for arg in args :
    if (os.path.isfile(arg)):
      f = file_reader.any_file(arg)
      if (f.file_type == "phil"):
        sources.append(f.file_object)
      elif (f.file_type == "hkl"):
        if (data_file_1 is None):
          data_file_1 = f
        elif (data_file_2 is None):
          data_file_2 = f
        else :
          raise Sorry("No more than two reflection files are supported.")
      else :
        raise Sorry("The file format for '%s' is not a supported input." % arg)
    else :
      try :
        sources.append(libtbx.phil.parse(arg))
      except RuntimeError as e :
        raise Sorry("Unrecognized argument '%s'.")
  params = master_phil.fetch(sources=sources).extract()
  if (params.data_file_1 is not None):
    data_file_1 = file_reader.any_file(params.data_file_1)
  if (params.data_file_2 is not None):
    data_file_2 = file_reader.any_file(params.data_file_2)
  if (None in [data_file_1, data_file_2]):
    raise Sorry("You must provide two datasets for this program to run.  If "+
      "they are different arrays within a specific file, you need to specify "+
      "the parameters explicitly, i.e. 'data_file_1=NAME data_file_2=NAME'.")
  symm_manager = iotbx.symmetry.manager()
  sg_err_1, uc_err_1 = symm_manager.process_reflections_file(data_file_1)
  sg_err_2, uc_err_2 = symm_manager.process_reflections_file(data_file_2)
  out_tmp = StringIO()
  symm_manager.show(out=out_tmp)
  if (sg_err_1) or (sg_err_2):
    raise Sorry(("Incompatible space groups in input files:\n%s\nAll files "+
      "must have the same point group (and ideally the same space group). ") %
      out_tmp.getvalue())
  elif (uc_err_1) or (uc_err_2):
    libtbx.call_back(message="warn",
      data=("Crystal symmetry mismatch:\n%s\nCalculations will continue "+
        "using the symmetry in the PDB file first reflection file, but the "+
        "maps should be treated with some suspicion.") % out_tmp.getvalue())
  obs_1 = patterson_map.extract_data(
    file_name=None,
    hkl_in=data_file_1,
    expected_labels=params.labels_1,
    out=out)
  obs_2 = patterson_map.extract_data(
    file_name=None,
    hkl_in=data_file_2,
    expected_labels=params.labels_2,
    out=out)
  assert (not None in [obs_1, obs_2])
  print("", file=out)
  print("Data array #1 (will use symmetry from here):", file=out)
  obs_1.show_summary(f=out, prefix="  ")
  print("", file=out)
  print("Data array #2:", file=out)
  obs_2.show_summary(f=out, prefix="  ")
  print("", file=out)
  obs_2 = obs_2.customized_copy(crystal_symmetry=obs_1)
  f_obs_1 = patterson_map.prepare_f_obs(obs_1, params).average_bijvoet_mates()
  f_obs_2 = patterson_map.prepare_f_obs(obs_2, params).average_bijvoet_mates()
  f_obs_1, f_obs_2 = f_obs_1.common_sets(other=f_obs_2)
  # XXX will this work as expected?
  print("Scaling second dataset...", file=out)
  f_obs_2 = f_obs_1.multiscale(other=f_obs_2)
  from scitbx.array_family import flex
  print("  max(F1) = %.2f  max(F2) = %.2f" % (flex.max(f_obs_1.data()),
    flex.max(f_obs_2.data())), file=out)
  delta_f = f_obs_1.customized_copy(data=f_obs_1.data()-f_obs_2.data())
  if (params.relative_diff_limit is not None):
    f_obs_1_non_zero = f_obs_1.select(f_obs_1.data() > 0)
    delta_f_non_zero = delta_f.common_set(other=f_obs_1_non_zero)
    relative_diff = delta_f_non_zero.customized_copy(
      data=delta_f_non_zero.data()/f_obs_1_non_zero.data())
    delta_f = delta_f_non_zero.select(
      relative_diff.data() < params.relative_diff_limit)
    n_discarded = len(delta_f_non_zero.data()) - len(delta_f.data())
    if (n_discarded > 0):
      print("Discarded %d reflections with excessive differences." % \
        n_discarded, file=out)
  print("  max(delta_F) = %.2f" % flex.max(delta_f.data()), file=out)
  map = patterson_map.calculate_patterson_map(data=delta_f, params=params,
    normalize=False)
  if (params.map_file_name is None):
    params.map_file_name = "iso_diff.ccp4"
  map.as_ccp4_map(
    file_name=params.map_file_name)
  print("Wrote %s" % params.map_file_name, file=out)

if (__name__ == "__main__"):
  run(sys.argv[1:])
