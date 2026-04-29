""" Export reflection file as scalepack unmerged"""
from __future__ import absolute_import, division, print_function
from libtbx.utils import Sorry
from libtbx import Auto
import os.path
import sys

master_phil_str = """
file_name = None
  .type = path
data_labels = None
  .type = str
space_group = None
  .type = space_group
ignore_merged = False
  .type = bool
batch_label = None
  .type = str
ignore_batch = False
  .type = bool
output {
  prefix = None
    .type = str
}
"""

def run(args, out=sys.stdout):
  import iotbx.phil
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil_string=master_phil_str,
    reflection_file_def="file_name",
    space_group_def="space_group")
  params = cmdline.work.extract()
  validate_params(params)
  hkl_in = cmdline.get_file(params.file_name)
  miller_arrays = hkl_in.file_object.as_miller_arrays(
    merge_equivalents=False)
  batch_numbers = None
  if (hkl_in.file_object.file_type() == "scalepack_no_merge_original_index"):
    if (not params.ignore_batch):
      batch_numbers = hkl_in.file_object.file_content().batch_numbers
  outputs = []
  for array in miller_arrays :
    labels = array.info().label_string()
    if (array.space_group() is None):
      if (params.space_group is None):
        raise Sorry("Space group needs to be explicitly specified.")
      else :
        array = array.customized_copy(
          space_group_info=params.space_group).set_info(array.info())
    elif (params.space_group is None):
      params.space_group = array.space_group_info()
    if (array.is_xray_intensity_array() and
        not array.is_unique_set_under_symmetry()):
      if (params.data_labels is None) or (labels == params.data_labels):
        outputs.append(array)
    elif (array.is_integer_array() and batch_numbers is None):
      if (not params.ignore_batch) and (params.batch_label is not None):
        if (labels == params.batch_label):
          batch_numbers = array
        elif (params.batch_label is Auto and
              "BATCH" in array.info().labels_string()):
          batch_numbers = array
  if (len(outputs) == 0):
    raise Sorry("No unmerged intensities found.")
  elif (params.space_group is None):
    raise Sorry("Space group needs to be explicitly specified.")
  for i_obs in outputs :
    print("Using intensities in %s" % i_obs.info().label_string(), file=out)
  if (batch_numbers is not None):
    if (type(batch_numbers).__name__ == "array"):
      print("Batch numbers will be taken from %s" % \
        batch_numbers.info().label_string(), file=out)
    else :
      print("Batch numbers taken from raw input file", file=out)
  if (params.output.prefix is None):
    params.output.prefix = os.path.splitext(
      os.path.basename(params.file_name))[0]
  output_file_names = []
  for i_obs in outputs :
    labels = i_obs.info().labels
    wavelength_id = None
    crystal_id = None
    for label in labels :
      if label.startswith("wavelength_id="):
        label2, w_id_str = label.split("=")
        wavelength_id = int(w_id_str)
      if label.startswith("crystal_id="):
        label2, c_id_str = label.split("=")
        crystal_id = int(c_id_str)
    file_base = params.output.prefix
    if (crystal_id is not None):
      file_base += "_c%d" % crystal_id
    if (wavelength_id is not None):
      file_base += "_w%d" % wavelength_id
    file_name = file_base + "_unmerged.sca"
    n = 1
    while file_name in output_file_names :
      file_name = file_base + "-" + str(n) + "_unmerged.sca"
      n += 1
    output_file_names.append(file_name)
    if (i_obs.space_group() is None):
      i_obs = i_obs.customized_copy(space_group_info=params.space_group)
    tmp_batch_numbers = batch_numbers
    if (batch_numbers is not None):
      mismatch = False
      if (type(batch_numbers).__name__ != "array"):
        if (len(batch_numbers) != len(i_obs.indices())):
          mismatch = True
      elif (batch_numbers.indices().all_eq(i_obs.indices())):
        mismatch = True
      if (mismatch):
        msg = ("The h,k,l indices for the batch number array (%s) " +
          "do not match the h,k,l indices for the intensity array %s.  ") % \
          (batch_numbers.info().label_string(),
           i_obs.info().label_string())
        if (not params.ignore_batch):
          raise Sorry(msg +
            "You can suppress this warning by specifying ignore_batch=True.")
        else :
          print(msg + "The batch numbers will not be output.", file=out)
          tmp_batch_numbers = None
    i_obs.export_as_scalepack_unmerged(
      file_name=file_name,
      batch_numbers=tmp_batch_numbers)
    print("Wrote %s" % file_name, file=out)
  return output_file_names

def validate_params(params):
  if (params.file_name is None):
    raise Sorry("Input file name not specified.")

if (__name__ == "__main__"):
  run(sys.argv[1:])
