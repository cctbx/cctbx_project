
# XXX this should probably be moved to iotbx/command_line at some point (once
# code for handling command-line arguments is added).  However, it is very
# similar to an existing program in the phenix tree (import_and_add_free.py),
# so there is some potential for confusion.

import libtbx.phil
from libtbx.utils import Sorry
import string
import os
import sys

master_phil = libtbx.phil.parse("""
import_data
  .short_caption = Import data and add/extend R-free flags
  .caption = This utility is used to combine new data with an existing set of \
    R-free flags, or create new flags if necessary.  All arrays in the data \
    file will be included to the output.  The reflection file \
    editor offers additional options for file manipulation.
  .style = box auto_align caption_img:icons/custom/phenix.reflection_file_editor.png
{
  data_file = None
    .type = path
    .short_caption = Data file
    .style = bold file_type:hkl input_file
  flags_file = None
    .type = path
    .short_caption = File with R-free flags
    .style = file_type:hkl input_file
  output_file = None
    .type = path
    .style = bold file_type:hkl new_file
  fraction_free = 0.05
    .type = float
    .short_caption = Fraction to flag for R-free
    .help = Typically between 0.05 and 0.1.  This option is ignored if an \
      existing set of flags is being extended.
  ignore_shells = False
    .type = bool
    .short_caption = Disable check for test sets in thin shells
    .style = noauto
}
""")

def run (args=(), params=None, out=None) :
  if (out is None) :
    out = sys.stdout
  validate_params(params)
  from iotbx import reflection_file_editor
  from iotbx import file_reader
  data_in = file_reader.any_file(params.import_data.data_file,
    force_type="hkl")
  data_in.check_file_type("hkl")
  miller_arrays = data_in.file_server.miller_arrays
  new_arrays = []
  labels = ["H","K","L"]
  warnings = []
  have_r_free = False
  for array in miller_arrays :
    if (reflection_file_editor.is_rfree_array(array, array.info())) :
      have_r_free = True
      if (params.import_data.flags_file is not None) :
        raise Sorry("The data file (%s) already contains R-free flags." %
          params.import_data.data_file)
    if (array.is_xray_reconstructed_amplitude_array()) :
      if ("F(+)" in labels) :
        labels.extend(["F_rec(+)", "SIGF_rec(+)", "F_rec(-)", "SIGF_rec(-)"])
      else :
        labels.extend(["F(+)", "SIGF(+)", "F(-)", "SIGF(-)"])
    else :
      labels.extend(array.info().labels)
    array = array.map_to_asu()
    if (not array.is_unique_set_under_symmetry()) :
      array = array.merge_equivalents().array()
    new_arrays.append(array)
  complete_set = reflection_file_editor.make_joined_set(
    new_arrays).complete_set()
  if (not have_r_free) :
    if (params.import_data.flags_file is not None) :
      flags_in = file_reader.any_file(params.import_data.flags_file,
        force_type="hkl")
      flags_in.check_file_type("hkl")
      flags_and_values = flags_in.file_server.get_r_free_flags(
        file_name=flags_in.file_name,
        label=None,
        test_flag_value=None,
        disable_suitability_test=False,
        parameter_scope=None,
        return_all_valid_arrays=True)
      if (len(flags_and_values) == 0) :
        raise Sorry("No R-free flags were found in the file %s." %
          params.import_data.flags_file)
      elif (len(flags_and_values) > 1) :
        raise Sorry(("Multiple valid sets of R-free flags were found in the "+
          "file %s.  Please use the reflection file editor to select a "+
          "single set of flags.") % params.import_data.flags_file)
      old_flags, test_flag_value = flags_and_values[0]
      labels.extend(old_flags.info().labels)
      old_flags = old_flags.map_to_asu().merge_equivalents().array()
      old_flags = old_flags.customized_copy(
        data=old_flags.data()==test_flag_value)
      missing_set = complete_set.lone_set(old_flags)
      n_missing = missing_set.indices().size()
      fraction_free = old_flags.data().count(True) / old_flags.data().size()
      if (n_missing != 0) :
        (n_bins, n_free, sse, accu) = reflection_file_editor.get_r_free_stats(
          miller_array=old_flags,
          test_flag_value=True)
        min_bins = int(old_flags.indices().size() * 0.005)
        if (n_bins < (n_free / 100)) or (sse > 0.005) or (n_bins < min_bins) :
          if (not params.import_data.ignore_shells) :
            raise Sorry(("The R-free flags in %s appear to have been "+
              "assigned in thin resolution shells.  PHENIX is unable to "+
              "properly extend flags created in this manner.  If you "+
              "prefer to ignore this check and want to create new flags using "+
              "random assignment, or if you think this message is in error, "+
              "you can use the reflection file editor instead. "+
              "(To view the distribution of R-free flags, click the "+
              "toolbar button \"Other tools\" and select \"Inspect R-free "+
              "flags\".)") % params.import_data.flags_file)
          else :
            print >> out, "WARNING: ignoring thin shells"
        if (n_missing <= 20) : # XXX hack
          from scitbx.array_family import flex
          missing_flags = missing_set.array(data=flex.bool(n_missing, False))
        else :
          missing_flags = missing_set.generate_r_free_flags(
            fraction_free=fraction_free,
            max_free=None,
            use_lattice_symmetry=True)
        new_flags = old_flags.concatenate(other=missing_flags)
      else :
        new_flags = old_flags
      new_arrays.append(new_flags)
  mtz_out = new_arrays[0].as_mtz_dataset(
    column_root_label="A")
  for i, array in enumerate(new_arrays[1:]) :
    mtz_out.add_miller_array(
      miller_array=array,
      column_root_label="%s" % string.uppercase[i+1])
  mtz_obj = mtz_out.mtz_object()
  for i, column in enumerate(mtz_obj.columns()) :
    column.set_label(labels[i])
  if (params.import_data.output_file is None) :
    base,ext = os.path.splitext(params.import_data.data_file)
    params.import_data.output_file = base + "_flags.mtz"
  mtz_obj.write(file_name=params.import_data.output_file)
  print >> out, "Data and flags written to %s" % params.import_data.output_file
  return params.import_data.output_file

def validate_params (params) :
  if (params.import_data.data_file is None) :
    raise Sorry("Please specify a data file.")
  elif (not os.path.isfile(params.import_data.data_file)) :
    raise Sorry("%s is not a recognizable file." %params.import_data.data_file)
  return True
