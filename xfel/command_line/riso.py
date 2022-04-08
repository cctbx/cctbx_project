from __future__ import absolute_import, division, print_function
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cctbx.riso
#
# $Id: riso.py idyoung $

import iotbx.phil
from dials.util.options import ArgumentParser
from iotbx import mtz
from cctbx.array_family import flex
from libtbx.str_utils import format_value
from xfel.cxi.cxi_cc import r1_factor, scale_factor
from cctbx.crystal import symmetry
from six.moves import zip

phil_scope = iotbx.phil.parse("""
input {
  data_1 = None
    .type = path
    .help = Reference dataset for Riso calculation.
  data_2 = None
    .type = path
    .help = Comparison dataset for Riso calculation.
  labels_1 = Iobs
    .type = str
    .help = Selected column in data_1.
  labels_2 = Iobs
    .type = str
    .help = Selected column in data_2.
}
d_min = None
  .type = float
  .help = High resolution cutoff for Riso calculation.
d_max = None
  .type = float
  .help = Low resolution cutoff for Riso calculation.
anomalous_flag = False
  .type = bool
  .help = Compare anomalous datasets.
output {
  n_bins = 20
    .type = int
    .help = Number of bins for the Riso calculation.
}
""")

def riso(data_1, data_2, params, show_tables=True):
  uniform = []

  # construct a list of intensity arrays to compare between the two datasets
  for item, label in zip(
    [data_1, data_2],
    [params.input.labels_1, params.input.labels_2]):
    for array in item.as_miller_arrays():
      this_label = array.info().labels[0]
      if this_label != label: continue
      # print this_label, array.observation_type()
      uniform.append(array.as_intensity_array())
  assert len(uniform) == 2, "Could not identify the two arrays to compare. "+\
    "Please check that columns %s and %s are available in the files provided."%\
    (params.labels_1, params.labels_2)

  # if anomalous data, generate Bijvoet mates for any arrays lacking them
  if params.anomalous_flag:
    for i in (0,1):
      if not uniform[i].anomalous_flag():
        uniform[i] = uniform[i].generate_bijvoet_mates()

  # reindex
  for i in (0,1):
    uniform[i] = uniform[i].change_basis("h,k,l").map_to_asu()

  assert uniform[0].space_group_info().symbol_and_number() == \
    uniform[1].space_group_info().symbol_and_number(),\
    "Incompatible space groups between the datasets provided."

  # copy the second array with the unit cell of the first
  d_min = max(params.d_min or 0, uniform[0].d_min(), uniform[1].d_min())
  d_max = min(params.d_max or 10000, 10000)
  common_set_1 = uniform[1].customized_copy(
    crystal_symmetry = symmetry(
      unit_cell=uniform[0].unit_cell(),
      space_group_info = uniform[0].space_group_info()),
    ).resolution_filter(d_min=d_min, d_max=d_max).map_to_asu()
  common_set_2 = uniform[0].common_set(common_set_1)
  common_set_1 = uniform[1].common_set(common_set_2)
  # set 1 intentionally repeated in case of low res missing reflections
  assert len(common_set_1.indices()) == len(common_set_2.indices())
  common = (common_set_1, common_set_2)
  print("%6d indices in common in the range %.2f-%.2f Angstroms"%\
    (common_set_1.size(),d_min, d_max))

  # bin for comparison
  for array in common:
    array.setup_binner(d_min=d_min, d_max=d_max,
      n_bins=params.output.n_bins)

  # calculate scale factor and Riso
  # XXX TODO: riso_scale_factor is not set up right yet
  riso_scale_factor = scale_factor(
    common_set_2, common_set_1,
    weights=flex.pow(common_set_2.sigmas(), -2),
    use_binning=True)
  riso_binned = r1_factor(
    common_set_2, common_set_1,
    scale_factor=riso_scale_factor,
    use_binning=True)
  riso_scale_factor_all = scale_factor(
    common_set_2, common_set_1,
    weights=flex.pow(common_set_2.sigmas(), -2),
    use_binning=False)
  riso_all = r1_factor(
    common_set_2, common_set_1,
    scale_factor=riso_scale_factor_all,
    use_binning=False)

  if show_tables:
    from libtbx import table_utils
    table_header = ["","","","R"]
    table_header2 = ["Bin","Resolution Range","Completeness","iso"]
    table_data = []
    table_data.append(table_header)
    table_data.append(table_header2)

    items = riso_binned.binner.range_used()
    cumulative_counts_given = 0
    cumulative_counts_complete = 0
    for bin in items:
      table_row = []
      table_row.append("%3d"%bin)
      table_row.append("%-13s"%riso_binned.binner.bin_legend(
        i_bin=bin,show_bin_number=False,show_bin_range=False,
        show_d_range=True, show_counts=False))
      table_row.append("%13s"%riso_binned.binner.bin_legend(
        i_bin=bin,show_bin_number=False,show_bin_range=False,
        show_d_range=False, show_counts=True))
      cumulative_counts_given += riso_binned.binner._counts_given[bin]
      cumulative_counts_complete += riso_binned.binner._counts_complete[bin]
      table_row.append("%.1f%%" % (100 * riso_binned.data[bin]))
      table_data.append(table_row)

    table_row = [format_value("%3s",   "All"),
                 format_value("%-13s", "                 "),
                 format_value("%13s",  "[%d/%d]"%(cumulative_counts_given,
                                                  cumulative_counts_complete)),
                 format_value("%.1f%%", 100 * riso_all)]
    table_data.append(table_row)
    print(table_utils.format(
      table_data, has_header=2, justify='center', delim=" "))
    print("Riso is the R1 factor between the two datasets supplied.")

  return riso_binned, riso_all

def run(args):
  import os
  if ("--help" in args) or ("-h" in args) or (len(args) == 0):
    print("""cctbx.riso: a command line script for calculating an R1 factor
    between two datasets.

    Example usage:

    cctbx.riso data_1=5kaf-sf.mtz data_2=5kai-sf.mtz \\
    labels_1=Iobs labels_2=Iobs
    """)
    return
  elif ("--config" in args) or ("-c" in args):
    iotbx.phil.parse(master_phil).show(attributes_level=2)
    return
  parser = ArgumentParser(phil=phil_scope)
  params, options = parser.parse_args(show_diff_phil=True)
  # XXX TODO: find out how to auto recognize .mtz files as data_X
  extracted_data = []
  for data, label in [(params.input.data_1, params.input.labels_1),
                      (params.input.data_2, params.input.labels_2)]:
    assert data is not None, "Please supply two mtz files."
    assert os.path.exists(data), "Cannot access data: %s"%data
    extracted = mtz.object(data)
    assert extracted.has_column(label), "%s does not contain column %s"%\
      (data, label)
    extracted_data.append(extracted)
  riso(extracted_data[0], extracted_data[1], params, show_tables=True)

if __name__ == "__main__":
  import sys
  result = run(sys.argv[1:])
  if not result:
    sys.exit(1)
