# LIBTBX_SET_DISPATCHER_NAME phenix.reindex

import libtbx.phil
from libtbx.utils import Sorry
import string
import os
import sys

master_phil = libtbx.phil.parse("""
reindex
  .caption = This utility will reindex the contents of an MTZ file \
    using a given change-of-basis operator.  For instance, to convert a \
    file with symmetry P22121 to P21212, the operator c,a,b could be used. \
    (If your reflections are in a format other than MTZ, you can convert them \
    in the reflection file editor, which also supports change-of-basis \
    operations.)
  .style = box auto_align caption_img:icons/custom/phenix.reflection_file_editor.png
{
  hkl_file = None
    .type = path
    .short_caption = MTZ file
    .style = bold file_type:hkl input_file
  output_file = None
    .type = path
    .short_caption = Output file
    .style = new_file file_type:hkl
  change_of_basis = None
    .type = str
    .short_caption = Change of basis
    .input_size=100
    .style = bold
}""")

def run (args=(), params=None, out=None) :
  if (out is None) :
    out = sys.stdout
  if (params is None) :
    if (len(args) == 0) :
      raise Usage("""
phenix.reindex data.mtz change_of_basis=<operator>

Change-of-basis operator: h,k,l or x,y,z or
                          to_reference_setting, to_primitive_setting,
                          to_niggli_cell, to_inverse_hand
""")
    import iotbx.phil
    cmdline = iotbx.phil.process_command_line_with_files(
      args=args,
      master_phil=master_phil,
      reflection_file_def="reindex.hkl_file")
    params = cmdline.work.extract()
  validate_params(params)
  cb_op = convert_operator(params.reindex.change_of_basis)
  from iotbx import file_reader
  hkl_in = file_reader.any_file(params.reindex.hkl_file)
  miller_arrays = hkl_in.file_server.miller_arrays
  new_arrays = []
  labels = ["H","K","L"]
  warnings = []
  for array in miller_arrays :
    if (array.is_xray_reconstructed_amplitude_array()) :
      if ("F(+)" in labels) :
        labels.extend(["F_rec(+)", "SIGF_rec(+)", "F_rec(-)", "SIGF_rec(-)"])
      else :
        labels.extend(["F(+)", "SIGF(+)", "F(-)", "SIGF(-)"])
    else :
      labels.extend(array.info().labels)
    array = array.change_basis(cb_op=cb_op)
    new_arrays.append(array)
  mtz_out = new_arrays[0].as_mtz_dataset(
    column_root_label="A")
  for i, array in enumerate(new_arrays[1:]) :
    mtz_out.add_miller_array(
      miller_array=array,
      column_root_label="%s" % string.uppercase[i+1])
  mtz_obj = mtz_out.mtz_object()
  for i, column in enumerate(mtz_obj.columns()) :
    column.set_label(labels[i])
  if (params.reindex.output_file is None) :
    base,ext = os.path.splitext(params.reindex.hkl_file)
    params.reindex.output_file = base + "_reindex.mtz"
  mtz_obj.write(file_name=params.reindex.output_file)
  print >> out, "Reindex reflections written to %s" % params.reindex.output_file
  return params.reindex.output_file

def convert_operator (change_of_basis) :
  from cctbx import sgtbx
  try :
    c_o_b = sgtbx.change_of_basis_op(change_of_basis)
  except RuntimeError, e :
    raise Sorry(str(e))
  else :
    return c_o_b

def validate_params (params) :
  if (params.reindex.hkl_file is None) :
    raise Sorry("Please specify a reflections file.")
  if (params.reindex.change_of_basis is None) :
    raise Sorry("Please specify a change-of-basis operator.")
  else :
    convert_operator(params.reindex.change_of_basis)
  return True

if (__name__ == "__main__") :
  run(args=sys.argv[1:])
