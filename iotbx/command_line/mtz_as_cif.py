"""Convert MTZ reflection file to mmCIF"""

from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.mtz_as_cif

import os
import sys
from cctbx.array_family import flex
from libtbx.utils import plural_s
import iotbx.phil
import iotbx.cif.model
from iotbx import reflection_file_utils

from iotbx.cif_mtz_data_labels import phenix_to_cif_labels_dict,\
  ccp4_to_cif_labels_dict
from six.moves import zip
# Probably we can align with what PDB choose to use
# http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v40.dic/Categories/refln.html
# comply with newer version:
# http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/refln.html
# http://www.ccp4.ac.uk/html/cif2mtz.html

def mtz_to_cif_label(mtz_to_cif_label_dict, mtz_label):
  if mtz_label.endswith("xray"):
    mtz_label = mtz_label[:-5]
  elif mtz_label.endswith("neutron"):
    mtz_label = mtz_label[:-8]
  elif mtz_label.endswith(("_X", "_N")):
    mtz_label = mtz_label[:-2]
  cif_label = mtz_to_cif_label_dict.get(mtz_label)
  if cif_label is None:
    # to catch e.g. IOBS_N(+), SIGIOBS_N(+), IOBS_N(-), SIGIOBS_N(-)
    mtz_label = mtz_label.replace("_N", "").replace("_X", "")
    cif_label = mtz_to_cif_label_dict.get(mtz_label)
  return cif_label

master_phil = iotbx.phil.parse("""
mtz_as_cif
  .short_caption = MTZ as mmCIF
  .caption = This program will convert reflections in MTZ format to mmCIF, \
    suitable for PDB deposition.  Note that phenix.refine can also write \
    mmCIF files directly if desired.
  .style = auto_align box \
    caption_img:icons/custom/phenix.reflection_file_editor.png
{
mtz_file = None
  .type = path
  .multiple = True
  .short_caption = MTZ file
  .style = file_type:mtz input_file bold
output_file = None
  .type = path
  .help = Optional output file name to override default
  .optional = True
  .help = 'Enter a .cif output name'
  .style = file_type:cif bold new_file
mtz_labels = None
  .help = Custom input labels for unknown MTZ columns
  .short_caption = Custom input labels
  .type = strings
cif_labels = None
  .help = Custom output labels for unknown mmCIF columns
  .short_caption = Custom output labels
  .type = strings
}
""")

def format_usage_message(log=sys.stdout):
  print("-"*79, file=log)
  msg = """\
phenix.mtz_as_cif: Convert mtz to CIF format.
The tool will automatically recognize most of the labels from Phenix and CCP4
and convert them to appropriate mmCIF format. If some labels are not recognized,
or another mmCIF labels are needed for them, one may use mtz_labels and
cif_labels parameters to provide pairs of mtz and cif labels. Note that provided
cif_labels must comply with mmCIF format. Number of mtz_labels should be equal
to number of cif_labels.

Usage: phenix.mtz_as_cif data.mtz [params.eff] [options ...]

Usage examples:
  phenix.mtz_as_cif data.mtz
  phenix.mtz_as_cif data.mtz output_file=custom.cif mtz_labels="FOBS SIGFOBS"
      cif_labels="_refln.custom1 _refln.custom2"
"""
  print(msg, file=log)
  print("-"*79, file=log)
  print(master_phil.show(), file=log)

def run(args, params=None, out=sys.stdout):
  from iotbx import file_reader
  work_params = params
  if (work_params is None):
    cmdline = iotbx.phil.process_command_line_with_files(
      args=args,
      master_phil=master_phil,
      reflection_file_def="mtz_as_cif.mtz_file")
    work_params = cmdline.work.extract()
  if (len(work_params.mtz_as_cif.mtz_file) == 0):
    format_usage_message(log = out)
    return
  work_params = work_params.mtz_as_cif
  mtz_objects = []
  for file_name in work_params.mtz_file :
    input_file = file_reader.any_file(file_name)
    input_file.check_file_type("hkl")
    if (input_file.file_object.file_type() != 'ccp4_mtz'):
      raise Sorry("Error reading '%s' - only MTZ files may be used as input."
        % file_name)
    mtz_objects.append(input_file.file_object.file_content())
  assert (len(mtz_objects) != 0)
  custom_cif_labels_dict = {}
  if work_params.mtz_labels is not None and work_params.cif_labels is not None:
    assert len(work_params.mtz_labels) == len(work_params.cif_labels)
    for mtz_label, cif_label in zip(work_params.mtz_labels, work_params.cif_labels):
      custom_cif_labels_dict.setdefault(mtz_label, cif_label)
  output_files = []
  for mtz_file_name, mtz_object in zip(work_params.mtz_file, mtz_objects):
    print("Converting %s" %mtz_file_name, file=out)
    cif_blocks = mtz_as_cif_blocks(
      mtz_object, custom_cif_labels_dict=custom_cif_labels_dict).cif_blocks

    prefix = os.path.splitext(os.path.basename(mtz_file_name))[0]
    output_file = work_params.output_file
    if output_file is None:
      output_file = prefix + ".reflections.cif"
    cif_model = iotbx.cif.model.cif()
    # This is gross... refactor some time the whole thing.
    for key in cif_blocks.keys():
      print(key)
      if key == "xray" and cif_blocks["xray"] is not None:
        cif_model[prefix] = cif_blocks["xray"].cif_block

      elif key == "neutron" and cif_blocks["neutron"] is not None:
        cif_model[prefix+"_neutron"] = cif_blocks["neutron"].cif_block

      elif cif_blocks[key] is not None:
        cif_model[prefix+"_"+key] = cif_blocks[key].cif_block
    with open(output_file, "w") as f:
      print("Writing data and map coefficients to CIF file:\n  %s" % \
        (f.name), file=out)
      print(cif_model, file=f)
      output_files.append(output_file)
  return output_files

class mtz_as_cif_blocks(object):

  def __init__(self, mtz_object, custom_cif_labels_dict=None, log=None,
      test_flag_value=None):

    self.cif_blocks = {
      'xray': None,
      'neutron': None
    }

    if log is None: log = sys.stdout

    miller_arrays = mtz_object.as_miller_arrays()

    miller_arrays_as_cif_block = None

    input_observations_xray = None
    input_observations_neutron = None
    r_free_xray = None
    r_free_neutron = None
    f_obs_filtered_xray = None
    f_obs_filtered_neutron = None

    mtz_to_cif_labels_dict = {}
    mtz_to_cif_labels_dict.update(phenix_to_cif_labels_dict)
    mtz_to_cif_labels_dict.update(ccp4_to_cif_labels_dict)
    if custom_cif_labels_dict is not None:
      mtz_to_cif_labels_dict.update(custom_cif_labels_dict)

    unknown_mtz_labels = []

    for array in miller_arrays:
      labels = array.info().labels
      label = labels[0]
      if reflection_file_utils.looks_like_r_free_flags_info(array.info()):
        if "(+)" in label:
          array = array.average_bijvoet_mates()
          labels = [label.replace("(+)", "")]
        if label.endswith(("neutron", "_N")):
          r_free_neutron = array
        else:
          r_free_xray = array
        continue # deal with these later
      elif label.startswith("F-obs-filtered"):
        if label.endswith(("neutron", "_N")):
          f_obs_filtered_neutron = array
        else:
          f_obs_filtered_xray = array
      elif label.startswith("F-obs") or label.startswith("I-obs"):
        if label.strip("(+)").endswith(("neutron", "_N")):
          input_observations_neutron = array
        else:
          input_observations_xray = array
      #elif label.startswith("R-free-flags"):
      column_names = []
      for mtz_label in labels:
        cif_label = mtz_to_cif_label(mtz_to_cif_labels_dict, mtz_label)
        column_names.append(cif_label)

      column_names = self.check_for_dano_and_convert(column_names, labels, array)
      if column_names.count(None) > 0:
        # I don't know what to do with this array
        for i, mtz_label in enumerate(labels):
          if column_names[i] is None:
            unknown_mtz_labels.append(mtz_label)
        continue
      assert column_names.count(None) == 0
      if labels[0].strip("(+)").endswith(("neutron", "_N")):
        data_type = "neutron"
      else:
        data_type = "xray"
      if column_names[0].startswith(("_refln.F_meas",
                                     "_refln.F_squared_meas",
                                     "_refln.pdbx_F_",
                                     "_refln.pdbx_I_")):
        if data_type == "neutron":
          input_observations_neutron = array
        else:
          input_observations_xray = array

      if self.cif_blocks.get(data_type) is None:
        self.cif_blocks[data_type] = iotbx.cif.miller_arrays_as_cif_block(
          array=array, column_names=column_names, format="mmcif")
      else:
        # check if it is taken already
        present = False
        for ln in column_names:
          if ln in self.cif_blocks[data_type].refln_loop.keys():
            present = True

        if present:
          labels_string = "_"
          for l in labels:
            labels_string += l+" "
          labels_string = labels_string.strip().replace(" ", "_")
          self.cif_blocks[data_type+labels_string] = iotbx.cif.miller_arrays_as_cif_block(
              array=array, column_names=column_names, format="mmcif")
        else:
          self.cif_blocks[data_type].add_miller_array(array, column_names=column_names)

    if len(unknown_mtz_labels):
      print("Warning: Unknown mtz label%s: %s" %(
        plural_s(len(unknown_mtz_labels))[1], ", ".join(unknown_mtz_labels)), file=log)
      print("  Use mtz_labels and cif_labels keywords to provide translation for custom labels.", file=log)

    data_types = ["xray"]
    if self.cif_blocks['neutron'] is not None:
      data_types.append("neutron")

    if input_observations_xray is None and f_obs_filtered_xray is not None:
      self.cif_blocks["xray"].add_miller_array(
        array=f_obs_filtered_xray,
        column_names=('_refln.F_meas_au','_refln.F_meas_sigma_au'))
    if input_observations_neutron is None and f_obs_filtered_neutron is not None:
      self.cif_blocks["neutron"].add_miller_array(
        array=f_obs_filtered_neutron,
        column_names=('_refln.F_meas_au','_refln.F_meas_sigma_au'))

    for data_type in data_types:
      if data_type == "xray":
        r_free = r_free_xray
        input_obs = input_observations_xray
        f_obs_filtered = f_obs_filtered_xray
        if (self.cif_blocks["xray"] is None and r_free_xray is not None and
            self.cif_blocks["neutron"] is not None and r_free_neutron is None):
          r_free_neutron = r_free_xray
      elif data_type == "neutron":
        r_free = r_free_neutron
        input_obs = input_observations_neutron
        f_obs_filtered = f_obs_filtered_neutron
      if self.cif_blocks[data_type] is not None and r_free is not None:
        self.cif_blocks[data_type].add_miller_array(
          array=r_free, column_name='_refln.pdbx_r_free_flag')

      if input_obs is None or r_free is None: continue
      # it may happen that there is an Rfree array but the values are all identical
      if (r_free.data().all_eq(r_free.data()[0])):
        refln_status = r_free.array(data=flex.std_string(r_free.size(), "o"))
        self.cif_blocks[data_type].add_miller_array(
          array=refln_status, column_name="_refln.status")
        continue
      if (test_flag_value is None):
        test_flag_value = reflection_file_utils.guess_r_free_flag_value(
          miller_array=r_free)
      assert (test_flag_value is not None)
      refln_status = r_free.array(data=flex.std_string(r_free.size(), "."))
      input_obs_non_anom = input_obs.average_bijvoet_mates()
      match = r_free.match_indices(input_obs_non_anom)
      refln_status.data().set_selected(match.pair_selection(0), "o")
      refln_status.data().set_selected(r_free.data() == test_flag_value, "f")
      if f_obs_filtered is not None:
        f_obs_filtered_non_anom = f_obs_filtered.average_bijvoet_mates()
        match = r_free.match_indices(f_obs_filtered_non_anom)
        refln_status.data().set_selected(match.single_selection(0), "<") # XXX
      self.cif_blocks[data_type].add_miller_array(
        array=refln_status, column_name="_refln.status")

  def check_for_dano_and_convert(self, column_names, labels, array):
    need_to_convert = False
    for l in labels:
      if l.lower().find("dano") >= 0:
        need_to_convert = True
        break
    if not need_to_convert:
      return column_names
    else:
      result = []
      #  assert array.anomalous_flag()  # No it is not anomalous (only one
      #    number per reflection...it is anomalous data though
      result = ["_refln.pdbx_anom_difference",
                "_refln.pdbx_anom_difference_sigma",]

      if labels[0].lower().find('i-obs') >= 0:
        raise Sorry("Cannot convert anomalous differences on intensity to CIF")
      return result


def validate_params(params):
  if (len(params.mtz_as_cif.mtz_file) == 0):
    raise Sorry("No MTZ file(s) specified!")
  return True

if __name__ == '__main__':
  run(sys.argv[1:])
