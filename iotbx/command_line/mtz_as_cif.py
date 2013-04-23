from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.mtz_as_cif

import os
import sys
from cctbx.array_family import flex
import iotbx.phil
import iotbx.cif.model

phenix_to_cif_labels_dict = {
  'FOBS': '_refln.F_meas_au',
  'SIGFOBS': '_refln.F_meas_sigma_au',
  'F-obs': '_refln.F_meas_au',
  'SIGF-obs': '_refln.F_meas_sigma_au',
  'F-obs(+)': '_refln.pdbx_F_plus',
  'SIGF-obs(+)': '_refln.pdbx_F_plus_sigma',
  'F-obs(-)': '_refln.pdbx_F_minus',
  'SIGF-obs(-)': '_refln.pdbx_F_minus_sigma',
  'I-obs': '_refln.F_squared_meas',
  'SIGI-obs': '_refln.F_squared_sigma',
  'I-obs(+)': '_refln.pdbx_I_plus',
  'SIGI-obs(+)': '_refln.pdbx_I_plus_sigma',
  'I-obs(-)': '_refln.pdbx_I_minus',
  'SIGI-obs(-)': '_refln.pdbx_I_minus_sigma',
  'R-free-flags': '_refln.phenix_R_free_flags',
  'HLA': '_refln.pdbx_HL_A_iso',
  'HLB': '_refln.pdbx_HL_B_iso',
  'HLC': '_refln.pdbx_HL_C_iso',
  'HLD': '_refln.pdbx_HL_D_iso',
  '2FOFCWT': '_refln.pdbx_FWT',
  'PH2FOFCWT': '_refln.pdbx_PHWT',
  'FOFCWT': '_refln.pdbx_DELFWT',
  'PHFOFCWT': '_refln.pdbx_DELPHWT',
  }

def mtz_to_cif_label(mtz_to_cif_label_dict, mtz_label):
  if mtz_label.endswith("xray"):
    mtz_label = mtz_label[:-5]
  elif mtz_label.endswith("neutron"):
    mtz_label = mtz_label[:-8]
  return mtz_to_cif_label_dict.get(mtz_label)

# Source: http://www.ccp4.ac.uk/html/cif2mtz.html
ccp4_to_cif_labels_dict = {
  'FREE': '_refln.status',
  'F': '_refln.F_meas_au',
  'SIGF': '_refln.F_meas_sigma_au',
  'FP': '_refln.F_meas_au',
  'SIGFP': '_refln.F_meas_sigma_au',
  'FC': '_refln.F_calc_au',
  'PHIC': '_refln.phase_calc',
  'PHIB': '_refln.phase_meas',
  'FOM': '_refln.fom',
  'I': '_refln.intensity_meas',
  'I': '_refln.F_squared_meas', # which I to prefer?
  'SIGI': '_refln.intensity_sigma',
  'SIGI': '_refln.F_squared_sigma', # which SIGI to prefer?
  'FPART': '_refln.F_part_au',
  'PHIP': '_refln.phase_part',
  'F(+)': '_refln.pdbx_F_plus',
  'SIGF(+)': '_refln.pdbx_F_plus_sigma',
  'F(-)': '_refln.pdbx_F_minus',
  'SIGF(-)': '_refln.pdbx_F_minus_sigma',
  'DP': '_refln.pdbx_anom_difference',
  'SIGDP': '_refln.pdbx_anom_difference_sigma',
  'I(+)': '_refln.pdbx_I_plus',
  'SIGI(+)': '_refln.pdbx_I_plus_sigma',
  'I(-)': '_refln.pdbx_I_minus',
  'SIGI(-)': '_refln.pdbx_I_minus_sigma',
  'HLA': '_refln.pdbx_HL_A_iso',
  'HLB': '_refln.pdbx_HL_B_iso',
  'HLC': '_refln.pdbx_HL_C_iso',
  'HLD': '_refln.pdbx_HL_D_iso',
}


master_phil = iotbx.phil.parse("""
mtz_labels = None
  .type = strings
cif_labels = None
  .type = strings
output_file = None
  .type = path
  .optional = True
  .help = 'Enter a .cif output name'

""")

def run(args):
  from iotbx import file_reader
  interpreter = master_phil.command_line_argument_interpreter()
  sources = []
  mtz_file_name = None
  mtz_object = None
  for arg in args:
    if os.path.isfile(arg):
      input_file = file_reader.any_file(arg)
      if input_file.file_type == "hkl":
        mtz_file_name = input_file.file_name
        assert input_file.file_object.file_type() == 'ccp4_mtz'
        assert mtz_object is None # only one file as input
        mtz_object = input_file.file_object.file_content()
      elif (input_file.file_type == "phil"):
        sources.append(input_file.file_object)
    else :
      arg_phil = interpreter.process(arg=arg)
      sources.append(arg_phil)
  work_phil = master_phil.fetch(sources=sources)
  work_params = work_phil.extract()

  custom_cif_labels_dict = {}
  if work_params.mtz_labels is not None and work_params.cif_labels is not None:
    assert len(work_params.mtz_labels) == len(work_params.cif_labels)
    for mtz_label, cif_label in zip(work_params.mtz_labels, work_params.cif_labels):
      custom_cif_labels_dict.setdefault(mtz_label, cif_label)

  cif_blocks = mtz_as_cif_blocks(
    mtz_object, custom_cif_labels_dict=custom_cif_labels_dict).cif_blocks

  prefix = os.path.splitext(os.path.basename(mtz_file_name))[0]
  output_file = work_params.output_file
  if output_file is None:
    output_file = prefix + ".reflections.cif"
  cif_model = iotbx.cif.model.cif()
  cif_model[prefix] = cif_blocks["xray"].cif_block
  if cif_blocks["neutron"] is not None:
    cif_model[prefix+"_neutron"] = self.cif_blocks["neutron"].cif_block
  with open(output_file, "wb") as f:
    print "Writing data and map coefficients to CIF file:\n  %s" % (f.name)
    print >> f, cif_model

class mtz_as_cif_blocks(object):

  def __init__(self, mtz_object, custom_cif_labels_dict=None, log=None):

    self.cif_blocks = {
      'xray': None,
      'neutron': None
    }

    if log is None: log = sys.stdout

    #mtz_object = mtz_dataset.mtz_object()
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

    for array in miller_arrays:
      labels = array.info().labels
      label = labels[0]
      if label.startswith("F-obs-filtered"):
        if label.endswith("neutron"):
          f_obs_filtered_neutron = array
        else:
          f_obs_filtered_xray = array
      elif label.startswith("F-obs") or label.startswith("I-obs"):
        if label.endswith("neutron"):
          input_observations_neutron = array
        else:
          input_observations_xray = array
      elif label.startswith("R-free-flags"):
        if "(+)" in label:
          array = array.average_bijvoet_mates()
          labels = [label.replace("(+)", "")]
        if label.endswith("neutron"):
          r_free_neutron = array
        else:
          r_free_xray = array
      column_names = []
      for mtz_label in labels:
        cif_label = mtz_to_cif_label(mtz_to_cif_labels_dict, mtz_label)
        column_names.append(cif_label)
      if column_names.count(None) == len(column_names):
        # I don't know what to do with this array
        continue
      assert column_names.count(None) == 0
      if labels[0].endswith("neutron"):
        data_type = "neutron"
      else:
        data_type = "xray"
      if column_names[0].startswith(("_refln.F_meas",
                                     "_refln.F_squared_meas",
                                     "_refln.pdbx_F_",
                                     "_refln.pdbx_I_")):
        input_observations_xray = array
      if self.cif_blocks.get(data_type) is None:
        self.cif_blocks[data_type] = iotbx.cif.miller_arrays_as_cif_block(
          array=array, column_names=column_names, format="mmcif")
      else:
        self.cif_blocks[data_type].add_miller_array(array, column_names=column_names)

    data_types = set(["xray"])
    if self.cif_blocks['neutron'] is not None:
      data_types.add("neutron")

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
      elif data_type == "neutron":
        r_free = r_free_neutron
        input_obs = input_observations_neutron
        f_obs_filtered = f_obs_filtered_neutron
      #assert input_obs is not None
      if input_obs is None or r_free is None: continue
      refln_status = r_free.array(data=flex.std_string(r_free.size(), "."))
      input_obs_non_anom = input_obs.average_bijvoet_mates()
      match = r_free.match_indices(input_obs_non_anom)
      refln_status.data().set_selected(match.pair_selection(0), "o")
      refln_status.data().set_selected(r_free.data() == 0, "f")
      if f_obs_filtered is not None:
        f_obs_filtered_non_anom = f_obs_filtered.average_bijvoet_mates()
        match = r_free.match_indices(f_obs_filtered_non_anom)
        refln_status.data().set_selected(match.single_selection(0), "<") # XXX
      self.cif_blocks[data_type].add_miller_array(
        array=refln_status, column_name="_refln.status")




if __name__ == '__main__':
  run(sys.argv[1:])
