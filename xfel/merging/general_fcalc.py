from __future__ import division
from __future__ import print_function
import mmtbx.command_line.fmodel
import mmtbx.utils
from iotbx import file_reader
import libtbx.phil.command_line
import sys
import math

def run (params) :
  # assume *.mtz has Iobs or IMEAN and was created in previous step
  # from mark1 algorithm
  if params.model[-4:]==".mtz":
    from iotbx import mtz
    data_SR = mtz.object(params.model)
    for array in data_SR.as_miller_arrays():
       this_label = array.info().label_string().lower()
       if True not in [this_label.find(tag)>=0 for tag in ["iobs","imean",params.scaling.mtz_column_F]]: continue
       i_array = array.as_intensity_array().change_basis(params.model_reindex_op).map_to_asu()
       c_array = i_array.complete_array(
                 d_min = params.d_min / math.pow(1 + params.unit_cell_length_tolerance, 1 / 3))
       # complete_array adds new Miller indices to complete the a.s.u., setting sigma=-1 to
       # indicate that the new entries have no data.  (sigma == -1) is tested to remove
       # undefined Miller indices from the scaling calculation.
       return c_array
    raise Exception("mtz did not contain expected label Iobs or IMEAN")

  pdb_in = file_reader.any_file(params.model, force_type="pdb")
  pdb_in.assert_file_type("pdb")
  xray_structure = pdb_in.file_object.xray_structure_simple()
  xray_structure.show_summary()
  phil2 = mmtbx.command_line.fmodel.fmodel_from_xray_structure_master_params
  params2 = phil2.extract()
  # adjust the cutoff of the generated intensities to assure that
  # statistics will be reported to the desired high-resolution limit
  # even if the observed unit cell differs slightly from the reference.
  params2.high_resolution = params.d_min / math.pow(
    1 + params.unit_cell_length_tolerance, 1 / 3)
  if params.d_max is not None:
    params2.low_resolution = params.d_max
  params2.output.type = "real"
  if (params.include_bulk_solvent) :
    params2.fmodel.k_sol = params.k_sol
    params2.fmodel.b_sol = params.b_sol
  f_model = mmtbx.utils.fmodel_from_xray_structure(
    xray_structure = xray_structure,
    f_obs          = None,
    add_sigmas     = True,
    params         = params2).f_model
  if not params.merge_anomalous:
    f_model = f_model.generate_bijvoet_mates()
  i_model = f_model.as_intensity_array().change_basis(params.model_reindex_op).map_to_asu()

  return i_model

def random_structure (params) :

  """We're going to do some very approximate stuff here.  Given a unit
   cell & SG, will put typical atomic contents in the unit cell & get
   structure factors.

  XXX This function is no longer called from either
  command_line/cxi_merge.py nor command_line/cxi_xmerge.py; it could
  probably be removed.
  """

  import random
  random.seed(0)
  from scitbx.array_family import flex
  flex.set_random_seed(0)
  from cctbx.development import random_structure

  uc_volume = params.target_unit_cell.volume()
  asu_volume = uc_volume / params.target_space_group.group().order_z()
  target_number_scatterers = int(asu_volume)//128 # Very approximate rule of thumb for proteins with ~50% solvent content
  element_unit = ['O']*19 + ['N']*18 + ['C']*62 + ['S']*1
  element_pallet = element_unit * (1 + ( target_number_scatterers//len(element_unit) ))
  assert len(element_pallet) >= target_number_scatterers
  # XXX Ersatz hard limit to prevent excessive execution time of
  # xray_structure() below.
  elements = element_pallet[:min(1000, target_number_scatterers)]

  xs = random_structure.xray_structure(
    space_group_info = params.target_space_group,
    unit_cell = params.target_unit_cell,
    elements=elements,
    min_distance=1.2)
  xs.show_summary()
  phil2 = mmtbx.command_line.fmodel.fmodel_from_xray_structure_master_params
  params2 = phil2.extract()
  # adjust the cutoff of the generated intensities to assure that
  # statistics will be reported to the desired high-resolution limit
  # even if the observed unit cell differs slightly from the reference.
  params2.high_resolution = params.d_min / math.pow(
    1 + params.unit_cell_length_tolerance, 1 / 3)
  params2.output.type = "real"
  if (params.include_bulk_solvent) :
    print("Sorry, can't include bulk solvent for randomly-generated sites.")
  f_model = mmtbx.utils.fmodel_from_xray_structure(
    xray_structure = xs,
    f_obs          = None,
    add_sigmas     = True,
    params         = params2).f_model
  if not params.merge_anomalous:
    f_model_possibly_anomalous = f_model.generate_bijvoet_mates()
  else:
    f_model_possibly_anomalous = f_model
  i_model = f_model_possibly_anomalous.as_intensity_array()

  if params.scaling.mtz_file is not None:
    f_fake = f_model.as_amplitude_array()
    # as the code that consumes the mtz f-obs expects non-anomalous data
    mtzdata = f_fake.as_mtz_dataset(column_root_label="f-obs")
    mtzdata.mtz_object().write(params.scaling.mtz_file)

  return i_model


if (__name__ == "__main__") :
  run(sys.argv[1:])
