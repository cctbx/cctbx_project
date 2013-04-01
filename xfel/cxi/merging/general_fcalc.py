from __future__ import division
import mmtbx.command_line.fmodel
import mmtbx.utils
from iotbx import file_reader
import libtbx.phil.command_line
import sys
import math

def run (params) :
  pdb_in = file_reader.any_file(params.model, force_type="pdb")
  pdb_in.assert_file_type("pdb")
  xray_structure = pdb_in.file_object.xray_structure_simple()
  xray_structure.show_summary()
  phil2 = mmtbx.command_line.fmodel.fmodel_from_xray_structure_master_params
  params2 = phil2.extract()
  # adjust the cutoff of the generated intensities to assure that
  # statistics will be reported to the desired high-resolution limit
  # even if the observed unit cell differs slightly from the reference.
  ISO_ALLOWANCE = 0.1 # isomorphous recip cell volume changes no more than 10%
  params2.high_resolution = params.d_min / math.pow( (1.+ISO_ALLOWANCE),(1./3.) )
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
  i_model = f_model.as_intensity_array()

  return i_model

def random_structure (params) :

  """We're going to do some very approximate stuff here. Given a unit cell & SG, will put
     typical atomic contents in the unit cell & get structure factors."""
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
  elements = element_pallet[:target_number_scatterers]

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
  ISO_ALLOWANCE = 0.1 # isomorphous recip cell volume changes no more than 10%
  params2.high_resolution = params.d_min / math.pow( (1.+ISO_ALLOWANCE),(1./3.) )
  params2.output.type = "real"
  if (params.include_bulk_solvent) :
    print "Sorry, can't include bulk solvent for randomly-generated sites."
  f_model = mmtbx.utils.fmodel_from_xray_structure(
    xray_structure = xs,
    f_obs          = None,
    add_sigmas     = True,
    params         = params2).f_model
  if not params.merge_anomalous:
    f_model = f_model.generate_bijvoet_mates()
  i_model = f_model.as_intensity_array()

  if params.scaling.mtz_file is not None:
    f_fake = f_model.as_amplitude_array()
    mtzdata = f_fake.as_mtz_dataset(column_root_label="f-obs")
    mtzdata.mtz_object().write(params.scaling.mtz_file)

  return i_model


if (__name__ == "__main__") :
  run(sys.argv[1:])
