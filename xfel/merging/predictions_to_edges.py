from __future__ import absolute_import, division, print_function
from six.moves import range
from dials.array_family import flex # implicit dependency
import libtbx.load_env # implicit dependency
from six.moves import cPickle as pickle
from cctbx.crystal import symmetry
from scitbx.matrix import col
from rstbx.bandpass import parameters_bp3, use_case_bp3
from libtbx import easy_pickle
from libtbx.phil import parse
from spotfinder.applications.xfel.cxi_phil import cxi_basic_start
from iotbx.detectors.npy import NpyImage
import math, os
from xfel.command_line.frame_unpickler import construct_reflection_table_and_experiment_list
from dials.algorithms.spot_prediction import StillsReflectionPredictor

def load_pickle(pickle_path):
  return pickle.load(open(pickle_path, "rb"))

class ImageInfo:
  def __init__(self, img_path, detector_phil=None):
    self.img_path = img_path
    self.img_data = load_pickle(img_path)
    self.img_size = self.img_data['SIZE2']
    self.pixel_size = self.img_data['PIXEL_SIZE']
    self.detector_phil=detector_phil
    self.tm = None
  def tiling_from_image(self):
    if self.tm is not None:
      return self.tm
    labelit_regression = libtbx.env.find_in_repositories(
          relative_path="labelit_regression",
          test=os.path.isdir)
    if detector_phil is None:
      detector_phil = os.path.join(labelit_regression, "xfel", "cxi-10.1.phil")
    detector_scope = parse(file_name=detector_phil)
    basic_scope = cxi_basic_start()
    new_scope = basic_scope.phil_scope.fetch(source=detector_scope)
    tiling_phil = new_scope.extract()
    img = NpyImage("dummy", source_data=img_data)
    img.readHeader(tiling_phil)
    self.tm = img.get_tile_manager(phil)
    return self.tm

def extend_predictions(pdata, int_pickle_path, image_info, dmin=1.5, dump=False, detector_phil=None):
  """
  Given a LABELIT format integration pickle, generate a new predictor for reflections
  extending to a higher resolution dmin matching the current unit cell, orientation,
  mosaicity and domain size.
  """
  # image_info is an instance of ImageInfo
  img_path = image_info.img_path
  img_size = image_info.img_size
  pixel_size = image_info.pixel_size

  # pdata is the integration pickle object
  ori = pdata['current_orientation'][0]
  ucell = ori.unit_cell()
  sg = pdata['pointgroup']
  cbop = pdata['current_cb_op_to_primitive'][0]
  xbeam = pdata['xbeam']
  ybeam = pdata['ybeam']
  wavelength = pdata['wavelength']

  if 'effective_tiling' in pdata.keys():
    tm_int = pdata['effective_tiling']
  else:
    tiling = image_info.tiling_from_image(detector_phil=detector_phil)
    tm_int = tiling.effective_tiling_as_flex_int(reapply_peripheral_margin=True,
                                                 encode_inactive_as_zeroes=True)

  xtal = symmetry(unit_cell=ucell, space_group="P1")
  indices = xtal.build_miller_set(anomalous_flag=True, d_min=dmin)
  params = parameters_bp3(indices=indices.indices(),
                          orientation=ori,
                          incident_beam=col((0.,0.,-1.)),
                          packed_tophat=col((1.,1.,0.)),
                          detector_normal=col((0.,0.,-1.)),
                          detector_fast=col((0.,1.,0.)), # if buggy, try changing sign
                          detector_slow=col((1.,0.,0.)), # if buggy, try changing sign
                          pixel_size=col((pixel_size,pixel_size,0.)),
                          pixel_offset=col((0.5,0.5,0.0)),
                          distance=pdata['distance'],
                          detector_origin=col((-ybeam, -xbeam, 0.))) # if buggy, try changing signs
  ucbp3 = use_case_bp3(parameters=params)

  # use the tiling manager above to construct the predictor parameters
  ucbp3.set_active_areas(tm_int)
  signal_penetration = 0.5 # from LG36 trial 94 params_2.phil
  ucbp3.set_sensor_model(thickness_mm=0.1, mu_rho=8.36644, signal_penetration=signal_penetration)
  # presume no subpixel corrections for now
  ucbp3.prescreen_indices(wavelength)
  ucbp3.set_orientation(ori)
  ucbp3.set_mosaicity(pdata['ML_half_mosaicity_deg'][0]*math.pi/180) # radians
  ucbp3.set_domain_size(pdata['ML_domain_size_ang'][0])
  bandpass = 1.E-3
  wave_hi = wavelength * (1.-(bandpass/2.))
  wave_lo = wavelength * (1.+(bandpass/2.))
  ucbp3.set_bandpass(wave_hi, wave_lo)
  ucbp3.picture_fast_slow()

  # the full set of predictable reflections can now be accessed
  predicted = ucbp3.selected_predictions_labelit_format()
  hkllist = ucbp3.selected_hkls()

  # construct the experiment list and other dials backbone to be able to write predictions
  frame = construct_reflection_table_and_experiment_list(int_pickle_path,
                                                         img_path,
                                                         pixel_size,
                                                         proceed_without_image=True)
  frame.assemble_experiments()
  frame.assemble_reflections()
  predictor = StillsReflectionPredictor(frame.experiment, dmin=1.5)
  Rcalc = flex.reflection_table.empty_standard(len(hkllist))
  Rcalc['miller_index'] = hkllist
  expt_xtal = frame.experiment.crystal
  predictor.for_reflection_table(Rcalc, expt_xtal.get_A())
  predicted = Rcalc['xyzcal.mm']

  # # apply the active area filter
  # from iotbx.detectors.active_area_filter import active_area_filter
  # active_areas = list(tm.effective_tiling_as_flex_int())
  # active_area_object = active_area_filter(active_areas)
  # aa_predicted, aa_hkllist = active_area_object(predicted, hkllist, 0.11)
  # extended_mapped_predictions = flex.vec2_double()
  # for i in range(len(aa_predicted)):
  # extended_mapped_predictions.append(aa_predicted[i][0:2])

  # return predictions without re-applying an active area filter
  newpreds = flex.vec2_double()
  for i in range(len(predicted)):
    newpreds.append((predicted[i][0]/pixel_size, img_size-predicted[i][1]/pixel_size))
  # finally, record new predictions as member data
  pdata['mapped_predictions_to_edge'] = newpreds
  pdata['indices_to_edge'] = hkllist

  if dump:
    newpath = int_pickle_path.split(".pickle")[0] + "_extended.pickle"
    easy_pickle.dump(newpath, pdata)

  else:
    return (hkllist, newpreds)

if __name__ == "__main__":
  import sys
  for path in sys.argv[1:]:
    pdata = load_pickle(path)
    extend_predictions(pdata, path, dmin=1.5, dump=True)
