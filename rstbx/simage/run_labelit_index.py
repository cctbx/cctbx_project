from __future__ import absolute_import, division, print_function
def get_spots_high_resolution(work_params, spots):
  from scitbx.array_family import flex
  import math
  dsx,dsy = work_params.detector.size
  dpx,dpy = work_params.detector.pixels
  dists = spots.ctr_mass_distances_from_direct_beam(
    detector_size=(dsx,dsy),
    detector_pixels=(dpx,dpy),
    xy_beam=(dsx/2-0.5,dsy/2-0.5))
  def resolution(distance_from_direct_beam):
    theta = 0.5 * math.atan2(
      distance_from_direct_beam,
      work_params.detector.distance)
    den = 2 * math.sin(theta)
    if (den == 0):
      return None
    return work_params.wavelength / den
  dists_max = flex.max(dists)
  result = resolution(dists_max)
  print("Highest-resolution of spots: %.3f" % result)
  print()
  return result

def process(work_params, spots, sampling_resolution_factor=0.5):
  import libtbx.load_env
  libtbx.env.require_module("labelit")
  spots_high_res = get_spots_high_resolution(
    work_params=work_params, spots=spots)
  if (spots_high_res is None):
    return
  uc = work_params.unit_cell
  uc_max_length = max(uc.parameters()[0:3])
  sampling = sampling_resolution_factor * spots_high_res / uc_max_length
  #
  from labelit.preferences import labelit_commands
  labelit_commands.model_refinement_minimum_N = 10
  labelit_commands.target_cell = uc
  print("labelit_commands.target_cell:", labelit_commands.target_cell)
  #
  from scitbx.array_family import flex
  raw_spot_input = flex.vec3_double()
  dsx,dsy = work_params.detector.size
  dpx,dpy = work_params.detector.pixels
  sopx,sopy = dsx/dpx, dsy/dpy
  for spot in spots: # XXX C++
    x, y = spot.ctr_mass_x(), spot.ctr_mass_y()
    raw_spot_input.append((x*sopx+0.5, y*sopy+0.5, 0.0))
  #
  from labelit.dptbx import AutoIndexEngine, Parameters
  from iotbx.detectors.context.endstation import EndStation
  ai = AutoIndexEngine(EndStation(), sampling)
  ai.setData(raw_spot_input)
  ai_params = Parameters(
    xbeam=dsx/2,
    ybeam=dsy/2,
    distance=work_params.detector.distance,
    twotheta=0.)
  ai.setBase(ai_params)
  ai.setWavelength(work_params.wavelength)
  ai.setMaxcell(1.25*max(uc.parameters()[0:3]))
  ai.setDeltaphi(0.0)
  f = ai.film_to_camera()
  c = ai.camera_to_xyz()
  #
  from labelit.dptbx.sampling import HemisphereSampler
  hem_samp = HemisphereSampler(
    max_grid=sampling,
    characteristic_grid=sampling,
    quick_grid=0.016) # all grid parameters in radians
  hem_samp.hemisphere(
    ai=ai, size=30, cutoff_divisor=4.) # never change these parameters
  #
  from labelit.dptbx.basis_choice import SelectBasisMetaprocedure
  pd = {}
  sbm = SelectBasisMetaprocedure(
    input_index_engine=ai,
    input_dictionary=pd,
    opt_rawframes=False,
    opt_target=True,
    reduce_target=False)
  ai.fixsign()
  return ai

def report_uc_cr(ai):
  miller_indices = ai.hklobserved()
  good_i_seqs = ai.get_observed_spot_i_seqs_good_only()
  assert miller_indices.size() == good_i_seqs.size()
  print("Number of spots indexed:", miller_indices.size())
  co = ai.getOrientation()
  print("Recovered unit cell:", co.unit_cell())
  print("Recovered crystal_rotation_matrix:")
  print(co.crystal_rotation_matrix())
  print()
  return good_i_seqs, miller_indices, co
