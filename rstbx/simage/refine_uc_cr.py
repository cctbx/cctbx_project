from __future__ import absolute_import, division, print_function
from six.moves import range
from six.moves import zip
class InfeasibleError(RuntimeError): pass

class refinery(object):

  def __init__(O,
        work_params,
        spots_xy0,
        miller_indices,
        unit_cell,
        crystal_rotation_uq):
    assert spots_xy0.size() == miller_indices.size()
    O.work_params = work_params
    O.spots_xy0 = spots_xy0
    O.miller_indices = miller_indices
    O.uq_scale = 100 # to balance gradients
    O.average_unit_cell = O.work_params.lattice_symmetry.group() \
      .average_unit_cell
    unit_cell = O.average_unit_cell(unit_cell)
    from scitbx.array_family import flex
    O.x = flex.double(
        unit_cell.parameters()
      + (crystal_rotation_uq * O.uq_scale).elems)
    assert O.x.size() == 10
    O.initial_functional = None
    import scitbx.lbfgs
    scitbx.lbfgs.run(target_evaluator=O)
    O.final_functional = O.compute_functional_and_gradients(
      functional_only=True)
    del O.average_unit_cell

  def compute_functional_and_gradients(O, functional_only=False):
    from cctbx import uctbx
    from scitbx import matrix
    def get_f():
      vals = tuple(O.x[:6])
      try:
        O.unit_cell = uctbx.unit_cell(vals)
      except RuntimeError as e:
        raise InfeasibleError(str(e))
      vals = matrix.col(O.x[6:]) / O.uq_scale
      try:
        vals = vals.normalize()
      except ZeroDivisionError as e:
        raise InfeasibleError(str(e))
      if (O.average_unit_cell is not None):
        O.unit_cell = O.average_unit_cell(O.unit_cell)
      O.crystal_rotation = vals.unit_quaternion_as_r3_rotation_matrix()
      from rstbx.simage import image_simple
      O.predicted_spots = image_simple(
        apply_detector_clipping=False,
        apply_proximity_filter=False,
        store_spots=True).compute(
          unit_cell=O.unit_cell,
          miller_indices=O.miller_indices,
          spot_intensity_factors=None,
          crystal_rotation_matrix=O.crystal_rotation,
          ewald_radius=1/O.work_params.wavelength,
          ewald_proximity=O.work_params.ewald_proximity,
          signal_max=O.work_params.signal_max,
          detector_distance=O.work_params.detector.distance,
          detector_size=O.work_params.detector.size,
          detector_pixels=O.work_params.detector.pixels,
          point_spread=O.work_params.point_spread,
          gaussian_falloff_scale=O.work_params.gaussian_falloff_scale).spots
      assert O.predicted_spots.size() == O.spots_xy0.size()
      return O.spots_xy0.rms_difference(O.predicted_spots)
    f = get_f()
    if (O.initial_functional is None):
      O.initial_functional = f
    if (functional_only):
      return f
    from scitbx.array_family import flex
    g = flex.double()
    g.reserve(10)
    eps = 1e-5
    for i in range(10):
      xi = O.x[i]
      O.x[i] = xi+eps
      f_eps = get_f()
      O.x[i] = xi
      g.append((f_eps-f)/eps)
    return f, g

  def outlier_removal(O, outlier_factor=3):
    if (O.spots_xy0.size() < 3): return None
    distances = (O.spots_xy0 - O.predicted_spots).dot()**0.5
    from scitbx.array_family import flex
    perm = flex.sort_permutation(distances, reverse=True)
    if (distances[perm[0]] > distances[perm[1]] * outlier_factor):
      return perm[1:]
    return None

  def show_summary(O):
    print("refinement target:")
    print("  initial: %.6g" % O.initial_functional)
    print("    final: %.6g" % O.final_functional)
    print("refined:")
    print(O.unit_cell)
    print(O.crystal_rotation)
    return O

  def show_distances(O):
    if (O.spots_xy0.size() == 0):
      return
    distances = (O.spots_xy0 - O.predicted_spots).dot()**0.5
    from scitbx.array_family import flex
    perm = flex.sort_permutation(distances, reverse=True)
    from itertools import count
    d0 = distances[perm[0]]
    for i,h,d in zip(
          count(),
          O.miller_indices.select(perm),
          distances.select(perm)):
      if (i >= 3 and (i >= 12 or d < d0*0.1)):
        j = perm.size() - i
        if (j > 1):
          print("... remaining %d distances not shown" % j)
          break
      print("%3d %3d %3d" % h, " %7.5f" % d)
    return O

def refine(
      work_params,
      spots,
      good_i_seqs,
      miller_indices,
      unit_cell,
      crystal_rotation):
  from scitbx.array_family import flex
  spots_xy0 = flex.vec3_double()
  for spot in spots.select(good_i_seqs):
    x,y = spot.ctr_mass_x()+0.5, spot.ctr_mass_y()+0.5
    spots_xy0.append((x,y,0))
  refined = refinery(
    work_params=work_params,
    spots_xy0=spots_xy0,
    miller_indices=miller_indices,
    unit_cell=unit_cell,
    crystal_rotation_uq=crystal_rotation
      .r3_rotation_matrix_as_unit_quaternion())
  refined.show_summary().show_distances()
  print()
  while True:
    remaining_sel = refined.outlier_removal()
    if (remaining_sel is None):
      break
    print("Removing one outlier and re-refining.")
    print()
    refined = refinery(
      work_params=refined.work_params,
      spots_xy0=refined.spots_xy0.select(remaining_sel),
      miller_indices=refined.miller_indices.select(remaining_sel),
      unit_cell=refined.unit_cell,
      crystal_rotation_uq=refined.crystal_rotation
        .r3_rotation_matrix_as_unit_quaternion())
    refined.show_summary().show_distances()
    print()
  return refined
