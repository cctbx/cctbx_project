import os
op = os.path
import sys
import time

_show_vm_info_time = time.time()

def show_vm_info(msg):
  print msg
  from libtbx import introspection
  introspection.virtual_memory_info().show(prefix="  ", show_max=True)
  global _show_vm_info_time
  t = time.time()
  print "  time since previous: %.2f seconds" % (t-_show_vm_info_time)
  _show_vm_info_time = t
  print
  sys.stdout.flush()

import libtbx

class image_model(libtbx.slots_getstate_setstate):

  __slots__ = [
    "pixels",
    "spot_positions",
    "spot_intensities",
    "miller_index_i_seqs",
    "unit_cell",
    "crystal_rotation",
    "partialities",
    "scale",
    "i_perm",
    "backup"]

  def __init__(O,
        pixels=None,
        spot_positions=None,
        spot_intensities=None,
        miller_index_i_seqs=None,
        unit_cell=None,
        crystal_rotation=None,
        partialities=None,
        scale=None,
        i_perm=None):
    O.pixels = pixels
    O.spot_positions = spot_positions
    O.spot_intensities = spot_intensities
    O.miller_index_i_seqs = miller_index_i_seqs
    O.unit_cell = unit_cell
    O.crystal_rotation = crystal_rotation
    O.partialities = partialities
    O.scale = scale
    O.i_perm = i_perm
    O.backup = None

  def make_backup(O):
    O.backup = image_model(
      spot_positions=O.spot_positions,
      spot_intensities=O.spot_intensities,
      miller_index_i_seqs=O.miller_index_i_seqs,
      unit_cell=O.unit_cell,
      crystal_rotation=O.crystal_rotation,
      partialities=O.partialities,
      scale=O.scale,
      i_perm=O.i_perm)

  def erase_spot_model(O):
    O.spot_positions = None
    O.spot_intensities = None
    O.miller_index_i_seqs = None
    O.unit_cell = None
    O.crystal_rotation = None
    O.partialities = None
    O.scale = None

  def reset_spot_model(O, other):
    if (other is None):
      O.erase_spot_model()
    else:
      O.spot_positions = other.spot_positions
      O.spot_intensities = other.spot_intensities
      O.miller_index_i_seqs = other.miller_index_i_seqs
      O.unit_cell = other.unit_cell
      O.crystal_rotation = other.crystal_rotation
      O.partialities = other.partialities
      O.scale = other.scale

  def reindex_in_place(O,
        reindexing_assistant=None,
        cb_op=None,
        miller_indices=None):
    assert [reindexing_assistant, cb_op].count(None) == 1
    assert (cb_op is None) == (miller_indices is None)
    if (reindexing_assistant is not None):
      assert O.i_perm is not None
      cb_op = reindexing_assistant.cb_ops[O.i_perm]
    if (not O.unit_cell.is_similar_to(
              other=O.unit_cell.change_basis(cb_op),
              relative_length_tolerance=1e-5,
              absolute_angle_tolerance=1e-3)):
      raise RuntimeError(
        "Unit cell is not compatible with reindexing operation.")
    if (reindexing_assistant is not None):
      assert O.i_perm is not None
      perm = reindexing_assistant.perms[O.i_perm]
      O.miller_index_i_seqs = perm.select(O.miller_index_i_seqs)
    else:
      mi_cb = cb_op.apply(miller_indices.select(O.miller_index_i_seqs))
      from cctbx import miller
      matches = miller.match_indices(miller_indices, mi_cb)
      assert matches.singles(1).size() == 0
      O.miller_index_i_seqs = matches.pairs().column(0)
    from scitbx.array_family import flex
    sort_perm = flex.sort_permutation(data=O.miller_index_i_seqs)
    O.miller_index_i_seqs = O.miller_index_i_seqs.select(sort_perm)
    O.spot_positions = O.spot_positions.select(sort_perm)
    O.spot_intensities = O.spot_intensities.select(sort_perm)
    from scitbx import matrix
    c_cart = matrix.sqr(O.unit_cell.matrix_cart(rot_mx=cb_op.c_inv().r()))
    O.crystal_rotation = (matrix.sqr(O.crystal_rotation) * c_cart).elems
    O.partialities = None
    O.i_perm = 0
    if (O.backup is not None):
      O.backup.i_perm = 0

  def reset_partialities(O, work_params, miller_indices):
    from rstbx.simage import image_simple
    O.partialities = image_simple(
      apply_detector_clipping=False,
      apply_proximity_filter=False,
      store_signals=True).compute(
        unit_cell=O.unit_cell,
        miller_indices=miller_indices.select(O.miller_index_i_seqs),
        spot_intensity_factors=None,
        crystal_rotation_matrix=O.crystal_rotation,
        ewald_radius=1/work_params.wavelength,
        ewald_proximity=work_params.ewald_proximity,
        signal_max=1,
        detector_distance=work_params.detector.distance,
        detector_size=work_params.detector.size,
        detector_pixels=work_params.detector.pixels,
        point_spread=work_params.point_spread,
        gaussian_falloff_scale=work_params.gaussian_falloff_scale).signals
    assert O.partialities.size() == O.miller_index_i_seqs.size()

  def usable(O, partiality_threshold):
    sel = O.partialities > partiality_threshold
    return libtbx.group_args(
      miis = O.miller_index_i_seqs.select(sel),
      esti = O.spot_intensities.select(sel) / O.partialities.select(sel))

  def extract_i_obs_est(O, work_params, miller_indices):
    assert O.partialities is not None
    usable = O.usable(work_params.usable_partiality_threshold)
    from cctbx import crystal
    return crystal.symmetry(
      unit_cell=work_params.unit_cell,
      space_group_symbol="P1").miller_set(
        indices=miller_indices.select(usable.miis),
        anomalous_flag=True).array(
          data=usable.esti)

class miller_image_map(libtbx.slots_getstate_setstate):

  __slots__ = ["miller_indices", "map"]

  def __init__(O, miller_indices):
    O.miller_indices = miller_indices
    O.map = [[] for i in xrange(O.miller_indices.size())]

  def enter(O, i_img, miller_index_i_seqs):
    map = O.map
    for ii_seq,i_seq in enumerate(miller_index_i_seqs):
      map[i_seq].append((i_img, ii_seq))

  def show_images_per_miller_index(O, first_block_size=20):
    print "Images per Miller index:"
    from libtbx import dict_with_default_0
    counts = dict_with_default_0()
    for iiis in O.map:
      counts[len(iiis)] += 1
    n_seq = O.miller_indices.size()
    have_break = False
    for n_imgs in sorted(counts.keys()):
      if (n_imgs > first_block_size and n_imgs < len(counts)-5):
        if (not have_break):
          have_break = True
          print "        ..."
      else:
        c = counts[n_imgs]
        print "  %6d %6d %8.6f" % (n_imgs, c, c/n_seq)
    print
    sys.stdout.flush()

def collect_estis(image_mdls_array, iiis, partiality_threshold):
  from scitbx.array_family import flex
  result = flex.double()
  for i_img,ii_seq in iiis:
    im = image_mdls_array[i_img]
    scale = im.scale
    if (scale != 0):
      signal = im.spot_intensities[ii_seq]
      if (signal == 0):
        result.append(0)
      else:
        part = im.partialities[ii_seq]
        if (part != 0 and part >= partiality_threshold):
          result.append(signal / part / scale)
  return result

class image_models(libtbx.slots_getstate_setstate):

  __slots__ = ["miller_indices", "array", "miller_image_map"]

  def __init__(O, miller_indices, array, miller_image_map=None):
    O.miller_indices = miller_indices
    O.array = array
    O.miller_image_map = miller_image_map

  def size(O):
    return len(O.array)

  def check_i_perm_vs_backup(O, reindexing_assistant):
    im0_i_perm = O.array[0].backup.i_perm
    for im in O.array:
      assert im.i_perm is not None
      assert im.i_perm == reindexing_assistant.i_j_inv_multiplication_table[
        im0_i_perm][
        im.backup.i_perm]

  def erase_spot_models(O):
    for im in O.array:
      im.erase_spot_model()

  def extract_scales(O):
    from scitbx.array_family import flex
    result = flex.double()
    for im in O.array:
      result.append(im.scale)
    return result

  def erase_scales(O):
    for im in O.array:
      im.scale = None

  def reset_scales(O, all_scales):
    for im,scale in zip(O.array, all_scales):
      im.scale = scale

  def iselection_entries_with_spot_model(O):
    from scitbx.array_family import flex
    result = flex.size_t()
    for i,im in enumerate(O.array):
      if (im.spot_positions is not None):
        result.append(i)
    return result

  def remove_all_entries_without_spot_model(O):
    remaining = []
    for im in O.array:
      if (im.spot_positions is not None):
        remaining.append(im)
    return image_models(miller_indices=O.miller_indices, array=remaining)

  def normalize_spot_intensities(O, target_mean):
    from scitbx.array_family import flex
    sum_si = 0
    num_si = 0
    for im in O.array:
      sum_si += flex.sum(im.spot_intensities)
      num_si += im.spot_intensities.size()
    if (sum_si != 0):
      global_scale = target_mean * num_si / sum_si
      for im in O.array:
        im.spot_intensities *= global_scale

  def reset_miller_image_map(O):
    O.miller_image_map = miller_image_map(miller_indices=O.miller_indices)
    for i_img,im in enumerate(O.array):
      O.miller_image_map.enter(
        i_img=i_img, miller_index_i_seqs=im.miller_index_i_seqs)

  def reset_partialities(O, work_params):
    for im in O.array:
      im.reset_partialities(work_params, O.miller_indices)

  def check_i_obs_vs_backup(O, work_params):
    print "Current i_obs vs. backup:"
    for im in O.array:
      im.backup.reset_partialities(work_params, O.miller_indices)
      b_obs = im.backup.extract_i_obs_est(work_params, O.miller_indices)
      im.reset_partialities(work_params, O.miller_indices)
      i_obs = im.extract_i_obs_est(work_params, O.miller_indices)
      max_common_size = -1
      max_cb_ci = None
      for s in work_params.lattice_symmetry.group():
        i_obs_cb = i_obs.change_basis(str(s))
        cb, ci = b_obs.common_sets(other=i_obs_cb)
        common_size = cb.indices().size()
        if (max_common_size < common_size):
          max_common_size = common_size
          max_cb_ci = cb, ci
      assert max_cb_ci is not None
      cb, ci = max_cb_ci
      from scitbx.array_family import flex
      num = flex.sum(cb.data()*ci.data())
      den = flex.sum_sq(cb.data())
      if (den == 0): scale = None
      else:          scale = num / den
      print " ", b_obs.indices().size(), i_obs.indices().size(), \
        cb.indices().size(), scale
    print

  def refinement_target(O, partiality_threshold):
    assert O.miller_image_map.map is not None
    from scitbx.array_family import flex
    result_num = 0
    result_den = 0
    for iiis in O.miller_image_map.map:
      estis = collect_estis(O.array, iiis, partiality_threshold)
      if (estis.size() < 2): continue
      i_obs_est = flex.mean(estis)
      result_num += flex.sum_sq(estis - i_obs_est)
      result_den += estis.size()
    return result_num / max(1, result_den)

  def extract_estimated_i_obs(O, partiality_threshold):
    from cctbx.array_family import flex
    indices = flex.miller_index()
    data = flex.double()
    mimmi = O.miller_image_map.miller_indices
    indices.reserve(mimmi.size())
    data.reserve(mimmi.size())
    for h,iiis in zip(mimmi, O.miller_image_map.map):
      estis = collect_estis(O.array, iiis, partiality_threshold)
      if (estis.size() != 0):
        indices.append(h)
        data.append(flex.mean(estis))
    return (indices, data)

  def write_to_mtz_files(O, common_unit_cell):
    from cctbx import crystal
    crystal_symmetry = crystal.symmetry(
      unit_cell=common_unit_cell,
      space_group_symbol="P1")
    def write_mtz(file_name, counts=None, miis=None):
      if (miis is None):
        isel = (counts != 0).iselection()
        data = counts.select(isel)
      else:
        isel = miis
        data = flex.size_t(isel.size(), 1)
      ma = crystal_symmetry.miller_set(
        indices=O.miller_indices.select(isel),
        anomalous_flag=True).array(data=data)
      ma.as_mtz_dataset(column_root_label="NOBS").mtz_object().write(
        file_name=file_name)
    n_indices = O.miller_indices.size()
    from scitbx.array_family import flex
    counts_all = flex.size_t(n_indices, 0)
    miis_0 = None
    for i_img,im in enumerate(O.array):
      miis = im.miller_index_i_seqs
      write_mtz(file_name="nobs_%03d.mtz" % i_img, miis=miis)
      counts_all.increment_and_track_up_from_zero(
        iselection=im.miller_index_i_seqs)
      if (miis_0 is None):
        miis_0 = miis
      else:
        counts_pair = flex.size_t(n_indices, 0)
        for isel in [miis_0, miis]:
          counts_pair.increment_and_track_up_from_zero(iselection=isel)
        write_mtz(file_name="nobs_000_%03d.mtz" % i_img, counts=counts_pair)
    write_mtz(file_name="nobs_all.mtz", counts=counts_all)

class refinement_target_eps(object):

  __slots__ = ["image_mdls", "partiality_threshold", "eps"]

  def __init__(O, image_mdls, partiality_threshold, eps):
    O.image_mdls = image_mdls
    O.partiality_threshold = partiality_threshold
    O.eps = eps

  def __call__(O, i_img):
    im = O.image_mdls.array[i_img]
    scale_orig = im.scale
    im.scale = scale_orig + O.eps
    result = O.image_mdls.refinement_target(O.partiality_threshold)
    im.scale = scale_orig
    return (i_img, result)

class refinery(object):

  def __init__(O, work_params, image_mdls):
    O.work_params = work_params
    O.image_mdls = image_mdls
    from scitbx.array_family import flex
    O.x = flex.double()
    O.x.reserve(O.image_mdls.size())
    for im in O.image_mdls.array:
      O.x.append(im.scale)
    O.initial_functional = None
    O.number_of_iterations = 0
    O.number_of_function_evaluations = 0
    import scitbx.lbfgs
    scitbx.lbfgs.run(
      target_evaluator=O,
      termination_params=scitbx.lbfgs.termination_parameters(
        max_iterations=O.work_params.refine_scales.max_iterations),
      exception_handling_params=scitbx.lbfgs.exception_handling_parameters(
        ignore_line_search_failed_step_at_lower_bound=True))
    O.image_mdls.reset_scales(all_scales=O.x)
    O.final_functional = O.image_mdls.refinement_target(
      partiality_threshold=O.work_params.usable_partiality_threshold)

  def compute_functional_and_gradients(O):
    O.image_mdls.reset_scales(all_scales=O.x)
    f = O.image_mdls.refinement_target(
      O.work_params.usable_partiality_threshold)
    if (O.initial_functional is None):
      O.initial_functional = f
    O.number_of_function_evaluations += 1
    n_mdls = O.x.size()
    from scitbx.array_family import flex
    g = flex.double()
    g.reserve(n_mdls)
    eps = O.work_params.refine_scales.finite_difference_eps
    if (not O.work_params.multiprocessing or n_mdls < 2):
      for im,x in zip(O.image_mdls.array, O.x):
        im.scale = x+eps
        f_eps = O.image_mdls.refinement_target(
          O.work_params.usable_partiality_threshold)
        im.scale = x
        g.append((f_eps-f)/eps)
    else:
      from libtbx import easy_mp
      mp_results = easy_mp.pool_map(
        fixed_func=refinement_target_eps(
          O.image_mdls, O.work_params.usable_partiality_threshold, eps),
        args=range(n_mdls),
        chunksize=1,
        log=sys.stdout)
      g.resize(n_mdls)
      for i,f_eps in mp_results:
        g[i] = (f_eps-f)/eps
    print "refine scale f, |g|: %.6g, %.6g" % (f, g.norm())
    sys.stdout.flush()
    return f, g

  def callback_after_step(O, minimizer):
    O.number_of_iterations += 1

  def show_summary(O):
    print "refinement target:"
    print "  initial: %.6g" % O.initial_functional
    print "    final: %.6g" % O.final_functional
    print "            iterations:", O.number_of_iterations
    print "  function evaluations:", O.number_of_function_evaluations
    print
    sys.stdout.flush()

# import before fork
from rstbx.simage import \
  run_spotfinder, \
  run_labelit_index, \
  refine_uc_cr, \
  integrate_crude

def index_and_integrate_one(work_params, image_mdls_miller_indices, pixels):
  spots = run_spotfinder.process(
    work_params=work_params, pixels=pixels, show_spots=False)
  if (spots.size() < work_params.min_number_of_spots_for_indexing):
    print "Insufficient number of spots for indexing."
    print
    sys.stdout.flush()
    return (spots.size(), None)
  ai = run_labelit_index.process(work_params=work_params, spots=spots)
  good_i_seqs, miller_indices, co = run_labelit_index.report_uc_cr(ai)
  refined = refine_uc_cr.refine(
    work_params=work_params,
    spots=spots,
    good_i_seqs=good_i_seqs,
    miller_indices=miller_indices,
    unit_cell=co.unit_cell(),
    crystal_rotation=co.crystal_rotation_matrix())
  predicted_spot_positions, \
  predicted_spot_miller_index_i_seqs = integrate_crude.predict_spot_positions(
    work_params=work_params,
    miller_indices=image_mdls_miller_indices,
    unit_cell=refined.unit_cell,
    crystal_rotation=refined.crystal_rotation)
  print "Number of predicted spot positions:", predicted_spot_positions.size()
  print
  spot_intensities = integrate_crude.collect_spot_intensities(
    pixels=pixels,
    spot_positions=predicted_spot_positions,
    point_spread_inner=work_params.point_spread,
    point_spread_outer=work_params.point_spread+4)
  sel = spot_intensities != 0
  return (
    spots.size(), image_model(
      spot_positions=predicted_spot_positions.select(sel),
      spot_intensities=spot_intensities.select(sel),
      miller_index_i_seqs=predicted_spot_miller_index_i_seqs.select(sel),
      unit_cell=refined.unit_cell,
      crystal_rotation=refined.crystal_rotation))

def index_and_integrate(work_params, image_mdls):
  n_mdls = image_mdls.size()
  if (not work_params.multiprocessing or n_mdls < 2):
    for im in image_mdls.array:
      n_spots, updated_im = index_and_integrate_one(
        work_params, image_mdls.miller_indices, im.pixels)
      im.reset_spot_model(other=updated_im)
  else:
    def mp_func(i_img):
      return index_and_integrate_one(
        work_params,
        image_mdls.miller_indices,
        image_mdls.array[i_img].pixels)
    from libtbx import easy_mp
    mp_results = easy_mp.pool_map(
      fixed_func=mp_func,
      args=range(n_mdls),
      chunksize=1,
      log=sys.stdout,
      buffer_stdout_stderr=True)
    print
    sys.stdout.flush()
    for i_img,(log,mp_result) in enumerate(mp_results):
      if (mp_result is None):
        print "ERROR index_and_integrate_one:"
        print "-"*80
        sys.stdout.write(log)
        print "-"*80
        print
      else:
        n_spots, updated_im = mp_result
        if (updated_im is None):
          uc = None
        else:
          uc = updated_im.unit_cell
        print "Refined unit cell %d (%d spots):" % (i_img, n_spots), uc
        image_mdls.array[i_img].reset_spot_model(other=updated_im)
      sys.stdout.flush()
    print
    if (work_params.show_refine_uc_cr):
      for _,(log,_) in enumerate(mp_results):
        print "v"*80
        sys.stdout.write(log)
        print "^"*80
        print
      sys.stdout.flush()

def check_refine_uc_cr(work_params, image_mdls,
      unit_cell_perturbation_factor=2,
      crystal_rotation_perturbation_angle=10):
  from cctbx import uctbx
  from scitbx.array_family import flex
  from scitbx import matrix
  for i_img,im in enumerate(image_mdls.array):
    print "Image number:", i_img
    mt = flex.mersenne_twister(seed=work_params.noise.random_seed+i_img)
    unit_cell = uctbx.unit_cell([
      v + unit_cell_perturbation_factor*(mt.random_double()-0.5)
        for v in im.unit_cell.parameters()])
    crystal_rotation = matrix.sqr(im.crystal_rotation) \
      * matrix.col(mt.random_double_point_on_sphere()) \
          .axis_and_angle_as_r3_rotation_matrix(
            angle=crystal_rotation_perturbation_angle, deg=True)
    refined = refine_uc_cr.refinery(
      work_params=work_params,
      spots_xy0=im.spot_positions,
      miller_indices=image_mdls.miller_indices.select(im.miller_index_i_seqs),
      unit_cell=unit_cell,
      crystal_rotation_uq=crystal_rotation
        .r3_rotation_matrix_as_unit_quaternion())
    refined.show_summary().show_distances()
    print

def build_images(work_params, i_calc, reindexing_assistant):
  result = []
  from create import add_noise
  from rstbx.simage import image_simple
  from cctbx.array_family import flex
  if (not work_params.apply_random_reindexing):
    i_calc_data_perms = [i_calc.data()]
  else:
    i_calc_data_perms = [i_calc.data().select(perm)
      for perm in reindexing_assistant.inv_perms]
  n_mdls = work_params.number_of_shots
  use_mp = (work_params.multiprocessing and n_mdls > 1)
  def build_one_image(i_img):
    mt = flex.mersenne_twister(seed=work_params.noise.random_seed+i_img)
    scale = int(work_params.signal_max*(0.1+0.9*mt.random_double()))
    crystal_rotation = mt.random_double_r3_rotation_matrix_arvo_1992()
    i_perm = mt.random_size_t() % len(i_calc_data_perms)
    image = image_simple(
      store_miller_index_i_seqs=True,
      store_spots=True,
      store_signals=True,
      set_pixels=True).compute(
        unit_cell=i_calc.unit_cell(),
        miller_indices=i_calc.indices(),
        spot_intensity_factors=i_calc_data_perms[i_perm],
        crystal_rotation_matrix=crystal_rotation,
        ewald_radius=1/work_params.wavelength,
        ewald_proximity=work_params.ewald_proximity,
        signal_max=scale,
        detector_distance=work_params.detector.distance,
        detector_size=work_params.detector.size,
        detector_pixels=work_params.detector.pixels,
        point_spread=work_params.point_spread,
        gaussian_falloff_scale=work_params.gaussian_falloff_scale)
    add_noise(work_params, pixels=image.pixels)
    if (not work_params.index_and_integrate):
      pixels = None
    else:
      pixels = image.pixels
    miller_index_i_seqs = image.miller_index_i_seqs
    if (use_mp):
      # to by-pass portable but slower pickling
      if (pixels is not None):
        assert pixels.is_0_based()
        assert not pixels.is_padded()
        assert pixels.all() == tuple(work_params.detector.pixels)
        pixels = pixels.copy_to_byte_str()
      miller_index_i_seqs = miller_index_i_seqs.copy_to_byte_str()
    return image_model(
      pixels=pixels,
      spot_positions=image.spots,
      spot_intensities=image.signals,
      unit_cell=i_calc.unit_cell(),
      crystal_rotation=crystal_rotation,
      miller_index_i_seqs=miller_index_i_seqs,
      scale=scale,
      i_perm=i_perm)
  if (not use_mp):
    for i_img in xrange(n_mdls):
      result.append(build_one_image(i_img))
  else:
    from libtbx import easy_mp
    result = easy_mp.pool_map(
      fixed_func=build_one_image,
      args=range(n_mdls),
      chunksize=1,
      log=sys.stdout)
    for im in result:
      if (im is None): raise RuntimeError("Failure building image.")
      if (im.pixels is not None):
        im.pixels = flex.int_from_byte_str(im.pixels)
        im.pixels.reshape(flex.grid(work_params.detector.pixels))
      im.miller_index_i_seqs = flex.size_t_from_byte_str(
        byte_str=im.miller_index_i_seqs)
  for im in result:
    im.make_backup()
  return image_models(miller_indices=i_calc.indices(), array=result)

class perm_rms_info(libtbx.slots_getstate_setstate):

  __slots__ = ["n", "scale", "rms"]

  def __init__(O, n, scale, rms):
    O.n = n
    O.scale = scale
    O.rms = rms

class perm_rms_list(libtbx.slots_getstate_setstate):

  __slots__ = ["array", "i_small", "score"]

  def __init__(O, array, i_small=None, score=None):
    O.array = array
    O.i_small = i_small
    O.score = score

  def set_score(O):
    if (len(O.array) == 1):
      O.i_small = 0
      _ = O.array[0]
      O.score = _.n / (1 + _.rms)
    else:
      from scitbx.array_family import flex
      rms_list = flex.double([_.rms for _ in O.array])
      sort_perm = flex.sort_permutation(rms_list)
      O.i_small = sort_perm[0]
      i_2nd = sort_perm[1]
      rms_min, rms_2nd = [rms_list[_] for _ in sort_perm[:2]]
      O.score = (rms_2nd - rms_min) * (O.array[O.i_small].n + O.array[i_2nd].n)

def build_usables(work_params, reindexing_assistant, image_mdls):
  from scitbx.array_family import flex
  usable_fractions = flex.double()
  usables = []
  for i_img,im in enumerate(image_mdls.array):
    usable = im.usable(
      partiality_threshold=work_params.usable_partiality_threshold)
    usable_fractions.append(usable.miis.size() / im.miller_index_i_seqs.size())
    miis_perms = []
    for perm in reindexing_assistant.inv_perms:
      m = perm.select(usable.miis)
      p = flex.sort_permutation(data=m)
      miis_perms.append((m.select(p), usable.esti.select(p)))
    usables.append(miis_perms)
  print "Usable fraction of estimated image intensities:"
  usable_fractions.min_max_mean().show(prefix="  ")
  print
  sys.stdout.flush()
  return usables

class i_perm_and_scale(object):

  __slots__ = ["i_perm", "scale"]

  def __init__(O, i_perm=None, scale=None):
    O.i_perm = i_perm
    O.scale = scale

class cluster_info(object):

  __slots__ = ["i_perm_and_scale_by_i_img", "miis_perms", "esti_perms"]

  def __init__(O,
        i_perm_and_scale_by_i_img=None,
        miis_perms=None,
        esti_perms=None):
    O.i_perm_and_scale_by_i_img = i_perm_and_scale_by_i_img
    O.miis_perms = miis_perms
    O.esti_perms = esti_perms

  def build_cluster_pair_info(O, other, work_params, reindexing_assistant):
    from scitbx.array_family import flex
    scale_max = work_params.scale_estimation_scale_max
    assert scale_max > 0
    scale_min = 1/scale_max
    miis_i, esti_i = O.miis_perms[0], O.esti_perms[0]
    result = []
    for j_perm in xrange(len(reindexing_assistant.cb_ops)):
      miis_j, esti_j = other.miis_perms[j_perm], other.esti_perms[j_perm]
      i_seqs, j_seqs = miis_i.intersection_i_seqs(other=miis_j)
      if (i_seqs.size() < 2):
        return None
      x = esti_i.select(i_seqs)
      y = esti_j.select(j_seqs)
      if (((x != 0) | (y != 0)).count(True) < 2):
        return None
      num = flex.sum(x*y)
      den = flex.sum_sq(x)
      if (num > den * scale_min and num < den * scale_max):
        scale = num / den
        rms = flex.mean_sq(x*scale-y)**0.5
        result.append(perm_rms_info(n=x.size(), scale=scale, rms=rms))
      else:
        return None
    result = perm_rms_list(array=result)
    result.set_score()
    return result

  def merge(O, other, pair_info, reindexing_assistant, image_mdls):
    # TODO: refine combined scales so that rms for entire cluster
    #       is minimal then compute esti
    miis_i, esti_i = O.miis_perms[0], O.esti_perms[0]
    j_perm = pair_info.i_small
    miis_j, esti_j = other.miis_perms[j_perm], other.esti_perms[j_perm]
    scale_j = pair_info.array[j_perm].scale
    mrg_miis = miis_i.concatenate(miis_j)
    mrg_esti = esti_i.concatenate(esti_j * (1/scale_j))
    from scitbx.array_family import flex
    sort_perm = flex.sort_permutation(mrg_miis)
    mrg_miis = mrg_miis.select(sort_perm)
    mrg_esti = mrg_esti.select(sort_perm)
    new_miis = flex.size_t()
    new_esti = flex.double()
    n = mrg_miis.size()
    i = 0
    while (i < n):
      new_miis.append(mrg_miis[i])
      if (i+1 == n or mrg_miis[i] != mrg_miis[i+1]):
        new_esti.append(mrg_esti[i])
        i += 1
      else:
        new_esti.append((mrg_esti[i] + mrg_esti[i+1]) / 2)
        i += 2
    for i_img,i_perm_and_scale_ in other.i_perm_and_scale_by_i_img.items():
      O.i_perm_and_scale_by_i_img[i_img] = i_perm_and_scale(
        i_perm=reindexing_assistant.i_inv_j_multiplication_table[
          j_perm][
          i_perm_and_scale_.i_perm],
        scale=scale_j*i_perm_and_scale_.scale)
    O.miis_perms = []
    O.esti_perms = []
    for perm in reindexing_assistant.inv_perms:
      m = perm.select(new_miis)
      p = flex.sort_permutation(data=m)
      O.miis_perms.append(m.select(p))
      O.esti_perms.append(new_esti.select(p))

def build_image_cluster(work_params, reindexing_assistant, image_mdls, usables):
  n_imgs = len(usables)
  clusters = []
  for i_img,miis_perms in enumerate(usables):
    clusters.append(cluster_info(
      i_perm_and_scale_by_i_img={i_img: i_perm_and_scale(0, 1)},
      miis_perms=[_ for _,__ in miis_perms],
      esti_perms=[_ for __,_ in miis_perms]))
  remaining = range(n_imgs)
  cluster_pairs = [{} for _ in xrange(n_imgs)]
  def process_cp(i_rem, j_rem):
    i_clu = remaining[i_rem]
    j_clu = remaining[j_rem]
    cp = clusters[i_clu].build_cluster_pair_info(
      other=clusters[j_clu],
      work_params=work_params,
      reindexing_assistant=reindexing_assistant)
    if (cp is not None):
      cluster_pairs[i_clu][j_clu] = cp
  while (len(remaining) != 1):
    if (len(remaining) == n_imgs):
      chunk_size = 3000 # ad-hoc
      if (not work_params.multiprocessing or n_imgs*(n_imgs-1) <= chunk_size):
        import time
        time_start = time.time()
        for i_rem in xrange(n_imgs):
          for j_rem in xrange(i_rem+1, n_imgs):
            process_cp(i_rem, j_rem)
        from libtbx.utils import show_wall_clock_time
        show_wall_clock_time(seconds=time.time()-time_start)
      else:
        def mp():
          ij_list = []
          for i_rem in xrange(n_imgs):
            for j_rem in xrange(i_rem+1, n_imgs):
              ij_list.append((i_rem,j_rem))
          n_chunks = len(ij_list) // chunk_size
          print "Number of chunks for computing cluster pairs:", n_chunks
          print
          def process_chunk(i_chunk):
            for j_chunk in xrange(chunk_size):
              i = i_chunk * chunk_size + j_chunk
              if (i == len(ij_list)):
                break
              i_rem, j_rem = ij_list[i]
              process_cp(i_rem, j_rem)
            return cluster_pairs
          from libtbx import easy_mp
          mp_results = easy_mp.pool_map(
            fixed_func=process_chunk,
            args=range(n_chunks),
            chunksize=1,
            log=sys.stdout)
          for cps in mp_results:
            for main,sub in zip(cluster_pairs,cps):
              main.update(sub)
        mp()
    else:
      for i_rem in xrange(max_j_rem):
        i_clu = remaining[i_rem]
        cps_i = cluster_pairs[i_clu]
        if (max_j_clu in cps_i):
          del cps_i[max_j_clu]
        if (i_rem < max_i_rem):
          if (max_i_clu in cps_i):
            del cps_i[max_i_clu]
          process_cp(i_rem, max_i_rem)
      for j_rem in xrange(max_i_rem+1, len(remaining)):
        process_cp(max_i_rem, j_rem)
    max_score = 0
    max_i_rem = None
    max_j_clu = None
    for i_rem,i_clu in enumerate(remaining):
      cps_i = cluster_pairs[i_clu]
      for j_clu,cp in cps_i.items():
        if (max_score < cp.score):
          max_score = cp.score
          max_i_rem = i_rem
          max_j_clu = j_clu
    if (max_i_rem is None):
      raise RuntimeError("Insufficient connectivity between images.")
    max_i_clu = remaining[max_i_rem]
    max_j_rem = remaining.index(max_j_clu)
    print "max_score:", max_score, (max_i_rem, max_j_rem)
    cp = cluster_pairs[max_i_clu][max_j_clu]
    clusters[max_i_clu].merge(
      other=clusters[max_j_clu],
      pair_info=cp,
      reindexing_assistant=reindexing_assistant,
      image_mdls=image_mdls)
    cluster_pairs[max_j_clu] = None
    clusters[max_j_clu] = None
    del remaining[max_j_rem]
  return clusters[remaining[0]]

def check_image_cluster(
      work_params,
      i_calc,
      reindexing_assistant,
      image_mdls,
      scales_input,
      cluster):
  from scitbx.array_family import flex
  for i_perm in xrange(len(cluster.miis_perms)):
    expected = i_calc.select(cluster.miis_perms[i_perm])
    reconstr = expected.customized_copy(data=cluster.esti_perms[i_perm])
    print "i_perm:", i_perm
    flex.linear_correlation(x=expected.data(), y=reconstr.data()).show_summary()
    r1 = expected.r1_factor(other=reconstr, scale_factor=libtbx.Auto)
    print "r1: %.5f" % r1
    print
  for i_img,i_perm_and_scale in cluster.i_perm_and_scale_by_i_img.items():
    im = image_mdls.array[i_img]
    im.i_perm = i_perm_and_scale.i_perm
    im.scale = i_perm_and_scale.scale
  if (not work_params.index_and_integrate):
    image_mdls.check_i_perm_vs_backup(reindexing_assistant)
  cluster_scales = image_mdls.extract_scales()
  print "input vs. cluster scales:"
  flex.linear_correlation(x=scales_input, y=cluster_scales).show_summary()
  print

def show_i_calc_reindexing_correlations(i_calc, reindexing_assistant):
  assert i_calc.indices().all_eq(reindexing_assistant.miller_indices)
  assert i_calc.space_group_info().type().number() == 1
  assert i_calc.anomalous_flag()
  from scitbx.array_family import flex
  data = i_calc.data()
  print "I-calc reindexing correlations:"
  for cb_op,inv_perm in zip(
        reindexing_assistant.cb_ops,
        reindexing_assistant.inv_perms):
    i_calc_cb = i_calc.change_basis(cb_op)
    assert i_calc_cb.indices().select(inv_perm).all_eq(i_calc.indices())
    data_perm = i_calc_cb.data().select(inv_perm)
    print "  %-12s  %8.5f" % (
      cb_op.c().r().as_hkl(),
      flex.linear_correlation(data, data_perm).coefficient())
  print

def process_core(work_params, i_calc, reindexing_assistant, image_mdls):
  show_i_calc_reindexing_correlations(i_calc, reindexing_assistant)
  input_im0_i_perm = image_mdls.array[0].backup.i_perm
  if (work_params.check_refine_uc_cr):
    check_refine_uc_cr(work_params, image_mdls)
  scales_input = image_mdls.extract_scales()
  image_mdls.erase_scales()
  if (work_params.index_and_integrate):
    image_mdls.erase_spot_models()
    index_and_integrate(work_params, image_mdls)
    show_vm_info("After index_and_integrate():")
    isel = image_mdls.iselection_entries_with_spot_model()
    print "Removing %d image models for which" \
      " indexing or integration failed." % (image_mdls.size() - isel.size())
    scales_input = scales_input.select(isel)
    image_mdls = image_mdls.remove_all_entries_without_spot_model()
    print
  image_mdls.normalize_spot_intensities(target_mean=100)
  image_mdls.check_i_obs_vs_backup(work_params)
  image_mdls.reset_miller_image_map()
  image_mdls.miller_image_map.show_images_per_miller_index()
  image_mdls.reset_partialities(work_params)
  if (work_params.pickle_image_models and work_params.index_and_integrate):
    from libtbx import easy_pickle
    file_name = "%s_image_mdls_index_and_integrate.pickle" % str(
      work_params.base36_timestamp)
    easy_pickle.dump(
      file_name=file_name,
      obj=(work_params, image_mdls, reindexing_assistant))
    show_vm_info("After %s:" % file_name)
  if (work_params.write_image_models_to_mtz_files):
    image_mdls.write_to_mtz_files(common_unit_cell=work_params.unit_cell)
    show_vm_info("After write_image_models_to_mtz_files:")
  usables = build_usables(work_params, reindexing_assistant, image_mdls)
  image_cluster = build_image_cluster(
    work_params, reindexing_assistant, image_mdls, usables)
  show_vm_info("After build_image_cluster():")
  check_image_cluster(
    work_params, i_calc, reindexing_assistant, image_mdls,
    scales_input, image_cluster)
  cluster_scales = image_mdls.extract_scales()
  for im in image_mdls.array:
    im.reindex_in_place(reindexing_assistant)
  image_mdls.reset_miller_image_map()
  image_mdls.miller_image_map.show_images_per_miller_index()
  image_mdls.reset_partialities(work_params)
  from scitbx.array_family import flex
  def show_correlation_of_scales(assert_perfect=False):
    expected = scales_input / scales_input[0]
    estimated = image_mdls.extract_scales()
    print "Correlation of expected and estimated scales:"
    flex.linear_correlation(expected, estimated).show_summary(prefix="  ")
    print
    sys.stdout.flush()
    if (assert_perfect):
      from libtbx.test_utils import approx_equal
      assert approx_equal(estimated, expected)
  show_correlation_of_scales(
    assert_perfect=not work_params.index_and_integrate)
  indices, data = image_mdls.extract_estimated_i_obs(
    work_params.usable_partiality_threshold)
  i_obs_cluster = i_calc.customized_copy(indices=indices, data=data)
  refined_scales = None
  if (work_params.refine_scales.max_iterations in [None, 0]):
    print "refinement target: %.6g" % image_mdls.refinement_target(
      work_params.usable_partiality_threshold)
    print
  else:
    refined = refinery(work_params, image_mdls)
    refined.show_summary()
    show_correlation_of_scales()
    refined_scales = image_mdls.extract_scales()
  indices, data = image_mdls.extract_estimated_i_obs(
    work_params.usable_partiality_threshold)
  i_obs_est = i_calc.customized_copy(indices=indices, data=data)
  from libtbx import easy_pickle
  from libtbx import group_args
  easy_pickle.dump(
    file_name="%s_solver_results.pickle" % work_params.base36_timestamp,
    obj=group_args(
      work_params=work_params,
      i_calc=i_calc,
      reindexing_assistant=reindexing_assistant,
      scales_input=scales_input,
      cluster_scales=cluster_scales,
      refined_scales=refined_scales,
      i_obs_cluster=i_obs_cluster,
      i_obs_est=i_obs_est))
  print "Input I-calc:"
  i_calc.show_comprehensive_summary(prefix="  ")
  print
  print "Estimated I-obs:"
  i_obs_est.show_comprehensive_summary(prefix="  ")
  print
  print "input_im0_i_perm:", input_im0_i_perm
  print
  for i_perm,cb_op in enumerate(reindexing_assistant.cb_ops):
    c, e = i_calc.change_basis(cb_op).common_sets(i_obs_est)
    if (c.indices().size() > 2):
      print "Correlation of input and estimated I-obs (i_perm=%d):" % i_perm
      print "  Number of points:", c.data().size()
      corr = flex.linear_correlation(c.data(), e.data())
      corr.show_summary(prefix="  ")
      if (not work_params.index_and_integrate and i_perm == input_im0_i_perm):
        from libtbx.test_utils import approx_equal
        assert approx_equal(corr.coefficient(), 1)
      print
  return True

def process(work_params, i_calc):
  from cctbx.miller import reindexing
  reindexing_assistant = reindexing.assistant(
    lattice_group=work_params.lattice_symmetry.group(),
    intensity_group=work_params.intensity_symmetry.group(),
    miller_indices=i_calc.indices())
  reindexing_assistant.show_summary()
  print
  image_mdls = build_images(work_params, i_calc, reindexing_assistant)
  show_vm_info("After build_images():")
  if (work_params.pickle_image_models):
    file_name = "%s_image_mdls.pickle" % work_params.base36_timestamp
    from libtbx import easy_pickle
    easy_pickle.dump(
      file_name=file_name,
      obj=(work_params, i_calc, reindexing_assistant, image_mdls))
    show_vm_info("After %s:" % file_name)
  process_core(work_params, i_calc, reindexing_assistant, image_mdls)

def run_with_pickle(file_name):
  from libtbx import easy_pickle
  work_params, i_calc, reindexing_assistant, image_mdls = easy_pickle.load(
    file_name)
  work_params.phil_master.format(work_params).show()
  print
  i_calc.show_comprehensive_summary()
  print
  reindexing_assistant.show_summary()
  print
  show_vm_info("After unpickling:")
  process_core(work_params, i_calc, reindexing_assistant, image_mdls)

def run_fresh(args):
  import run_spotfinder
  work_params = run_spotfinder.process_args(
    args=args,
    extra_phil_str="""\
number_of_shots = 10
  .type = int
usable_partiality_threshold = 0.1
  .type = float
scale_estimation_scale_max = 1e3
  .type = float
min_number_of_spots_for_indexing = 16
  .type = int
sample_random_seeds = None
  .type = int
check_refine_uc_cr = False
  .type = bool
index_and_integrate = False
  .type = bool
show_refine_uc_cr = False
  .type = bool
apply_random_reindexing = False
  .type = bool
multiprocessing = False
  .type = bool
refine_scales {
  max_iterations = None
    .type = int
  finite_difference_eps = 1e-4
    .type = float
}
pickle_image_models = False
  .type = bool
write_image_models_to_mtz_files = False
  .type = bool
""")
  from create import build_i_calc
  i_calc, _ = build_i_calc(work_params)
  i_calc.show_comprehensive_summary()
  print
  show_vm_info("After build_i_calc:")
  if (work_params.sample_random_seeds is None):
    process(work_params, i_calc)
  else:
    for work_params.noise.random_seed in xrange(
          work_params.sample_random_seeds):
      process(work_params, i_calc)
  show_vm_info("Final:")

def run(args):
  import libtbx.utils
  libtbx.utils.show_times_at_exit()
  if (len(args) == 1):
    file_name = args[0]
    if (file_name.endswith(".pickle") and op.isfile(file_name)):
      return run_with_pickle(file_name)
  return run_fresh(args)
