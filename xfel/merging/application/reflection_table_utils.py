from __future__ import absolute_import, division, print_function
import numpy as np
from six.moves import range
from dials.array_family import flex
import hashlib
import math
from xfel.merging.application.utils.outlier import is_outlier

class reflection_table_utils(object):

  @staticmethod
  def get_next_hkl_reflection_table(reflections):
    '''Generate asu hkl slices from an asu hkl-sorted reflection table'''
    if reflections.size() == 0:
      yield reflections

    i_begin = 0
    hkl_ref = reflections['miller_index_asymmetric'][0]
    for i in range(reflections.size()):
      hkl = reflections['miller_index_asymmetric'][i]
      if hkl == hkl_ref:
        continue
      else:
        yield reflections[i_begin:i]
        i_begin = i
        hkl_ref = hkl

    yield reflections[i_begin:i+1]

  @staticmethod
  def is_even(identifier):
    return ord(hashlib.md5(identifier.encode()).hexdigest()[-1]) % 2 == 0

  @staticmethod
  def select_odd_experiment_reflections(reflections, col='id'):
    'Select reflections from experiments with odd identifiers.'
    id_map = reflections.experiment_identifiers()
    reverse_map = {id_map[id_]: id_ for id_ in id_map.keys()}
    sel = flex.bool(len(reflections), False)
    for identifier in reverse_map:
      if not reflection_table_utils.is_even(identifier):
        sel |= reflections[col] == reverse_map[identifier]
    reflections["is_odd_experiment"] = sel  # store this for later use, NOTE this is un-prunable if expanded_bookkeeping=True
    return reflections.select(sel)

  @staticmethod
  def select_even_experiment_reflections(reflections, col='id'):
    'Select reflections from experiments with even identifiers'
    id_map = reflections.experiment_identifiers()
    reverse_map = {id_map[id_]: id_ for id_ in id_map.keys()}
    sel = flex.bool(len(reflections), False)
    for identifier in reverse_map:
      if reflection_table_utils.is_even(identifier):
        sel |= reflections[col] == reverse_map[identifier]
    return reflections.select(sel)

  @staticmethod
  def merged_reflection_table():
    '''Create a reflection table for storing merged HKLs'''
    table = flex.reflection_table()
    table['miller_index'] = flex.miller_index()
    table['intensity'] = flex.double()
    table['sigma'] = flex.double()
    table['multiplicity'] = flex.int()
    return table

  @staticmethod
  def merge_reflections(reflections, min_multiplicity, nameprefix=None, thresh=None):
    '''Merge intensities of multiply-measured symmetry-reduced HKLs. The input reflection table must be sorted by symmetry-reduced HKLs.'''
    merged_reflections = reflection_table_utils.merged_reflection_table()
    for i_refls,refls in enumerate(reflection_table_utils.get_next_hkl_reflection_table(reflections=reflections)):
      if refls.size() == 0:
        break # unless the input "reflections" list is empty, generated "refls" lists cannot be empty

      hkl = refls[0]['miller_index_asymmetric']
      # This assert is timeconsuming when using a small number of cores
      #assert not (hkl in merged_reflections['miller_index']) # i.e. assert that the input reflection table came in sorted

      refls = refls.select(refls['intensity.sum.variance'] > 0.0)

      if refls.size() >= min_multiplicity:
        weighted_intensity_array = refls['intensity.sum.value'] / refls['intensity.sum.variance']
        weights_array = flex.double(refls.size(), 1.0) / refls['intensity.sum.variance']

        weighted_mean_intensity = flex.sum(weighted_intensity_array) / flex.sum(weights_array)
        standard_error_of_weighted_mean_intensity = 1.0/math.sqrt(flex.sum(weights_array))

        if thresh is not None:
          vals = refls["intensity.sum.value"].as_numpy_array()
          good = ~is_outlier(vals, thresh)
          good_vals = vals[good]
          weighted_mean_intensity = np.mean(good_vals)
          vals_var = refls["intensity.sum.variance"].as_numpy_array()
          num_good = good.sum()
          standard_error_of_weighted_mean_intensity = np.sqrt(vals_var[good].sum())/num_good
        merged_reflections.append(
                                  {'miller_index' : hkl,
                                  'intensity' : weighted_mean_intensity,
                                  'sigma' : standard_error_of_weighted_mean_intensity,
                                  'multiplicity' : refls.size()})
    return merged_reflections

  @staticmethod
  def prune_reflection_table_keys(reflections, keys_to_delete=None, keys_to_keep=None,
                                  keys_to_ignore=None):
    '''Remove reflection table keys: either inclusive or exclusive, columns in keys_to_ignore will always remain'''
    # These columns were created by the merging application, and we want to retain them
    if keys_to_delete is not None:
      keys_to_delete = [k for k in keys_to_delete if k not in keys_to_ignore]
    if keys_to_keep is not None:
      keys_to_keep += [k for k in keys_to_ignore if k not in keys_to_keep]

    if len(reflections) != 0:
      all_keys = list()
      for key in reflections.keys():
        all_keys.append(key)
      if keys_to_delete != None:
        for key in keys_to_delete:
          if key in all_keys:
            del reflections[key]
      elif keys_to_keep != None:
        for key in all_keys:
          if not key in keys_to_keep:
            del reflections[key]
    return reflections

  @staticmethod
  def get_next_reflection_table_slice(reflections, n_slices, reflection_table_stub):
    '''Generate an exact number of slices from a reflection table. Make slices as even as possible. If not enough reflections, generate empty tables'''
    assert n_slices >= 0

    if n_slices == 1:
      yield reflections
    else:
      import math

      generated_slices = 0
      count = len(reflections)

      if count > 0:
        # how many non-empty slices should we generate and with what stride?
        nonempty_slices = min(count, n_slices)
        stride = int(math.ceil(count / nonempty_slices))

        # generate all non-empty slices
        for i in range(0, count, stride):
          generated_slices += 1
          i2 = i + stride
          if generated_slices == nonempty_slices:
            i2 = count
          yield reflections[i:i2]

      # generate some empty slices if necessary
      empty_slices = max(0, n_slices - generated_slices)
      for i in range(empty_slices):
        yield reflection_table_stub(reflections)

  @staticmethod
  def iterate_experiments_and_load_imagesets(params, experiments, reflections = None):

    # Re-generate the image sets using their format classes so we can read the raw data
    current_imageset = None
    current_imageset_path = None

    if self.params.mp.psana2_mode:
      paths_ = list(set([p for iset in experiments.imagesets() for p in iset.paths()]))
      paths = list(set([p for plist in self.mpi_helper.comm.allgather(paths_) for p in plist]))

      for path in paths:
        current_imageset_path = paths[0]
        current_imageset = ImageSetFactory.make_imageset([current_imageset_path])
      
        for expt_id, expt in enumerate(experiments):
          assert len(expt.imageset.paths()) == 1 and len(expt.imageset) == 1
          if expt.imageset.paths()[0] != path: continue
          
          idx = expt.imageset.indices()[0]
          expt.imageset = current_imageset[idx:idx+1]

          if reflections:
            refls = reflections.select(reflections['id'] == expt_id)
            refls['id'] = flex.int(len(refls), 0)
            idents = refls.experiment_identifiers()
            del idents[expt_id]
            idents[0] = expt.identifier
          else:
            refls = None

          yield expt, refls

    else:
      for expt_id, expt in enumerate(experiments):
        assert len(expt.imageset.paths()) == 1 and len(expt.imageset) == 1

        refls = reflections.select(reflections['id'] == expt_id)

        if expt.imageset.paths()[0] != current_imageset_path:
          current_imageset_path = expt.imageset.paths()[0]
          current_imageset = ImageSetFactory.make_imageset(expt.imageset.paths())

        idx = expt.imageset.indices()[0]
        expt.imageset = current_imageset[idx:idx+1]
        idents = refls.experiment_identifiers()
        del idents[expt_id]
        idents[0] = expt.identifier
        refls['id'] = flex.int(len(refls), 0)

        yield expt, refls

