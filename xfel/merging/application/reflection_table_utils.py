from __future__ import absolute_import, division, print_function
from six.moves import range
from dials.array_family import flex
import math
import pickle
import numpy as np
from scipy.stats import exponnorm
from libtbx.mpi4py import MPI


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
  def select_odd_experiment_reflections(reflections):
    'Select reflections from experiments with odd ids. An experiment id must be a string representing a hexadecimal number'
    exp_index_map = {exp_uid: i for i, exp_uid in enumerate(set(reflections["exp_id"]))}
    sel = [exp_index_map[exp_id] % 2 == 1 for exp_id in reflections["exp_id"]]
    sel = flex.bool(sel)
    reflections["is_odd_experiment"] = sel  # store this for later use, NOTE this is un-prunable if expanded_bookkeeping=True
    return reflections.select(sel)

  @staticmethod
  def select_even_experiment_reflections(reflections):
    'Select reflections from experiments with even ids. An experiment id must be a string representing a hexadecimal number'
    exp_index_map = {exp_uid: i for i, exp_uid in enumerate(set(reflections["exp_id"]))}
    sel = [exp_index_map[exp_id] % 2 == 0 for exp_id in reflections["exp_id"]]
    sel = flex.bool(sel)
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
  def merge_reflections(reflections, min_multiplicity, debug=True):
    '''Merge intensities of multiply-measured symmetry-reduced HKLs. The input reflection table must be sorted by symmetry-reduced HKLs.'''
    merged_reflections = reflection_table_utils.merged_reflection_table()
    for refls in reflection_table_utils.get_next_hkl_reflection_table(reflections=reflections):
      if refls.size() == 0:
        break # unless the input "reflections" list is empty, generated "refls" lists cannot be empty

      hkl = refls[0]['miller_index_asymmetric']
      # This assert is timeconsuming when using a small number of cores
      #assert not (hkl in merged_reflections['miller_index']) # i.e. assert that the input reflection table came in sorted

      refls = refls.select(refls['intensity.sum.variance'] > 0.0)
      sel_low = refls['intensity.sum.value'] >= np.percentile(refls['intensity.sum.value'], 35)
      sel_high = refls['intensity.sum.value'] <= np.percentile(refls['intensity.sum.value'], 85)
      sel = sel_low & sel_high
      #refls = refls.select(sel)

      if refls.size() >= min_multiplicity:
        weighted_intensity_array = refls['intensity.sum.value'] / refls['intensity.sum.variance']
        weights_array = flex.double(refls.size(), 1.0) / refls['intensity.sum.variance']

        weighted_mean_intensity = flex.sum(weighted_intensity_array) / flex.sum(weights_array)
        log_i = math.log(max(weighted_mean_intensity, 1), 10)
        log_i = max(log_i, 1)
        pct = 17 + 15*log_i
        pct = min(pct, 72)
        iarray = refls['intensity.sum.value']
        fitpars = exponnorm.fit(iarray)
        try:
          weighted_mean_intensity = exponnorm.ppf(0.95, *fitpars)
        except:
          print(fitpars)
        #weighted_mean_intensity = np.percentile(refls['intensity.sum.value'], pct)
        standard_error_of_weighted_mean_intensity = 1.0/math.sqrt(flex.sum(weights_array))
        if debug:
          #pass
          _h,_k,_l = hkl
          outname = 'pickles2/refl.%s.%s.%s.pkl' % (_h, _k, _l)
          with open(outname, 'wb') as f: pickle.dump(refls, f)


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
