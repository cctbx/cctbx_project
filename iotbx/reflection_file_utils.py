"""Tools for extracting information from reflection files
"""
from __future__ import absolute_import, division, print_function
from iotbx import reflection_file_reader
from cctbx import miller
from cctbx.array_family import flex
import libtbx.path
from libtbx.str_utils import show_string
from libtbx.utils import Sorry, null_out
import libtbx.phil
from itertools import count
import math
import sys, os
from functools import cmp_to_key
from six.moves import range
from six.moves import zip

class Sorry_No_array_of_the_required_type(Sorry):
  __orig_module__ = __module__
  __module__ = "exceptions"

class Sorry_Not_a_suitable_array(Sorry):
  __orig_module__ = __module__
  __module__ = "exceptions"

def find_labels(search_labels, info_string):
  for search_label in search_labels:
    if (info_string.find(search_label) < 0):
      return False
  return True

class label_table(object):

  def __init__(self, miller_arrays, err=None):
    self.miller_arrays = miller_arrays
    if (err is None): self.err = sys.stderr
    else: self.err = err
    self.info_strings = []
    self.info_labels = []
    for p_array,miller_array in zip(count(1), miller_arrays):
      info = miller_array.info()
      if (info is not None):
        self.info_strings.append(str(info))
      else:
        self.info_strings.append(str(p_array))
      self.info_labels.append(getattr(info, "labels", None))

  def scores(self, label=None, labels=None):
    assert [label, labels].count(None) == 1
    if (labels is None):
      labels = [label]
    else:
      assert len(labels) > 0
    if(type(labels)==type("")): labels = [labels]
    result = []
    labels_lower = [lbl.lower() for lbl in labels]
    for info_string,info_labels in zip(self.info_strings, self.info_labels):
      if (not find_labels(
               search_labels=labels_lower,
               info_string=info_string.lower())):
        result.append(0)
      elif (not find_labels(
                  search_labels=labels,
                  info_string=info_string)):
        result.append(1)
      else:
        n_exact_matches = 0
        if (info_labels is not None):
          for info_label in info_labels:
            if (info_label in labels):
              n_exact_matches += 1
        result.append(2 + n_exact_matches)
    return result

  def show_possible_choices(self,
        f=None,
        scores=None,
        minimum_score=None,
        parameter_name=None):
    if (f is None): f = self.err
    print("Possible choices:", file=f)
    if (scores is None):
      for info_string in self.info_strings:
        print(" ", info_string, file=f)
    else:
      for info_string,score in zip(self.info_strings, scores):
        if (score >= minimum_score):
          print(" ", info_string, file=f)
    print(file=f)
    if (parameter_name is None): hint = ""
    else: hint = "use %s\nto " % parameter_name
    print("Please %sspecify an unambiguous substring of the target label." % hint, file=f)
    print(file=f)

  def match_data_label(self, label, command_line_switch, f=None):
    if (f is None): f = self.err
    assert label is not None
    scores = self.scores(label=label)
    selected_array = None
    for high_score in range(max(scores),0,-1):
      if (scores.count(high_score) > 0):
        if (scores.count(high_score) > 1):
          print(file=f)
          print("Ambiguous %s=%s" % (command_line_switch, label), file=f)
          print(file=f)
          self.show_possible_choices(
            f=f, scores=scores, minimum_score=high_score)
          return None
        return self.miller_arrays[scores.index(high_score)]
    print(file=f)
    print("Unknown %s=%s" % (command_line_switch, label), file=f)
    print(file=f)
    self.show_possible_choices(f=f)
    return None

  def select_array(self, label, command_line_switch, f=None):
    if (f is None): f = self.err
    if (len(self.miller_arrays) == 0):
      print(file=f)
      print("No reflection arrays available.", file=f)
      print(file=f)
      return None
    if (len(self.miller_arrays) == 1):
      return self.miller_arrays[0]
    if (label is None):
      print(file=f)
      s = command_line_switch
      print("Please use %s to select a reflection array." % s, file=f)
      print("For example: %s=%s" % (s, show_string(str(
        self.miller_arrays[1].info()).split(":")[-1])), file=f)
      print(file=f)
      self.show_possible_choices(f=f)
      return None
    return self.match_data_label(
      label=label,
      command_line_switch=command_line_switch)

def presumably_from_mtz_FQDQY(miller_array):
  lbls = miller_array.info().labels
  return (lbls is not None and len(lbls) == 5)

def get_amplitude_scores(miller_arrays, strict=False):
  result = []
  for miller_array in miller_arrays:
    score = 0
    if (miller_array.is_complex_array()) and (not strict):
      score = 1
    elif (miller_array.is_real_array()):
      if (miller_array.is_xray_reconstructed_amplitude_array()):
        if (presumably_from_mtz_FQDQY(miller_array)):
          score = 3
        else:
          score = 2
      elif (miller_array.is_xray_amplitude_array()):
        score = 5
      elif (miller_array.is_xray_intensity_array()) and (not strict):
        score = 4
    result.append(score)
  return result

def get_phase_scores(miller_arrays):
  result = []
  for miller_array in miller_arrays:
    score = 0
    if (   miller_array.is_complex_array()
        or miller_array.is_hendrickson_lattman_array()):
      score = 4
    elif (miller_array.is_real_array()):
      if (miller_array.is_xray_reconstructed_amplitude_array()):
        pass
      elif (miller_array.is_xray_amplitude_array()):
        pass
      elif (miller_array.is_xray_intensity_array()):
        pass
      elif (miller_array.data().size() == 0):
        pass
      else:
        m = flex.mean(flex.abs(miller_array.data()))
        if (m < 5):
          score = 2
        elif (m < 500):
          score = 3
        else:
          score = 1
    result.append(score)
  return result

def get_xray_data_scores(miller_arrays, ignore_all_zeros,
    prefer_anomalous=None, prefer_amplitudes=None):
  anomalous_bonus = 4
  intensity_bonus = 2
  amplitude_bonus = 0
  if (prefer_amplitudes):
    intensity_bonus = 0
    amplitude_bonus = 2
  result = []
  for miller_array in miller_arrays:
    if (not miller_array.is_real_array()):
      result.append(0)
    else:
      score = None
      if (miller_array.data().all_eq(0)):
        if (ignore_all_zeros):
          score = 0
        else:
          score = 1
      elif (miller_array.is_xray_intensity_array()):
        score = 8 + intensity_bonus
        if (prefer_anomalous is not None):
          if (((prefer_anomalous) and (miller_array.anomalous_flag())) or
              ((not prefer_anomalous) and (not miller_array.anomalous_flag()))):
            score += anomalous_bonus
          else :
            score -= anomalous_bonus
      elif (miller_array.is_xray_amplitude_array()):
        if (miller_array.is_xray_reconstructed_amplitude_array()):
          if (presumably_from_mtz_FQDQY(miller_array)):
            if (prefer_anomalous):
              score = 8 + amplitude_bonus
            else :
              score = 6 + amplitude_bonus
          else:
            score = 4
        else:
          score = 6 + amplitude_bonus
          if (prefer_anomalous is not None):
            if (((prefer_anomalous) and (miller_array.anomalous_flag())) or
                ((not prefer_anomalous) and
                 (not miller_array.anomalous_flag()))):
              score += anomalous_bonus
            else :
              score -= anomalous_bonus
      else:
        score = 2
      assert score is not None
      if (    miller_array.sigmas() is not None
          and isinstance(miller_array.sigmas(), flex.double)
          and miller_array.sigmas_are_sensible()
          and score != 0):
        score += 1
      if (prefer_anomalous is not None) and (miller_array.anomalous_flag()):
        if (prefer_anomalous):
          score += 1
        else :
          score -= 1
      result.append(score)
  return result

def looks_like_r_free_flags_info(array_info):
  if (not isinstance(array_info, miller.array_info)): return False
  if (len(array_info.labels) > 2): return False
  label = array_info.labels[-1].lower()
  for word in ["free", "test", "cross", "status", "flag"]:
    if (label.find(word) >= 0): return True
  return False

class get_r_free_flags_score(object):

  def __init__(self, test_flag_value, n, n_free, miller_array_info):
    if (test_flag_value is not None or n_free < n*0.50):
      self.reversed = False
    else:
      self.reversed = True
      n_free = n - n_free
    self.flag_score = 0
    if (min(1000,n*0.01) < n_free < n*0.35):
      if (   looks_like_r_free_flags_info(miller_array_info)
          or min(2000,n*0.04) < n_free < n*0.20):
        self.flag_score = 3
      else:
        self.flag_score = 2

class get_r_free_flags_scores(object):

  def __init__(self, miller_arrays, test_flag_value):
    self.scores = []
    self.test_flag_values = []
    for miller_array in miller_arrays:
      flag_score = 0
      effective_test_flag_value = None
      data = miller_array.data()
      if (miller_array.is_bool_array()):
        trial_test_flag_value = (
          test_flag_value is None or bool(test_flag_value))
        n_free = data.count(trial_test_flag_value)
        scoring = get_r_free_flags_score(
          test_flag_value=test_flag_value,
          n=data.size(),
          n_free=n_free,
          miller_array_info=miller_array.info())
        if (scoring.flag_score != 0):
          flag_score = scoring.flag_score
          if (scoring.reversed):
            trial_test_flag_value = not trial_test_flag_value
          effective_test_flag_value = trial_test_flag_value
      elif (miller_array.is_integer_array()):
        try: counts = data.counts(max_keys=200)
        except RuntimeError: pass
        else:
          c_keys = list(counts.keys())
          c_values = list(counts.values())
          if (   test_flag_value is None
              or test_flag_value in c_keys):
            if (counts.size() == 2):
              if (test_flag_value is None):
                if (c_values[1] < c_values[0]):
                  i_free = 1
                else:
                  i_free = 0
              elif (test_flag_value == c_keys[0]):
                i_free = 0
              else:
                i_free = 1
              scoring = get_r_free_flags_score(
                test_flag_value=test_flag_value,
                n=data.size(),
                n_free=c_values[i_free],
                miller_array_info=miller_array.info())
              if (scoring.flag_score != 0):
                flag_score = scoring.flag_score
                if (scoring.reversed): i_free = 1-i_free
                effective_test_flag_value = c_keys[i_free]
            elif (counts.size() >= 3):
              c_keys_min = min(c_keys)
              c_keys_max = max(c_keys)
              if (((c_keys_max - c_keys_min) < data.size()) and
                  (c_keys == list(range(c_keys_min, c_keys_max+1)))):
                # XXX 0.55 may be too close a margin - the routine to export
                # R-free flags for CCP4 seems to get this wrong frequently.
                if (min(c_values) > max(c_values)*0.55):
                  if (looks_like_r_free_flags_info(miller_array.info())):
                    flag_score = 3
                  else:
                    flag_score = 2
                elif (looks_like_r_free_flags_info(miller_array.info())):
                  flag_score = 1
                if (flag_score != 0):
                  if (test_flag_value is None):
                    c_keys.sort()
                    # XXX this appears to be a special case (from Axel); I'm
                    # erring on the side of consistency with CNS here
                    if ((c_keys == [-1, 0, 1]) and
                        (counts[0] >= data.size()*0.55)):
                      effective_test_flag_value = 1
                    else :
                      effective_test_flag_value = min(c_keys)
                  else:
                    effective_test_flag_value = test_flag_value
              else : # XXX gross fix to avoid memory leak for corrupted flags
                n_binary = counts.get(0, 0) + counts.get(1, 0)
                if (n_binary >= data.size()*0.95):
                  if (looks_like_r_free_flags_info(miller_array.info())):
                    flag_score = 2
                  else :
                    flag_score = 1
                if (flag_score != 0):
                  if (test_flag_value is None):
                    if (counts[0] >= data.size()*0.55):
                      effective_test_flag_value = 1
                    else :
                      effective_test_flag_value = min(c_keys)
                  else :
                    effective_test_flag_value = test_flag_value
      elif miller_array.is_string_array():
        trial_test_flag_value = "f"
        n_free = data.count(trial_test_flag_value)
        scoring = get_r_free_flags_score(
          test_flag_value=test_flag_value,
          n=data.size(),
          n_free=n_free,
          miller_array_info=miller_array.info())
        if (scoring.flag_score != 0):
          flag_score = scoring.flag_score
          if (scoring.reversed):
            trial_test_flag_value = not trial_test_flag_value
          effective_test_flag_value = trial_test_flag_value
      self.scores.append(flag_score)
      self.test_flag_values.append(effective_test_flag_value)
    assert len(self.scores) == len(miller_arrays)
    assert len(self.test_flag_values) == len(miller_arrays)

def get_experimental_phases_scores(miller_arrays, ignore_all_zeros):
  result = []
  for miller_array in miller_arrays:
    if (miller_array.is_hendrickson_lattman_array()):
      if (not miller_array.data().all_eq((0,0,0,0))):
        result.append(2)
      elif (not ignore_all_zeros):
        result.append(1)
      else:
        result.append(0)
    else:
      result.append(0)
  return result

def sort_arrays_by_score(miller_arrays, array_scores, minimum_score):
  assert len(miller_arrays) == len(array_scores)
  i = 0
  scored_arrays = []
  for array in miller_arrays :
    scored_arrays.append( (array, array_scores[i]) )
    i += 1
  def cmp_fn(x,y):
    return y[1] - x[1]
  scored_arrays.sort(key=cmp_to_key(cmp_fn))
  valid_arrays = []
  for (array, score) in scored_arrays :
    if score >= minimum_score :
      valid_arrays.append(array)
  return valid_arrays

def select_array(
      parameter_name,
      labels,
      miller_arrays,
      data_scores,
      err,
      error_message_no_array,
      error_message_not_a_suitable_array,
      error_message_multiple_equally_suitable,
      raise_no_array=True):
  if (labels is not None): assert parameter_name is not None
  if (len(miller_arrays) == 0):
    raise Sorry_No_array_of_the_required_type(
      "No reflection arrays available.")
  if (data_scores is not None):
    assert max(data_scores) >= 0
  else:
    data_scores = [1]*len(miller_arrays)
  lbl_tab = label_table(miller_arrays=miller_arrays, err=err)
  if (labels is None):
    label_scores = None
  else:
    label_scores = lbl_tab.scores(labels=labels)
  if (label_scores is not None and max(label_scores) == 0):
    error = "No matching array: %s=%s" % (parameter_name, " ".join(labels))
    print("\n" + error + "\n", file=err)
    if (max(data_scores) > 0):
      lbl_tab.show_possible_choices(
        scores=data_scores,
        minimum_score=1,
        parameter_name=parameter_name)
    raise Sorry(error)
  if (max(data_scores) == 0):
    if (label_scores is None):
      if(raise_no_array):
        print("\n" + error_message_no_array + "\n", file=err)
        raise Sorry_No_array_of_the_required_type(error_message_no_array)
      else: return
    error = "%s%s=%s" % (
      error_message_not_a_suitable_array, parameter_name, " ".join(labels))
    print("\n" + error + "\n", file=err)
    raise Sorry_Not_a_suitable_array(error)
  if (label_scores is None):
    combined_scores = data_scores
  else:
    n = max(data_scores) + 1
    combined_scores = []
    for label_score,data_score in zip(label_scores, data_scores):
      combined_scores.append(label_score*n+data_score)
  i = combined_scores.index(max(combined_scores))
  if (combined_scores.count(combined_scores[i]) > 1):
    error = error_message_multiple_equally_suitable
    print("\n" + error + "\n", file=err)
    lbl_tab.show_possible_choices(
      scores=combined_scores,
      minimum_score=max(combined_scores),
      parameter_name=parameter_name)
    raise Sorry(error)
  return i

class reflection_file_server(object):

  def __init__(self,
        crystal_symmetry=None,
        force_symmetry=None,
        reflection_files=None,
        miller_arrays=None,
        err=None):
    self.crystal_symmetry = crystal_symmetry
    self.force_symmetry = force_symmetry
    if (err is None): self.err = sys.stderr
    else: self.err = err
    self.miller_arrays = []
    if (reflection_files is not None):
      for reflection_file in reflection_files:
        self.miller_arrays.extend(reflection_file.as_miller_arrays(
          crystal_symmetry=self.crystal_symmetry,
          force_symmetry=self.force_symmetry))
    if (miller_arrays is not None):
      self.miller_arrays.extend(miller_arrays)
    self.file_name_miller_arrays = {}
    for miller_array in self.miller_arrays:
      self.file_name_miller_arrays.setdefault(
        libtbx.path.canonical_path(
          miller_array.info().source), []).append(miller_array)

  def update_crystal_symmetry(self, crystal_symmetry):
    self.crystal_symmetry = crystal_symmetry
    for i, ma in enumerate(self.miller_arrays):
      info = ma.info()
      self.miller_arrays[i] = ma.customized_copy(
        crystal_symmetry = self.crystal_symmetry)
      self.miller_arrays[i].set_info(info)

  def get_miller_arrays(self, file_name):
    if (file_name is None): return self.miller_arrays
    canonical_file_name = libtbx.path.canonical_path(file_name)
    result = self.file_name_miller_arrays.get(canonical_file_name, None)
    if (result is None and hasattr(os.path, "samefile")):
      for tabulated_file_name in self.file_name_miller_arrays.keys():
        if (os.path.samefile(canonical_file_name, tabulated_file_name)):
          result = self.file_name_miller_arrays[canonical_file_name] \
                 = self.file_name_miller_arrays[tabulated_file_name]
          break
    if (result is None):
      reflection_file = reflection_file_reader.any_reflection_file(
        file_name=file_name)
      if (reflection_file.file_type() is None):
        self.file_name_miller_arrays[canonical_file_name] = None
      else:
        result = self.file_name_miller_arrays[canonical_file_name] \
               = reflection_file.as_miller_arrays(
                   crystal_symmetry=self.crystal_symmetry,
                   force_symmetry=self.force_symmetry)
    if (result is None):
      raise Sorry("No reflection data in file: %s" % file_name)
    return result

  def get_miller_array(self, labels, file_name=None):
    if (file_name is None):
      miller_arrays = self.miller_arrays
    else :
      canonical_file_name = libtbx.path.canonical_path(file_name)
      miller_arrays = self.file_name_miller_arrays[canonical_file_name]
    for array in miller_arrays :
      if (isinstance(labels, str)):
        if (array.info().label_string() == labels):
          return array
      else :
        assert (isinstance(labels, list))
        if (array.info().labels == labels):
          return array
    return None

  def get_amplitudes(self,
        file_name,
        labels,
        convert_to_amplitudes_if_necessary,
        parameter_scope,
        parameter_name,
        return_all_valid_arrays=False,
        minimum_score=1,
        strict=False):
    miller_arrays = self.get_miller_arrays(file_name=file_name)
    data_scores = get_amplitude_scores(miller_arrays=miller_arrays,
      strict=strict)
    if (parameter_scope is not None):
      parameter_name = parameter_scope + "." + parameter_name
    if return_all_valid_arrays :
      return sort_arrays_by_score(miller_arrays, data_scores, minimum_score)
    i = select_array(
      parameter_name=parameter_name,
      labels=labels,
      miller_arrays=miller_arrays,
      data_scores=data_scores,
      err=self.err,
      error_message_no_array
        ="No array of amplitudes found.",
      error_message_not_a_suitable_array
        ="Not a suitable array of amplitudes: ",
      error_message_multiple_equally_suitable
        ="Multiple equally suitable arrays of amplitudes found.")
    result = miller_arrays[i]
    if (convert_to_amplitudes_if_necessary):
      info = result.info()
      if (info is None):
        info_labels = None
      else:
        info_labels = info.labels
      if (result.is_complex_array()):
        result = result.amplitudes()
        if (info_labels is not None):
          result.set_info(info.customized_copy(labels=info_labels[:1]))
        else:
          result.set_info(info=info)
      elif (result.is_xray_intensity_array()):
        result = result.as_amplitude_array()
        if (info_labels is not None):
          result.set_info(info.customized_copy(
            labels=info_labels[:1]+["as_amplitude_array"]))
        else:
          result.set_info(info=info)
    return result

  def get_phases_deg(self,
        file_name,
        labels,
        convert_to_phases_if_necessary,
        original_phase_units,
        parameter_scope,
        parameter_name,
        return_all_valid_arrays=False,
        minimum_score=1):
    assert original_phase_units in [None, "deg", "rad"]
    miller_arrays = self.get_miller_arrays(file_name=file_name)
    data_scores = get_phase_scores(miller_arrays=miller_arrays)
    if (parameter_scope is not None):
      parameter_name = parameter_scope + "." + parameter_name
    if return_all_valid_arrays :
      return sort_arrays_by_score(miller_arrays, data_scores, minimum_score)
    i = select_array(
      parameter_name=parameter_name,
      labels=labels,
      miller_arrays=miller_arrays,
      data_scores=data_scores,
      err=self.err,
      error_message_no_array
        ="No array of phases found.",
      error_message_not_a_suitable_array
        ="Not a suitable array of phases: ",
      error_message_multiple_equally_suitable
        ="Multiple equally suitable arrays of phases found.")
    result = miller_arrays[i]
    info = result.info()
    if (info is None):
      info_labels = None
    else:
      info_labels = info.labels
    if (convert_to_phases_if_necessary):
      if (result.is_complex_array()):
        result = result.phases(deg=True)
        if (info_labels is not None and len(info_labels) == 2):
          result.set_info(info.customized_copy(labels=[info_labels[1]]))
        else:
          result.set_info(info=info)
      elif (result.is_hendrickson_lattman_array()):
        result = result.phase_integrals().phases(deg=True)
        if (info_labels is not None):
          result.set_info(info.customized_copy(
            labels=info_labels+["converted_to_centroid_phases"]))
        else:
          result.set_info(info=info)
    elif (    not result.is_complex_array()
          and original_phase_units == "rad"):
      result = result.customized_copy(data=result.data()*(180/math.pi))
      if (info_labels is not None):
        result.set_info(info.customized_copy(
          labels=info_labels+["converted_to_deg"]))
      else:
        result.set_info(info=info)
    return result

  def get_xray_data(self,
        file_name,
        labels,
        ignore_all_zeros,
        parameter_scope,
        parameter_name="labels",
        return_all_valid_arrays=False,
        minimum_score=1,
        prefer_anomalous=None,
        prefer_amplitudes=None):
    miller_arrays = self.get_miller_arrays(file_name=file_name)
    data_scores = get_xray_data_scores(
      miller_arrays=miller_arrays,
      ignore_all_zeros=ignore_all_zeros,
      prefer_anomalous=prefer_anomalous,
      prefer_amplitudes=prefer_amplitudes)
    if return_all_valid_arrays :
      return sort_arrays_by_score(miller_arrays, data_scores, minimum_score)
    # Recognize phenix.refine file and do the "right thing". May be too ad hoc..
    # XXX can we just check the first label instead?
    new_miller_arrays = miller_arrays
    new_data_scores = data_scores
    if(labels is None):
      new_miller_arrays = []
      new_data_scores = []
      # This allows still use "filtered" if no other option is available
      for j in range(len(miller_arrays)):
        ma = miller_arrays[j]
        ds = data_scores[j]
        if((ma.info().labels in [
           ['F-obs-filtered', 'SIGF-obs-filtered'],
           ['F-obs-filtered(+)', 'SIGF-obs-filtered(+)', 'F-obs-filtered(-)',
            'SIGF-obs-filtered(-)',],
            ]) ):
          data_scores[j] = ds-1
      #
      for ma, ds in zip(miller_arrays, data_scores):
        if(not ((ma.info().labels in [
           ['F-model', 'PHIF-model'],
           ['F-model(+)', 'PHIF-model(+)', 'F-model(-)', 'PHIF-model(-)'] ]) or
           isinstance(ma.data(), flex.complex_double))):
          new_miller_arrays.append(ma)
          new_data_scores.append(ds)
    #
    parameter_name = parameter_name.strip()
    if(len(parameter_name)==0):
      parameter_name_ = parameter_scope
    else:
      parameter_name_ = parameter_scope+"."+parameter_name
    i = select_array(
      parameter_name=parameter_name_,
      labels=labels,
      miller_arrays=new_miller_arrays,
      data_scores=new_data_scores,
      err=self.err,
      error_message_no_array
        ="No array of observed xray data found.",
      error_message_not_a_suitable_array
        ="Not a suitable array of observed xray data: ",
      error_message_multiple_equally_suitable
        ="Multiple equally suitable arrays of observed xray data found.")
    return new_miller_arrays[i]

  def get_r_free_flags(self,
        file_name,
        label,
        test_flag_value,
        disable_suitability_test,
        parameter_scope,
        return_all_valid_arrays=False,
        minimum_score=1):
    miller_arrays = self.get_miller_arrays(file_name=file_name)
    if (disable_suitability_test):
      if (label is None or test_flag_value is None):
        raise Sorry((
          "%s=True: Suitability test for R-free flags can only be disabled"
          " if both %s and %s are defined.") % (
            parameter_scope+".disable_suitability_test",
            parameter_scope+".label",
            parameter_scope+".test_flag_value"))
      elif return_all_valid_arrays :
        raise Sorry("return_all_valid_arrays=True: Suitability test can not "+
          "be disabled in this mode.")
      flag_scores = None
      data_scores = None
    else:
      flag_scores = get_r_free_flags_scores(
        miller_arrays=miller_arrays,
        test_flag_value=test_flag_value)
      data_scores = flag_scores.scores
    if return_all_valid_arrays : # used in PHENIX GUI
      test_flag_values = flag_scores.test_flag_values
      scored_arrays = []
      for i, array in enumerate(miller_arrays):
        scored_arrays.append( (array, test_flag_values[i], data_scores[i]) )
      def cmp_fn(x,y):
        return y[2] - x[2]
      scored_arrays.sort(key=cmp_to_key(cmp_fn))
      valid_arrays_and_flags = []
      for (array, flag_value, score) in scored_arrays :
        if score >= minimum_score :
          if array.is_string_array():
            array, flag_value = cif_status_flags_as_int_r_free_flags(
              array, flag_value)
          valid_arrays_and_flags.append((array, flag_value))
      return valid_arrays_and_flags
    if (label is None): labels = None
    else: labels=[label]
    try:
      i = select_array(
        parameter_name=parameter_scope,
        labels=labels,
        miller_arrays=miller_arrays,
        data_scores=data_scores,
        err=self.err,
        error_message_no_array
          ="No array of R-free flags found.\n\n"
          +"For manual selection define:\n"
          +"  %s.test_flag_value\n"%parameter_scope
          +"  %s.disable_suitability_test=True"%parameter_scope,
        error_message_not_a_suitable_array
          ="Not a suitable array of R-free flags: ",
        error_message_multiple_equally_suitable
          ="Multiple equally suitable arrays of R-free flags found.")
    except Sorry_Not_a_suitable_array as e:
      raise Sorry_Not_a_suitable_array(
        str(e) + "\nTo override the suitability test define:"
               + " %s.disable_suitability_test=True" % parameter_scope)
    miller_array = miller_arrays[i]
    if data_scores is not None:
      test_flag_value = flag_scores.test_flag_values[i]
    if miller_array.is_string_array():
      miller_array, test_flag_value = cif_status_flags_as_int_r_free_flags(
        miller_array, test_flag_value)
    return miller_array, test_flag_value

  def get_experimental_phases(self,
        file_name,
        labels,
        ignore_all_zeros,
        parameter_scope,
        raise_no_array=True,
        parameter_name="labels",
        return_all_valid_arrays=False,
        minimum_score=1):
    miller_arrays = self.get_miller_arrays(file_name=file_name)
    data_scores = get_experimental_phases_scores(
      miller_arrays=miller_arrays,
      ignore_all_zeros=ignore_all_zeros)
    if return_all_valid_arrays : # used in PHENIX GUI
      return sort_arrays_by_score(miller_arrays, data_scores, minimum_score)
    i = select_array(
      parameter_name=parameter_scope+"."+parameter_name,
      labels=labels,
      miller_arrays=miller_arrays,
      data_scores=data_scores,
      err=self.err,
      raise_no_array=raise_no_array,
      error_message_no_array
        ="No array of experimental phases found.",
      error_message_not_a_suitable_array
        ="Not a suitable array of experimental phases: ",
      error_message_multiple_equally_suitable
        ="Multiple equally suitable arrays of experimental phases found.")
    if i is None: return None
    return miller_arrays[i]

def cif_status_flags_as_int_r_free_flags(miller_array, test_flag_value):
  assert test_flag_value == 'f'
  if miller_array.is_string_array():
    selection = (miller_array.data() == 'o') | (miller_array.data() == 'f')
    info = miller_array.info()
    miller_array = miller_array.select(selection)
    data = flex.int(miller_array.size(), 1)
    data.set_selected(miller_array.data() == 'f', 0)
    miller_array = miller_array.array(data=data)
    miller_array.set_info(info)
    test_flag_value = 0
  return miller_array, test_flag_value

def guess_r_free_flag_value(miller_array, test_flag_value=None):
  flag_scores = get_r_free_flags_scores(
    miller_arrays=[miller_array],
    test_flag_value=test_flag_value)
  return flag_scores.test_flag_values[0]

def construct_output_file_name(input_file_names,
                               user_file_name,
                               file_type_label,
                               file_extension,
                               extension_seperator="."):
  if (user_file_name == "."):
    if (len(input_file_names) > 1):
      raise Sorry(
        "Ambiguous name for output %s file (more than one input file)."
          % file_type_label)
    user_file_name = os.path.basename(input_file_names[0])
  if (not user_file_name.lower().endswith(file_extension)):
    user_file_name += extension_seperator + file_extension
  if sys.platform == "win32":
    if (os.path.isfile(user_file_name)
        and user_file_name == input_file_names[0]):
      user_file_name += extension_seperator + file_extension
  else:
    if (os.path.isfile(user_file_name)
        and os.path.samefile(user_file_name, input_file_names[0])):
      user_file_name += extension_seperator + file_extension
  return user_file_name

def make_joined_set(miller_arrays):
  if(len(miller_arrays)==0): return None
  cs0 = miller_arrays[0].crystal_symmetry()
  for ma in miller_arrays:
    if([ma.crystal_symmetry().unit_cell(), cs0.unit_cell()].count(None)>0):
      return None
    if(not ma.crystal_symmetry().is_similar_symmetry(cs0)): return None
  from cctbx import miller
  master_set = miller.set(
    crystal_symmetry=miller_arrays[0].crystal_symmetry(),
    indices=miller_arrays[0].indices(),
    anomalous_flag=False)
  master_indices = miller_arrays[0].indices().deep_copy()
  for array in miller_arrays[1:] :
    current_indices = array.indices()
    missing_isel = miller.match_indices(master_indices,
      current_indices).singles(1)
    missing_indices = current_indices.select(missing_isel)
    master_indices.extend(missing_indices)
  master_set = miller.set(
    crystal_symmetry=miller_arrays[0].crystal_symmetry(),
    indices=master_indices,
    anomalous_flag=False)
  return \
    master_set.map_to_asu().unique_under_symmetry().remove_systematic_absences()

def extract_miller_array_from_file(file_name, label=None, type=None, log=None):
  if(log is None): log = sys.stdout
  assert type in ["complex", "real", None]
  result = None
  miller_arrays = reflection_file_reader.any_reflection_file(file_name =
    file_name).as_miller_arrays()
  def get_flag(ma):
    return (type == "complex" and ma.is_complex_array()) or \
           (type == "real"    and ma.is_real_array()) or \
           type is None
  print("  Available suitable arrays:", file=log)
  suitable_arrays = []
  suitable_labels = []
  for ma in miller_arrays:
    if(get_flag(ma=ma)):
      print("    ", ma.info().label_string(), file=log)
      suitable_arrays.append(ma)
      suitable_labels.append(ma.info().label_string())
  print(file=log)
  if(  len(suitable_arrays) == 0): raise Sorry("No suitable arrays.")
  elif(len(suitable_arrays) == 1): result = suitable_arrays[0]
  elif(len(suitable_arrays) >  1):
    if(label is None):
      msg='''Multiple choices available. No map coefficients array selected.

  See choices listed above and use "label=" to select one.
  Example: label="2FOFCWT,PH2FOFCWT"'''
      raise Sorry(msg)
    else:
      for ma in miller_arrays:
        if(get_flag(ma=ma) and (ma.info().label_string() == label)):
          print("  Selected:", ma.info().label_string(), file=log)
          result = ma
  return result

class process_raw_data(object):
  """
  Automation wrapper - prepares single-wavelength experimental data (and
  optional R-free flags and experimental phases) for any future step in the
  structure determination process.  Used in Phenix for ligand pipeline and
  automated re-refinement of PDB entries.
  """
  __slots__ = [
    "f_obs",
    "r_free_flags",
    "test_flag_value",
    "phases",
    "_generate_new",
  ]
  def __init__(self,
      obs,
      r_free_flags,
      test_flag_value,
      phases=None,
      d_min=None,
      d_max=None,
      r_free_flags_params=None,
      merge_anomalous=False,
      log=sys.stdout,
      verbose=True):
    assert (log is not None) and (obs is not None)
    if (r_free_flags_params is None):
      from cctbx.r_free_utils import generate_r_free_params_str
      r_free_flags_params = libtbx.phil.parse(
        generate_r_free_params_str).extract()
    obs_info = obs.info()
    r_free_flags_info = phases_info = None
    sg = obs.space_group_info()
    obs = obs.map_to_asu().merge_equivalents().array()
    obs = obs.eliminate_sys_absent(log=log)
    obs = obs.resolution_filter(d_min=d_min, d_max=d_max)
    if (obs.is_xray_intensity_array()):
      from cctbx import french_wilson
      if (verbose):
        fw_out = log
      else :
        fw_out = null_out()
      obs = french_wilson.french_wilson_scale(
        miller_array=obs,
        params=None,
        log=fw_out)
    assert (obs is not None)
    merged_obs = obs.average_bijvoet_mates()
    if (merged_obs.completeness() < 0.9):
      print("""
  WARNING: data are incomplete (%.1f%% of possible reflections measured to
  %.2fA).  This may cause problems if you plan to use the maps for building
  and/or ligand fitting!
    """ % (100*merged_obs.completeness(), merged_obs.d_min()), file=log)
    # XXX this is kind of a hack (the reconstructed arrays break some of my
    # assumptions about labels)
    if (merge_anomalous):
      obs = obs.average_bijvoet_mates()
    if (r_free_flags is not None):
      r_free_flags_info = r_free_flags.info()
      format = "cns"
      if (test_flag_value == 0):
        format = "ccp4"
      elif (test_flag_value == -1):
        format = "shelx"
      if (r_free_flags.anomalous_flag()):
        r_free_flags = r_free_flags.average_bijvoet_mates()
      is_compatible_symmetry = False
      obs_pg = obs.space_group().build_derived_point_group()
      flags_pg = r_free_flags.space_group().build_derived_point_group()
      if (obs_pg.type().number() == flags_pg.type().number()):
        is_compatible_symmetry = True
      else :
        pass # TODO unit cell comparison?
      if (is_compatible_symmetry):
        r_free_flags = r_free_flags.map_to_asu().merge_equivalents().array()
        r_free_flags = r_free_flags.eliminate_sys_absent(log=log)
        if (format == "cns"):
          r_free_flags = r_free_flags.customized_copy(
            crystal_symmetry=obs.crystal_symmetry(),
            data=(r_free_flags.data() == test_flag_value))
          test_flag_value = True
        obs_tmp = obs.deep_copy()
        if (obs.anomalous_flag()):
          obs_tmp = obs.average_bijvoet_mates()
        r_free_flags = r_free_flags.common_set(other=obs_tmp)
        n_r_free = r_free_flags.indices().size()
        n_obs = obs_tmp.indices().size()
        if ((test_flag_value is None) or
            (r_free_flags.data().all_eq(r_free_flags.data()[0]))):
          print("""
  WARNING: uniform R-free flags detected; a new test set will be generated,
  but this will bias the refinement statistics.
""", file=log)
          r_free_flags = None
        elif (n_r_free != n_obs):
          missing_set = obs_tmp.lone_set(other=r_free_flags)
          n_missing = missing_set.indices().size()
          if (n_missing > 0):
            print("""
  WARNING: R-free flags are incomplete relative to experimental
  data (%d vs. %d reflections).  The flags will be extended to
  complete the set, but we recommend supplying flags that are already
  generated to the maximum expected resolution.
""" % (n_r_free, n_obs), file=log)
            if (n_missing < 20) : # FIXME
              if (format == "cns"):
                missing_flags = missing_set.array(data=flex.bool(n_missing,
                  False))
              else :
                missing_flags = missing_set.array(data=flex.int(n_missing, 1))
            else :
              missing_flags = missing_set.generate_r_free_flags(
                fraction=(r_free_flags.data().count(test_flag_value)/n_r_free),
                max_free=None,
                use_lattice_symmetry=True,
                format=format)
            r_free_flags = r_free_flags.concatenate(other=missing_flags)
        if (r_free_flags is not None):
          assert (r_free_flags.indices().size() == obs_tmp.indices().size())
      else :
        print("""
    NOTE: incompatible symmetry between the data and the R-free flags:
         Data  : %s  %s
         Flags : %s  %s
       A new test set will be generated.
""" % (str(obs.space_group_info()),
          " ".join([ "%g" % x for x in obs.unit_cell().parameters() ]),
          str(r_free_flags.space_group_info()),
          " ".join(["%g" % x for x in r_free_flags.unit_cell().parameters()])), file=log)
    else :
      print("""
 WARNING: R-free flags not supplied.  This may bias the refinement if the
     structures are very nearly isomorphous!
""", file=log)
    self._generate_new = False
    if (r_free_flags is None):
      r_free_flags = obs.generate_r_free_flags(
        fraction=r_free_flags_params.fraction,
        max_free=r_free_flags_params.max_free,
        use_lattice_symmetry=r_free_flags_params.use_lattice_symmetry,
        use_dataman_shells=r_free_flags_params.use_dataman_shells,
        n_shells=r_free_flags_params.n_shells,
        format="ccp4")
      test_flag_value = 0
      self._generate_new = True
    if (r_free_flags.anomalous_flag()):
      r_free_flags = r_free_flags.average_bijvoet_mates()
    if (phases is not None):
      phases_info = phases.info()
      phases = phases.map_to_asu().resolution_filter(d_min=d_min, d_max=d_max)
    assert (obs.is_xray_amplitude_array())
    self.f_obs = obs.set_info(obs_info)
    self.r_free_flags = r_free_flags.set_info(r_free_flags_info)
    self.test_flag_value = test_flag_value
    self.phases = None
    if (phases is not None):
      self.phases = phases.set_info(phases_info)

  def data_labels(self):
    if (self.f_obs.is_xray_reconstructed_amplitude_array()):
      return "F,SIGF,DANO,SIGDANO,ISYM"
    elif (self.f_obs.anomalous_flag()):
      if (self.f_obs.sigmas() is not None):
        return "F(+),SIGF(+),F(-),SIGF(-)"
      else :
        return "F(+),F(-)"
    elif (self.f_obs.sigmas() is not None):
      return "F,SIGF"
    else :
      return "F"

  def r_free_flags_label(self):
    return "FreeR_flag"

  def r_free_flags_as_boolean_array(self):
    flags = self.r_free_flags.customized_copy(
      data=self.r_free_flags.data()==self.test_flag_value)
    if (self.f_obs.anomalous_flag()) and (not flags.anomalous_flag()):
      flags = flags.generate_bijvoet_mates()
    return flags

  def data_and_flags(self):
    return self.f_obs.common_sets(other=self.r_free_flags_as_boolean_array())

  def phase_labels(self):
    if (self.phases is not None):
      return "HLA,HLB,HLC,HLD"
    return None

  def n_obs(self):
    return self.f_obs.data().size()

  def fraction_free(self):
    return (self.r_free_flags.data().count(self.test_flag_value) /
            self.r_free_flags.data().size())

  def flags_are_new(self):
    return self._generate_new

  def write_mtz_file(self, file_name,
      title=None,
      wavelength=None,
      single_dataset=True):
    mtz_data = self.f_obs.as_mtz_dataset(
      column_root_label="F",
      wavelength=wavelength)
    if (self.f_obs.anomalous_flag()) and (not single_dataset):
      mtz_data.add_miller_array(
        miller_array=self.f_obs.average_bijvoet_mates(),
        column_root_label="F")
    if (self.phases is not None):
      mtz_data.add_miller_array(self.phases,
        column_root_label="HL")
    mtz_data.add_miller_array(self.r_free_flags,
      column_root_label="FreeR_flag")
    mtz_data.mtz_object().write(file_name)

def change_space_group(file_name, space_group_info):
  """
  Update the space group in an MTZ file, writing it in place.
  """
  mtz_in = reflection_file_reader.any_reflection_file(file_name)
  assert (mtz_in.file_type() == "ccp4_mtz")
  mtz_object = mtz_in.file_content()
  mtz_new = mtz_object.set_space_group_info(space_group_info)
  mtz_new.write(file_name)

def load_f_obs_and_r_free(file_name, anomalous_flag=False, phases=False):
  """
  Automation wrapper for reading in MTZ files generated by the process_raw_data
  class.
  """
  mtz_in = reflection_file_reader.any_reflection_file(file_name)
  assert (mtz_in.file_type() == "ccp4_mtz")
  file_server = reflection_file_server(
    crystal_symmetry=None,
    force_symmetry=True,
    reflection_files=[mtz_in],
    err=sys.stderr)
  f_obs = f_obs_anom = r_free = None
  r_free, test_flag_value = file_server.get_r_free_flags(
     file_name=file_name,
     label="FreeR_flag",
     test_flag_value=None,
     disable_suitability_test=False,
     parameter_scope="")
  r_free_info = r_free.info()
  assert (test_flag_value is not None)
  r_free = r_free.customized_copy(data=r_free.data()==test_flag_value)
  f_obs_info = f_obs_anom_info = None
  for array in file_server.miller_arrays :
    label = array.info().label_string()
    if (array.is_xray_amplitude_array()) and (array.anomalous_flag()):
      f_obs_anom = array
      f_obs_anom_info = f_obs_anom.info()
    elif (label == "F,SIGF"):
      f_obs = array
      f_obs_info = array.info()
  if (f_obs is None) and (f_obs_anom is not None):
    f_obs = f_obs_anom.average_bijvoet_mates()
  f_obs = f_obs.eliminate_sys_absent()
  # XXX this may still be necessary
  f_obs = f_obs.common_set(other=r_free)
  r_free = r_free.common_set(other=f_obs)
  assert (not None in [f_obs, r_free])
  if (f_obs_anom is not None) and (anomalous_flag):
    f_obs_anom = f_obs_anom.eliminate_sys_absent()
    r_free = r_free.generate_bijvoet_mates()
    f_obs = f_obs_anom.common_set(other=r_free).set_info(f_obs_anom_info)
    r_free = r_free.common_set(other=f_obs).set_info(r_free_info)
    return f_obs, r_free
  return f_obs.set_info(f_obs_info), r_free.set_info(r_free_info)
