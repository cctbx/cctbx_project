from iotbx import reflection_file_reader
from cctbx import miller
from cctbx.array_family import flex
import libtbx.path
from libtbx.str_utils import show_string
from libtbx.utils import Sorry
from itertools import count
import math
import sys, os

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
    print >> f, "Possible choices:"
    if (scores is None):
      for info_string in self.info_strings:
        print >> f, " ", info_string
    else:
      for info_string,score in zip(self.info_strings, scores):
        if (score >= minimum_score):
          print >> f, " ", info_string
    print >> f
    if (parameter_name is None): hint = ""
    else: hint = "use %s\nto " % parameter_name
    print >> f, \
      "Please %sspecify an unambiguous substring of the target label." % hint
    print >> f

  def match_data_label(self, label, command_line_switch, f=None):
    if (f is None): f = self.err
    assert label is not None
    scores = self.scores(label=label)
    selected_array = None
    for high_score in xrange(max(scores),0,-1):
      if (scores.count(high_score) > 0):
        if (scores.count(high_score) > 1):
          print >> f
          print >> f, "Ambiguous %s=%s" % (command_line_switch, label)
          print >> f
          self.show_possible_choices(
            f=f, scores=scores, minimum_score=high_score)
          return None
        return self.miller_arrays[scores.index(high_score)]
    print >> f
    print >> f, "Unknown %s=%s" % (command_line_switch, label)
    print >> f
    self.show_possible_choices(f=f)
    return None

  def select_array(self, label, command_line_switch, f=None):
    if (f is None): f = self.err
    if (len(self.miller_arrays) == 0):
      print >> f
      print >> f, "No reflection arrays available."
      print >> f
      return None
    if (len(self.miller_arrays) == 1):
      return self.miller_arrays[0]
    if (label is None):
      print >> f
      s = command_line_switch
      print >> f, "Please use %s to select a reflection array." % s
      print >> f, "For example: %s=%s" % (s, show_string(str(
        self.miller_arrays[1].info()).split(":")[-1]))
      print >> f
      self.show_possible_choices(f=f)
      return None
    return self.match_data_label(
      label=label,
      command_line_switch=command_line_switch)

def get_amplitude_scores(miller_arrays):
  result = []
  for miller_array in miller_arrays:
    score = 0
    if (miller_array.is_complex_array()):
      score = 1
    elif (miller_array.is_real_array()):
      if (miller_array.is_xray_amplitude_array()):
        score = 4
      elif (miller_array.is_xray_intensity_array()):
        score = 3
      elif (miller_array.is_xray_reconstructed_amplitude_array()):
        score = 2
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
      if (miller_array.is_xray_amplitude_array()):
        pass
      elif (miller_array.is_xray_intensity_array()):
        pass
      elif (miller_array.is_xray_reconstructed_amplitude_array()):
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

def get_xray_data_scores(miller_arrays, ignore_all_zeros):
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
        score = 8
      elif (miller_array.is_xray_amplitude_array()):
        if (miller_array.is_xray_reconstructed_amplitude_array()):
          score = 4
        else:
          score = 6
      else:
        score = 2
      assert score is not None
      if (    miller_array.sigmas() is not None
          and isinstance(miller_array.sigmas(), flex.double)
          and miller_array.sigmas_are_sensible()
          and score != 0):
        score += 1
      result.append(score)
  return result

def looks_like_r_free_flags_info(array_info):
  if (not isinstance(array_info, miller.array_info)): return False
  if (len(array_info.labels) > 2): return False
  label = array_info.labels[-1].lower()
  for word in ["free", "test", "cross", "status"]:
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
          c_keys = counts.keys()
          c_values = counts.values()
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
            elif (counts.size() > 3):
              if (c_keys == range(min(c_keys), max(c_keys)+1)):
                if (min(c_values) > max(c_values)*0.55):
                  if (looks_like_r_free_flags_info(miller_array.info())):
                    flag_score = 3
                  else:
                    flag_score = 2
                elif (looks_like_r_free_flags_info(miller_array.info())):
                  flag_score = 1
                if (flag_score != 0):
                  if (test_flag_value is None):
                    effective_test_flag_value = min(c_keys)
                  else:
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

def sort_arrays_by_score (miller_arrays, array_scores, minimum_score) :
  assert len(miller_arrays) == len(array_scores)
  i = 0
  scored_arrays = []
  for array in miller_arrays :
    scored_arrays.append( (array, array_scores[i]) )
    i += 1
  scored_arrays.sort(lambda x, y: y[1] - x[1])
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
      error_message_multiple_equally_suitable):
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
    print >> err, "\n" + error + "\n"
    if (max(data_scores) > 0):
      lbl_tab.show_possible_choices(
        scores=data_scores,
        minimum_score=1,
        parameter_name=parameter_name)
    raise Sorry(error)
  if (max(data_scores) == 0):
    if (label_scores is None):
      print >> err, "\n" + error_message_no_array + "\n"
      raise Sorry_No_array_of_the_required_type(error_message_no_array)
    error = "%s%s=%s" % (
      error_message_not_a_suitable_array, parameter_name, " ".join(labels))
    print >> err, "\n" + error + "\n"
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
    print >> err, "\n" + error + "\n"
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
    self.file_name_miller_arrays = {}
    for miller_array in self.miller_arrays:
      self.file_name_miller_arrays.setdefault(
        libtbx.path.canonical_path(
          miller_array.info().source), []).append(miller_array)

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

  def get_amplitudes(self,
        file_name,
        labels,
        convert_to_amplitudes_if_necessary,
        parameter_scope,
        parameter_name,
        return_all_valid_arrays=False,
        minimum_score=1):
    miller_arrays = self.get_miller_arrays(file_name=file_name)
    data_scores = get_amplitude_scores(miller_arrays=miller_arrays)
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
        minimum_score=1) :
    miller_arrays = self.get_miller_arrays(file_name=file_name)
    data_scores = get_xray_data_scores(
      miller_arrays=miller_arrays,
      ignore_all_zeros=ignore_all_zeros)
    if return_all_valid_arrays :
      return sort_arrays_by_score(miller_arrays, data_scores, minimum_score)
    # Recognize phenix.refine file and do the "right thing". May be too ad hoc..
    new_miller_arrays = []
    new_data_scores = []
    for ma, ds in zip(miller_arrays, data_scores):
      if(not (ma.info().labels == ['F-obs-filtered', 'SIGF-obs-filtered'] or
         ma.info().labels == ['F-model', 'PHIF-model'] or
         isinstance(ma.data(), flex.complex_double))):
        new_miller_arrays.append(ma)
        new_data_scores.append(ds)
    #
    i = select_array(
      parameter_name=parameter_scope+"."+parameter_name,
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
        minimum_score=1) :
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
      for i, array in enumerate(miller_arrays) :
        scored_arrays.append( (array, test_flag_values[i], data_scores[i]) )
      scored_arrays.sort(lambda x, y: y[2] - x[2])
      valid_arrays_and_flags = []
      for (array, flag_value, score) in scored_arrays :
        if score >= minimum_score :
          valid_arrays_and_flags.append((array, flag_value))
      return valid_arrays_and_flags
    if (label is None): labels = None
    else: labels=[label]
    try:
      i = select_array(
        parameter_name=parameter_scope+".label",
        labels=labels,
        miller_arrays=miller_arrays,
        data_scores=data_scores,
        err=self.err,
        error_message_no_array
          ="No array of R-free flags found.\n\n"
          +"For manual selection define:\n"
          +"  %s.label\n"%parameter_scope
          +"  %s.test_flag_value\n"%parameter_scope
          +"  %s.disable_suitability_test=True"%parameter_scope,
        error_message_not_a_suitable_array
          ="Not a suitable array of R-free flags: ",
        error_message_multiple_equally_suitable
          ="Multiple equally suitable arrays of R-free flags found.")
    except Sorry_Not_a_suitable_array, e:
      raise Sorry_Not_a_suitable_array(
        str(e) + "\nTo override the suitability test define:"
               + " %s.disable_suitability_test=True" % parameter_scope)
    if (data_scores is None):
      return miller_arrays[i], test_flag_value
    return miller_arrays[i], flag_scores.test_flag_values[i]

  def get_experimental_phases(self,
        file_name,
        labels,
        ignore_all_zeros,
        parameter_scope,
        parameter_name="labels",
        return_all_valid_arrays=False,
        minimum_score=1) :
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
      error_message_no_array
        ="No array of experimental phases found.",
      error_message_not_a_suitable_array
        ="Not a suitable array of experimental phases: ",
      error_message_multiple_equally_suitable
        ="Multiple equally suitable arrays of experimental phases found.")
    return miller_arrays[i]

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
  if (    os.path.isfile(user_file_name)
      and os.path.samefile(user_file_name, input_file_names[0])):
    user_file_name += extension_seperator + file_extension
  return user_file_name
