from libtbx.itertbx import count
from libtbx.utils import UserError
import sys, os

def find_labels(search_labels, candidate_labels):
  for search_label in search_labels:
    for candidate_label in candidate_labels:
      if (candidate_label.find(search_label) < 0):
        return False
  return True

class label_table:

  def __init__(self, miller_arrays):
    self.miller_arrays = miller_arrays
    self.labels = []
    for p_array,miller_array in zip(count(1), miller_arrays):
      if (miller_array.info() is not None):
        self.labels.append(str(miller_array.info()))
      else:
        self.labels.append(str(p_array))

  def show_possible_choices(self, f=None, scores=None, high_score=None):
    if (f is None): f = sys.stdout
    print >> f, "Possible choices:"
    if (scores is None):
      for p_array,lbl in zip(count(1), self.labels):
        print >> f, "  %d:" % p_array, lbl
    else:
      for p_array,lbl,score in zip(count(1), self.labels, scores):
        if (score == high_score):
          print >> f, "  %d:" % p_array, lbl
    print >> f
    print >> f, "Please specify a number or an unambiguous substring of the",
    print >> f, "target data label."
    print >> f

  def match_data_label(self,
        label=None,
        labels=None,
        f=None,
        command_line_switch=None,
        parameter_name=None):
    assert label is None or labels is None
    assert label is not None or labels is not None
    if (f is None): f = sys.stdout
    assert command_line_switch is None or parameter_name is None
    if (labels is None):
      labels = [label]
    else:
      assert len(labels) > 0
      label = labels[0]
    if (command_line_switch is None):
      if (parameter_name is None):
        label_identifier = "--label"
      else:
        label_identifier = parameter_name
    else:
      label_identifier = parameter_name
    if (len(labels) == 1):
      try: i = int(label)-1
      except: pass
      else:
        if (0 <= i < len(self.miller_arrays)):
          return self.miller_arrays[i]
    scores = []
    labels_lower = [lbl.lower() for lbl in labels]
    for lbl in self.labels:
      if (not find_labels(
               search_labels=labels_lower, candidate_labels=[lbl.lower()])):
        scores.append(0)
      elif (not find_labels(search_labels=labels, candidate_labels=[lbl])):
        scores.append(1)
      else:
        scores.append(2)
    selected_array = None
    for high_score in [2,1]:
      if (scores.count(high_score) > 0):
        if (scores.count(high_score) > 1):
          print >> f
          print >> f, "Ambiguous %s=%s" % (label_identifier, " ".join(labels))
          print >> f
          self.show_possible_choices(f=f, scores=scores, high_score=high_score)
          return None
        return self.miller_arrays[scores.index(high_score)]
    print >> f
    print >> f, "Unknown %s=%s" % (label_identifier, " ".join(labels))
    print >> f
    self.show_possible_choices(f=f)
    return None

def construct_output_file_name(input_file_names,
                               user_file_name,
                               file_type_label,
                               file_extension,
                               extension_seperator="."):
  if (user_file_name == "."):
    if (len(input_file_names) > 1):
      raise UserError(
        "Ambiguous name for output %s file (more than one input file)."
          % file_type_label)
    user_file_name = os.path.basename(input_file_names[0])
  if (not user_file_name.lower().endswith(file_extension)):
    user_file_name += extension_seperator + file_extension
  if (    os.path.isfile(user_file_name)
      and os.path.samefile(user_file_name, input_file_names[0])):
    user_file_name += extension_seperator + file_extension
  return user_file_name
