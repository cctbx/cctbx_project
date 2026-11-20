"""
Utility for matching a dials-format non-image pickle file to an experiment list and producing a cctbx-format dictionary
from the pair.
"""
from __future__ import absolute_import, division, print_function

import os
from dials.util.options import Importer, flatten_reflections, flatten_experiments
from serialtbx.util.construct_frame import ConstructFrame

def find_json(pickle, pickle_ext=None, json_ext=None):
  """find the matching json file for a given dials-format non-image pickle file"""
  name = os.path.basename(pickle).split(".refl")[0]
  dirname = os.path.dirname(pickle)
  kind = None # the reflection table's suffix, e.g. _integrated
  if pickle_ext is not None:
    if pickle_ext == "":
      base = name
    else:
      base = name.split(pickle_ext)[0]
      kind = pickle_ext
  elif name.endswith("_indexed"):
    base = name.split("_indexed")[0]
    kind = "_indexed"
  elif name.endswith("_integrated"):
    base = name.split("_integrated")[0]
    kind = "_integrated"
  else:
    base = name
  if json_ext is not None:
    return os.path.join(dirname, base + json_ext + ".expt")
  # dials.stills_process writes each reflection table next to an experiment
  # list with the same suffix. Older versions paired _integrated.refl and
  # _indexed.refl with _refined.expt, so that is tried next.
  candidates = []
  if kind:
    candidates.append(kind + ".expt")
  candidates.extend(["_refined.expt", "_indexed.expt"])
  for ext in candidates:
    json = os.path.join(dirname, base + ext)
    if os.path.exists(json):
      return json
  return None

class read_dials_pickle(object):
  """given a dials-format non-image pickle file with matching json file, return a dictionary pickle"""
  def __init__(self, path, json=None, pickle_ext=None, json_ext=None):
    if json is None:
      json = find_json(path, pickle_ext, json_ext)
    if json is None:
      importer = Importer([path], read_experiments=False, read_reflections=True, check_format=False)
      print("unable to find experiment list")
      self.experiments = None
    else:
      importer = Importer([path, json], read_experiments=True, read_reflections=True, check_format=False)
      try:
        self.experiments = flatten_experiments(importer.experiments)[0]
      except IndexError:
        print("unable to read experiment list")
        self.experiments = None
    try:
      self.reflections = flatten_reflections(importer.reflections)[0]
    except IndexError:
      print("unable to read reflection table")
      self.reflections = None
  def make_pickle(self):
    if (self.experiments is None) or (self.reflections is None):
      self.dictionary = None
    else:
      contents = ConstructFrame(self.reflections, self.experiments)
      self.dictionary = contents.make_frame()
