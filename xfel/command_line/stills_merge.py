# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# LIBTBX_SET_DISPATCHER_NAME stills.merge
#
# $Id: stills_merge.py 20545 2014-08-25 22:22:15Z idyoung $

from __future__ import division
from __future__ import print_function
from six.moves import range

import iotbx.phil
from cctbx.array_family import flex
from cctbx.sgtbx.bravais_types import bravais_lattice
from libtbx.utils import Sorry
from dials.array_family import flex
from xfel.command_line import cxi_merge
from xfel.command_line import frame_extractor
from xfel.command_line.cxi_merge import OutlierCellError, WrongBravaisError
import os
import math
import sys
op = os.path

def eval_ending (file_name):
  ordered_endings_mapping = [
    ("refined_experiments.json", "integrated.pickle"),
    ("experiments.json", "integrated.pickle"),
    ]
  dir_name = os.path.dirname(file_name)
  basename = os.path.basename(file_name)
  for pair in ordered_endings_mapping:
    if basename.endswith(pair[0]):
      refl_name = os.path.join(dir_name, basename.split(pair[0])[0] + pair[1])
      fragment = basename.split(pair[0])[0]
      digit = "0" # if fragment contains no digit
      for i in range(len(fragment)):
        if fragment[-i-1].isdigit():
          digit = fragment[-i-1]
          break
      return (digit, file_name, refl_name)
  return None

def get_observations (data_dirs,data_subset):
  print("Step 1.  Get a list of all files")
  file_names = []
  for dir_name in data_dirs :
    if not os.path.isdir(dir_name):
      continue
    for file_name in os.listdir(dir_name):
      if eval_ending(file_name) is not None: # only accept jsons
        if data_subset == 0 or \
          (data_subset == 1 and int(eval_ending(file_name)[0][-1]) % 2 == 1 or \
          (data_subset == 2 and int(eval_ending(file_name)[0][-1]) % 2 == 0)):
          file_names.append(os.path.join(dir_name, file_name))
  print("Number of frames found:", len(file_names))
  return file_names

cxi_merge.get_observations = get_observations

def load_result (file_name,
                 ref_bravais_type,
                 reference_cell,
                 params,
                 reindex_op,
                 out) :
  # Pull relevant information from integrated.pickle and refined_experiments.json
  # files to construct the equivalent of a single integration pickle (frame).
  try:
    frame = frame_extractor.ConstructFrameFromFiles(eval_ending(file_name)[2], eval_ending(file_name)[1]).make_frame()
  except Exception:
    return None

  # If @p file_name cannot be read, the load_result() function returns
  # @c None.

  print("Step 2.  Load frame obj and filter on lattice & cell with",reindex_op)
  """
  Take a frame with all expected contents of an integration pickle, confirm
  that it contains the appropriate data, and check the lattice type and unit
  cell against the reference settings - if rejected, raises an exception
  (for tracking statistics).
  """
  # Ignore frames with no integrated reflections.
  obj = frame
  if ("observations" not in obj) :
    return None

  if reindex_op == "h,k,l":
    pass
  else:
    obj["observations"][0].apply_change_of_basis(reindex_op)
    pass

  result_array = obj["observations"][0]
  unit_cell = result_array.unit_cell()
  sg_info = result_array.space_group_info()
  print("", file=out)
  print("-" * 80, file=out)
  print(file_name, file=out)
  print(sg_info, file=out)
  print(unit_cell, file=out)

  #Check for pixel size (at this point we are assuming we have square pixels, all experiments described in one
  #refined_experiments.json file use the same detector, and all panels on the detector have the same pixel size)

  if params.pixel_size is not None:
    pixel_size = params.pixel_size
  elif "pixel_size" in obj:
    pixel_size = obj["pixel_size"]
  else:
    raise Sorry("Cannot find pixel size. Specify appropriate pixel size in mm for your detector in phil file.")

  #Calculate displacements based on pixel size
  assert obj['mapped_predictions'][0].size() == obj["observations"][0].size()
  mm_predictions = pixel_size*(obj['mapped_predictions'][0])
  mm_displacements = flex.vec3_double()
  cos_two_polar_angle = flex.double()
  for pred in mm_predictions:
    mm_displacements.append((pred[0]-obj["xbeam"],pred[1]-obj["ybeam"],0.0))
    cos_two_polar_angle.append( math.cos( 2. * math.atan2(pred[1]-obj["ybeam"],pred[0]-obj["xbeam"]) ) )
  obj["cos_two_polar_angle"] = cos_two_polar_angle
  #then convert to polar angle and compute polarization correction

  if (not bravais_lattice(sg_info.type().number()) == ref_bravais_type) :
    raise WrongBravaisError("Skipping cell in different Bravais type (%s)" %
      str(sg_info))
  if (not unit_cell.is_similar_to(
      other=reference_cell,
      relative_length_tolerance=params.unit_cell_length_tolerance,
      absolute_angle_tolerance=params.unit_cell_angle_tolerance)) :
    raise OutlierCellError(
      "Skipping cell with outlier dimensions (%g %g %g %g %g %g" %
      unit_cell.parameters())
  print("Integrated data:", file=out)
  result_array.show_summary(f=out, prefix="  ")
  # XXX don't force reference setting here, it will be done later, after the
  # original unit cell is recorded
  return obj

cxi_merge.load_result = load_result


if (__name__ == "__main__"):
  show_plots = False
  if ("--plots" in sys.argv) :
    sys.argv.remove("--plots")
    show_plots = True
  result = cxi_merge.run(args=sys.argv[1:])
  if result is None:
    sys.exit(1)
  if (show_plots) :
    try :
      result.plots.show_all_pyplot()
      from wxtbx.command_line import loggraph
      loggraph.run([result.loggraph_file])
    except Exception as e :
      print("Can't display plots")
      print("You should be able to view them by running this command:")
      print("  wxtbx.loggraph %s" % result.loggraph_file)
      raise e
