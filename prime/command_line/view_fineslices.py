# LIBTBX_SET_DISPATCHER_NAME prime.view_fineslices
"""
Author      : Uervirojnangkoorn, M.
Desc        : Read pickle files (hopefully fine-sliced, and named after the
              order of the slices. Grab selected reflections and plot the
              intensity.
"""
from __future__ import absolute_import, division, print_function
import sys
from six.moves import cPickle as pickle
from dials.array_family import flex
import matplotlib.pyplot as plt
import numpy as np
from prime.postrefine.mod_input import read_frame, read_pickles
from six.moves import range
from six.moves import zip

def read_input(args):
  data = []
  data_sweep = ''
  for i in range(len(args)):
    pair=args[i].split('=')
    if len(pair) == 2:
      if pair[0]=='data':
        data.append(pair[1])
      elif pair[0]=='data_sweep':
        data_sweep = pair[1]
  if len(data) == 0:
    print("Please input all parameters")
    exit()
  return data, data_sweep

if (__name__ == "__main__"):
  #Read input parameters and frames (pickle files)
  data_fine, data_sweep = read_input(args = sys.argv[1:])
  pixel_size_mm = 0.079346
  #read all the pickle files from fine-slice data
  frames_fine = read_pickles(data_fine)
  #get a sample reflections
  sample_no = 0
  obs_fine_sample = None
  for i in range(2000):
    frame = frames_fine[i]
    pickle_fine = read_frame(frame)
    obs_fine = pickle_fine["observations"][0]
    obs_fine = obs_fine.select(obs_fine.data() > 0)
    if len(obs_fine.data())> 5:
      print(frame)
      for index, d, I, sigI in zip(obs_fine.indices(), obs_fine.d_spacings().data(),\
                                   obs_fine.data(), obs_fine.sigmas()):
        print(index, d, I, sigI)
      obs_fine_sample = obs_fine.deep_copy()
      #break
    sample_no += 1
  """
  #find matching indices in the next consecutive frames
  n_frames_sel = 10
  colors = np.random.rand(n_frames_sel)
  i = 0
  for frame in frames_fine[sample_no:sample_no+n_frames_sel]:
    pickle_fine = read_frame(frame)
    obs_fine = pickle_fine["observations"][0]
    mm_predictions = pixel_size_mm * (pickle_fine['mapped_predictions'][0])
    xbeam = pickle_fine["xbeam"]
    ybeam = pickle_fine["ybeam"]
    spot_pred_x_mm = flex.double([pred[0]-xbeam for pred in mm_predictions])
    spot_pred_y_mm = flex.double([pred[1]-ybeam for pred in mm_predictions])
    #plt.subplot(2,5,i+1)
    plt.scatter(spot_pred_x_mm, spot_pred_y_mm, s=10, c=colors)
    #plt.title(frame)
    plt.grid(True)
    obs_fine = obs_fine.select(obs_fine.data() > 0)
    obs_fine_sample_match, obs_fine_match = obs_fine_sample.common_sets(obs_fine)
    print frame
    for index, d, I, sigI in zip(obs_fine_match.indices(), obs_fine_match.d_spacings().data(),\
                                 obs_fine_match.data(), obs_fine_match.sigmas()):
      print index, d, I, sigI
    i += 1
  plt.show()
  """
