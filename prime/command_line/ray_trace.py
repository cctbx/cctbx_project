# LIBTBX_SET_DISPATCHER_NAME prime.ray_trace
"""
Author      : Uervirojnangkoorn, M.
Created     : 11/1/2015
Description : read integration pickles and view systemetic absences and beam X, Y position
"""
from __future__ import absolute_import, division, print_function

from six.moves import cPickle as pickle
from cctbx.array_family import flex
from cctbx.uctbx import unit_cell
from scitbx.matrix import sqr
import sys, math
import numpy as np
import matplotlib.pyplot as plt
from prime.postrefine.mod_leastsqr import good_unit_cell
from prime.postrefine.mod_partiality import partiality_handler
from prime.postrefine.mod_mx import mx_handler
from prime.postrefine.mod_leastsqr import good_unit_cell
from prime.postrefine.mod_input import read_frame, read_pickles
from six.moves import range

def read_input(args):
  if len(args) == 0:
    print("prime.ray_trace: read integration pickles and ray trace the reflections")
    print("usage: prime.ray_trace data=integrated.lst")
    exit()
  data = []
  hklrefin = None
  pixel_size_mm = None
  target_unit_cell = None
  d_min = 0
  d_max = 99
  for i in range(len(args)):
    pair=args[i].split('=')
    if pair[0]=='data':
      data.append(pair[1])
    if pair[0]=='hklrefin':
      hklrefin = pair[1]
    if pair[0]=='pixel_size_mm':
      pixel_size_mm = float(pair[1])
    if pair[0]=='target_unit_cell':
      try:
        tuc = pair[1].split(',')
        a,b,c,alpha,beta,gamma = (float(tuc[0]), float(tuc[1]), float(tuc[2]), \
          float(tuc[3]), float(tuc[4]), float(tuc[5]))
        target_unit_cell = unit_cell((a,b,c,alpha,beta,gamma))
      except Exception:
        pass
    if pair[0]=='d_min':
      d_min = float(pair[1])
    if pair[0]=='d_max':
      d_max = float(pair[1])
  if len(data)==0:
    print("Please provide data path. (eg. data=/path/to/pickle/)")
    exit()
  if pixel_size_mm is None:
    print("Please specify pixel size (eg. pixel_size_mm=0.079346)")
    exit()
  return data, hklrefin, pixel_size_mm, target_unit_cell, d_min, d_max


if (__name__ == "__main__"):
  uc_tol = 3
  ry, rz, re, rotx, roty = (0, 0, 0.008, 0, 0)
  flag_beam_divergence = False
  lambda_template = flex.double(range(-50,50,1))/1000
  #0 .read input parameters and frames (pickle files)
  data, hklrefin, pixel_size_mm, target_unit_cell, \
    d_min, d_max = read_input(args = sys.argv[1:])
  frame_files = read_pickles(data)
  for pickle_filename in frame_files:
    observations_pickle = read_frame(pickle_filename)
    pickle_filename_arr = pickle_filename.split('/')
    pickle_filename_only = pickle_filename_arr[len(pickle_filename_arr)-1]
    mxh = mx_handler()
    flag_hklisoin_found, miller_array_iso = mxh.get_miller_array_from_reflection_file(hklrefin)
    observations = observations_pickle["observations"][0]
    #check if the uc is good
    flag_good_unit_cell = good_unit_cell(observations.unit_cell().parameters(), None, uc_tol, target_unit_cell=target_unit_cell)
    #update lambda_set
    lambda_set = lambda_template + observations_pickle["wavelength"]
    crystal_init_orientation = observations_pickle["current_orientation"][0]
    detector_distance_mm = observations_pickle['distance']
    mm_predictions = pixel_size_mm*(observations_pickle['mapped_predictions'][0])
    xbeam = observations_pickle["xbeam"]
    ybeam = observations_pickle["ybeam"]
    alpha_angle = flex.double([math.atan(abs(pred[0]-xbeam)/abs(pred[1]-ybeam)) \
                                   for pred in mm_predictions])
    spot_pred_x_mm = flex.double([pred[0]-xbeam for pred in mm_predictions])
    spot_pred_y_mm = flex.double([pred[1]-ybeam for pred in mm_predictions])
    #resoultion filter
    i_sel_res = observations.resolution_filter_selection(d_min=d_min, d_max=d_max)
    observations = observations.select(i_sel_res)
    alpha_angle = alpha_angle.select(i_sel_res)
    spot_pred_x_mm = spot_pred_x_mm.select(i_sel_res)
    spot_pred_y_mm = spot_pred_y_mm.select(i_sel_res)
    #collect offset spots (left and right)
    off_left_x_mm = flex.double()
    off_left_y_mm = flex.double()
    off_right_x_mm = flex.double()
    off_right_y_mm = flex.double()
    off_mid_x_mm = flex.double()
    off_mid_y_mm = flex.double()
    if flag_good_unit_cell and len(spot_pred_x_mm) > 0:
      #calculate rh set for lambda_set
      rh_array = np.zeros([len(observations.data()), len(lambda_set)])
      p_array = np.zeros([len(observations.data()), len(lambda_set)])
      ph = partiality_handler()
      r0 = ph.calc_spot_radius(sqr(crystal_init_orientation.reciprocal_matrix()),
                                observations.indices(), observations_pickle["wavelength"])
      for lambda_now, i_lambda in zip(lambda_set, range(len(lambda_set))):
        two_theta = observations.two_theta(wavelength=lambda_now).data()
        partiality, delta_xy, rs_set, rh_set = ph.calc_partiality_anisotropy_set(\
                                                crystal_init_orientation.unit_cell(),
                                                rotx, roty, observations.indices(),
                                                ry, rz, r0, re, 0, two_theta, alpha_angle,
                                                lambda_now, crystal_init_orientation,
                                                spot_pred_x_mm, spot_pred_y_mm,
                                                detector_distance_mm, "Lorentzian",
                                                flag_beam_divergence)
        rh_array[:, i_lambda] = list(flex.abs(rh_set))
        p_array[:, i_lambda] = list(partiality)
      #find minimum distance for all reflections
      assigned_lambda = np.argmin(rh_array, axis=1)
      #assigned_lambda = np.argmax(p_array, axis=1)
      for i in range(len(observations.data())):
        if assigned_lambda[i] == 0:
          off_left_x_mm.append(spot_pred_x_mm[i])
          off_left_y_mm.append(spot_pred_y_mm[i])
        elif assigned_lambda[i] == len(lambda_template)-1:
          off_right_x_mm.append(spot_pred_x_mm[i])
          off_right_y_mm.append(spot_pred_y_mm[i])
        else:
          off_mid_x_mm.append(spot_pred_x_mm[i])
          off_mid_y_mm.append(spot_pred_y_mm[i])
        if assigned_lambda[i] == 0 or assigned_lambda[i] == len(lambda_template)-1:
          dummy = 0
          # uncomment below to see best wavelength for each spot
          """
          print p_array[i,:], assigned_lambda[i]
          plt.plot(lambda_set, p_array[i,:])
          plt.title('Best wavelength=%8.6f (i=%2.0f) res=%6.2f'%(lambda_set[assigned_lambda[i]], assigned_lambda[i], observations.d_spacings().data()[i]))
          plt.show()
          """
      # plot which part of the spectrum this reflection is close to
      plt.scatter(off_mid_x_mm, off_mid_y_mm, s=10, c='g', marker='o')
      plt.scatter(off_left_x_mm, off_left_y_mm, s=20, c='b', marker='<')
      plt.scatter(off_right_x_mm, off_right_y_mm, s=20, c='r', marker='>')
      plt.title("Ray Trace Result o=right on, <=high energy, > low energy");
      plt.xlim([-100,100])
      plt.ylim([-100,100])
      plt.grid(True)
      plt.show()
      n, bins, patches = plt.hist(assigned_lambda, 30, normed=0, facecolor='b', alpha=0.5)
      plt.show()
