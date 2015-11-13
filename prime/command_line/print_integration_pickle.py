from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME prime.print_integration_pickle
"""
Author      : Uervirojnangkoorn, M.
Created     : 11/1/2015
Description : read integration pickles and view systemetic absences and beam X, Y position
"""

import cPickle as pickle
from cctbx.array_family import flex
from iotbx import reflection_file_reader
import sys, os, math
from scitbx.matrix import sqr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

def get_miller_array_from_mtz(mtz_filename):
  flag_hklisoin_found = False
  miller_array_iso = None
  if mtz_filename is not None:
    flag_hklisoin_found = True
    reflection_file_iso = reflection_file_reader.any_reflection_file(mtz_filename)
    miller_arrays_iso=reflection_file_iso.as_miller_arrays()
    is_found_iso_as_intensity_array = False
    is_found_iso_as_amplitude_array = False
    for miller_array in miller_arrays_iso:
      if miller_array.is_xray_intensity_array():
        miller_array_iso = miller_array.deep_copy()
        is_found_iso_as_intensity_array = True
        break
      elif miller_array.is_xray_amplitude_array():
        is_found_iso_as_amplitude_array = True
        miller_array_converted_to_intensity = miller_array.as_intensity_array()
    if is_found_iso_as_intensity_array == False:
      if is_found_iso_as_amplitude_array:
        miller_array_iso = miller_array_converted_to_intensity.deep_copy()
      else:
        flag_hklisoin_found = False
  return flag_hklisoin_found, miller_array_iso


def read_pickles(data):
  frame_files = []
  for p in data:
    if os.path.isdir(p) == False:
      #check if list-of-pickle text file is given
      pickle_list_file = open(p,'r')
      pickle_list = pickle_list_file.read().split("\n")
      for pickle_filename in pickle_list:
        if os.path.isfile(pickle_filename):
          frame_files.append(pickle_filename)
    else:
      for pickle_filename in os.listdir(p):
        if pickle_filename.endswith('.pickle'):
          frame_files.append(p+'/'+pickle_filename)

  #check if pickle_dir is given in input file instead of from cmd arguments.
  if len(frame_files)==0:
    print 'No pickle files found.'
    exit()

  return frame_files

def read_input(args):
  if len(args) == 0:
    print "prime.print_integration_pickle: for viewing systematic absences and beam xy position from integration pickles."
    print "usage: prime.print_integration_pickle data=integrated.lst pixel_size_mm=0.079346 check_sys_absent=True target_space_group=P212121"
    exit()
  data = []
  hklrefin = None
  pixel_size_mm = None
  check_sys_absent = False
  target_space_group = None
  target_anomalous_flag = False
  for i in range(len(args)):
    pair=args[i].split('=')
    if pair[0]=='data':
      data.append(pair[1])
    if pair[0]=='hklrefin':
      hklrefin = pair[1]
    if pair[0]=='pixel_size_mm':
      pixel_size_mm = float(pair[1])
    if pair[0]=='check_sys_absent':
      check_sys_absent = bool(pair[1])
    if pair[0]=='target_space_group':
      target_space_group = pair[1]
    if pair[0]=='target_anomalous_flag':
      target_anomalous_flag = bool(pair[1])

  if len(data)==0:
    print "Please provide data path. (eg. data=/path/to/pickle/)"
    exit()
  if check_sys_absent:
    if target_space_group is None:
      print "Please provide target space group if you want to check systematic absence (eg. target_space_group=P212121)"
      exit()
  if pixel_size_mm is None:
    print "Please specify pixel size (eg. pixel_size_mm=0.079346)"
    exit()

  return data, hklrefin, pixel_size_mm, check_sys_absent, target_space_group, target_anomalous_flag


if (__name__ == "__main__"):

  #0 .read input parameters and frames (pickle files)
  data, hklrefin, pixel_size_mm, check_sys_absent_input, target_space_group, target_anomalous_flag = read_input(args = sys.argv[1:])
  frame_files = read_pickles(data)

  xbeam_set = flex.double()
  ybeam_set = flex.double()
  sys_abs_set = flex.double()
  for pickle_filename in frame_files:
    check_sys_absent = check_sys_absent_input
    observations_pickle = pickle.load(open(pickle_filename,"rb"))
    pickle_filename_arr = pickle_filename.split('/')
    pickle_filename_only = pickle_filename_arr[len(pickle_filename_arr)-1]

    flag_hklisoin_found, miller_array_iso = get_miller_array_from_mtz(hklrefin)

    observations = observations_pickle["observations"][0]

    #apply constrain using the crystal system
    if check_sys_absent:
      try:
        from cctbx.crystal import symmetry
        miller_set = symmetry(
            unit_cell=observations.unit_cell().parameters(),
            space_group_symbol=target_space_group
          ).build_miller_set(
            anomalous_flag=target_anomalous_flag,
            d_min=observations.d_min())

        observations = observations.customized_copy(anomalous_flag=target_anomalous_flag,
                          crystal_symmetry=miller_set.crystal_symmetry())
      except Exception:
        print 'Cannot apply target space group: observed space group=', observations.space_group_info()
        check_sys_absent = False

    #calculate partiality
    wavelength = observations_pickle["wavelength"]
    crystal_init_orientation = observations_pickle["current_orientation"][0]
    detector_distance_mm = observations_pickle['distance']
    mm_predictions = pixel_size_mm*(observations_pickle['mapped_predictions'][0])
    xbeam = observations_pickle["xbeam"]
    ybeam = observations_pickle["ybeam"]
    xbeam_set.append(xbeam)
    ybeam_set.append(ybeam)
    alpha_angle = flex.double([math.atan(abs(pred[0]-xbeam)/abs(pred[1]-ybeam)) \
                                   for pred in mm_predictions])
    spot_pred_x_mm = flex.double([pred[0]-xbeam for pred in mm_predictions])
    spot_pred_y_mm = flex.double([pred[1]-ybeam for pred in mm_predictions])

    from prime.postrefine.mod_leastsqr import calc_spot_radius
    r0 = calc_spot_radius(sqr(crystal_init_orientation.reciprocal_matrix()),
                                          observations.indices(), wavelength)

    from prime.postrefine.mod_leastsqr import calc_partiality_anisotropy_set
    two_theta = observations.two_theta(wavelength=wavelength).data()
    sin_theta_over_lambda_sq = observations.two_theta(wavelength=wavelength).sin_theta_over_lambda_sq().data()
    ry, rz, re, rotx, roty = (0, 0, 0.002, 0, 0)
    flag_beam_divergence = False

    partiality_init, delta_xy_init, rs_init, rh_init = calc_partiality_anisotropy_set(crystal_init_orientation.unit_cell(),
                                                          rotx, roty, observations.indices(),
                                                          ry, rz, r0, re, two_theta, alpha_angle, wavelength,
                                                          crystal_init_orientation, spot_pred_x_mm, spot_pred_y_mm,
                                                          detector_distance_mm, "Lorentzian",
                                                          flag_beam_divergence)
    I_full = observations.data()/ partiality_init
    sigI_full = observations.sigmas()/ partiality_init
    observations_full = observations.customized_copy(data=I_full, sigmas=sigI_full)

    #calculate R and cc with reference
    cc_iso, cc_full_iso = (0, 0)
    observations_asu = observations.map_to_asu()
    if flag_hklisoin_found:
      corr = observations_asu.correlation(miller_array_iso, use_binning=False, assert_is_similar_symmetry=False)
      cc_iso = corr.coefficient()

      observations_full_asu = observations_full.map_to_asu()
      corr = observations_full_asu.correlation(miller_array_iso, use_binning=False, assert_is_similar_symmetry=False)
      cc_full_iso = corr.coefficient()

    a, b, c, alpha, beta, gamma = observations.unit_cell().parameters()
    txt_out_head= '{0:25} {1:5.2f} {2:6.0f} {3:6.2f} {4:6.2f} {5:6.2f} {6:8.6f} {7:6.2f} {8:6.2f} {9:6.2f} {10:6.2f} {11:6.2f} {12:6.2f} {13:6.2f} {17:6.2f} {14:6.4f} {15:6.4f} {16:6.2f} {18:8.2f} {19:15}'.format(pickle_filename_only, observations.d_min(), len(observations.data()), detector_distance_mm, xbeam, ybeam, wavelength, a, b, c, alpha, beta, gamma, cc_iso, observations_pickle["residual"], observations_pickle["ML_half_mosaicity_deg"][0], observations_pickle["mosaicity"], cc_full_iso, observations_pickle["ML_domain_size_ang"][0], observations.space_group_info())


    if check_sys_absent:
      cn_refl = 0
      for sys_absent_flag, miller_index_ori, miller_index_asu, I, sigI in zip(observations.sys_absent_flags(), observations.indices(), observations_asu.indices(), observations.data(), observations.sigmas()):
        if sys_absent_flag[1]:
          if True:
            txt_out = txt_out_head + '{0:3} {1:3} {2:3} {3:3} {4:3} {5:3} {6:8.2f} {7:8.2f} {8:6.2f}'.format(miller_index_ori[0], miller_index_ori[1], miller_index_ori[2], miller_index_asu[0], miller_index_asu[1], miller_index_asu[2], I, sigI, I/sigI)
            print txt_out
            cn_refl +=1
            sys_abs_set.append(I/sigI)

  #plot beamxy
  xbeam_mean = flex.mean(xbeam_set)
  xbeam_std = np.std(xbeam_set)
  ybeam_mean = flex.mean(ybeam_set)
  ybeam_std = np.std(ybeam_set)

  beam_thres = 3.0
  xbeam_filtered_set = flex.double()
  ybeam_filtered_set = flex.double()
  frame_filtered_set = []
  txt_out = ''
  for pickle_filename, xbeam, ybeam in zip(frame_files, xbeam_set, ybeam_set):
    if abs(xbeam - xbeam_mean)/xbeam_std < beam_thres and \
      abs(ybeam - ybeam_mean)/ybeam_std < beam_thres:
        xbeam_filtered_set.append(xbeam)
        ybeam_filtered_set.append(ybeam)
        frame_filtered_set.append(pickle_filename)
        txt_out += pickle_filename + '\n'
  print 'Xbeam mean=%8.4f std=%6.4f'%(xbeam_mean, xbeam_std)
  print 'Ybeam mean=%8.4f std=%6.4f'%(ybeam_mean, ybeam_std)
  print 'N_frames = %6.0f After filter = %6.0f'%(len(frame_files), len(frame_filtered_set))

  #write out filtered beamxy pickle files
  f = open('integration_pickle_beam_filtered.lst', 'w')
  f.write(txt_out)
  f.close()

  plt.plot(xbeam_set,ybeam_set,"r.")
  plt.axes().set_aspect("equal")
  plt.title('Raw data')
  plt.grid(True)
  plt.show()

  plt.plot(xbeam_filtered_set,ybeam_filtered_set,"r.")
  plt.axes().set_aspect("equal")
  plt.title('After beamxy position filter')
  plt.grid(True)
  plt.show()

  #plot I/sigI histogram for systematic absences
  if len(sys_abs_set) > 0:
    x = sys_abs_set.as_numpy_array()
    mu = np.mean(x)
    med = np.median(x)
    sigma = np.std(x)
    num_bins = 20
    n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
    y = mlab.normpdf(bins, mu, sigma)
    plt.plot(bins, y, 'r--')
    plt.ylabel('Frequencies')
    plt.title('I/sigI distribution of systematic absences\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))
    plt.show()
