from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME prime.print_integration_pickle
"""
Author      : Uervirojnangkoorn, M.
Created     : 11/1/2015
Description : read integration pickles and view systemetic absences and beam X, Y position
"""

import cPickle as pickle
from cctbx.array_family import flex
import sys, os, math
from scitbx.matrix import sqr
import numpy as np
import matplotlib.pyplot as plt
from cctbx.uctbx import unit_cell
from prime.postrefine.mod_leastsqr import good_unit_cell
from prime.postrefine.mod_mx import mx_handler
from cctbx import statistics
import iotbx.phil
from libtbx.utils import Usage, Sorry
from prime.postrefine.mod_input import read_pickles

master_phil = iotbx.phil.parse("""
data = None
  .type = path
  .multiple = True
  .help = Directory containing integrated data in pickle format.  Repeat to \
    specify additional directories.
hklisoin = None
  .type = path
  .help = Mtz file for the calculation of CCiso
pixel_size_mm = None
  .type = float
  .help = Pixel size in mm. (MAR = 0.079346)
  .optional = False
target_unit_cell = None
  .type = unit_cell
  .help = Target unit-cell parameters are used to discard outlier cells.
target_space_group = None
  .type = str
  .help = Target space group.
target_anomalous_flag = False
  .type = bool
  .help = Target anomalous flag (False = Not anomalous data)
flag_plot = False
  .type = bool
  .help = Normal plots.
d_min = 0.1
  .type = float
  .help = Minimum resolution.
d_max = 99
  .type = float
  .help = Maximum resolution.
n_residues = None
  .type = int
  .help = No. of amino acid residues.
frame_accept_min_cc = 0.25
  .type = float
  .help = CC cut-off for the rejection of frames for reporting only.
beam_accept_offset_mm = 0.5
  .type = float
  .help = Beam center cut-off in mm. (calculated from mean values).
uc_tolerance = 5
  .type = float
  .help = Unit-cell tolerance in percent.
""")
txt_help = """**************************************************************************************************

prime.print_integration_pickle

For viewing systematic absences and beam xy position from integration pickles.

Usage: prime.print_integration_pickle parameter.phil

With this command, you can specify all parameters required by prime in your parameter.phil file.
To obtain the template of these parameters, you can perform a dry run (simply run prime.print_integration_pickle).
You can then change the values of the parameters.

For feedback, please contact monarin@stanford.edu.

**************************************************************************************************

List of available parameters:
"""

def calc_wilson(observations_full, n_residues):
  """
  Caculate isotropic Wilson G and B-factors
  """
  G,B  = (0,0)
  if n_residues:
    mxh = mx_handler()
    asu_contents = mxh.get_asu_contents(n_residues)
    try:
      observations_as_f = observations_full.as_amplitude_array()
      binner_template_asu = observations_as_f.setup_binner(auto_binning=True)
      wp = statistics.wilson_plot(observations_as_f, asu_contents, e_statistics=True)
      G = wp.wilson_intensity_scale_factor
      B = wp.wilson_b
    except Exception, e:
      print "Warning (Wilson Plot Calculation)", e
  return G, B

if (__name__ == "__main__"):
  #setup phil parameters
  user_phil = []
  if sys.argv[1:] == None:
    master_phil.show()
    raise Usage("Use the above list of parameters to generate your input file (.phil). For more information, run prime.print_integration_pickle -h.")
  else:
    for arg in sys.argv[1:]:
      if os.path.isfile(arg):
        user_phil.append(iotbx.phil.parse(open(arg).read()))
      elif (os.path.isdir(arg)) :
        user_phil.append(iotbx.phil.parse("""data=\"%s\"""" % arg))
      else :
        if arg == '--help' or arg == '-h':
          print txt_help
          master_phil.show(attributes_level=1)
          raise Usage("Run prime.print_integration_pickle to generate a list of initial parameters.")
        else:
          try:
            user_phil.append(iotbx.phil.parse(arg))
          except RuntimeError, e :
            raise Sorry("Unrecognized argument '%s' (error: %s)" % (arg, str(e)))
  working_phil = master_phil.fetch(sources=user_phil)
  params = working_phil.extract()
  if not params.data:
    raise Usage("Error: Data is required. Please specify path to your data folder (data=/path/to/integration/results).")
  if params.target_unit_cell:
    if not params.target_space_group: raise Usage("Error: Please specify target space group (target_space_group=xx).")
  if not params.pixel_size_mm:
    raise Usage("Error: Please specify pixel size in mm. (pixel_size_mm=xx).")
  #0 .read input parameters and frames (pickle files)
  frame_files = read_pickles(params.data)
  xbeam_set,ybeam_set,sys_abs_all,cc_bin_low_set = (flex.double(),flex.double(),flex.double(),flex.double())
  uc_a, uc_b, uc_c, dd_mm, wavelength_set = (flex.double(),flex.double(),flex.double(),flex.double(),flex.double())
  sys_abs_set,cc_bins_set,d_bins_set,oodsqr_bins_set,flag_good_unit_cell_set = ([],[],[],[],[])
  print 'Summary of integration pickles:'
  print '(image file, min. res., max. res, beamx, beamy, n_refl, cciso, <cciso_bin>, a, b, c, mosaicity, residual, dd_mm, wavelength, good_cell?, G, B)'
  txt_out_sys_abs = '\nSummary of systematic absences:\n'
  txt_out_sys_abs += 'Image File, Resolution, Original Index, Asym. Index, I, sigI, I/sigI\n'
  for pickle_filename in frame_files:
    observations_pickle = pickle.load(open(pickle_filename,"rb"))
    pickle_filename_only = os.path.basename(pickle_filename)
    mxh = mx_handler()
    flag_hklisoin_found, miller_array_iso = mxh.get_miller_array_from_reflection_file(params.hklisoin)
    observations = observations_pickle["observations"][0]
    swap_dict = {'test_xx1.pickle':0, \
      'tes_xx2.pickle':0}
    if pickle_filename_only in swap_dict:
      from cctbx import sgtbx
      cb_op = sgtbx.change_of_basis_op('a,c,b')
      observations = observations.change_basis(cb_op)
    if params.target_unit_cell:
      try:
        from cctbx.crystal import symmetry
        miller_set = symmetry(
            unit_cell=observations.unit_cell().parameters(),
            space_group_symbol=params.target_space_group
          ).build_miller_set(
            anomalous_flag=params.target_anomalous_flag,
            d_min=observations.d_min())
        observations = observations.customized_copy(anomalous_flag=params.target_anomalous_flag,
                          crystal_symmetry=miller_set.crystal_symmetry())
      except Exception:
        print 'Cannot apply target space group: observed space group=', observations.space_group_info()
    #check if the uc is good
    if params.target_unit_cell is None:
      flag_good_unit_cell = True
    else:
      flag_good_unit_cell = good_unit_cell(observations.unit_cell().parameters(), None, params.uc_tolerance, target_unit_cell=params.target_unit_cell)
    flag_good_unit_cell_set.append(flag_good_unit_cell)
    #calculate partiality
    wavelength = observations_pickle["wavelength"]
    crystal_init_orientation = observations_pickle["current_orientation"][0]
    detector_distance_mm = observations_pickle['distance']
    mm_predictions = params.pixel_size_mm*(observations_pickle['mapped_predictions'][0])
    xbeam = observations_pickle["xbeam"]
    ybeam = observations_pickle["ybeam"]
    xbeam_set.append(xbeam)
    ybeam_set.append(ybeam)
    alpha_angle = flex.double([math.atan(abs(pred[0]-xbeam)/abs(pred[1]-ybeam)) \
                                   for pred in mm_predictions])
    spot_pred_x_mm = flex.double([pred[0]-xbeam for pred in mm_predictions])
    spot_pred_y_mm = flex.double([pred[1]-ybeam for pred in mm_predictions])
    #calculate mean unit cell
    if flag_good_unit_cell:
      uc_a.append(observations.unit_cell().parameters()[0])
      uc_b.append(observations.unit_cell().parameters()[1])
      uc_c.append(observations.unit_cell().parameters()[2])
      dd_mm.append(detector_distance_mm)
      wavelength_set.append(wavelength)
    #resoultion filter
    i_sel_res = observations.resolution_filter_selection(d_min=params.d_min, d_max=params.d_max)
    observations = observations.select(i_sel_res)
    alpha_angle = alpha_angle.select(i_sel_res)
    spot_pred_x_mm = spot_pred_x_mm.select(i_sel_res)
    spot_pred_y_mm = spot_pred_y_mm.select(i_sel_res)
    #sort by resolution
    perm = observations.sort_permutation(by_value="resolution", reverse=True)
    observations = observations.select(perm)
    alpha_angle = alpha_angle.select(perm)
    spot_pred_x_mm = spot_pred_x_mm.select(perm)
    spot_pred_y_mm = spot_pred_y_mm.select(perm)
    from prime.postrefine.mod_partiality import partiality_handler
    ph = partiality_handler()
    r0 = ph.calc_spot_radius(sqr(crystal_init_orientation.reciprocal_matrix()),
                                          observations.indices(), wavelength)
    two_theta = observations.two_theta(wavelength=wavelength).data()
    sin_theta_over_lambda_sq = observations.two_theta(wavelength=wavelength).sin_theta_over_lambda_sq().data()
    ry, rz, re, nu, rotx, roty = (0, 0, 0.003,0.5, 0, 0)
    flag_beam_divergence = False
    partiality_init, delta_xy_init, rs_init, rh_init = ph.calc_partiality_anisotropy_set(crystal_init_orientation.unit_cell(),
                                                          rotx, roty, observations.indices(),
                                                          ry, rz, r0, re, nu, two_theta, alpha_angle, wavelength,
                                                          crystal_init_orientation, spot_pred_x_mm, spot_pred_y_mm,
                                                          detector_distance_mm, "Lorentzian",
                                                          flag_beam_divergence)
    I_full = observations.data()/ partiality_init
    sigI_full = observations.sigmas()/ partiality_init
    observations_full = observations.customized_copy(data=I_full, sigmas=sigI_full)
    wilson_G, wilson_B = calc_wilson(observations_full, params.n_residues)
    #calculate R and cc with reference
    cc_iso, cc_full_iso, cc_bin_low, cc_bin_med = (0, 0, 0, 0)
    observations_asu = observations.map_to_asu()
    observations_full_asu = observations_full.map_to_asu()
    cc_bins = flex.double()
    oodsqr_bins = flex.double()
    d_bins = flex.double()
    if flag_hklisoin_found:
      #build observations dict
      obs_dict = {}
      for mi_asu, I, sigI, I_full, sigI_full in zip(observations_asu.indices(), \
        observations_asu.data(), observations_asu.sigmas(), \
        observations_full_asu.data(), observations_full_asu.sigmas()):
        obs_dict[mi_asu]=(I, sigI, I_full, sigI_full)
      I_match = flex.double()
      I_full_match = flex.double()
      I_iso_match = flex.double()
      d_match = flex.double()
      oodsqr_match = flex.double()
      for mi_asu, d, I_iso in zip(miller_array_iso.indices(), miller_array_iso.d_spacings().data(), miller_array_iso.data()):
        if mi_asu in obs_dict:
          I, sigI, I_full, sigI_full = obs_dict[mi_asu]
          I_match.append(I)
          I_full_match.append(I_full)
          I_iso_match.append(I_iso)
          oodsqr_match.append(1/(d**2))
          d_match.append(d)
      #calculate correlation
      try:
        cc_iso = np.corrcoef(I_iso_match, I_match)[0,1]
        cc_full_iso = np.corrcoef(I_iso_match, I_full_match)[0,1]
      except Exception:
        pass
      #scatter plot
      if params.flag_plot and len(frame_files)==1:
        plt.subplot(211)
        #plt.scatter(oodsqr_match, I_iso_match, s=10, marker='o', c='r')
        plt.plot(oodsqr_match, I_iso_match)
        plt.title('Reference intensity')
        plt.xlabel('1/d^2')
        plt.ylabel('I_ref')
        plt.subplot(212)
        #plt.scatter(oodsqr_match, I_match, s=10, marker='o', c='r')
        plt.plot(oodsqr_match, I_match)
        plt.title('Observed intensity CC=%.4g'%(cc_iso))
        plt.xlabel('1/d^2')
        plt.ylabel('I_obs')
        plt.show()
      #scatter bin plot
      n_bins = 10
      n_refl = int(round(len(I_match)/n_bins))
      if len(I_match) > 0:
        for i_bin in range(n_bins):
          i_start = i_bin*n_refl
          if i_bin < n_bins - 1:
            i_end = (i_bin*n_refl) + n_refl
          else:
            i_end = -1
          I_iso_bin = I_iso_match[i_start:i_end]
          I_bin = I_match[i_start:i_end]
          d_bin = d_match[i_start:i_end]
          try:
            cc_bin = np.corrcoef(I_iso_bin, I_bin)[0,1]
          except Exception:
            cc_bin = 0
          cc_bins.append(cc_bin)
          try:
            min_d_bin = np.min(d_bin)
          except Exception:
            min_d_bin = 1
          d_bins.append(min_d_bin)
          oodsqr_bins.append(1/(min_d_bin**2))
          if i_bin == 0:
            cc_bin_low = cc_bin
          if i_bin == 5:
            cc_bin_med = cc_bin
          if params.flag_plot and len(frame_files)==1:
            plt.subplot(2,5,i_bin+1)
            plt.scatter(I_iso_bin, I_bin,s=10, marker='o', c='r')
            plt.title('Bin %2.0f (%6.2f-%6.2f A) CC=%6.2f'%(i_bin+1, np.max(d_bin), np.min(d_bin), cc_bin))
            if i_bin == 0:
              plt.xlabel('I_ref')
              plt.ylabel('I_obs')
        if params.flag_plot and len(frame_files)==1:
          plt.show()
      #print full detail if given a single file
      if len(frame_files) == 1:
        print 'Crystal orientation'
        print crystal_init_orientation.crystal_rotation_matrix()
        print 'Direct matrix'
        print crystal_init_orientation.direct_matrix()
    a, b, c, alpha, beta, gamma = observations.unit_cell().parameters()
    txt_out_head= '{0:40} {1:5.2f} {2:5.2f} {3:5.2f} {4:5.2f} {5:5.0f} {6:6.2f} {7:6.2f} {8:6.2f} {9:6.2f} {10:6.2f} {11:6.2f} {12:6.2f} {13:6.2f} {14:6.2f} {15:6.2f} {16:6.2f} {20:6.4f} {17:5} {18:6.2f} {19:6.2f}'.format(pickle_filename_only, observations.d_min(), np.max(observations.d_spacings().data()), xbeam, ybeam, len(observations.data()), cc_iso, np.mean(cc_bins), a, b, c, alpha, beta, gamma, observations_pickle["mosaicity"], observations_pickle["residual"], detector_distance_mm, flag_good_unit_cell, wilson_G, wilson_B, wavelength)
    print txt_out_head
    cc_bin_low_set.append(cc_iso)
    cc_bins_set.append(cc_bins)
    d_bins_set.append(d_bins)
    oodsqr_bins_set.append(oodsqr_bins)
    sys_abs_lst = flex.double()
    if params.target_unit_cell:
      cn_refl = 0
      for sys_absent_flag, d_spacing, miller_index_ori, miller_index_asu, I, sigI in zip(observations.sys_absent_flags(), observations.d_spacings().data(), observations.indices(), observations_asu.indices(), observations.data(), observations.sigmas()):
        if sys_absent_flag[1]:
          txt_out_sys_abs += '{9:40} {10:6.2f} {0:3},{1:3},{2:3} {3:3},{4:3},{5:3} {6:8.2f} {7:8.2f} {8:6.2f}\n'.format(miller_index_ori[0], miller_index_ori[1], miller_index_ori[2], miller_index_asu[0], miller_index_asu[1], miller_index_asu[2], I, sigI, I/sigI, pickle_filename_only, d_spacing)
          cn_refl +=1
          sys_abs_lst.append(I/sigI)
          sys_abs_all.append(I/sigI)
    sys_abs_set.append(sys_abs_lst)
  #report systematic absences
  print txt_out_sys_abs
  #collect beamxy
  xbeam_mean = flex.mean(xbeam_set)
  xbeam_std = np.std(xbeam_set)
  ybeam_mean = flex.mean(ybeam_set)
  ybeam_std = np.std(ybeam_set)
  xbeam_filtered_set = flex.double()
  ybeam_filtered_set = flex.double()
  frame_filtered_set = []
  sys_abs_all_filtered = flex.double()
  txt_out = ''
  txt_out_mix = ''
  txt_out_uc = ''
  txt_out_report_beam_filter = 'Images with beam center displaced > %6.2f mm.:\n'%(params.beam_accept_offset_mm)
  txt_out_report_cc_filter = 'Images with cc < %6.2f:\n'%(params.frame_accept_min_cc)
  from scitbx.matrix import col
  for pickle_filename, xbeam, ybeam, sys_abs_lst, cc_bin_low, flag_good_unit_cell in \
    zip(frame_files, xbeam_set, \
    ybeam_set, sys_abs_set, cc_bin_low_set, flag_good_unit_cell_set):
    pickle_filename_arr = pickle_filename.split('/')
    pickle_filename_only = pickle_filename_arr[len(pickle_filename_arr)-1]
    pred_xy = col((xbeam, ybeam))
    calc_xy = col((xbeam_mean, ybeam_mean))
    diff_xy = pred_xy - calc_xy
    txt_out_report_tmp = '{0:80} {1:6.2f} {2:6.2f} {3:6.2f} {4:6.4f}\n'.format(pickle_filename_only, xbeam, ybeam, cc_bin_low, diff_xy.length())
    if diff_xy.length() < params.beam_accept_offset_mm:
      xbeam_filtered_set.append(xbeam)
      ybeam_filtered_set.append(ybeam)
      frame_filtered_set.append(pickle_filename)
      txt_out += pickle_filename + '\n'
      sys_abs_all_filtered.extend(sys_abs_lst)
      if cc_bin_low > params.frame_accept_min_cc and flag_good_unit_cell:
        txt_out_mix += pickle_filename + '\n'
      else:
        txt_out_report_cc_filter += txt_out_report_tmp
      if flag_good_unit_cell:
        txt_out_uc += pickle_filename + '\n'
    else:
      txt_out_report_beam_filter += txt_out_report_tmp
  print
  print 'CC mean=%6.2f median=%6.2f std=%6.2f'%(flex.mean(cc_bin_low_set), np.median(cc_bin_low_set), np.std(cc_bin_low_set))
  print 'Xbeam mean=%8.4f std=%6.4f'%(xbeam_mean, xbeam_std)
  print 'Ybeam mean=%8.4f std=%6.4f'%(ybeam_mean, ybeam_std)
  if uc_a:
    print 'UC mean a=%8.4f (%8.4f) b=%8.4f (%8.4f) c=%9.4f (%8.4f)'%(flex.mean(uc_a), np.std(uc_a), flex.mean(uc_b), np.std(uc_b), flex.mean(uc_c), np.std(uc_c))
  print 'Detector distance mean=%8.4f (%8.4f)'%(flex.mean(dd_mm), np.std(dd_mm))
  print 'Wavelength mean=%8.4f (%8.4f)'%(flex.mean(wavelength_set), np.std(wavelength_set))
  print 'No. of frames: All = %6.0f Beam outliers = %6.0f CC filter=%6.0f'%(len(frame_files), len(frame_files) - (len(txt_out.split('\n'))-1), len(frame_files) - (len(txt_out_mix.split('\n'))-1))
  print
  print 'Reporting outliers (image name, xbeam, ybeam, cciso, delta_xy)'
  print txt_out_report_beam_filter
  if params.hklisoin: print txt_out_report_cc_filter
  #write out filtered beamxy pickle files
  f = open('integration_pickle_beam_filter.lst', 'w')
  f.write(txt_out)
  f.close()
  #write out mix filter pickle files
  f = open('integration_mix_filter.lst', 'w')
  f.write(txt_out_mix)
  f.close()
  #write out filtered beamxy and uc pickle files
  f = open('integration_pickle_beam_uc_filter.lst', 'w')
  f.write(txt_out_uc)
  f.close()
  if params.flag_plot:
    plt.subplot(211)
    plt.plot(xbeam_set,ybeam_set,"r.")
    plt.xlim([xbeam_mean-2, xbeam_mean+2])
    plt.ylim([ybeam_mean-2, ybeam_mean+2])
    #plt.axes().set_aspect("equal")
    plt.title('Raw data')
    plt.grid(True)
    plt.subplot(212)
    plt.plot(xbeam_filtered_set,ybeam_filtered_set,"r.")
    plt.xlim([xbeam_mean-2, xbeam_mean+2])
    plt.ylim([ybeam_mean-2, ybeam_mean+2])
    #plt.axes().set_aspect("equal")
    plt.title('After beamxy position filter')
    plt.grid(True)
    plt.show()
    #plot I/sigI histogram for systematic absences
    if len(sys_abs_set) > 0:
      plt.subplot(211)
      sys_abs_sel = sys_abs_all.select(flex.abs(sys_abs_all)>2)
      x = sys_abs_sel.as_numpy_array()
      mu = np.mean(x)
      med = np.median(x)
      sigma = np.std(x)
      num_bins = 25
      n, bins, patches = plt.hist(x, num_bins, normed=False, facecolor='blue', alpha=0.5)
      plt.ylabel('Frequencies')
      plt.grid()
      plt.title('Systematic absences with |I/sigI| > 2.0\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))
      plt.subplot(212)
      sys_abs_all_filtered_sel = sys_abs_all_filtered.select(flex.abs(sys_abs_all_filtered)>2)
      x = sys_abs_all_filtered_sel.as_numpy_array()
      mu = np.mean(x)
      med = np.median(x)
      sigma = np.std(x)
      num_bins = 25
      n, bins, patches = plt.hist(x, num_bins, normed=False, facecolor='blue', alpha=0.5)
      plt.ylabel('Frequencies')
      plt.grid()
      plt.title('Systematic absences with |I/sigI| > 2.0 (After BeamXY filter)\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))
      plt.show()
    cn_i = 0
    for cc_bins, d_bins, oodsqrt_bins, cc_bin_low, pickle_filename in zip(cc_bins_set, d_bins_set, oodsqr_bins_set, cc_bin_low_set, frame_files):
      pickle_filename_arr = pickle_filename.split('/')
      pickle_filename_only = pickle_filename_arr[len(pickle_filename_arr)-1]
      if cc_bin_low < params.frame_accept_min_cc:
        plt.subplot(2,1,1)
        plt.plot(oodsqrt_bins, cc_bins)
        plt.title('CC by resolutions for cc < %6.2f'%(params.frame_accept_min_cc))
        plt.xlabel('1/d^2')
        plt.ylim([0,1])
        plt.grid(True)
      else:
        plt.subplot(2,1,2)
        plt.plot(oodsqrt_bins, cc_bins)
        plt.title('CC by resolutions for cc > %6.2f'%(params.frame_accept_min_cc))
        plt.xlabel('1/d^2')
        plt.ylim([0,1])
        plt.grid(True)
        cn_i += 1
    plt.show()
  #print twin operators
  from prime.index_ambiguity.mod_indexing_ambiguity import indamb_handler
  idah = indamb_handler()
  print "\nReporting all images with alternative indexing choices:"
  for pickle_filename in frame_files:
    observations_pickle = pickle.load(open(pickle_filename,"rb"))
    observations = observations_pickle["observations"][0]
    operators = idah.generate_twin_operators(observations, flag_all=True)
    ops_hkl = [op.operator.r().as_hkl() for op in operators]
    print os.path.basename(pickle_filename), '%6.1f %6.1f %6.1f %6.1f %6.1f %6.1f'%observations.unit_cell().parameters(), ops_hkl
