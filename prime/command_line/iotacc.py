# LIBTBX_SET_DISPATCHER_NAME prime.iotacc
'''
Author      : Uervirojnangkoorn, M.
Created     : 11/25/2014
Description : iotacc selects iota integration results base on CC with ref. set.
'''
from __future__ import absolute_import, division, print_function
import os, sys
import numpy as np
import math
from iotbx import reflection_file_reader
from cctbx.array_family import flex
from cctbx import miller
from six.moves import cPickle as pickle
from cctbx.crystal import symmetry
from scitbx.matrix import sqr
import shutil
from libtbx.easy_mp import pool_map
from prime.postrefine.mod_input import read_frame

def read_input(args):
  from prime.postrefine.mod_input import process_input
  iparams, txt_out_input = process_input(args)
  return iparams, txt_out_input

def read_pickles(data):
  frame_files = []
  for file_name in os.listdir(data):
    if file_name.endswith('.pickle'):
      frame_files.append(data+'/'+file_name)

  return frame_files

def get_Eoc_corrected_observations(pickle_filename, iparams):
  '''
  Read pickle file name and output Eoc corrected original miller array object.
  '''
  observations_pickle = read_frame(pickle_filename)
  observations_original = observations_pickle["observations"][0]
  wavelength = observations_pickle["wavelength"]
  crystal_init_orientation = observations_pickle["current_orientation"][0]
  uc_params = observations_original.unit_cell().parameters()
  detector_distance_mm = observations_pickle['distance']
  mm_predictions = iparams.pixel_size_mm*(observations_pickle['mapped_predictions'][0])
  xbeam = observations_pickle["xbeam"]
  ybeam = observations_pickle["ybeam"]
  alpha_angle = flex.double([math.atan(abs(pred[0]-xbeam)/abs(pred[1]-ybeam)) \
                                   for pred in mm_predictions])
  spot_pred_x_mm = flex.double([pred[0]-xbeam for pred in mm_predictions])
  spot_pred_y_mm = flex.double([pred[1]-ybeam for pred in mm_predictions])

  #filter resolution
  i_sel_res = observations_original.resolution_filter_selection(\
        d_max=iparams.iotacc.d_max, d_min=iparams.iotacc.d_min)
  observations_original = observations_original.customized_copy(\
        indices=observations_original.indices().select(i_sel_res), \
        data=observations_original.data().select(i_sel_res), \
        sigmas=observations_original.sigmas().select(i_sel_res))
  alpha_angle = alpha_angle.select(i_sel_res)
  spot_pred_x_mm = spot_pred_x_mm.select(i_sel_res)
  spot_pred_y_mm = spot_pred_y_mm.select(i_sel_res)


  #Filter weak
  i_sel = (observations_original.data()/observations_original.sigmas()) > iparams.iotacc.sigma_min
  observations_original = observations_original.customized_copy(\
        indices=observations_original.indices().select(i_sel),\
        data=observations_original.data().select(i_sel),\
        sigmas=observations_original.sigmas().select(i_sel))
  alpha_angle = alpha_angle.select(i_sel)
  spot_pred_x_mm = spot_pred_x_mm.select(i_sel)
  spot_pred_y_mm = spot_pred_y_mm.select(i_sel)

  #calculate spot radius
  from prime.postrefine.mod_leastsqr import calc_spot_radius
  spot_radius = calc_spot_radius(sqr(crystal_init_orientation.reciprocal_matrix()),
                                          observations_original.indices(), wavelength)

  G = 1
  B = 0
  from prime.postrefine.mod_leastsqr import calc_partiality_anisotropy_set
  two_theta = observations_original.two_theta(wavelength=wavelength).data()
  sin_theta_over_lambda_sq = observations_original.two_theta(wavelength=wavelength).sin_theta_over_lambda_sq().data()
  ry = 0
  rz = 0
  re = iparams.gamma_e
  rotx = 0
  roty = 0

  #calc partiality
  partiality_init, delta_xy_init, rs_init, rh_init = calc_partiality_anisotropy_set(crystal_init_orientation.unit_cell(), rotx, roty, observations_original.indices(), ry, rz, spot_radius, re, two_theta, alpha_angle, wavelength, crystal_init_orientation, spot_pred_x_mm, spot_pred_y_mm, detector_distance_mm, iparams.partiality_model, iparams.flag_beam_divergence)

  from prime.postrefine.mod_leastsqr import calc_full_refl
  I_full = calc_full_refl(observations_original.data(), sin_theta_over_lambda_sq,
                            G, B, partiality_init, rs_init, iparams.flag_volume_correction)
  sigI_full = calc_full_refl(observations_original.sigmas(), sin_theta_over_lambda_sq,
                            G, B, partiality_init, rs_init, iparams.flag_volume_correction)
  observations_Eoc_corrected = observations_original.customized_copy(data=I_full, sigmas=sigI_full)

  return observations_Eoc_corrected, observations_original

def get_observations_ref(hklrefin):
  '''
  Read and return ref. observations.
  '''
  reflection_file_iso = reflection_file_reader.any_reflection_file(hklrefin)
  miller_arrays_iso=reflection_file_iso.as_miller_arrays()
  is_found_iso_as_intensity_array = False
  is_found_iso_as_amplitude_array = False
  miller_array_iso = None
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

  return miller_array_iso

def calc_cc(miller_array_ref, miller_array_obs):
  '''
  Calculate cc between matched intensities.
  '''
  matches_ref = miller.match_multi_indices(
                miller_indices_unique=miller_array_ref.indices(),
                miller_indices=miller_array_obs.indices())

  I_ref = flex.double([miller_array_ref.data()[pair[0]] for pair in matches_ref.pairs()])
  I_o = flex.double([miller_array_obs.data()[pair[1]] for pair in matches_ref.pairs()])

  cc = 0
  n_refl = 0
  if len(matches_ref.pairs()) > 0 :
    cc = np.corrcoef(I_o, I_ref)[0,1]
    if math.isnan(cc):
      cc = 0
    n_refl = len(matches_ref.pairs())

  return cc, n_refl

def good_unit_cell(uc_params, target_unit_cell, uc_tol):
  '''
  check unit cell
  '''
  flag_good_uc = False
  if (abs(uc_params[0]-target_unit_cell[0]) \
      <= (uc_tol*target_unit_cell[0]/100) \
                and abs(uc_params[1]-target_unit_cell[1]) \
                <= (uc_tol*target_unit_cell[1]/100) \
                and abs(uc_params[2]-target_unit_cell[2]) \
                <= (uc_tol*target_unit_cell[2]/100) \
                and abs(uc_params[3]-target_unit_cell[3]) \
                <= (uc_tol*target_unit_cell[3]/100) \
                and abs(uc_params[4]-target_unit_cell[4]) \
                <= (uc_tol*target_unit_cell[4]/100) \
                and abs(uc_params[5]-target_unit_cell[5]) \
                <= (uc_tol*target_unit_cell[5]/100)):
    flag_good_uc = True
  return flag_good_uc

def select_best_by_cc_mproc(shot_no, cryst_id, iparams):
  if iparams.iotacc.set_id is not None:
    data_dir = iparams.data[0] + '/' +iparams.iotacc.set_id
  else:
    data_dir = iparams.data[0]

  cryst_id_shot_no = cryst_id+'_0_'+str(shot_no).zfill(iparams.iotacc.LEN_SHOT_SEQ)
  #cryst_id_shot_no = cryst_id+str(shot_no).zfill(LEN_SHOT_SEQ)
  cc_list = flex.double()
  n_refl_list = flex.int()
  pickle_filename_list = []
  flag_success = False
  d_min, d_max = (iparams.iotacc.d_min, iparams.iotacc.d_max)
  txt_out = 'ID'.center(3)+'shot'.center(6)+'res.'.center(7)+'n_refl'.center(7)+ \
    'cc_pt'.center(7)+'cc_fl'.center(7)+'n_ref'.center(7)+'cc_pt'.center(7)+'cc_fl'.center(7)+'n_iso'.center(7) + \
    'a'.center(7)+'b'.center(7)+'c'.center(7)+'alp'.center(7)+'gam'.center(7)+'bet'.center(7)+'file name'.center(20)+'\n'
  for tmp_shot_dir in os.listdir(data_dir+'/'+cryst_id):
    if tmp_shot_dir.find(cryst_id_shot_no) > 0:
      data = data_dir+'/'+cryst_id+'/'+tmp_shot_dir
      #data = data_dir+'/'+tmp_shot_dir
      frame_files = []
      if os.path.isdir(data):
        frame_files = read_pickles(data)
      if len(frame_files)==0:
        print(data, ' - no pickle file found.')
      else:
        for pickle_filename in frame_files:
          pickle_filename_split = pickle_filename.split('/')
          pickle_filename_only = pickle_filename_split[len(pickle_filename_split)-1]
          pickle_filename_only_split = pickle_filename_only.split('_')

          observations_Eoc_corrected, observations_original = get_Eoc_corrected_observations(pickle_filename, iparams)

          #check unitcell
          uc_params = observations_Eoc_corrected.unit_cell().parameters()
          cc_eoc, n_refl_eoc, cc_ori, n_refl_ori, cc_iso_eoc, n_refl_iso_eoc, cc_iso_ori, n_refl_iso_ori = (0,0,0,0,0,0,0,0)
          a,b,c,alpha,beta,gamma = uc_params
          if good_unit_cell(uc_params, iparams.target_unit_cell.parameters(), iparams.iotacc.uc_tolerance):
            try:
              miller_set = symmetry(
                  unit_cell=observations_Eoc_corrected.unit_cell().parameters(),
                  space_group_symbol=iparams.target_space_group
                ).build_miller_set(
                  anomalous_flag=iparams.target_anomalous_flag,
                  d_min=d_min)

              observations_Eoc_corrected = observations_Eoc_corrected.customized_copy(anomalous_flag=iparams.target_anomalous_flag,
                              crystal_symmetry=miller_set.crystal_symmetry())

              miller_set = symmetry(
                  unit_cell=observations_original.unit_cell().parameters(),
                  space_group_symbol=iparams.target_space_group
                ).build_miller_set(
                  anomalous_flag=iparams.target_anomalous_flag,
                  d_min=d_min)

              observations_original = observations_original.customized_copy(anomalous_flag=iparams.target_anomalous_flag,
                                  crystal_symmetry=miller_set.crystal_symmetry())

              #convert to asymmetric unit index
              observations_Eoc_corrected_asu = observations_Eoc_corrected.map_to_asu()
              observations_original_asu = observations_original.map_to_asu()

              #calculate CCori CCEoc
              miller_array_ref = get_observations_ref(iparams.hklrefin)
              cc_eoc, n_refl_eoc = calc_cc(miller_array_ref, observations_Eoc_corrected_asu)
              cc_ori, n_refl_ori = calc_cc(miller_array_ref, observations_original_asu)

              miller_array_iso = get_observations_ref(iparams.hklisoin)
              cc_iso_eoc, n_refl_iso_eoc = calc_cc(miller_array_iso, observations_Eoc_corrected_asu)
              cc_iso_ori, n_refl_iso_ori = calc_cc(miller_array_iso, observations_original_asu)
            except Exception:
              dummy = 0
              #print cryst_id, shot_no, sel_type.ljust(15,' '), 'mismatch spacegroup %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f'%(a,b,c,alpha,beta,gamma)

          if cc_eoc > 0:
            pickle_filename_list.append(pickle_filename)
            cc_list.append(cc_eoc)
            n_refl_list.append(n_refl_eoc)

          txt_out += cryst_id + '  ' + str(shot_no).zfill(iparams.iotacc.LEN_SHOT_SEQ) + \
                 ' %6.2f %6.0f %6.2f %6.2f %6.0f %6.2f %6.2f %6.0f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f'%(\
                  observations_original.d_min(), len(observations_original.data()), cc_ori, cc_eoc, n_refl_eoc, \
                  cc_iso_ori, cc_iso_eoc, n_refl_iso_eoc, a,b,c,alpha,beta,gamma) + ' ' + pickle_filename_only + '\n'

  #select and copy best_by_cc shot
  if len(pickle_filename_list) > 0:
    #sort by cc and select best n results
    cc_thres = flex.max(cc_list) - (flex.max(cc_list)*iparams.iotacc.percent_top_cc/100)
    i_sel = cc_list > cc_thres

    cc_list_sel = cc_list.select(i_sel)
    n_refl_list_sel = n_refl_list.select(i_sel)
    i_seq_sel = flex.int(range(len(cc_list))).select(i_sel)
    pickle_filename_list_sel = [pickle_filename_list[ii_sel] for ii_sel in i_seq_sel]

    #sort by n_refl and select the first
    i_sorted = flex.sort_permutation(n_refl_list_sel, reverse=True)
    n_refl_sel = n_refl_list_sel[i_sorted[0]]
    cc_sel = cc_list_sel[i_sorted[0]]
    pickle_filename_sel = pickle_filename_list_sel[i_sorted[0]]

    pickle_filename_sel_split = pickle_filename_sel.split('/')
    pickle_filename_sel_only = pickle_filename_sel_split[len(pickle_filename_sel_split)-1]
    pickle_filename_out = iparams.run_no + '/' + iparams.iotacc.set_id + '_' + pickle_filename_sel_only

    #shutil.copyfile(pickle_filename_sel, pickle_filename_out)
    txt_out += 'Select --> CC_thres=%6.2f Best CC=%6.2f N_refl=%4.0f '%(cc_thres, cc_sel, n_refl_sel) + pickle_filename_out + '\n'
    flag_success = True

  else:
    txt_out += 'No image with CC > 0\n'

  if flag_success:
    return pickle_filename_sel, txt_out
  else:
    return None


if __name__=="__main__":
  iparams, txt_out_input = read_input(sys.argv[:1])
  print(txt_out_input)

  if iparams.iotacc.set_id is not None:
    data_dir = iparams.data[0] + '/' +iparams.iotacc.set_id
  else:
    data_dir = iparams.data[0]

  cryst_ids = []
  for cryst_id in os.listdir(data_dir):
    cryst_ids.append(cryst_id)

  shot_nos = range(iparams.iotacc.n_shots)
  run_no = iparams.run_no

  if os.path.exists(run_no):
    shutil.rmtree(run_no)

  os.makedirs(run_no)

  txt_out_log = 'iotacc\n'
  txt_out_log += txt_out_input
  txt_out_pickle_filename_sel = ''
  for cryst_id in cryst_ids:
    def select_best_by_cc_mproc_wrapper(arg):
      return select_best_by_cc_mproc(arg, cryst_id, iparams)

    select_best_by_cc_results = pool_map(
            args=shot_nos,
            func=select_best_by_cc_mproc_wrapper,
            processes=iparams.n_processors)


    for result in select_best_by_cc_results:
      if result is not None:
        pickle_filename_sel, txt_out = result
        txt_out_log += txt_out
        txt_out_pickle_filename_sel += pickle_filename_sel + '\n'
        print(txt_out)

  f = open(run_no+'/pickle_selected.txt', 'w')
  f.write(txt_out_pickle_filename_sel)
  f.close()

  f = open(run_no+'/log.txt', 'w')
  f.write(txt_out_log)
  f.close()
