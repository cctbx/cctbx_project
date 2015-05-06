from __future__ import division
from scitbx.matrix import sqr
from cctbx.uctbx import unit_cell
from cctbx import miller
from cctbx import crystal
from cctbx.array_family import flex
from iotbx import mtz
from iotbx import reflection_file_reader
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from libtbx.utils import Sorry
from cctbx.uctbx import unit_cell

class intensities_scaler(object):
  '''
  Author      : Uervirojnangkoorn, M.
  Created     : 7/13/2014
  Merge equivalent reflections and report intensity and refinement statistics.
  '''

  def __init__(self):
    '''
    Constructor
    '''
    self.CONST_SE_MIN_WEIGHT = 0.17
    self.CONST_SE_MAX_WEIGHT = 1.0
    self.CONST_SIG_I_FACTOR = 1.5

  def calc_avg_I_cpp(self, group_no, group_id_list, miller_index, miller_indices_ori, I, sigI, G, B,
                     p_set, rs_set, wavelength_set, sin_theta_over_lambda_sq, SE, avg_mode,
                     iparams, pickle_filename_set):
    from prime import Average_Mode, averaging_engine
    if avg_mode == 'average': avg_mode_cpp = Average_Mode.Average
    elif avg_mode == 'weighted': avg_mode_cpp = Average_Mode.Weighted
    elif avg_mode == 'final': avg_mode_cpp = Average_Mode.Final
    else: raise Sorry("Bad averaging mode selected: %s"%avg_mode)

    sigma_max = iparams.sigma_rejection

    engine = averaging_engine(group_no, group_id_list, miller_index, miller_indices_ori, I, sigI, G, B,p_set, rs_set, wavelength_set, sin_theta_over_lambda_sq, SE, pickle_filename_set)
    engine.avg_mode = avg_mode_cpp
    engine.sigma_max = sigma_max
    engine.flag_volume_correction = iparams.flag_volume_correction
    engine.n_rejection_cycle = iparams.n_rejection_cycle
    engine.flag_output_verbose = iparams.flag_output_verbose

    results = engine.calc_avg_I()

    return results.miller_index, results.I_avg, results.sigI_avg, (results.r_meas_w_top, results.r_meas_w_btm, results.r_meas_top, results.r_meas_btm, results.multiplicity), (results.I_avg_even, results.I_avg_odd, results.I_avg_even_h, results.I_avg_odd_h, results.I_avg_even_k, results.I_avg_odd_k, results.I_avg_even_l, results.I_avg_odd_l), results.txt_obs_out, results.txt_reject_out


  def calc_mean_unit_cell(self, results):
    a_all = flex.double()
    b_all = flex.double()
    c_all = flex.double()
    alpha_all = flex.double()
    beta_all = flex.double()
    gamma_all = flex.double()
    for pres in results:
      if pres is not None:
        a_all.append(pres.uc_params[0])
        b_all.append(pres.uc_params[1])
        c_all.append(pres.uc_params[2])
        alpha_all.append(pres.uc_params[3])
        beta_all.append(pres.uc_params[4])
        gamma_all.append(pres.uc_params[5])

    uc_mean = flex.double([np.mean(a_all), np.mean(b_all), np.mean(c_all), np.mean(alpha_all), np.mean(beta_all), np.mean(gamma_all)])
    uc_med = flex.double([np.median(a_all), np.median(b_all), np.median(c_all), np.median(alpha_all), np.median(beta_all), np.median(gamma_all)])
    uc_std = flex.double([np.std(a_all), np.std(b_all), np.std(c_all), np.std(alpha_all), np.std(beta_all), np.std(gamma_all)])

    return uc_mean, uc_med, uc_std

  def calc_mean_postref_parameters(self, results):
    G_all = flex.double()
    B_all = flex.double()
    rotx_all = flex.double()
    roty_all = flex.double()
    ry_all = flex.double()
    rz_all = flex.double()
    re_all = flex.double()
    r0_all = flex.double()
    R_final_all = flex.double()
    R_xy_final_all = flex.double()
    SE_all = flex.double()
    for pres in results:
      if pres is not None:
        if not math.isnan(pres.G):
          G_all.append(pres.G)
        if not math.isnan(pres.B):
          B_all.append(pres.B)
        if not math.isnan(pres.rotx):
          rotx_all.append(pres.rotx)
        if not math.isnan(pres.roty):
          roty_all.append(pres.roty)
        if not math.isnan(pres.ry):
          ry_all.append(pres.ry)
        if not math.isnan(pres.rz):
          rz_all.append(pres.rz)
        if not math.isnan(pres.re):
          re_all.append(pres.re)
        if not math.isnan(pres.r0):
          r0_all.append(pres.r0)
        if not math.isnan(pres.R_final):
          R_final_all.append(pres.R_final)
        if not math.isnan(pres.R_xy_final):
          R_xy_final_all.append(pres.R_xy_final)
        if not math.isnan(pres.SE):
          SE_all.append(pres.SE)

    pr_params_mean = flex.double([np.mean(G_all), np.mean(B_all),
                                  np.mean(flex.abs(ry_all)), np.mean(flex.abs(rz_all)),
                                  np.mean(flex.abs(re_all)), np.mean(flex.abs(r0_all)),
                                  np.mean(flex.abs(rotx_all)), np.mean(flex.abs(roty_all)),
                                  np.mean(R_final_all), np.mean(R_xy_final_all),
                                  np.mean(SE_all)])
    pr_params_med = flex.double([np.median(G_all), np.median(B_all),
                                  np.median(flex.abs(ry_all)), np.median(flex.abs(rz_all)),
                                  np.median(flex.abs(re_all)), np.median(flex.abs(r0_all)),
                                  np.median(flex.abs(rotx_all)), np.median(flex.abs(roty_all)),
                                  np.median(R_final_all), np.median(R_xy_final_all),
                                  np.median(SE_all)])
    pr_params_std = flex.double([np.std(G_all), np.std(B_all),
                                  np.std(flex.abs(ry_all)), np.std(flex.abs(rz_all)),
                                  np.std(flex.abs(re_all)), np.std(flex.abs(r0_all)),
                                  np.std(flex.abs(rotx_all)), np.std(flex.abs(roty_all)),
                                  np.std(R_final_all), np.std(R_xy_final_all),
                                  np.std(SE_all)])

    return pr_params_mean, pr_params_med, pr_params_std

  def prepare_output(self, results, iparams, avg_mode):
    if avg_mode == 'average':
      cc_thres = 0
    else:
      cc_thres = iparams.frame_accept_min_cc

    std_filter = iparams.sigma_rejection

    if avg_mode == 'final':
      target_anomalous_flag = iparams.target_anomalous_flag
    else:
      target_anomalous_flag = False

    pr_params_mean, pr_params_med, pr_params_std = self.calc_mean_postref_parameters(results)
    G_mean, B_mean, ry_mean, rz_mean, re_mean, r0_mean, rotx_mean, roty_mean, R_mean, R_xy_mean, SE_mean = pr_params_mean
    G_med, B_med, ry_med, rz_med, re_med, r0_med, rotx_med, roty_med, R_med, R_xy_med, SE_med = pr_params_med
    G_std, B_std, ry_std, rz_std, re_std, r0_std, rotx_std, roty_std, R_std, R_xy_std, SE_std = pr_params_std

    #prepare data for merging
    miller_indices_all = flex.miller_index()
    miller_indices_ori_all = flex.miller_index()
    I_all = flex.double()
    sigI_all = flex.double()
    G_all = flex.double()
    B_all = flex.double()
    p_all = flex.double()
    rx_all = flex.double()
    rs_all = flex.double()
    rh_all = flex.double()
    SE_all = flex.double()
    sin_sq_all = flex.double()
    wavelength_all = flex.double()
    R_init_all = flex.double()
    R_final_all = flex.double()
    R_xy_init_all = flex.double()
    R_xy_final_all = flex.double()
    pickle_filename_all = []
    filtered_results = []
    cn_good_frame, cn_bad_frame_SE, cn_bad_frame_uc, cn_bad_frame_cc, cn_bad_frame_G, cn_bad_frame_re = (0,0,0,0,0,0)
    i_seq = flex.int()
    for pres in results:
      if pres is not None:
        pickle_filepath = pres.pickle_filename.split('/')
        img_filename = pickle_filepath[len(pickle_filepath)-1]
        flag_pres_ok = True

        #check SE, CC, UC, G, B, gamma_e
        if math.isnan(pres.G):
          flag_pres_ok = False

        if math.isnan(pres.SE) or np.isinf(pres.SE):
          flag_pres_ok = False

        if flag_pres_ok and SE_med > 0:
          if abs(pres.SE-SE_med)/SE_std > std_filter:
            flag_pres_ok = False
            cn_bad_frame_SE += 1
            print pres.frame_no, img_filename, ' discarded (SE = %6.2f)'%(pres.SE)

        if flag_pres_ok and pres.CC_final < cc_thres:
          flag_pres_ok = False
          cn_bad_frame_cc += 1
          print pres.frame_no, img_filename, ' discarded (CC = %6.2f)'%(pres.CC_final)

        if flag_pres_ok:
          if G_std > 0:
            if abs(pres.G-G_med)/G_std > std_filter:
              flag_pres_ok = False
              cn_bad_frame_G += 1
              print pres.frame_no, img_filename, ' discarded (G = %6.2f)'%(pres.G)

        if flag_pres_ok:
          if re_std > 0:
            if abs(pres.re-re_med)/(math.sqrt(re_med)) > std_filter:
              flag_pres_ok = False
              cn_bad_frame_re += 1
              print pres.frame_no, img_filename, ' discarded (gamma_e = %6.2f)'%(pres.re)

        from mod_leastsqr import good_unit_cell
        if flag_pres_ok and not good_unit_cell(pres.uc_params, iparams, iparams.merge.uc_tolerance):
          flag_pres_ok = False
          cn_bad_frame_uc += 1
          print pres.frame_no, img_filename, ' discarded - unit-cell exceeds the limits (%6.2f %6.2f %6.2f %5.2f %5.2f %5.2f)'%(pres.uc_params[0], pres.uc_params[1], pres.uc_params[2], pres.uc_params[3], pres.uc_params[4], pres.uc_params[5])

        if flag_pres_ok:
          cn_good_frame += 1
          sin_theta_over_lambda_sq = pres.observations.two_theta(wavelength=pres.wavelength).sin_theta_over_lambda_sq().data()
          filtered_results.append(pres)
          R_init_all.append(pres.R_init)
          R_final_all.append(pres.R_final)
          R_xy_init_all.append(pres.R_xy_init)
          R_xy_final_all.append(pres.R_xy_final)

          miller_indices_all.extend(pres.observations.indices())
          miller_indices_ori_all.extend(pres.observations_original.indices())
          I_all.extend(pres.observations.data())
          sigI_all.extend(pres.observations.sigmas())
          G_all.extend(flex.double([pres.G]*len(pres.observations.data())))
          B_all.extend(flex.double([pres.B]*len(pres.observations.data())))
          p_all.extend(pres.partiality)
          rs_all.extend(pres.rs_set)
          rh_all.extend(pres.rh_set)
          sin_sq_all.extend(sin_theta_over_lambda_sq)
          SE_all.extend(flex.double([pres.SE]*len(pres.observations.data())))
          wavelength_all.extend(flex.double([pres.wavelength]*len(pres.observations.data())))
          pickle_filename_all += [pres.pickle_filename for i in range(len(pres.observations.data()))]
          i_seq.extend(flex.int([i for i in range(len(i_seq), len(i_seq)+len(pres.observations.data()))]))
          print pres.frame_no+1, img_filename, ' merged'


    #plot stats
    self.plot_stats(filtered_results, iparams)

    #calculate average unit cell
    uc_mean, uc_med, uc_std = self.calc_mean_unit_cell(filtered_results)
    unit_cell_mean = unit_cell((uc_mean[0], uc_mean[1], uc_mean[2], uc_mean[3], uc_mean[4], uc_mean[5]))

    #recalculate stats for pr parameters
    pr_params_mean, pr_params_med, pr_params_std = self.calc_mean_postref_parameters(filtered_results)
    G_mean, B_mean, ry_mean, rz_mean, re_mean, r0_mean, rotx_mean, roty_mean, R_mean, R_xy_mean, SE_mean = pr_params_mean
    G_med, B_med, ry_med, rz_med, re_med, r0_med, rotx_med, roty_med, R_med, R_xy_med, SE_med = pr_params_med
    G_std, B_std, ry_std, rz_std, re_std, r0_std, rotx_std, roty_std, R_std, R_xy_std, SE_std = pr_params_std

    #from all observations merge them
    crystal_symmetry = crystal.symmetry(
        unit_cell=(uc_mean[0], uc_mean[1], uc_mean[2], uc_mean[3], uc_mean[4], uc_mean[5]),
        space_group_symbol=iparams.target_space_group)
    miller_set_all=miller.set(
                crystal_symmetry=crystal_symmetry,
                indices=miller_indices_all,
                anomalous_flag=target_anomalous_flag)
    miller_array_all = miller_set_all.array(
              data=I_all,
              sigmas=sigI_all).set_observation_type_xray_intensity()

    #sort reflections according to asymmetric-unit symmetry hkl
    perm = miller_array_all.sort_permutation(by_value="packed_indices")
    miller_indices_all_sort = miller_array_all.indices().select(perm)
    miller_indices_ori_all_sort = miller_indices_ori_all.select(perm)
    I_obs_all_sort = miller_array_all.data().select(perm)
    sigI_obs_all_sort = miller_array_all.sigmas().select(perm)
    G_all_sort = G_all.select(perm)
    B_all_sort = B_all.select(perm)
    p_all_sort = p_all.select(perm)
    rs_all_sort = rs_all.select(perm)
    wavelength_all_sort = wavelength_all.select(perm)
    sin_sq_all_sort = sin_sq_all.select(perm)
    SE_all_sort = SE_all.select(perm)
    i_seq_sort = i_seq.select(perm)
    pickle_filename_all_sort = [pickle_filename_all[i] for i in i_seq_sort]

    miller_array_uniq = miller_array_all.merge_equivalents().array().complete_array(d_min=iparams.merge.d_min, d_max=iparams.merge.d_max)
    matches_uniq = miller.match_multi_indices(
                  miller_indices_unique=miller_array_uniq.indices(),
                  miller_indices=miller_indices_all_sort)

    pair_0 = flex.int([pair[0] for pair in matches_uniq.pairs()])
    pair_1 = flex.int([pair[1] for pair in matches_uniq.pairs()])
    group_id_list = flex.int([pair_0[pair_1[i]] for i in range(len(matches_uniq.pairs()))])

    from collections import Counter
    tally = Counter()
    for elem in group_id_list:
      tally[elem] += 1
    cn_group = len(tally)

    txt_out = 'Summary of refinement and merging\n'
    txt_out += ' No. good frames:          %12.0f\n'%(cn_good_frame)
    txt_out += ' No. bad cc frames:        %12.0f\n'%(cn_bad_frame_cc)
    txt_out += ' No. bad G frames) :       %12.0f\n'%(cn_bad_frame_G)
    txt_out += ' No. bad unit cell frames: %12.0f\n'%(cn_bad_frame_uc)
    txt_out += ' No. bad gamma_e frames:   %12.0f\n'%(cn_bad_frame_re)
    txt_out += ' No. bad SE:               %12.0f\n'%(cn_bad_frame_SE)
    txt_out += ' No. observations:         %12.0f\n'%(len(I_obs_all_sort))
    txt_out += 'Mean target value (BEFORE: Mean Median (Std.))\n'
    txt_out += ' post-refinement:          %12.2f %12.2f (%9.2f)\n'%(np.mean(R_init_all), np.median(R_init_all), np.std(R_init_all))
    txt_out += ' (x,y) restraints:         %12.2f %12.2f (%9.2f)\n'%(np.mean(R_xy_init_all), np.median(R_xy_init_all), np.std(R_xy_init_all))
    txt_out += 'Mean target value (AFTER: Mean Median (Std.))\n'
    txt_out += ' post-refinement:          %12.2f %12.2f (%9.2f)\n'%(np.mean(R_final_all), np.median(R_final_all), np.std(R_final_all))
    txt_out += ' (x,y) restraints:         %12.2f %12.2f (%9.2f)\n'%(np.mean(R_xy_final_all), np.median(R_xy_final_all), np.std(R_xy_final_all))
    txt_out += ' SE:                       %12.2f %12.2f (%9.2f)\n'%(SE_mean, SE_med, SE_std)
    txt_out += ' G:                        %12.2f %12.2f (%9.2f)\n'%(G_mean, G_med, G_std)
    txt_out += ' B:                        %12.2f %12.2f (%9.2f)\n'%(B_mean, B_med, B_std)
    txt_out += ' Rot.x:                    %12.2f %12.2f (%9.2f)\n'%(rotx_mean*180/math.pi, rotx_med*180/math.pi, rotx_std*180/math.pi)
    txt_out += ' Rot.y:                    %12.2f %12.2f (%9.2f)\n'%(roty_mean*180/math.pi, roty_med*180/math.pi, roty_std*180/math.pi)
    txt_out += ' gamma_y:                  %12.5f %12.5f (%9.5f)\n'%(ry_mean, ry_med, ry_std)
    txt_out += ' gamma_z:                  %12.5f %12.5f (%9.5f)\n'%(rz_mean, rz_med, rz_std)
    txt_out += ' gamma_0:                  %12.5f %12.5f (%9.5f)\n'%(r0_mean, r0_med, r0_std)
    txt_out += ' gamma_e:                  %12.5f %12.5f (%9.5f)\n'%(re_mean, re_med, re_std)
    txt_out += ' unit cell\n'
    txt_out += '   a:                      %12.2f %12.2f (%9.2f)\n'%(uc_mean[0], uc_med[0], uc_std[0])
    txt_out += '   b:                      %12.2f %12.2f (%9.2f)\n'%(uc_mean[1], uc_med[1], uc_std[1])
    txt_out += '   c:                      %12.2f %12.2f (%9.2f)\n'%(uc_mean[2], uc_med[2], uc_std[2])
    txt_out += '   alpha:                  %12.2f %12.2f (%9.2f)\n'%(uc_mean[3], uc_med[3], uc_std[3])
    txt_out += '   beta:                   %12.2f %12.2f (%9.2f)\n'%(uc_mean[4], uc_med[4], uc_std[4])
    txt_out += '   gamma:                  %12.2f %12.2f (%9.2f)\n'%(uc_mean[5], uc_med[5], uc_std[5])
    txt_out += '* (standard deviation)\n'

    return cn_group, group_id_list, miller_indices_all_sort, miller_indices_ori_all_sort, \
           I_obs_all_sort, sigI_obs_all_sort,G_all_sort, B_all_sort, \
           p_all_sort, rs_all_sort, wavelength_all_sort, sin_sq_all_sort, SE_all_sort, uc_mean, \
           np.mean(wavelength_all), pickle_filename_all_sort, txt_out


  def write_output(self, miller_indices_merge, I_merge, sigI_merge, stat_all,
                   I_two_halves_tuple, iparams, uc_mean, wavelength_mean,
                   output_mtz_file_prefix, avg_mode):

    if avg_mode == 'final':
      target_anomalous_flag = iparams.target_anomalous_flag
    else:
      target_anomalous_flag = False

    #extract I_even and I_odd pair
    I_even, I_odd, I_even_h, I_odd_h, I_even_k, I_odd_k, I_even_l, I_odd_l = I_two_halves_tuple

    #output mtz file and report binning stat
    miller_set_merge = crystal.symmetry(
          unit_cell=unit_cell((uc_mean[0],uc_mean[1],uc_mean[2],uc_mean[3],uc_mean[4],uc_mean[5])),
          space_group_symbol=iparams.target_space_group
        ).build_miller_set(
          anomalous_flag=target_anomalous_flag,
          d_min=iparams.merge.d_min)

    miller_array_merge = miller_set_merge.array().customized_copy(indices=miller_indices_merge,
              data=I_merge, sigmas=sigI_merge, anomalous_flag=target_anomalous_flag).map_to_asu().set_observation_type_xray_intensity()
    miller_array_complete = miller_set_merge.array()
    fake_data = flex.double([1.0]*len(miller_array_complete.indices()))
    miller_array_template_asu = miller_array_complete.customized_copy(data=fake_data, \
              sigmas=fake_data).resolution_filter(d_min=iparams.merge.d_min, \
              d_max=iparams.merge.d_max)

    #do another resolution filter here
    i_sel_res = miller_array_merge.resolution_filter_selection(d_min=iparams.merge.d_min, d_max=iparams.merge.d_max)
    miller_indices_merge = miller_array_merge.indices().select(i_sel_res)
    I_merge = I_merge.select(i_sel_res)
    sigI_merge = sigI_merge.select(i_sel_res)
    miller_array_merge = miller_array_merge.customized_copy(indices=miller_indices_merge,
        data=I_merge,
        sigmas=sigI_merge)
    i_seq=flex.int(range(0,len(i_sel_res)))
    i_sel_seq=i_seq.select(i_sel_res == True)
    stat_all_tmp = [stat_all[j] for j in i_seq.select(i_sel_res == True)]
    stat_all = stat_all_tmp
    I_even = I_even.select(i_sel_res)
    I_odd = I_odd.select(i_sel_res)
    I_even_h = I_even_h.select(i_sel_res)
    I_odd_h = I_odd_h.select(i_sel_res)
    I_even_k = I_even_k.select(i_sel_res)
    I_odd_k = I_odd_k.select(i_sel_res)
    I_even_l = I_even_l.select(i_sel_res)
    I_odd_l = I_odd_l.select(i_sel_res)

    print 'N_refl after resolution filter', len(miller_indices_merge)

    #remove outliers
    I_bin_sigma_filter = 10
    n_bin_sigma_filter = 200
    n_rejection_iterations = iparams.n_rejection_cycle

    for i_rejection in range(n_rejection_iterations):
      binner_merge = miller_array_merge.setup_binner(n_bins=n_bin_sigma_filter)
      binner_merge_indices = binner_merge.bin_indices()
      i_seq = flex.int([j for j in range(len(miller_array_merge.data()))])
      i_seq_sel = flex.int()
      print 'Outlier rejection cycle %3.0f'%(i_rejection+1)
      print 'N_refl before outlier rejection', len(miller_array_merge.indices())
      #print 'Bin  Nrefl  Nrefl_rejected      <I>       std(I)'
      for i_bin in range(1,n_bin_sigma_filter+1):
        i_binner = (binner_merge_indices == i_bin)
        if len(miller_array_merge.data().select(i_binner)) > 0:
          I_obs_bin = miller_array_merge.data().select(i_binner)
          i_seq_bin = i_seq.select(i_binner)
          med_I_bin = np.median(I_obs_bin)
          i_filter = flex.abs((I_obs_bin - med_I_bin)/np.std(I_obs_bin)) > I_bin_sigma_filter
          i_seq_bin_filter = i_seq_bin.select(i_filter)
          i_seq_sel.extend(i_seq_bin_filter)
          #print '%2.0f %6.0f %6.0f %14.2f %12.2f'%(i_bin, len(I_obs_bin), len(i_seq_bin_filter), med_I_bin, math.sqrt(med_I_bin)*2.5)

      flag_sel = flex.bool([True]*len(miller_array_merge.data()))
      for i_i_seq_sel in i_seq_sel:
        flag_sel[i_i_seq_sel] = False

      miller_array_merge = miller_array_merge.customized_copy(indices=miller_array_merge.indices().select(flag_sel),
                                                              data=miller_array_merge.data().select(flag_sel),
                                                              sigmas=miller_array_merge.sigmas().select(flag_sel))

      stat_all_tmp = [stat_all[j] for j in i_seq.select(flag_sel == True)]
      stat_all = stat_all_tmp
      I_even = I_even.select(flag_sel)
      I_odd = I_odd.select(flag_sel)
      I_even_h = I_even_h.select(flag_sel)
      I_odd_h = I_odd_h.select(flag_sel)
      I_even_k = I_even_k.select(flag_sel)
      I_odd_k = I_odd_k.select(flag_sel)
      I_even_l = I_even_l.select(flag_sel)
      I_odd_l = I_odd_l.select(flag_sel)

      print 'N_refl after outlier rejection', len(miller_array_merge.indices())

    #get iso if given
    flag_hklisoin_found = False
    if iparams.hklisoin is not None:
      flag_hklisoin_found = True
      from iotbx import reflection_file_reader
      reflection_file_iso = reflection_file_reader.any_reflection_file(iparams.hklisoin)
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


    #write output files
    if output_mtz_file_prefix != '':
      #write as mtz file
      miller_array_merge_unique = miller_array_merge.merge_equivalents().array()
      mtz_dataset_merge = miller_array_merge_unique.as_mtz_dataset(column_root_label="IOBS")
      mtz_dataset_merge.mtz_object().write(file_name=output_mtz_file_prefix+'_merge.mtz')
      #write as cns file
      f_cns = open(output_mtz_file_prefix+'_merge.hkl', 'w')
      miller_array_merge_unique.export_as_cns_hkl(file_object=f_cns)
      f_cns.close()


    #calculate total cc_anom for two halves
    cc_anom_acentric, cc_anom_centric, nrefl_anom_acentric, nrefl_anom_centric = (0,0,0,0)
    if miller_array_merge.anomalous_flag():
      miller_array_merge_even = miller_array_merge.customized_copy(data = I_even)
      miller_array_merge_odd = miller_array_merge.customized_copy(data = I_odd)
      ma_anom_dif_even = miller_array_merge_even.anomalous_differences()
      ma_anom_dif_odd = miller_array_merge_odd.anomalous_differences()
      ma_anom_centric_flags = ma_anom_dif_even.centric_flags()
      anom_dif_even = ma_anom_dif_even.data()
      anom_dif_odd = ma_anom_dif_odd.data()
      anom_dif_centric_flags = ma_anom_centric_flags.data()
      i_acentric = (anom_dif_centric_flags == False)
      i_centric = (anom_dif_centric_flags == True)

      anom_dif_even_acentric = anom_dif_even.select(i_acentric)
      anom_dif_even_centric = anom_dif_even.select(i_centric)
      anom_dif_odd_acentric = anom_dif_odd.select(i_acentric)
      anom_dif_odd_centric = anom_dif_odd.select(i_centric)
      mat_anom_acentric = np.corrcoef(anom_dif_even_acentric, anom_dif_odd_acentric)
      mat_anom_centric = np.corrcoef(anom_dif_even_centric, anom_dif_odd_centric)
      if len(mat_anom_acentric) > 0:
        cc_anom_acentric = mat_anom_acentric[0,1]
      if len(mat_anom_centric) > 0:
        cc_anom_centric = mat_anom_centric[0,1]
      nrefl_anom_acentric = len(anom_dif_even_acentric)
      nrefl_anom_centric = len(anom_dif_even_centric)

      if iparams.flag_plot:
        plt.subplot(211)
        plt.scatter(anom_dif_even_acentric, anom_dif_odd_acentric,s=10, marker='x', c='r')
        plt.title('CCanoma=%5.2f N_refl=%6.0f'%(cc_anom_acentric, nrefl_anom_acentric))
        plt.xlabel('delta_I_even')
        plt.ylabel('delta_I_odd')
        plt.subplot(212)
        plt.scatter(anom_dif_even_centric, anom_dif_odd_centric,s=10, marker='x', c='r')
        plt.title('CCanomc=%5.2f N_refl=%6.0f'%(cc_anom_centric, nrefl_anom_centric))
        plt.xlabel('delta_I_even')
        plt.ylabel('delta_I_odd')
        plt.show()

    #select single cone reflections on the three crystal axes
    fraction_percent = iparams.percent_cone_fraction
    miller_array_merge_astar = miller_array_merge.remove_cone(fraction_percent, axis_point_2=(1,0,0), negate=True)
    miller_array_merge_bstar = miller_array_merge.remove_cone(fraction_percent, axis_point_2=(0,1,0), negate=True)
    miller_array_merge_cstar = miller_array_merge.remove_cone(fraction_percent, axis_point_2=(0,0,1), negate=True)

    #report binning stats
    binner_template_asu = miller_array_template_asu.setup_binner(n_bins=iparams.n_bins)
    binner_template_asu_indices = binner_template_asu.bin_indices()

    txt_out = '\n'
    txt_out += 'Summary for '+output_mtz_file_prefix+'_merge.mtz\n'
    txt_out += 'Bin Resolution Range     Completeness      <N_obs> |Rmerge  Rsplit   CC1/2   N_ind |CCiso   N_ind|CCanoma  N_ind CCanomc  N_ind| <I/sigI>   <I>\n'
    txt_out += '--------------------------------------------------------------------------------------------------------------------------------------------------\n'
    sum_r_meas_w_top, sum_r_meas_w_btm, sum_r_meas_top, sum_r_meas_btm, sum_refl_obs, sum_refl_complete, n_refl_obs_total  = (0,0,0,0,0,0,0)
    cc12_list = []
    avgI_star_list = []
    cc12_n_refl_list = []
    avgI_list = []
    avgI_n_refl_list = []
    I_even_astar = flex.double()
    I_odd_astar = flex.double()
    I_even_bstar = flex.double()
    I_odd_bstar = flex.double()
    I_even_cstar = flex.double()
    I_odd_cstar = flex.double()
    for i in range(1,iparams.n_bins+1):
      i_binner = (binner_template_asu_indices == i)
      miller_indices_bin = miller_array_template_asu.indices().select(i_binner)

      matches_template = miller.match_multi_indices(
                  miller_indices_unique=miller_indices_bin,
                  miller_indices=miller_array_merge.indices())

      I_bin = flex.double([miller_array_merge.data()[pair[1]] for pair in matches_template.pairs()])
      sigI_bin = flex.double([miller_array_merge.sigmas()[pair[1]] for pair in matches_template.pairs()])
      miller_indices_obs_bin = flex.miller_index([miller_array_merge.indices()[pair[1]] for pair in matches_template.pairs()])

      #caculate CCanom for the two halves.
      cc_anom_bin_acentric, cc_anom_bin_centric, nrefl_anom_bin_acentric, nrefl_anom_bin_centric = (0,0,0,0)
      if miller_array_merge.anomalous_flag():
        matches_anom_dif_even = miller.match_multi_indices(
                    miller_indices_unique=miller_indices_bin,
                    miller_indices=ma_anom_dif_even.indices())
        anom_dif_even = flex.double([ma_anom_dif_even.data()[pair[1]] for pair in matches_anom_dif_even.pairs()])
        anom_dif_centric_flags = flex.bool([ma_anom_centric_flags.data()[pair[1]] for pair in matches_anom_dif_even.pairs()])

        matches_anom_dif_odd = miller.match_multi_indices(
                    miller_indices_unique=miller_indices_bin,
                    miller_indices=ma_anom_dif_odd.indices())
        anom_dif_odd = flex.double([ma_anom_dif_odd.data()[pair[1]] for pair in matches_anom_dif_odd.pairs()])
        i_acentric = (anom_dif_centric_flags == False)
        i_centric = (anom_dif_centric_flags == True)

        anom_dif_even_acentric = anom_dif_even.select(i_acentric)
        anom_dif_even_centric = anom_dif_even.select(i_centric)
        anom_dif_odd_acentric = anom_dif_odd.select(i_acentric)
        anom_dif_odd_centric = anom_dif_odd.select(i_centric)

        mat_anom_acentric = np.corrcoef(anom_dif_even_acentric, anom_dif_odd_acentric)
        mat_anom_centric = np.corrcoef(anom_dif_even_centric, anom_dif_odd_centric)
        if len(mat_anom_acentric) > 0:
          cc_anom_bin_acentric = mat_anom_acentric[0,1]
        if len(mat_anom_centric) > 0:
          cc_anom_bin_centric = mat_anom_centric[0,1]
        nrefl_anom_bin_acentric = len(anom_dif_even_acentric)
        nrefl_anom_bin_centric = len(anom_dif_even_centric)

      #prepare the calculation of these parameters
      mean_i_over_sigi_bin, multiplicity_bin, r_meas_w_bin, r_meas_bin, n_refl_cc12_bin, r_split_bin  = (0,0,0,0,0,0)
      cc12_bin, cc12_bin_astar, cc12_bin_bstar, cc12_bin_cstar = (0,0,0,0)
      avgI_bin, avgI_bin_h, avgI_bin_k, avgI_bin_l = (0,0,0,0)
      avgI_bin_astar, avgI_bin_bstar, avgI_bin_cstar = (0,0,0)
      if len(I_bin) > 0:
        #calculate <I>
        avgI_bin = flex.mean(I_bin)
        mean_i_over_sigi_bin = flex.mean(I_bin/sigI_bin)

        #calculation of Rmeas
        stat_bin = [stat_all[pair[1]] for pair in matches_template.pairs()]
        sum_r_meas_w_top_bin, sum_r_meas_w_btm_bin, sum_r_meas_top_bin, sum_r_meas_btm_bin, sum_mul_bin = (0,0,0,0,0)
        for stat in stat_bin:
          r_meas_w_top, r_meas_w_btm, r_meas_top, r_meas_btm, mul = stat

          sum_r_meas_w_top_bin += r_meas_w_top
          sum_r_meas_w_btm_bin += r_meas_w_btm
          sum_r_meas_top_bin += r_meas_top
          sum_r_meas_btm_bin += r_meas_btm
          sum_mul_bin += mul
          sum_r_meas_w_top += r_meas_w_top
          sum_r_meas_w_btm += r_meas_w_btm
          sum_r_meas_top += r_meas_top
          sum_r_meas_btm += r_meas_btm
          n_refl_obs_total += mul

        multiplicity_bin = sum_mul_bin/len(I_bin)
        if sum_r_meas_w_btm_bin > 0:
          r_meas_w_bin = sum_r_meas_w_top_bin/ sum_r_meas_w_btm_bin
        else:
          r_meas_w_bin = float('Inf')

        if sum_r_meas_btm_bin > 0:
          r_meas_bin = sum_r_meas_top_bin/ sum_r_meas_btm_bin
        else:
          r_meas_bin = float('Inf')

        #filter I_even and I_odd for this bin and select only >0 values
        I_even_bin = flex.double([I_even[pair[1]] for pair in matches_template.pairs()])
        I_odd_bin = flex.double([I_odd[pair[1]] for pair in matches_template.pairs()])

        i_even_filter_sel = (I_even_bin > 0)
        I_even_bin_sel = I_even_bin.select(i_even_filter_sel)
        I_odd_bin_sel = I_odd_bin.select(i_even_filter_sel)
        miller_indices_bin_sel = miller_indices_obs_bin.select(i_even_filter_sel)
        n_refl_halfset_bin = len(I_even_bin_sel)

        #calculation of CC1/2 and <I> on the three crystal axes (a*, b*, c*) ---------------
        if n_refl_halfset_bin > 0:
          cc12_bin = np.corrcoef(I_even_bin_sel, I_odd_bin_sel)[0,1]
          r_split_bin = (1/math.sqrt(2))*(flex.sum(flex.abs(I_even_bin_sel - I_odd_bin_sel))/(flex.sum(I_even_bin_sel + I_odd_bin_sel)*0.5))

        I_even_bin_astar = flex.double()
        I_odd_bin_astar = flex.double()
        I_even_bin_bstar = flex.double()
        I_odd_bin_bstar = flex.double()
        I_even_bin_cstar = flex.double()
        I_odd_bin_cstar = flex.double()
        try:
          matches_astar = miller.match_multi_indices(
                      miller_indices_unique=miller_array_merge_astar.indices(),
                      miller_indices=miller_indices_bin_sel)
          I_even_bin_astar = flex.double([I_even_bin_sel[pair[1]] for pair in matches_astar.pairs()])
          I_odd_bin_astar = flex.double([I_odd_bin_sel[pair[1]] for pair in matches_astar.pairs()])

          if len(I_even_bin_astar) > 0:
            cc12_bin_astar = np.corrcoef(I_even_bin_astar, I_odd_bin_astar)[0,1]
            avgI_bin_astar = flex.mean(I_even_bin_astar)
          I_even_astar.extend(I_even_bin_astar)
          I_odd_astar.extend(I_odd_bin_astar)
        except Exception:
          dummy = 1

        try:
          matches_bstar = miller.match_multi_indices(
                      miller_indices_unique=miller_array_merge_bstar.indices(),
                      miller_indices=miller_indices_bin_sel)
          I_even_bin_bstar = flex.double([I_even_bin_sel[pair[1]] for pair in matches_bstar.pairs()])
          I_odd_bin_bstar = flex.double([I_odd_bin_sel[pair[1]] for pair in matches_bstar.pairs()])

          if len(I_even_bin_bstar) > 0:
            cc12_bin_bstar = np.corrcoef(I_even_bin_bstar, I_odd_bin_bstar)[0,1]
            avgI_bin_bstar = flex.mean(I_even_bin_bstar)
          I_even_bstar.extend(I_even_bin_bstar)
          I_odd_bstar.extend(I_odd_bin_bstar)
        except Exception:
          dummy = 1

        try:
          matches_cstar = miller.match_multi_indices(
                      miller_indices_unique=miller_array_merge_cstar.indices(),
                      miller_indices=miller_indices_bin_sel)
          I_even_bin_cstar = flex.double([I_even_bin_sel[pair[1]] for pair in matches_cstar.pairs()])
          I_odd_bin_cstar = flex.double([I_odd_bin_sel[pair[1]] for pair in matches_cstar.pairs()])

          if len(I_even_bin_cstar) > 0:
            cc12_bin_cstar = np.corrcoef(I_even_bin_cstar, I_odd_bin_cstar)[0,1]
            avgI_bin_cstar = flex.mean(I_even_bin_cstar)
          I_even_cstar.extend(I_even_bin_cstar)
          I_odd_cstar.extend(I_odd_bin_cstar)
        except Exception:
          dummy = 1

        #caculation of <I> on the three lab axes (h,k,l) ---------------------------------
        I_even_h_bin = flex.double([I_even_h[pair[1]] for pair in matches_template.pairs()])
        I_odd_h_bin = flex.double([I_odd_h[pair[1]] for pair in matches_template.pairs()])
        I_even_k_bin = flex.double([I_even_k[pair[1]] for pair in matches_template.pairs()])
        I_odd_k_bin = flex.double([I_odd_k[pair[1]] for pair in matches_template.pairs()])
        I_even_l_bin = flex.double([I_even_l[pair[1]] for pair in matches_template.pairs()])
        I_odd_l_bin = flex.double([I_odd_l[pair[1]] for pair in matches_template.pairs()])

        if len(I_even_h_bin.select(I_even_h_bin > 0)) > 0:
          avgI_bin_h = flex.mean(I_even_h_bin.select(I_even_h_bin > 0))
        if len(I_even_k_bin.select(I_even_k_bin > 0)) > 0:
          avgI_bin_k = flex.mean(I_even_k_bin.select(I_even_k_bin > 0))
        if len(I_even_l_bin.select(I_even_l_bin > 0)) > 0:
          avgI_bin_l = flex.mean(I_even_l_bin.select(I_even_l_bin > 0))


      #collect CC1/2 and <I> on the crystal axes
      cc12_list.append([cc12_bin, cc12_bin_astar, cc12_bin_bstar, cc12_bin_cstar])
      avgI_star_list.append([avgI_bin, avgI_bin_astar, avgI_bin_bstar, avgI_bin_cstar])
      cc12_n_refl_list.append([n_refl_halfset_bin, len(I_even_bin_astar), len(I_even_bin_bstar), len(I_even_bin_cstar)])

      #collect <I> on the lab axes
      avgI_list.append([avgI_bin, avgI_bin_h, avgI_bin_k, avgI_bin_l])
      avgI_n_refl_list.append([len(I_bin), len(I_even_h_bin.select(I_even_h_bin > 0)), len(I_even_k_bin.select(I_even_k_bin > 0)), len(I_even_l_bin.select(I_even_l_bin > 0))])
      completeness = len(miller_indices_obs_bin)/len(miller_indices_bin)
      sum_refl_obs += len(miller_indices_obs_bin)
      sum_refl_complete += len(miller_indices_bin)

      #calculate CCiso
      cc_iso_bin = 0
      r_iso_bin = 0
      n_refl_cciso_bin = 0

      if flag_hklisoin_found:
        matches_iso = miller.match_multi_indices(
                  miller_indices_unique=miller_array_iso.indices(),
                  miller_indices=miller_indices_obs_bin)

        I_iso = flex.double([miller_array_iso.data()[pair[0]] for pair in matches_iso.pairs()])
        I_merge_match_iso = flex.double([I_bin[pair[1]] for pair in matches_iso.pairs()])
        sigI_merge_match_iso = flex.double([sigI_bin[pair[1]] for pair in matches_iso.pairs()])
        n_refl_cciso_bin = len(matches_iso.pairs())
        if len(matches_iso.pairs()) > 0 :
          cc_iso_bin = np.corrcoef(I_merge_match_iso, I_iso)[0,1]

          #calculate r_iso (need to find slope and intercept)
          from scipy import stats
          slope, intercept, r_value, p_value, std_err = stats.linregress(I_iso, I_merge_match_iso)
          r_iso_bin = (np.sum(((I_merge_match_iso - ((slope*I_iso)+intercept))/sigI_merge_match_iso)**2) * \
            math.sqrt(n_refl_cciso_bin/(max(1, n_refl_cciso_bin-1))))/ \
            flex.sum((I_merge_match_iso/sigI_merge_match_iso)**2)


      txt_out += '%02d %7.2f - %7.2f %5.1f %6.0f / %6.0f %7.2f %7.2f %7.2f %7.2f %6.0f %7.2f %6.0f %7.2f %6.0f %7.2f %6.0f %8.2f %10.2f' \
          %(i, binner_template_asu.bin_d_range(i)[0], binner_template_asu.bin_d_range(i)[1], completeness*100, \
          len(miller_indices_obs_bin), len(miller_indices_bin),\
          multiplicity_bin, r_meas_bin*100, r_split_bin*100, cc12_bin*100, n_refl_halfset_bin, cc_iso_bin*100, n_refl_cciso_bin, \
          cc_anom_bin_acentric, nrefl_anom_bin_acentric, \
          cc_anom_bin_centric, nrefl_anom_bin_acentric, \
          mean_i_over_sigi_bin, np.mean(I_bin))
      txt_out += '\n'


    #calculate CCiso
    cc_iso = 0
    n_refl_iso = 0
    if flag_hklisoin_found:
      matches_iso = miller.match_multi_indices(
                miller_indices_unique=miller_array_iso.indices(),
                miller_indices=miller_array_merge.indices())

      I_iso = flex.double([miller_array_iso.data()[pair[0]] for pair in matches_iso.pairs()])
      I_merge_match_iso = flex.double([miller_array_merge.data()[pair[1]] for pair in matches_iso.pairs()])
      sigI_merge_match_iso = flex.double([miller_array_merge.sigmas()[pair[1]] for pair in matches_iso.pairs()])

      if len(matches_iso.pairs()) > 0 :
        cc_iso = np.corrcoef(I_merge_match_iso, I_iso)[0,1]
        n_refl_iso = len(matches_iso.pairs())
        from scipy import stats
        slope, intercept, r_value, p_value, std_err = stats.linregress(I_iso, I_merge_match_iso)
        r_iso = (np.sum(((I_merge_match_iso - ((slope*I_iso)+intercept))/sigI_merge_match_iso)**2) * \
            math.sqrt(n_refl_cciso_bin/(n_refl_cciso_bin-1)))/ \
            flex.sum((I_merge_match_iso/sigI_merge_match_iso)**2)
      if iparams.flag_plot:
        plt.scatter(I_iso, I_merge_match_iso,s=10, marker='x', c='r')
        plt.title('CC=%.4g'%(cc_iso))
        plt.xlabel('I_ref')
        plt.ylabel('I_obs')
        plt.show()

    #calculate cc12
    i_even_filter_sel = (I_even > 0)
    cc12, cc12_astar, cc12_bstar, cc12_cstar, r_split, n_refl_half_total  = (0,0,0,0,0,0)
    try:
      I_even_filter = I_even.select(i_even_filter_sel)
      I_odd_filter = I_odd.select(i_even_filter_sel)
      cc12 = np.corrcoef(I_even_filter, I_odd_filter)[0,1]
      cc12_astar = np.corrcoef(I_even_astar, I_odd_astar)[0,1]
      cc12_bstar = np.corrcoef(I_even_bstar, I_odd_bstar)[0,1]
      cc12_cstar = np.corrcoef(I_even_cstar, I_odd_cstar)[0,1]
      r_split = (1/math.sqrt(2))*(flex.sum(flex.abs(I_even_filter - I_odd_filter))/(flex.sum(I_even_filter + I_odd_filter)*0.5))
      n_refl_half_total = len(I_even_filter)
    except Exception:
      dummy = 0

    #calculate Qmeas and Qw
    if sum_r_meas_w_btm > 0:
      r_meas_w = sum_r_meas_w_top/sum_r_meas_w_btm
    else:
      r_meas_w = float('Inf')

    if sum_r_meas_btm > 0:
      r_meas = sum_r_meas_top/sum_r_meas_btm
    else:
      r_meas = float('Inf')

    txt_out += '--------------------------------------------------------------------------------------------------------------------------------------------------\n'
    txt_out += '        TOTAL        %5.1f %6.0f / %6.0f %7.2f %7.2f %7.2f %7.2f %6.0f %7.2f %6.0f %7.2f %6.0f %7.2f %6.0f %8.2f %10.2f\n' \
    %((sum_refl_obs/sum_refl_complete)*100, sum_refl_obs, \
     sum_refl_complete, n_refl_obs_total/sum_refl_obs, \
     r_meas*100, r_split*100, cc12*100, len(I_even.select(i_even_filter_sel)), \
     cc_iso*100, n_refl_iso, \
     cc_anom_acentric, nrefl_anom_acentric, \
     cc_anom_centric, nrefl_anom_centric, \
     np.mean(miller_array_merge.data()/miller_array_merge.sigmas()), np.mean(miller_array_merge.data()))
    txt_out += '--------------------------------------------------------------------------------------------------------------------------------------------------\n'
    txt_out += '\n'

    #output CC1/2 on the three crystal axes
    txt_out += 'Summary of CC1/2 on three crystal axes\n'
    txt_out += 'Bin Resolution Range                CC1/2                                <I>                               N_refl           \n'
    txt_out += '                        All      a*      b*      c*  |   All         a*        b*        c*      | All      a*      b*     c* \n'
    txt_out += '-------------------------------------------------------------------------------------------------------------------------------\n'
    cn12_n_sum, cc12_n_astar_sum, cc12_n_bstar_sum, cc12_n_cstar_sum = (0,0,0,0)
    for i in range(1,iparams.n_bins+1):
      i_binner = (binner_template_asu_indices == i)
      _cc12, _cc12_astar, _cc12_bstar, _cc12_cstar = cc12_list[i-1]
      _cc12_n, _cc12_n_astar, _cc12_n_bstar, _cc12_n_cstar = cc12_n_refl_list[i-1]
      _avgI, _avgI_astar, _avgI_bstar, _avgI_cstar= avgI_star_list[i-1]
      cn12_n_sum += _cc12_n
      cc12_n_astar_sum += _cc12_n_astar
      cc12_n_bstar_sum += _cc12_n_bstar
      cc12_n_cstar_sum += _cc12_n_cstar
      txt_out += '%02d %7.2f - %7.2f %7.2f %7.2f %7.2f %7.2f %10.1f %10.1f %10.1f %10.1f %6.0f %6.0f %6.0f %6.0f\n' \
          %(i, binner_template_asu.bin_d_range(i)[0], binner_template_asu.bin_d_range(i)[1], \
          _cc12*100, _cc12_astar*100, _cc12_bstar*100, _cc12_cstar*100, \
          _avgI, _avgI_astar, _avgI_bstar, _avgI_cstar, \
          _cc12_n, _cc12_n_astar, _cc12_n_bstar, _cc12_n_cstar)

    txt_out += '-------------------------------------------------------------------------------------------------------------------------------\n'
    txt_out += '        TOTAL        %7.2f %7.2f %7.2f %7.2f %10.1f %10.1f %10.1f %10.1f %6.0f %6.0f %6.0f %6.0f\n' \
          %(cc12*100, cc12_astar*100, cc12_bstar*100, cc12_cstar*100, \
          np.mean(miller_array_merge.data()), np.mean(I_even_astar), np.mean(I_even_bstar), np.mean(I_even_cstar), \
          cn12_n_sum, cc12_n_astar_sum, cc12_n_bstar_sum, cc12_n_cstar_sum)
    txt_out += '-------------------------------------------------------------------------------------------------------------------------------\n'
    txt_out += '\n'
    """
    #output <I> on the three lab coordinates
    txt_out += 'Summary of <I> on the three lab coordinates\n'
    txt_out += 'Bin Resolution Range                  <I>                             N_refl        \n'
    txt_out += '                        All       h=0       k=0      l=0  | All    h=0    k=0    l=0\n'
    txt_out += '------------------------------------------------------------------------------------\n'
    avgI_n_sum, avgI_n_h_sum, avgI_n_k_sum, avgI_n_l_sum = (0,0,0,0)
    for i in range(1,iparams.n_bins+1):
      i_binner = (binner_template_asu_indices == i)
      _avgI_n, _avgI_n_h, _avgI_n_k, _avgI_n_l = avgI_n_refl_list[i-1]
      _avgI, _avgI_h, _avgI_k, _avgI_l = avgI_list[i-1]
      avgI_n_sum += _avgI_n
      avgI_n_h_sum += _avgI_n_h
      avgI_n_k_sum += _avgI_n_k
      avgI_n_l_sum += _avgI_n_l
      txt_out += '%02d %7.2f - %7.2f %8.2f %8.2f %8.2f %8.2f %6.0f %6.0f %6.0f %6.0f\n' \
          %(i, binner_template_asu.bin_d_range(i)[0], binner_template_asu.bin_d_range(i)[1], \
          _avgI, _avgI_h, _avgI_k, _avgI_l,
          _avgI_n, _avgI_n_h, _avgI_n_k, _avgI_n_l)

    txt_out += '------------------------------------------------------------------------------------\n'
    txt_out += '        TOTAL        %8.2f %8.2f %8.2f %8.2f %6.0f %6.0f %6.0f %6.0f\n' \
          %(np.mean(miller_array_merge.data()), np.mean(I_even_h.select(I_even_h>0)), np.mean(I_even_k.select(I_even_k>0)), np.mean(I_even_l.select(I_even_l>0)), \
          avgI_n_sum, avgI_n_h_sum, avgI_n_k_sum, avgI_n_l_sum)
    txt_out += '------------------------------------------------------------------------------------\n'
    txt_out += '\n'
    """

    return miller_array_merge, txt_out

  def plot_stats(self, results, iparams):
    #retrieve stats from results and plot them
    G_frame = flex.double()
    B_frame = flex.double()
    rotx_frame = flex.double()
    roty_frame = flex.double()
    ry_frame = flex.double()
    rz_frame = flex.double()
    re_frame = flex.double()
    uc_a_frame = flex.double()
    uc_b_frame = flex.double()
    uc_c_frame = flex.double()
    uc_alpha_frame = flex.double()
    uc_beta_frame = flex.double()
    uc_gamma_frame = flex.double()
    SE_all = flex.double()
    R_sq_all = flex.double()
    cc_all = flex.double()

    for pres in results:
      G_frame.append(pres.G)
      B_frame.append(pres.B)
      rotx_frame.append(pres.rotx*180/math.pi)
      roty_frame.append(pres.roty*180/math.pi)
      ry_frame.append(pres.ry)
      rz_frame.append(pres.rz)
      re_frame.append(pres.re)
      uc_a_frame.append(pres.uc_params[0])
      uc_b_frame.append(pres.uc_params[1])
      uc_c_frame.append(pres.uc_params[2])
      uc_alpha_frame.append(pres.uc_params[3])
      uc_beta_frame.append(pres.uc_params[4])
      uc_gamma_frame.append(pres.uc_params[5])
      SE_all.append(pres.SE)
      R_sq_all.append(pres.R_sq)
      cc_all.append(pres.CC_final)


    if iparams.flag_plot:
      plt.subplot(121)
      x = SE_all.as_numpy_array()
      mu = np.mean(x)
      med = np.median(x)
      sigma = np.std(x)
      num_bins = 10
      n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
      y = mlab.normpdf(bins, mu, sigma)
      plt.plot(bins, y, 'r--')
      plt.ylabel('Frequencies')
      plt.title('SE distribution\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))

      plt.subplot(122)
      x = cc_all.as_numpy_array()
      mu = np.mean(x)
      med = np.median(x)
      sigma = np.std(x)
      num_bins = 10
      n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
      y = mlab.normpdf(bins, mu, sigma)
      plt.plot(bins, y, 'r--')
      plt.ylabel('Frequencies')
      plt.title('CC distribution\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))
      plt.show()


      plt.subplot(241)
      x = G_frame.as_numpy_array()
      mu = np.mean(x)
      med = np.median(x)
      sigma = np.std(x)
      num_bins = 10
      n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
      y = mlab.normpdf(bins, mu, sigma)
      plt.plot(bins, y, 'r--')
      plt.ylabel('Frequencies')
      plt.title('G distribution\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))

      plt.subplot(242)
      x = B_frame.as_numpy_array()
      mu = np.mean(x)
      med = np.median(x)
      sigma = np.std(x)
      num_bins = 10
      n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
      y = mlab.normpdf(bins, mu, sigma)
      plt.plot(bins, y, 'r--')
      plt.ylabel('Frequencies')
      plt.title('B distribution\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))

      plt.subplot(243)
      x = rotx_frame.as_numpy_array()
      mu = np.mean(x)
      med = np.median(x)
      sigma = np.std(x)
      num_bins = 10
      n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
      y = mlab.normpdf(bins, mu, sigma)
      plt.plot(bins, y, 'r--')
      plt.ylabel('Frequencies')
      plt.title('Delta rot_x distribution\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))

      plt.subplot(244)
      x = roty_frame.as_numpy_array()
      mu = np.mean(x)
      med = np.median(x)
      sigma = np.std(x)
      num_bins = 10
      n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
      y = mlab.normpdf(bins, mu, sigma)
      plt.plot(bins, y, 'r--')
      plt.ylabel('Frequencies')
      plt.title('Delta rot_y distribution\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))

      plt.subplot(245)
      x = ry_frame.as_numpy_array()
      mu = np.mean(x)
      med = np.median(x)
      sigma = np.std(x)
      num_bins = 10
      n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
      y = mlab.normpdf(bins, mu, sigma)
      plt.plot(bins, y, 'r--')
      plt.ylabel('Frequencies')
      plt.title('ry distribution\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))

      plt.subplot(246)
      x = rz_frame.as_numpy_array()
      mu = np.mean(x)
      med = np.median(x)
      sigma = np.std(x)
      num_bins = 10
      n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
      y = mlab.normpdf(bins, mu, sigma)
      plt.plot(bins, y, 'r--')
      plt.ylabel('Frequencies')
      plt.title('rz distribution\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))

      plt.subplot(247)
      x = re_frame.as_numpy_array()
      mu = np.mean(x)
      med = np.median(x)
      sigma = np.std(x)
      num_bins = 10
      n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
      y = mlab.normpdf(bins, mu, sigma)
      plt.plot(bins, y, 'r--')
      plt.ylabel('Frequencies')
      plt.title('re distribution\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))

      plt.show()

      plt.subplot(231)
      x = uc_a_frame.as_numpy_array()
      mu = np.mean(x)
      med = np.median(x)
      sigma = np.std(x)
      num_bins = 10
      n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
      y = mlab.normpdf(bins, mu, sigma)
      plt.plot(bins, y, 'r--')
      plt.ylabel('Frequencies')
      plt.title('a distribution\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))

      plt.subplot(232)
      x = uc_b_frame.as_numpy_array()
      mu = np.mean(x)
      med = np.median(x)
      sigma = np.std(x)
      num_bins = 10
      n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
      y = mlab.normpdf(bins, mu, sigma)
      plt.plot(bins, y, 'r--')
      plt.ylabel('Frequencies')
      plt.title('b distribution\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))

      plt.subplot(233)
      x = uc_c_frame.as_numpy_array()
      mu = np.mean(x)
      med = np.median(x)
      sigma = np.std(x)
      num_bins = 10
      n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
      y = mlab.normpdf(bins, mu, sigma)
      plt.plot(bins, y, 'r--')
      plt.ylabel('Frequencies')
      plt.title('c distribution\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))

      plt.subplot(234)
      x = uc_alpha_frame.as_numpy_array()
      mu = np.mean(x)
      med = np.median(x)
      sigma = np.std(x)
      num_bins = 10
      n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
      y = mlab.normpdf(bins, mu, sigma)
      plt.plot(bins, y, 'r--')
      plt.ylabel('Frequencies')
      plt.title('alpha distribution\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))

      plt.subplot(235)
      x = uc_beta_frame.as_numpy_array()
      mu = np.mean(x)
      med = np.median(x)
      sigma = np.std(x)
      num_bins = 10
      n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
      y = mlab.normpdf(bins, mu, sigma)
      plt.plot(bins, y, 'r--')
      plt.ylabel('Frequencies')
      plt.title('beta distribution\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))

      plt.subplot(236)
      x = uc_gamma_frame.as_numpy_array()
      mu = np.mean(x)
      med = np.median(x)
      sigma = np.std(x)
      num_bins = 10
      n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
      y = mlab.normpdf(bins, mu, sigma)
      plt.plot(bins, y, 'r--')
      plt.ylabel('Frequencies')
      plt.title('gamma distribution\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))
      plt.show()


class basis_handler(object):
  '''
  classdocs
  '''

  def __init__(self):
    '''
    Constructor
    '''

  def calc_direct_space_matrix(self, my_unit_cell, rotation_matrix):

    #calculate the conversion matrix (from fractional to cartesian coordinates
    frac2cart_matrix = my_unit_cell.orthogonalization_matrix()
    frac2cart_matrix = sqr(frac2cart_matrix)

    #calculate direct_space matrix
    direct_space_matrix = frac2cart_matrix.transpose()*rotation_matrix

    return direct_space_matrix


class svd_handler(object):
  '''
  Singular value decomposion
  Solve linear equations with best fit basis
  '''


  # Input: expects Nx3 matrix of points
  # Returns R,t
  # R = 3x3 rotation matrix
  # t = 3x1 column vector

  def __init__(self):
    '''
    Constructor
    '''

  def rigid_transform_3D(self, A, B):

      assert len(A) == len(B)

      N = A.shape[0]; # total points

      centroid_A = np.mean(A, axis=0)
      centroid_B = np.mean(B, axis=0)

      # centre the points
      AA = A - np.tile(centroid_A, (N, 1))
      BB = B - np.tile(centroid_B, (N, 1))

      # dot is matrix multiplication for array
      H = np.transpose(AA) * BB

      U, S, Vt = np.linalg.svd(H)

      R = Vt.T * U.T

      # special reflection case
      if np.linalg.det(R) < 0:
         #print "Reflection detected"
         Vt[2,:] *= -1
         R = Vt.T * U.T

      t = -R*centroid_A.T + centroid_B.T


      return R, t

