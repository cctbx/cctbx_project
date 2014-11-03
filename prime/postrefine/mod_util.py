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
    self.CONST_SIG_I_FACTOR = 5.0

  def calc_avg_I(self, group_no, miller_index, miller_indices_ori, I, sigI, G, B,
                     p_set, rs_set, wavelength_set, sin_theta_over_lambda_sq, SE, avg_mode,
                     iparams, pickle_filename_set):


    from mod_leastsqr import calc_full_refl
    I_full = calc_full_refl(I, sin_theta_over_lambda_sq,
                            G, B, p_set, rs_set, iparams.flag_volume_correction)

    txt_obs_out = 'Reflection: '+str(miller_index[0])+','+str(miller_index[1])+','+str(miller_index[2])+'\n'
    txt_obs_out += 'meanI    medI  sigI_est sigI_true delta_sigI   n_refl\n'

    median_I = np.median(I_full)
    mean_I = np.mean(I_full)
    try:
      std_I_est = math.sqrt(median_I)
    except Exception:
      print 'Error <I>:', median_I
      return None
    std_I = np.std(I_full)

    txt_obs_out += '%6.2f %6.2f %8.2f %8.2f %8.2f %8.0f\n'%(mean_I, median_I, std_I_est, std_I, abs(std_I_est - std_I), len(I))

    #filter out outliers
    if avg_mode == 'average':
      sigma_max = 99
    else:
      sigma_max = iparams.sigma_rejection

    txt_reject_out = ''
    if len(I_full) > 2:
      for i_rejection in range(iparams.n_rejection_cycle):
        median_I = np.median(I_full)
        mean_I = np.mean(I_full)
        try:
          std_I_est = math.sqrt(median_I)
        except Exception:
          print 'Error <I>:', median_I
          return None

        std_I = np.std(I_full)

        I_full_as_sigma = (I_full - median_I)/ std_I

        i_seq = flex.int([i for i in range(len(I_full))])
        i_sel_inv = (flex.abs(I_full_as_sigma) > sigma_max)
        i_seq_sel_inv = i_seq.select(i_sel_inv)
        for _i in i_seq_sel_inv:
          txt_reject_out += pickle_filename_set[_i] + '%3.0f %3.0f %3.0f %10.2f %10.2f\n'%(\
            miller_indices_ori[_i][0], miller_indices_ori[_i][1], miller_indices_ori[_i][2], I[_i], sigI[_i])

        i_sel = (flex.abs(I_full_as_sigma) <= sigma_max)
        I_full = I_full.select(i_sel)
        SE = SE.select(i_sel)
        p_set = p_set.select(i_sel)
        rs_set = rs_set.select(i_sel)
        wavelength_set = wavelength_set.select(i_sel)
        G = G.select(i_sel)
        B = B.select(i_sel)
        I = I.select(i_sel)
        sigI = sigI.select(i_sel)
        miller_indices_ori = miller_indices_ori.select(i_sel)

        txt_obs_out += '%6.2f %6.2f %8.2f %8.2f %8.2f %8.0f\n'%(mean_I, median_I, std_I_est, std_I, abs(std_I_est - std_I), len(I))

        if (len(I) <=3) or (std_I < std_I_est):
          break

    if len(I_full) == 0:
      return None
    #normalize the SE
    max_w = self.CONST_SE_MAX_WEIGHT
    min_w = math.sqrt(self.CONST_SE_MIN_WEIGHT)
    if len(SE) == 1 or ((flex.min(SE)-flex.max(SE)) == 0) or avg_mode == 'average':
      SE_norm = flex.double([min_w+((max_w - min_w)/2)]*len(SE))
      SE_std_norm = flex.double([1.0+((self.CONST_SIG_I_FACTOR - 1.0)/2)]*len(SE))
    else:
      m = (max_w - min_w)/(flex.min(SE)-flex.max(SE))
      b = max_w - (m*flex.min(SE))
      SE_norm = (m*SE) + b
      m_std = (self.CONST_SIG_I_FACTOR - 1.0)/(flex.min(SE)-flex.max(SE))
      b_std = self.CONST_SIG_I_FACTOR - (m_std*flex.min(SE))
      SE_std_norm = (m_std*SE) + b_std


    I_avg = flex.sum(SE_norm * I_full)/flex.sum(SE_norm)

    #test calculation of sigI
    sigI_full = flex.sqrt(I_full) * SE_std_norm
    sigI_avg = np.mean(sigI_full)

    #Rmeas, Rmeas_w, multiplicity
    multiplicity = len(I_full)
    if multiplicity == 1:
      r_meas_w_top = 0
      r_meas_w_btm = 0
      r_meas_top = 0
      r_meas_btm = 0
      r_meas = 0
      r_meas_w = 0
    else:
      n_obs = multiplicity
      r_meas_w_top = flex.sum(((I_full - I_avg)*SE_norm)**2)*math.sqrt(n_obs/(n_obs-1))
      r_meas_w_btm = flex.sum((I_full*SE_norm)**2)
      r_meas_top = flex.sum(flex.abs((I_full - I_avg)*SE_norm))*math.sqrt(n_obs/(n_obs-1))
      r_meas_btm = flex.sum(flex.abs(I_full*SE_norm))
      r_meas = r_meas_top/r_meas_btm
      r_meas_w = r_meas_w_top/r_meas_w_btm

    #for calculattion of cc1/2
    #sepearte the observations into two groups
    if multiplicity == 1:
      I_avg_even = 0
      I_avg_odd = 0
    else:
      i_even = range(0,len(I_full),2)
      i_odd = range(1,len(I_full),2)
      I_even = I_full.select(i_even)
      SE_norm_even = SE_norm.select(i_even)
      I_odd = I_full.select(i_odd)
      SE_norm_odd = SE_norm.select(i_odd)
      if len(i_even) > len(i_odd):
        I_odd.append(I_even[len(I_even)-1])
        SE_norm_odd.append(SE_norm_even[len(I_even)-1])

      if avg_mode == 'weighted' or avg_mode == 'final':
        I_avg_even = np.sum(SE_norm_even * I_even)/np.sum(SE_norm_even)
        I_avg_odd = np.sum(SE_norm_odd * I_odd)/np.sum(SE_norm_odd)
      elif avg_mode== 'average':
        I_avg_even = np.mean(I_even)
        I_avg_odd = np.mean(I_odd)


    #calculate mosaic spread
    mosaic_radian_set = 2 * (rs_set) * (wavelength_set)

    txt_obs_out += '    I_o        sigI_o    G      B     Eoc      rs    lambda rocking(deg) W     Wsig     I_full     sigI_full\n'
    for i_o, sigi_o, g, b, eoc, rs, wavelength, mosaic_radian, se_norm, se_std_norm, i_full, sigi_full  in \
      zip(I, sigI, G, B, p_set, rs_set, wavelength_set, mosaic_radian_set, SE_norm, SE_std_norm, I_full, sigI_full):
      txt_obs_out += '%10.2f %10.2f %6.2f %6.2f %6.2f %8.5f %8.5f %8.5f %6.2f %6.2f %10.2f %10.2f\n'%(\
        i_o, sigi_o, 1/g, b, eoc, rs, wavelength, mosaic_radian*180/math.pi, se_norm, se_std_norm, i_full, sigi_full)
    txt_obs_out += 'Merged I, sigI: %6.2f, %6.2f\n'%(I_avg, sigI_avg)
    txt_obs_out += 'Rmeas: %6.2f Qw: %6.2f\n'%(r_meas, r_meas_w)
    txt_obs_out += 'No. total observed: %4.0f No. after rejection: %4.0f\n'%(len(sin_theta_over_lambda_sq), len(I_full))
    txt_obs_out += 'List of rejected observations:\n'
    txt_obs_out += txt_reject_out

    return miller_index, I_avg, sigI_avg, (r_meas_w_top, r_meas_w_btm, r_meas_top, r_meas_btm, multiplicity), I_avg_even, I_avg_odd, txt_obs_out, txt_reject_out

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
    uc_std = flex.double([np.std(a_all), np.std(b_all), np.std(c_all), np.std(alpha_all), np.std(beta_all), np.std(gamma_all)])

    return uc_mean, uc_std

  def calc_mean_postref_parameters(self, results):
    G_all = flex.double()
    B_all = flex.double()
    rotx_all = flex.double()
    roty_all = flex.double()
    ry_all = flex.double()
    rz_all = flex.double()
    re_all = flex.double()
    r0_all = flex.double()
    for pres in results:
      if pres is not None:
        G_all.append(pres.G)
        B_all.append(pres.B)
        rotx_all.append(pres.rotx)
        roty_all.append(pres.roty)
        ry_all.append(pres.ry)
        rz_all.append(pres.rz)
        re_all.append(pres.re)
        r0_all.append(pres.spot_radius)

    pr_params_mean = flex.double([np.mean(G_all), np.mean(B_all),
                                  np.mean(flex.abs(ry_all)), np.mean(flex.abs(rz_all)),
                                  np.mean(re_all), np.mean(r0_all),
                                  np.mean(flex.abs(rotx_all)), np.mean(flex.abs(roty_all))])
    pr_params_std = flex.double([np.std(G_all), np.std(B_all),
                                  np.std(flex.abs(ry_all)), np.std(flex.abs(rz_all)),
                                  np.std(re_all), np.std(r0_all),
                                  np.std(flex.abs(rotx_all)), np.std(flex.abs(roty_all))])

    return pr_params_mean, pr_params_std

  def prepare_output(self, results, iparams, avg_mode):
    if avg_mode == 'average':
      cc_thres = 0
      std_G_filter = 99
    else:
      cc_thres = iparams.frame_accept_min_cc
      std_G_filter = iparams.sigma_rejection

    #calculate distribution of R (target post-refinement and x,y restraints).
    R_final_all = flex.double()
    R_xy_final_all = flex.double()
    for pres in results:
      if pres is not None:
        R_final_all.append(pres.R_final)
        R_xy_final_all.append(pres.R_xy_final)
    mean_R_final = np.median(R_final_all)
    mean_R_xy_final = np.median(R_xy_final_all)

    #calculate mean G for filtering
    G_stats = flex.double()
    for pres in results:
      if pres is not None:
        if pres.R_final <= (mean_R_final * 6):
          if pres.R_xy_final <= (mean_R_xy_final * 6):
            if pres.CC_final >= cc_thres:
              G_stats.extend(flex.double([pres.G]*len(pres.observations.data())))
    mean_G_stats = np.median(G_stats)
    std_G_stats = math.sqrt(mean_G_stats)

    #prepare data for merging
    miller_indices_all = flex.miller_index()
    miller_indices_ori_all = flex.miller_index()
    I_all = flex.double()
    sigI_all = flex.double()
    G_all = flex.double()
    B_all = flex.double()
    k_all = flex.double()
    p_all = flex.double()
    rs_all = flex.double()
    rh_all = flex.double()
    SE_all = flex.double()
    sin_sq_all = flex.double()
    wavelength_all = flex.double()
    cn_good_frame = 0
    cn_bad_frame_uc = 0
    cn_bad_frame_cc = 0
    cn_bad_R = 0
    cn_bad_R_xy = 0
    cn_bad_frame_G = 0
    R_init_all = flex.double()
    R_final_all = flex.double()
    R_xy_init_all = flex.double()
    R_xy_final_all = flex.double()
    pickle_filename_all = []
    filtered_results = []
    i_seq = flex.int()
    for pres in results:
      if pres is not None:
        pickle_filepath = pres.pickle_filename.split('/')
        img_filename = pickle_filepath[len(pickle_filepath)-1]
        if pres.R_final <= (mean_R_final * 6):
          if pres.R_xy_final <= (mean_R_xy_final * 6):
            if pres.CC_final >= cc_thres:
              #check unit-cell
              if (abs(pres.uc_params[0]-iparams.target_unit_cell.parameters()[0]) <= (iparams.merge.uc_tolerance*iparams.target_unit_cell.parameters()[0]/100) \
                  and abs(pres.uc_params[1]-iparams.target_unit_cell.parameters()[1]) <= (iparams.merge.uc_tolerance*iparams.target_unit_cell.parameters()[1]/100) \
                  and abs(pres.uc_params[2]-iparams.target_unit_cell.parameters()[2]) <= (iparams.merge.uc_tolerance*iparams.target_unit_cell.parameters()[2]/100) \
                  and abs(pres.uc_params[3]-iparams.target_unit_cell.parameters()[3]) <= (iparams.merge.uc_tolerance*iparams.target_unit_cell.parameters()[3]/100) \
                  and abs(pres.uc_params[4]-iparams.target_unit_cell.parameters()[4]) <= (iparams.merge.uc_tolerance*iparams.target_unit_cell.parameters()[4]/100) \
                  and abs(pres.uc_params[5]-iparams.target_unit_cell.parameters()[5]) <= (iparams.merge.uc_tolerance*iparams.target_unit_cell.parameters()[5]/100)):

                if abs((pres.G - mean_G_stats)/std_G_stats) < std_G_filter:
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
                  print pres.frame_no, img_filename, ' merged'
                else:
                  print pres.frame_no, img_filename, ' discarded - G too high/low (G=%5.2f%%)'%(pres.G)
                  cn_bad_frame_G += 1
              else:
                print pres.frame_no, img_filename, ' discarded - unit-cell exceeds the limits (%6.2f %6.2f %6.2f %5.2f %5.2f %5.2f)'%(pres.uc_params[0], pres.uc_params[1], pres.uc_params[2], pres.uc_params[3], pres.uc_params[4], pres.uc_params[5])
                cn_bad_frame_uc += 1
            else:
              print pres.frame_no, img_filename, ' discarded - C.C. too low (C.C.=%5.2f%%)'%(pres.CC_final*100)
              cn_bad_frame_cc += 1
          else:
            print pres.frame_no, img_filename, ' discarded - target R (x,y) too high (R=%5.2f%%)'%(pres.R_xy_final*100)
            cn_bad_R_xy +=1
        else:
          print pres.frame_no, img_filename, ' discarded - target R too high (R=%5.2f%%)'%(pres.R_final*100)
          cn_bad_R +=1


    #plot stats
    self.plot_stats(filtered_results, iparams)


    #calculate average unit cell
    uc_mean, uc_std = self.calc_mean_unit_cell(filtered_results)
    unit_cell_mean = unit_cell((uc_mean[0], uc_mean[1], uc_mean[2], uc_mean[3], uc_mean[4], uc_mean[5]))

    pr_params_mean, pr_params_std = self.calc_mean_postref_parameters(filtered_results)

    #from all observations merge them
    crystal_symmetry = crystal.symmetry(
        unit_cell=(uc_mean[0], uc_mean[1], uc_mean[2], uc_mean[3], uc_mean[4], uc_mean[5]),
        space_group_symbol=iparams.target_space_group)
    miller_set_all=miller.set(
                crystal_symmetry=crystal_symmetry,
                indices=miller_indices_all,
                anomalous_flag=iparams.target_anomalous_flag)
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
    txt_out += ' No. bad G frames) :       %12.0f (median=%6.2f std.=%6.2f)\n'%(cn_bad_frame_G, mean_G_stats, std_G_stats)
    txt_out += ' No. bad unit cell frames: %12.0f\n'%(cn_bad_frame_uc)
    txt_out += ' No. bad target R:         %12.0f\n'%(cn_bad_R)
    txt_out += ' No. bad target R(x,y) :   %12.0f\n'%(cn_bad_R_xy)
    txt_out += ' No. observations:         %12.0f\n'%(len(I_obs_all_sort))
    txt_out += 'Mean target value (BEFORE: Sum Mean (Std.))\n'
    txt_out += ' post-refinement:          %12.2f %12.2f (%9.2f)\n'%(np.sum(R_init_all), np.mean(R_init_all), np.std(R_init_all))
    txt_out += ' (x,y) restraints:         %12.2f %12.2f (%9.2f)\n'%(np.sum(R_xy_init_all), np.mean(R_xy_init_all), np.std(R_xy_init_all))
    txt_out += 'Mean target value (AFTER: Sum Mean (Std.))\n'
    txt_out += ' post-refinement:          %12.2f %12.2f (%9.2f)\n'%(np.sum(R_final_all), np.mean(R_final_all), np.std(R_final_all))
    txt_out += ' (x,y) restraints:         %12.2f %12.2f (%9.2f)\n'%(np.sum(R_xy_final_all), np.mean(R_xy_final_all), np.std(R_xy_final_all))
    txt_out += ' G:                        %12.2f (%7.2f)\n'%(pr_params_mean[0], pr_params_std[0])
    txt_out += ' B:                        %12.2f (%7.2f)\n'%(pr_params_mean[1], pr_params_std[1])
    txt_out += ' Rot.x:                    %12.2f (%7.2f)\n'%(pr_params_mean[6]*180/math.pi, pr_params_std[6]*180/math.pi)
    txt_out += ' Rot.y:                    %12.2f (%7.2f)\n'%(pr_params_mean[7]*180/math.pi, pr_params_std[7]*180/math.pi)
    txt_out += ' gamma_y:                  %12.5f (%7.5f)\n'%(pr_params_mean[2], pr_params_std[2])
    txt_out += ' gamma_z:                  %12.5f (%7.5f)\n'%(pr_params_mean[3], pr_params_std[3])
    txt_out += ' gamma_e:                  %12.5f (%7.5f)\n'%(pr_params_mean[4], pr_params_std[4])
    txt_out += ' unit cell\n'
    txt_out += '   a:                        %12.5f (%7.5f)\n'%(uc_mean[0], uc_std[0])
    txt_out += '   b:                        %12.5f (%7.5f)\n'%(uc_mean[1], uc_std[1])
    txt_out += '   c:                        %12.5f (%7.5f)\n'%(uc_mean[2], uc_std[2])
    txt_out += '   alpha:                    %12.5f (%7.5f)\n'%(uc_mean[3], uc_std[3])
    txt_out += '   beta:                     %12.5f (%7.5f)\n'%(uc_mean[4], uc_std[4])
    txt_out += '   gamma:                    %12.5f (%7.5f)\n'%(uc_mean[5], uc_std[5])
    txt_out += '* (standard deviation)\n'

    return cn_group, group_id_list, miller_indices_all_sort, miller_indices_ori_all_sort, \
           I_obs_all_sort, sigI_obs_all_sort,G_all_sort, B_all_sort, \
           p_all_sort, rs_all_sort, wavelength_all_sort, sin_sq_all_sort, SE_all_sort, uc_mean, \
           np.median(wavelength_all), pickle_filename_all_sort, txt_out


  def write_output(self, miller_indices_merge, I_merge, sigI_merge, stat_all,
                   I_even, I_odd, iparams, uc_mean, wavelength_mean,
                   output_mtz_file_prefix, avg_mode):

    #output mtz file and report binning stat
    crystal_symmetry = crystal.symmetry(
        unit_cell=(uc_mean[0], uc_mean[1], uc_mean[2], uc_mean[3], uc_mean[4], uc_mean[5]),
        space_group_symbol=iparams.target_space_group)
    miller_set_merge=miller.set(
              crystal_symmetry=crystal_symmetry,
              indices=miller_indices_merge,
              anomalous_flag=iparams.target_anomalous_flag)
    miller_array_merge = miller_set_merge.array(data=I_merge,
              sigmas=sigI_merge).set_observation_type_xray_intensity()

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


    #remove outliers
    if avg_mode == 'average':
      sigma_filter = 99
    else:
      sigma_filter = iparams.sigma_rejection
    binner_merge = miller_array_merge.setup_binner(n_bins=iparams.n_bins)
    binner_merge_indices = binner_merge.bin_indices()
    miller_indices_merge_filter = flex.miller_index()
    I_merge_filter = flex.double()
    sigI_merge_filter = flex.double()
    I_even_filter = flex.double()
    I_odd_filter = flex.double()
    stat_filter = []
    i_seq = flex.int([j for j in range(len(binner_merge_indices))])
    for i in range(1,iparams.n_bins+1):
      i_binner = (binner_merge_indices == i)
      if len(miller_array_merge.data().select(i_binner)) > 0:
        I_obs_bin = miller_array_merge.data().select(i_binner)
        sigI_obs_bin = miller_array_merge.sigmas().select(i_binner)
        miller_indices_bin = miller_array_merge.indices().select(i_binner)
        stat_bin = [stat_all[j] for j in i_seq.select(i_binner)]
        I_even_bin = I_even.select(i_binner)
        I_odd_bin = I_odd.select(i_binner)

        i_filter = flex.abs((I_obs_bin - np.median(I_obs_bin))/np.median(I_obs_bin)) < sigma_filter
        I_obs_bin_filter = I_obs_bin.select(i_filter)
        sigI_obs_bin_filter = sigI_obs_bin.select(i_filter)
        miller_indices_bin_filter = miller_indices_bin.select(i_filter)
        i_seq_bin = flex.int([j for j in range(len(i_filter))])
        stat_bin_filter = [stat_bin[j] for j in i_seq_bin.select(i_filter)]
        I_even_bin_filter = I_even_bin.select(i_filter)
        I_odd_bin_filter = I_odd_bin.select(i_filter)

        for i_obs, sigi_obs, miller_index, stat, i_even, i_odd in zip(I_obs_bin_filter, sigI_obs_bin_filter,
            miller_indices_bin_filter, stat_bin_filter, I_even_bin_filter, I_odd_bin_filter):
          I_merge_filter.append(i_obs)
          sigI_merge_filter.append(sigi_obs)
          miller_indices_merge_filter.append(miller_index)
          stat_filter.append(stat)
          I_even_filter.append(i_even)
          I_odd_filter.append(i_odd)

    miller_array_merge = miller_array_merge.customized_copy(indices=miller_indices_merge_filter,
                                                            data=I_merge_filter,
                                                            sigmas=sigI_merge_filter)
    stat_all = stat_filter
    I_even = I_even_filter
    I_odd = I_odd_filter

    if output_mtz_file_prefix != '':
      #write as mtz file
      mtz_dataset_merge = miller_array_merge.as_mtz_dataset(column_root_label="IOBS")
      mtz_dataset_merge.mtz_object().write(file_name=output_mtz_file_prefix+'_merge.mtz')
      #write as cns file
      f_cns = open(output_mtz_file_prefix+'_merge.hkl', 'w')
      miller_array_merge.export_as_cns_hkl(file_object=f_cns)
      f_cns.close()

    #report binning stats
    miller_array_template_asu = miller_array_merge.complete_set().resolution_filter(
      d_min=iparams.merge.d_min, d_max=iparams.merge.d_max)
    binner_template_asu = miller_array_template_asu.setup_binner(n_bins=iparams.n_bins)
    binner_template_asu_indices = binner_template_asu.bin_indices()
    avg_I_by_bin = flex.double()
    one_dsqr_by_bin = flex.double()

    csv_out = ""
    csv_out +='Bin, Low, High, Completeness, <N_obs>, Rmeas, Qw, CC1/2, N_ind, CCiso, N_ind, <I/sigI>\n'
    txt_out = '\n'
    txt_out += 'Summary for '+output_mtz_file_prefix+'_merge.mtz\n'
    txt_out += 'Bin Resolution Range     Completeness      <N_obs>  |CC1/2   N_ind |CCiso   N_ind| <I/sigI>   <I>\n'
    txt_out += '---------------------------------------------------------------------------------------------------\n'
    sum_r_meas_w_top = 0
    sum_r_meas_w_btm = 0
    sum_r_meas_top = 0
    sum_r_meas_btm = 0
    sum_refl_obs = 0
    sum_refl_complete = 0

    n_refl_obs_total = 0
    for i in range(1,iparams.n_bins+1):
      i_binner = (binner_template_asu_indices == i)
      miller_indices_bin = miller_array_template_asu.indices().select(i_binner)

      matches_template = miller.match_multi_indices(
                  miller_indices_unique=miller_indices_bin,
                  miller_indices=miller_array_merge.indices())

      I_bin = flex.double([miller_array_merge.data()[pair[1]] for pair in matches_template.pairs()])
      sigI_bin = flex.double([miller_array_merge.sigmas()[pair[1]] for pair in matches_template.pairs()])
      miller_indices_obs_bin = flex.miller_index([miller_array_merge.indices()[pair[1]] for pair in matches_template.pairs()])

      avg_I_by_bin.append(np.mean(I_bin))
      one_dsqr_by_bin.append(1/binner_template_asu.bin_d_range(i)[1]**2)
      if len(I_bin) == 0:
        mean_i_over_sigi_bin = 0
        multiplicity_bin = 0
        r_meas_w_bin = 0
        r_meas_bin = 0
        cc12 = 0
      else:
        mean_i_over_sigi_bin = flex.mean(I_bin/sigI_bin)
        stat_bin = [stat_all[pair[1]] for pair in matches_template.pairs()]
        sum_r_meas_w_top_bin = 0
        sum_r_meas_w_btm_bin = 0
        sum_r_meas_top_bin = 0
        sum_r_meas_btm_bin = 0
        sum_mul_bin = 0
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

        I_even_bin = flex.double([I_even[pair[1]] for pair in matches_template.pairs()])
        I_odd_bin = flex.double([I_odd[pair[1]] for pair in matches_template.pairs()])
        #for cc1/2, use only non-zero I (zero when there is only one observation)
        i_even_filter_sel = (I_even_bin > 0)
        n_refl_cc12_bin = len(I_even_bin.select(i_even_filter_sel))
        cc12_bin = 0
        if n_refl_cc12_bin > 0:
          cc12_bin = np.corrcoef(I_even_bin.select(i_even_filter_sel), I_odd_bin.select(i_even_filter_sel))[0,1]

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
            math.sqrt(n_refl_cciso_bin/(n_refl_cciso_bin-1)))/ \
            flex.sum((I_merge_match_iso/sigI_merge_match_iso)**2)


      txt_out += '%02d %7.2f - %7.2f %5.1f %6.0f / %6.0f %7.2f %7.2f %6.0f %7.2f %6.0f %8.2f %10.2f' \
          %(i, binner_template_asu.bin_d_range(i)[0], binner_template_asu.bin_d_range(i)[1], completeness*100, \
          len(miller_indices_obs_bin), len(miller_indices_bin),\
          multiplicity_bin, cc12_bin*100, n_refl_cc12_bin, cc_iso_bin*100, n_refl_cciso_bin, mean_i_over_sigi_bin, np.mean(I_bin))
      txt_out += '\n'


    if iparams.flag_plot:
      x_axis = one_dsqr_by_bin
      plt.plot(x_axis, avg_I_by_bin, linestyle='-', linewidth=2.0, c='b', label='<I>')
      plt.title('Intensity plot by resolutions')
      plt.xlabel('1/d^2')
      plt.ylabel('<I>')
      plt.show()
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
    cc12 = np.corrcoef(I_even.select(i_even_filter_sel), I_odd.select(i_even_filter_sel))[0,1]

    #calculate Qmeas and Qw
    if sum_r_meas_w_btm > 0:
      r_meas_w = sum_r_meas_w_top/sum_r_meas_w_btm
    else:
      r_meas_w = float('Inf')

    if sum_r_meas_btm > 0:
      r_meas = sum_r_meas_top/sum_r_meas_btm
    else:
      r_meas = float('Inf')

    txt_out += '---------------------------------------------------------------------------------------------------\n'
    txt_out += '        TOTAL        %5.1f %6.0f / %6.0f %7.2f %7.2f %6.0f %7.2f %6.0f %8.2f %10.2f\n' \
    %((sum_refl_obs/sum_refl_complete)*100, sum_refl_obs, \
     sum_refl_complete, n_refl_obs_total/sum_refl_obs, \
     cc12*100, len(I_even.select(i_even_filter_sel)), cc_iso*100, \
     n_refl_iso, np.mean(miller_array_merge.data()/miller_array_merge.sigmas()), np.mean(miller_array_merge.data()))
    txt_out += '---------------------------------------------------------------------------------------------------\n'

    return miller_array_merge, txt_out, csv_out

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
