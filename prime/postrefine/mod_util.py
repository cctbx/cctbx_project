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

  def calc_avg_I(self, group_no, miller_index, I, sigI, G, B,
                     p_set, rs_set, sin_theta_over_lambda_sq, SE, avg_mode, iparams):

    from mod_leastsqr import calc_full_refl
    I_full = calc_full_refl(I, sin_theta_over_lambda_sq,
                            G, B, p_set, rs_set, iparams.flag_volume_correction)
    sigI_full = calc_full_refl(sigI, sin_theta_over_lambda_sq,
                               G, B, p_set, rs_set, iparams.flag_volume_correction)

    #filter outliers iteratively (Read, 1999)
    sigma_max = 3
    n_iters = len(I_full)*2
    for i_iter in range(n_iters):
      if len(I_full) > 2:
        import random
        i_seq = flex.int([i for i in range(len(I_full))])
        i_r = random.randint(0,len(I_full)-1)
        i_seq_sel = (i_seq != i_r)
        i_full_sel = I_full.select(i_seq_sel)
        if abs((I_full[i_r] - np.median(i_full_sel))/np.std(i_full_sel)) > sigma_max:
          #discard this observation
          I_full = I_full.select(i_seq_sel)
          sigI_full = sigI_full.select(i_seq_sel)
          SE = SE.select(i_seq_sel)

      else:
        break

    #normalize the SE
    max_w = 1.0
    min_w = 0.6
    if len(SE) == 1 or ((flex.min(SE)-flex.max(SE)) == 0):
      SE_norm = flex.double([min_w+((max_w - min_w)/2)]*len(SE))
    else:
      m = (max_w - min_w)/(flex.min(SE)-flex.max(SE))
      b = max_w - (m*flex.min(SE))
      SE_norm = (m*SE) + b

    if avg_mode == 'weighted':
      I_avg = np.sum(SE_norm * I_full)/np.sum(SE_norm)
      sigI_avg = np.sum(SE_norm * sigI_full)/np.sum(SE_norm)
    elif avg_mode== 'average':
      I_avg = np.mean(I_full)
      sigI_avg = np.mean(sigI_full)

    #Rmeas, Rmeas_w, multiplicity
    multiplicity = len(I_full)
    if multiplicity == 1:
      r_meas_w_top = 0
      r_meas_w_btm = 0
      r_meas_top = 0
      r_meas_btm = 0
    else:
      n_obs = multiplicity
      r_meas_w_top = flex.sum(((I_full - I_avg)*SE_norm)**2)*math.sqrt(n_obs/(n_obs-1))
      r_meas_w_btm = flex.sum((I_full*SE_norm)**2)
      r_meas_top = flex.sum((I_full - I_avg)**2)*math.sqrt(n_obs/(n_obs-1))
      r_meas_btm = flex.sum((I_full)**2)


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

      if avg_mode == 'weighted':
        I_avg_even = np.sum(SE_norm_even * I_even)/np.sum(SE_norm_even)
        I_avg_odd = np.sum(SE_norm_odd * I_odd)/np.sum(SE_norm_odd)
      elif avg_mode== 'average':
        I_avg_even = np.mean(I_even)
        I_avg_odd = np.mean(I_odd)

    return miller_index, I_avg, sigI_avg, (r_meas_w_top, r_meas_w_btm, r_meas_top, r_meas_btm, multiplicity), I_avg_even, I_avg_odd

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

    uc_mean = flex.double([np.median(a_all), np.median(b_all), np.median(c_all), np.median(alpha_all), np.median(beta_all), np.median(gamma_all)])
    uc_std = flex.double([np.std(a_all), np.std(b_all), np.std(c_all), np.std(alpha_all), np.std(beta_all), np.std(gamma_all)])

    return uc_mean, uc_std

  def calc_mean_postref_parameters(self, results):
    G_all = flex.double()
    B_all = flex.double()
    ry_all = flex.double()
    rz_all = flex.double()
    re_all = flex.double()
    for pres in results:
      if pres is not None:
        G_all.append(pres.G)
        B_all.append(pres.B)
        ry_all.append(pres.ry)
        rz_all.append(pres.rz)
        re_all.append(pres.re)

    pr_params_mean = flex.double([np.median(G_all), np.median(B_all),
                                  np.median(ry_all), np.median(rz_all), np.median(re_all)])
    pr_params_std = flex.double([np.std(G_all), np.std(B_all),
                                  np.std(ry_all), np.std(rz_all), np.std(re_all)])

    return pr_params_mean, pr_params_std

  def prepare_output(self, results, iparams, avg_mode):
    if avg_mode == 'average':
      cc_thres = 0
    else:
      cc_thres = iparams.frame_accept_min_cc

    #calculate distribution of R (target post-refinement and x,y restraints).
    R_final_all = flex.double()
    R_xy_final_all = flex.double()
    for pres in results:
      if pres is not None:
        R_final_all.append(pres.R_final)
        R_xy_final_all.append(pres.R_xy_final)
    mean_R_final = np.median(R_final_all)
    mean_R_xy_final = np.median(R_xy_final_all)


    #prepare data for merging
    miller_indices_all = flex.miller_index()
    I_all = flex.double()
    sigI_all = flex.double()
    G_all = flex.double()
    B_all = flex.double()
    k_all = flex.double()
    p_all = flex.double()
    rs_all = flex.double()
    SE_all = flex.double()
    sin_sq_all = flex.double()
    wavelength_all = flex.double()
    cn_good_frame = 0
    cn_bad_frame_uc = 0
    cn_bad_frame_cc = 0
    cn_bad_R = 0
    cn_bad_R_xy = 0
    R_init_all = flex.double()
    R_final_all = flex.double()
    R_xy_init_all = flex.double()
    R_xy_final_all = flex.double()
    filtered_results = []
    for pres in results:
      if pres is not None:
        img_filename = pres.pickle_filename
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

                cn_good_frame += 1
                sin_theta_over_lambda_sq = pres.observations.two_theta(wavelength=pres.wavelength).sin_theta_over_lambda_sq().data()
                filtered_results.append(pres)
                R_init_all.append(pres.R_init)
                R_final_all.append(pres.R_final)
                R_xy_init_all.append(pres.R_xy_init)
                R_xy_final_all.append(pres.R_xy_final)

                miller_indices_all.extend(pres.observations.indices())
                I_all.extend(pres.observations.data())
                sigI_all.extend(pres.observations.sigmas())
                G_all.extend(flex.double([pres.G]*len(pres.observations.data())))
                B_all.extend(flex.double([pres.B]*len(pres.observations.data())))
                p_all.extend(pres.partiality)
                rs_all.extend(pres.rs_set)
                sin_sq_all.extend(sin_theta_over_lambda_sq)
                SE_all.extend(flex.double([pres.SE]*len(pres.observations.data())))
                wavelength_all.extend(flex.double([pres.wavelength]*len(pres.observations.data())))

                print pres.frame_no, img_filename, ' merged'
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
    I_obs_all_sort = miller_array_all.data().select(perm)
    sigI_obs_all_sort = miller_array_all.sigmas().select(perm)
    G_all_sort = G_all.select(perm)
    B_all_sort = B_all.select(perm)
    p_all_sort = p_all.select(perm)
    rs_all_sort = rs_all.select(perm)
    sin_sq_all_sort = sin_sq_all.select(perm)
    SE_all_sort = SE_all.select(perm)

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
    txt_out += ' No. bad unit cell frames: %12.0f\n'%(cn_bad_frame_uc)
    txt_out += ' No. bad target R:         %12.0f\n'%(cn_bad_R)
    txt_out += ' No. bad target R(x,y) :   %12.0f\n'%(cn_bad_R_xy)
    txt_out += ' No. observations:         %12.0f\n'%(len(I_obs_all_sort))
    txt_out += 'Mean target value (BEFORE)\n'
    txt_out += ' post-refinement:          %12.2f (%7.2f)\n'%(np.median(R_init_all),
                                                                     np.std(R_init_all))
    txt_out += ' (x,y) restraints:         %12.2f (%7.2f)\n'%(np.median(R_xy_init_all),
                                                                     np.std(R_xy_init_all))
    txt_out += 'Mean target value (AFTER)\n'
    txt_out += ' post-refinement:          %12.2f (%7.2f)\n'%(np.median(R_final_all),
                                                                     np.std(R_final_all))
    txt_out += ' (x,y) restraints:         %12.2f (%7.2f)\n'%(np.median(R_xy_final_all),
                                                                     np.std(R_xy_final_all))
    txt_out += ' G:                        %12.2f (%7.2f)\n'%(pr_params_mean[0], pr_params_std[0])
    txt_out += ' B:                        %12.2f (%7.2f)\n'%(pr_params_mean[1], pr_params_std[1])
    txt_out += ' gamma_y:                  %12.5f (%7.2f)\n'%(pr_params_mean[2], pr_params_std[2])
    txt_out += ' gamma_z:                  %12.5f (%7.2f)\n'%(pr_params_mean[3], pr_params_std[3])
    txt_out += ' gamma_e:                  %12.5f (%7.2f)\n'%(pr_params_mean[4], pr_params_std[4])
    txt_out += ' unit cell:                %5.2f(%5.2f) %5.2f(%5.2f) %5.2f(%5.2f)\n' \
    %(uc_mean[0], uc_std[0], uc_mean[1], uc_std[1], uc_mean[2], uc_std[2])
    txt_out += '                           %5.2f(%5.2f) %5.2f(%5.2f) %5.2f(%5.2f)\n' \
    %(uc_mean[3], uc_std[3], uc_mean[4], uc_std[4], uc_mean[5], uc_std[5])
    txt_out += '* (standard deviation)\n'

    return cn_group, group_id_list, miller_indices_all_sort, \
           I_obs_all_sort, sigI_obs_all_sort,G_all_sort, B_all_sort, \
           p_all_sort, rs_all_sort, sin_sq_all_sort, SE_all_sort, uc_mean, \
           np.median(wavelength_all), txt_out


  def write_output(self, miller_indices_merge, I_merge, sigI_merge, stat_all,
                   I_even, I_odd, iparams, uc_mean, wavelength_mean, output_mtz_file_prefix):

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
    from mod_outlier import outlier_handler
    olh = outlier_handler()
    i_sel = olh.good_i_flags(miller_array_merge, iparams, flag_show_summary=True)

    miller_array_merge = miller_array_merge.customized_copy(
      indices=miller_indices_merge.select(i_sel),
      data=I_merge.select(i_sel),
      sigmas=sigI_merge.select(i_sel))
    i_seq = flex.int([i for i in range(len(i_sel))])
    stat_sel = [stat_all[i] for i in i_seq.select(i_sel)]
    stat_all = stat_sel
    I_even = I_even.select(i_sel)
    I_odd = I_odd.select(i_sel)

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
    csv_out +='Bin, Low, High, Completeness, <N_obs>, Qmeas, Qw, CC1/2, N_ind, CCiso, N_ind, <I/sigI>\n'
    txt_out = '\n'
    txt_out += 'Summary for '+output_mtz_file_prefix+'_merge.mtz\n'
    txt_out += 'Bin Resolution Range     Completeness      <N_obs>  |Qmeas    Qw     CC1/2   N_ind |CCiso   N_ind| <I/sigI>   <I>\n'
    txt_out += '--------------------------------------------------------------------------------------------------------------------\n'
    sum_r_meas_w_top = 0
    sum_r_meas_w_btm = 0
    sum_r_meas_top = 0
    sum_r_meas_btm = 0
    sum_refl_obs = 0
    sum_refl_complete = 0
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
      multiplicity_all = flex.double()
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
      multiplicity_all.append(multiplicity_bin)
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


      txt_out += '%02d %7.2f - %7.2f %5.1f %6.0f / %6.0f %7.2f %7.2f %7.2f %7.2f %6.0f %7.2f %6.0f %8.2f %10.2f' \
          %(i, binner_template_asu.bin_d_range(i)[0], binner_template_asu.bin_d_range(i)[1], completeness*100, \
          len(miller_indices_obs_bin), len(miller_indices_bin),\
          multiplicity_bin, r_meas_bin*100, r_meas_w_bin*100, cc12_bin*100, n_refl_cc12_bin, cc_iso_bin*100, n_refl_cciso_bin, mean_i_over_sigi_bin, np.mean(I_bin))
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

    txt_out += '--------------------------------------------------------------------------------------------------------------------\n'
    txt_out += '        TOTAL        %5.1f %6.0f / %6.0f %7.2f %7.2f %7.2f %7.2f %6.0f %7.2f %6.0f %8.2f %10.2f\n' \
    %((sum_refl_obs/sum_refl_complete)*100, sum_refl_obs, \
     sum_refl_complete, np.mean(multiplicity_all), \
     r_meas*100, r_meas_w*100, cc12*100, len(I_even.select(i_even_filter_sel)), cc_iso*100, \
     n_refl_iso, np.mean(miller_array_merge.data()/miller_array_merge.sigmas()), np.mean(miller_array_merge.data()))
    txt_out += '--------------------------------------------------------------------------------------------------------------------\n'

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
