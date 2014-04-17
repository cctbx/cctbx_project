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
  classdocs
  '''

  def __init__(self):
    '''
    Constructor
    '''

  def calc_average_I_sigI(self, I, sigI, G, B, p, SE_I, sin_theta_over_lambda_sq, avg_mode, SE):
    I_full = I/(G * flex.exp(-2*B*sin_theta_over_lambda_sq) * p)
    sigI_full = sigI/(G * flex.exp(-2*B*sin_theta_over_lambda_sq) * p)

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
      I_avg = flex.sum(SE_norm * I_full)/flex.sum(SE_norm)
      sigI_avg = flex.sum(SE_norm * sigI_full)/flex.sum(SE_norm)
    elif avg_mode== 'average':
      I_avg = flex.mean(I_full)
      sigI_avg = flex.mean(sigI_full)

    #Rmeas, Rmeas_w, multiplicity
    multiplicity = len(I)
    if multiplicity == 1:
      n_obs = 1.25
    else:
      n_obs = multiplicity

    r_meas_w_top = flex.sum(((I_full - I_avg)*SE_norm)**2)*math.sqrt(n_obs/(n_obs-1))
    r_meas_w_btm = flex.sum((I_full*SE_norm)**2)
    r_meas_top = flex.sum((I_full - I_avg)**2)*math.sqrt(n_obs/(n_obs-1))
    r_meas_btm = flex.sum((I_full)**2)


    #for calculattion of cc1/2
    #sepearte the observations into two groups
    i_even = range(0,len(I_full),2)
    i_odd = range(1,len(I_full),2)
    I_even = I_full.select(i_even)
    sigI_even = sigI_full.select(i_even)
    SE_norm_even = SE_norm.select(i_even)
    I_odd = I_full.select(i_odd)
    sigI_odd = sigI_full.select(i_odd)
    SE_norm_odd = SE_norm.select(i_odd)
    if len(i_even) > len(i_odd):
      I_odd.append(I_even[len(I_even)-1])
      sigI_odd.append(sigI_even[len(I_even)-1])
      SE_norm_odd.append(SE_norm_even[len(I_even)-1])

    if avg_mode == 'weighted':
      I_avg_even = flex.sum(SE_norm_even * I_even)/flex.sum(SE_norm_even)
      I_avg_odd = flex.sum(SE_norm_odd * I_odd)/flex.sum(SE_norm_odd)
    elif avg_mode== 'average':
      I_avg_even = flex.mean(I_even)
      I_avg_odd = flex.mean(I_odd)

    return I_avg, sigI_avg, (r_meas_w_top, r_meas_w_btm, r_meas_top, r_meas_btm, multiplicity), I_avg_even, I_avg_odd

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

    return uc_mean

  def output_mtz_files(self, results, iph, output_mtz_file_prefix, avg_mode):
    partiality_filter = 0.1
    sigma_filter = 8
    uc_len_tol = 3.5
    uc_angle_tol = 2

    #prepare data for merging
    miller_indices_all = flex.miller_index()
    I_all = flex.double()
    sigI_all = flex.double()
    G_all = flex.double()
    B_all = flex.double()
    k_all = flex.double()
    p_all = flex.double()
    SE_I_all = flex.double()
    SE_all = flex.double()
    sin_sq_all = flex.double()
    cn_good_frame = 0

    for pres in results:
      if pres is not None:
        #check unit-cell
        if (abs(pres.uc_params[0]-iph.target_unit_cell[0]) <= uc_len_tol and abs(pres.uc_params[1]-iph.target_unit_cell[1]) <= uc_len_tol \
        and abs(pres.uc_params[2]-iph.target_unit_cell[2]) <= uc_len_tol and abs(pres.uc_params[3]-iph.target_unit_cell[3]) <= uc_angle_tol \
        and abs(pres.uc_params[4]-iph.target_unit_cell[4]) <= uc_angle_tol and abs(pres.uc_params[5]-iph.target_unit_cell[5]) <= uc_angle_tol):
          cn_good_frame += 1
          sin_theta_over_lambda_sq = pres.observations.two_theta(wavelength=pres.wavelength).sin_theta_over_lambda_sq().data()
          for miller_index, i_obs, sigi_obs, p, se_i, sin_sq in zip(
              pres.observations.indices(), pres.observations.data(),
              pres.observations.sigmas(), pres.partiality, pres.SE_I, sin_theta_over_lambda_sq):

            miller_indices_all.append(miller_index)
            I_all.append(i_obs)
            sigI_all.append(sigi_obs)
            G_all.append(pres.G)
            B_all.append(pres.B)
            p_all.append(p)
            SE_I_all.append(se_i)
            sin_sq_all.append(sin_sq)
            SE_all.append(pres.stats[0])
        else:
          print pres.frame_no, ' unit-cell parameters exceed the limits', list(pres.uc_params)

    #plot stats
    self.plot_stats(results, iph)

    #calculate average unit cell
    uc_mean = self.calc_mean_unit_cell(results)

    #from all observations merge them, using mean
    crystal_symmetry = crystal.symmetry(
        unit_cell=iph.target_unit_cell,
        space_group_symbol=iph.target_space_group)
    miller_set_all=miller.set(
                crystal_symmetry=crystal_symmetry,
                indices=miller_indices_all,
                anomalous_flag=iph.target_anomalous_flag)
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
    SE_I_all_sort = SE_I_all.select(perm)
    sin_sq_all_sort = sin_sq_all.select(perm)
    SE_all_sort = SE_all.select(perm)

    refl_now = 0
    miller_indices_merge = flex.miller_index()
    I_merge = flex.double()
    sigI_merge = flex.double()
    stat_all = []
    I_even = flex.double()
    I_odd = flex.double()
    while refl_now < len(I_obs_all_sort)-1:
      miller_index_group = miller_indices_all_sort[refl_now]
      I_obs_group = flex.double()
      sigI_obs_group = flex.double()
      G_group = flex.double()
      B_group = flex.double()
      p_group = flex.double()
      SE_I_group = flex.double()
      sin_sq_group = flex.double()
      SE_group = flex.double()
      for i in range(refl_now, len(I_obs_all_sort)):
        if miller_indices_all_sort[i][0] == miller_index_group[0] and \
            miller_indices_all_sort[i][1] == miller_index_group[1] and \
            miller_indices_all_sort[i][2] == miller_index_group[2]:

          #select only reflections with higher partiality
          if p_all_sort[i] >= partiality_filter:
            I_obs_group.append(I_obs_all_sort[i])
            sigI_obs_group.append(sigI_obs_all_sort[i])
            G_group.append(G_all_sort[i])
            B_group.append(B_all_sort[i])
            p_group.append(p_all_sort[i])
            SE_I_group.append(SE_I_all_sort[i])
            sin_sq_group.append(sin_sq_all_sort[i])
            SE_group.append(SE_all_sort[i])
          if i == (len(I_obs_all_sort) - 1):
            refl_now = i
            break
        else:
          refl_now = i
          break

      if len(I_obs_group) > 0:
        I_avg, sigI_avg, stat, I_avg_even, I_avg_odd = self.calc_average_I_sigI(I_obs_group, sigI_obs_group,
            G_group, B_group, p_group, SE_I_group, sin_sq_group, avg_mode, SE_group)

        if math.isnan(stat[0]) or math.isinf(stat[0]) or math.isnan(stat[1]) or math.isinf(stat[1]):
          print miller_index_group, ' not merged (Qw=%.4g/%.4g)'%(stat[0],stat[1])
        else:
          miller_indices_merge.append(miller_index_group)
          I_merge.append(I_avg)
          sigI_merge.append(sigI_avg)
          stat_all.append(stat)
          I_even.append(I_avg_even)
          I_odd.append(I_avg_odd)


    #output mtz file and report binning stat
    miller_set_merge=miller.set(
              crystal_symmetry=crystal_symmetry,
              indices=miller_indices_merge,
              anomalous_flag=iph.target_anomalous_flag)
    miller_array_merge = miller_set_merge.array(data=I_merge,
              sigmas=sigI_merge).set_observation_type_xray_intensity()

    #remove outliers
    binner_merge = miller_array_merge.setup_binner(n_bins=iph.n_bins)
    binner_merge_indices = binner_merge.bin_indices()
    miller_indices_merge_filter = flex.miller_index()
    I_merge_filter = flex.double()
    sigI_merge_filter = flex.double()
    I_even_filter = flex.double()
    I_odd_filter = flex.double()
    stat_filter = []
    i_seq = flex.int([j for j in range(len(binner_merge_indices))])
    for i in range(1,iph.n_bins+1):
      i_binner = (binner_merge_indices == i)
      if len(miller_array_merge.data().select(i_binner)) > 0:
        I_obs_bin = miller_array_merge.data().select(i_binner)
        sigI_obs_bin = miller_array_merge.sigmas().select(i_binner)
        miller_indices_bin = miller_array_merge.indices().select(i_binner)
        stat_bin = [stat_all[j] for j in i_seq.select(i_binner)]
        I_even_bin = I_even.select(i_binner)
        I_odd_bin = I_odd.select(i_binner)

        mean_I_obs_bin = np.mean(I_obs_bin)
        std_I_obs_bin = np.std(I_obs_bin)
        i_filter = flex.abs((I_obs_bin - mean_I_obs_bin)/std_I_obs_bin) < sigma_filter
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

    miller_set_merge=miller.set(
              crystal_symmetry=crystal_symmetry,
              indices=miller_indices_merge_filter,
              anomalous_flag=iph.target_anomalous_flag)
    miller_array_merge = miller_set_merge.array(data=I_merge_filter,
              sigmas=sigI_merge_filter).set_observation_type_xray_intensity()

    if output_mtz_file_prefix != '':
      mtz_dataset_merge = miller_array_merge.as_mtz_dataset(column_root_label="IOBS")
      mtz_dataset_merge.mtz_object().write(file_name=output_mtz_file_prefix+'_merge.mtz')

    #report binning stats
    binner_merge = miller_array_merge.setup_binner(n_bins=iph.n_bins)
    binner_merge_indices = binner_merge.bin_indices()
    completeness_merge = miller_array_merge.completeness(use_binning=True)

    txt_out = '\n'
    txt_out += 'Summary for '+output_mtz_file_prefix+'_merge.mtz\n'
    txt_out += 'Bin Resolution Range     Completeness    <N_obs>  Qmeas    Qw    CC1/2   CCiso  <I/sigI>\n'
    txt_out += '----------------------------------------------------------------------------------------\n'
    sum_bin_completeness = 0
    sum_r_meas_w_top = 0
    sum_r_meas_w_btm = 0
    sum_r_meas_top = 0
    sum_r_meas_btm = 0
    i_seq = flex.int([j for j in range(len(binner_merge_indices))])
    for i in range(1,iph.n_bins+1):
      i_binner = (binner_merge_indices == i)
      if len(miller_array_merge.data().select(i_binner)) == 0:
        mean_i_over_sigi_bin = 0
        multiplicity_bin = 0
        r_meas_w_bin = 0
        r_meas_bin = 0
        cc12 = 0
      else:
        mean_i_over_sigi_bin = flex.mean(miller_array_merge.data().select(i_binner)/miller_array_merge.sigmas().select(i_binner))
        stat_bin = [stat_filter[j] for j in i_seq.select(i_binner)]
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

        multiplicity_bin = sum_mul_bin/completeness_merge.binner.counts_given()[i]
        r_meas_w_bin = sum_r_meas_w_top_bin/ sum_r_meas_w_btm_bin
        r_meas_bin = sum_r_meas_top_bin/ sum_r_meas_btm_bin

        cc12_bin = np.corrcoef(I_even_filter.select(i_binner), I_odd_filter.select(i_binner))[0,1]

      completeness = completeness_merge.data[i]
      sum_bin_completeness += completeness

      #calculate CCiso
      cc_iso_bin = 0
      if iph.file_name_iso_mtz != '':
        matches_iso = miller.match_multi_indices(
                  miller_indices_unique=iph.miller_array_iso.indices(),
                  miller_indices=miller_array_merge.indices().select(i_binner))

        I_iso = flex.double([iph.miller_array_iso.data()[pair[0]] for pair in matches_iso.pairs()])
        I_merge_match_iso = flex.double([miller_array_merge.data().select(i_binner)[pair[1]] for pair in matches_iso.pairs()])
        if len(matches_iso.pairs()) > 0 :
          cc_iso_bin = np.corrcoef(I_merge_match_iso, I_iso)[0,1]

        if iph.flag_plot:
          plt.scatter(I_iso, I_merge_match_iso,s=10, marker='x', c='r')
          plt.title('bin %3.0f CC=%.4g meanI=%.4g std=%.4g sqrt_meanI=%.4g mul=%.4g'%(i, cc_iso_bin, np.mean(I_merge_match_iso), np.std(I_merge_match_iso), math.sqrt(np.mean(I_merge_match_iso)), math.sqrt(np.mean(I_merge_match_iso))*2.5))
          plt.xlabel('I_ref')
          plt.ylabel('I_obs')
          plt.show()

      txt_out += '%02d %7.2f - %7.2f %5.1f %6.0f / %6.0f %5.2f %7.2f %7.2f %7.2f %7.2f %7.2f' \
          %(i, binner_merge.bin_d_range(i)[0], binner_merge.bin_d_range(i)[1], completeness*100, \
          completeness_merge.binner.counts_given()[i], completeness_merge.binner.counts_complete()[i],\
          multiplicity_bin, r_meas_bin*100, r_meas_w_bin*100, cc12_bin, cc_iso_bin, mean_i_over_sigi_bin)
      txt_out += '\n'

    #calculate CCiso
    cc_iso = 0
    if iph.file_name_iso_mtz != '':
      matches_iso = miller.match_multi_indices(
                miller_indices_unique=iph.miller_array_iso.indices(),
                miller_indices=miller_array_merge.indices())

      I_iso = flex.double([iph.miller_array_iso.data()[pair[0]] for pair in matches_iso.pairs()])
      I_merge_match_iso = flex.double([miller_array_merge.data()[pair[1]] for pair in matches_iso.pairs()])
      if len(matches_iso.pairs()) > 0 :
        cc_iso = np.corrcoef(I_merge_match_iso, I_iso)[0,1]
      if iph.flag_plot:
        plt.scatter(I_iso, I_merge_match_iso,s=10, marker='x', c='r')
        plt.title('CC=%.4g'%(cc_iso))
        plt.xlabel('I_ref')
        plt.ylabel('I_obs')
        plt.show()

    #calculate cc12
    cc12 = np.corrcoef(I_even_filter, I_odd_filter)[0,1]

    txt_out += '----------------------------------------------------------------------------------------\n'
    txt_out += '        TOTAL        %5.1f %6.0f / %6.0f %5.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n' \
    %((sum_bin_completeness/iph.n_bins)*100, np.sum(completeness_merge.binner.counts_given()), \
     np.sum(completeness_merge.binner.counts_complete()), len(miller_indices_all)/len(miller_array_merge.data()), \
     (sum_r_meas_top/sum_r_meas_btm)*100, (sum_r_meas_w_top/sum_r_meas_w_btm)*100, cc12, cc_iso, \
     np.mean(miller_array_merge.data()/miller_array_merge.sigmas()))
    txt_out += '----------------------------------------------------------------------------------------\n'
    txt_out += 'No. of total observed reflections: %9.0f from %5.0f frames' %(len(miller_indices_all), cn_good_frame)
    txt_out += '\n'
    txt_out += 'Average unit-cell parameters: (%6.2f, %6.2f, %6.2f %6.2f, %6.2f, %6.2f)'%(uc_mean[0], uc_mean[1], uc_mean[2], uc_mean[3], uc_mean[4], uc_mean[5])
    txt_out += '\n'
    print txt_out


    return miller_array_merge, txt_out

  def plot_stats(self, results, iph):
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
    sum_var_I_p_all = flex.double()
    sum_var_k_all = flex.double()
    sum_var_p_all = flex.double()
    SE_I_all = flex.double()

    for pres in results:
      if pres is not None:
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
        SE_all.append(pres.stats[0])
        R_sq_all.append(pres.stats[1])
        cc_all.append(pres.stats[2])
        sum_var_I_p_all.append(np.median(pres.var_I_p))
        sum_var_k_all.append(np.median(pres.var_k))
        sum_var_p_all.append(np.median(pres.var_p))

        for se_i in pres.SE_I:
          SE_I_all.append(se_i)


    if iph.flag_plot:
      plt.subplot(231)
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

      plt.subplot(232)
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

      plt.subplot(233)
      x = SE_I_all.as_numpy_array()
      mu = np.mean(x)
      med = np.median(x)
      sigma = np.std(x)
      num_bins = 10
      n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
      y = mlab.normpdf(bins, mu, sigma)
      plt.plot(bins, y, 'r--')
      plt.ylabel('Frequencies')
      plt.title('SE I all observations distribution\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))

      plt.subplot(234)
      x = sum_var_I_p_all.as_numpy_array()
      mu = np.mean(x)
      med = np.median(x)
      sigma = np.std(x)
      num_bins = 10
      n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
      y = mlab.normpdf(bins, mu, sigma)
      plt.plot(bins, y, 'r--')
      plt.ylabel('Frequencies')
      plt.title('SE I distribution\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))

      plt.subplot(235)
      x = sum_var_k_all.as_numpy_array()
      mu = np.mean(x)
      med = np.median(x)
      sigma = np.std(x)
      num_bins = 10
      n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
      y = mlab.normpdf(bins, mu, sigma)
      plt.plot(bins, y, 'r--')
      plt.ylabel('Frequencies')
      plt.title('SE G distribution\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))

      plt.subplot(236)
      x = sum_var_p_all.as_numpy_array()
      mu = np.mean(x)
      med = np.median(x)
      sigma = np.std(x)
      num_bins = 10
      n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
      y = mlab.normpdf(bins, mu, sigma)
      plt.plot(bins, y, 'r--')
      plt.ylabel('Frequencies')
      plt.title('SE p distribution\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))

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

class input_handler(object):
  '''
  handle reading txt file
  '''

  def __init__(self):
    '''
    Constructor
    '''

  def read_input(self, file_name_input):

    self.run_no = ''
    self.title = ''
    self.frame_start = 0
    self.frame_end = 0
    self.flag_plot = False
    self.flag_polar = False
    self.file_name_iso_mtz = ''
    self.file_name_ref_mtz = ''
    self.file_name_in_energy = ''
    self.file_name_in_img = ''
    self.pickle_dir = ''
    self.d_min = 0
    self.d_min_merge = 0
    self.d_max = 45
    self.sigma_max = 15
    self.sigma_max_merge=99
    self.flag_reference_a_matrix = False
    self.target_unit_cell = ''
    self.target_space_group = ''
    self.target_anomalous_flag = False
    self.file_name_pdb = ''
    self.n_postref_cycle = 1
    self.miller_array_iso = None
    self.n_bins=25
    self.flag_on_screen_output=True
    self.q_w_merge = 1.0
    self.pixel_size_mm = 0


    file_input = open(file_name_input, 'r')
    data_input = file_input.read().split('\n')

    for line_input in data_input:
      pair=line_input.split('=')
      if len(pair) == 2:
        param_name = pair[0].strip()
        param_val = pair[1].strip()
        if param_name=='run_no':
          self.run_no=param_val
        elif param_name=='title':
          self.title=param_val
        elif param_name=='frame_start':
          self.frame_start=int(param_val)
        elif param_name=='frame_end':
          self.frame_end=int(param_val)
        elif param_name=='flag_plot':
          if param_val=='True':
            self.flag_plot=True
        elif param_name=='flag_polar':
          if param_val=='True':
            self.flag_polar=True
        elif param_name=='flag_reference_a_matrix':
          if param_val=='True':
            self.flag_reference_a_matrix=True
        elif param_name=='hklisoin':
          self.file_name_iso_mtz=param_val
        elif param_name=='hklrefin':
          self.file_name_ref_mtz=param_val
        elif param_name=='energyin':
          self.file_name_in_energy=param_val
        elif param_name=='imagein':
          self.file_name_in_img=param_val
        elif param_name=='pdbin':
          self.file_name_pdb=param_val
        elif param_name=='pickle_dir':
          self.pickle_dir=param_val
        elif param_name=='d_min':
          self.d_min=float(param_val)
        elif param_name=='d_min_merge':
          self.d_min_merge=float(param_val)
        elif param_name=='d_max':
          self.d_max=float(param_val)
        elif param_name=='sigma_max':
          self.sigma_max=float(param_val)
        elif param_name=='sigma_max_merge':
          self.sigma_max_merge=float(param_val)
        elif param_name=='target_unit_cell':
          tmp_uc = param_val.split(',')
          if len(tmp_uc) != 6:
            print 'Parameter: target_unit_cell has wrong format (usage: target_unit_cell= a,b,c,alpha,beta,gamma)'
            exit()
          else:
            self.target_unit_cell = (float(tmp_uc[0]), float(tmp_uc[1]), float(tmp_uc[2]), \
              float(tmp_uc[3]), float(tmp_uc[4]), float(tmp_uc[5]))
        elif param_name=='target_space_group':
          self.target_space_group=param_val
        elif param_name=='target_anomalous_flag':
          if param_val=='True':
            self.target_anomalous_flag=True
        elif param_name=='n_postref_cycle':
          self.n_postref_cycle=int(param_val)
        elif param_name=='n_bins':
          self.n_bins=int(param_val)
        elif param_name=='flag_on_screen_output':
          if param_val=='False':
            self.flag_on_screen_output=False
        elif param_name=='q_w_merge':
          self.q_w_merge=float(param_val)
        elif param_name=='pixel_size_mm':
          self.pixel_size_mm=float(param_val)

    if self.frame_end == 0:
      print 'Parameter: frame_end - please specifiy at least one frame (usage: frame_end=1)'
      exit()

    if self.pickle_dir == '':
      print 'Parameter: pickle_dir - please specify file path to pickle files (usage: pickle_dir=/path/to/pickles)'
      exit()

    if self.target_space_group == '':
      print 'Parameter: target_space_group - please specify space_group (usage: target_space_group=SGSYMBOL)'
      exit()

    if self.flag_polar and self.file_name_iso_mtz =='':
      print 'Conflict of parameters: you turned flag_polar on, please also input isomorphous-reference mtz file (usage: hklisoin = 1jw8-sf-asu.mtz)'
      exit()

    if self.pixel_size_mm == 0:
      print 'pixel size (in mm) is required (usage: pixel_size_mm = 0.08)'
      exit()

    #fetch isomorphous structure
    if self.file_name_iso_mtz != '':
      reflection_file_iso = reflection_file_reader.any_reflection_file(self.file_name_iso_mtz)
      miller_arrays_iso=reflection_file_iso.as_miller_arrays()
      is_found_iso_as_intensity_array = False
      is_found_iso_as_amplitude_array = False
      for miller_array_iso in miller_arrays_iso:
        if miller_array_iso.is_xray_intensity_array():
          self.miller_array_iso = miller_array_iso
          is_found_iso_as_intensity_array = True
          break
        elif miller_array_iso.is_xray_amplitude_array():
          is_found_iso_as_amplitude_array = True
          miller_array_iso_converted_to_intensity = miller_array_iso.as_intensity_array()
      if is_found_iso_as_intensity_array == False:
        if is_found_iso_as_amplitude_array:
          print 'Found amplitude array, convert it to intensity array'
          self.miller_array_iso = miller_array_iso_converted_to_intensity
        else:
          print 'Cannot find intensity array in the isomorphous-reference mtz'
          exit()

      self.miller_array_iso = self.miller_array_iso.expand_to_p1().generate_bijvoet_mates()

    self.txt_out = ''
    self.txt_out += 'Input parameters\n'
    self.txt_out += 'run_no '+str(self.run_no)+'\n'
    self.txt_out += 'title '+str(self.title)+'\n'
    self.txt_out += 'frame '+str(self.frame_start)+'-'+str(self.frame_end)+'\n'
    self.txt_out += 'flag_plot '+str(self.flag_plot)+'\n'
    self.txt_out += 'flag_polar '+str(self.flag_polar)+'\n'
    self.txt_out += 'flag_reference_a_matrix '+str(self.flag_reference_a_matrix)+'\n'
    self.txt_out += 'hklisoin '+str(self.file_name_iso_mtz)+'\n'
    self.txt_out += 'hklrefin '+str(self.file_name_ref_mtz)+'\n'
    self.txt_out += 'energyin '+str(self.file_name_in_energy)+'\n'
    self.txt_out += 'imagein '+str(self.file_name_in_img)+'\n'
    self.txt_out += 'pdbin '+str(self.file_name_pdb)+'\n'
    self.txt_out += 'picke_dir '+str(self.pickle_dir)+'\n'
    self.txt_out += 'd_min '+str(self.d_min)+'\n'
    self.txt_out += 'd_min_merge '+str(self.d_min_merge)+'\n'
    self.txt_out += 'd_max '+str(self.d_max)+'\n'
    self.txt_out += 'sigma_max '+str(self.sigma_max)+'\n'
    self.txt_out += 'sigma_max_merge '+str(self.sigma_max_merge)+'\n'
    self.txt_out += 'q_w_merge '+str(self.q_w_merge)+'\n'
    self.txt_out += 'target_unit_cell '+str(self.target_unit_cell)+'\n'
    self.txt_out += 'target_space_group '+str(self.target_space_group)+'\n'
    self.txt_out += 'target_anomalous_flag '+str(self.target_anomalous_flag)+'\n'

    print self.txt_out

class file_handler(object):
  '''
  handle reading txt file
  '''

  def __init__(self):
    '''
    Constructor
    '''

  def get_imgname_from_pickle_filename(self, file_name_in_img, pickle_filename):

    img_filename = ''

    file_img = open(file_name_in_img, 'r')
    data_img = file_img.read().split('\n')
    for line_img in data_img:
      data_img = line_img.split(' ')
      if pickle_filename.find(data_img[0]) > 0:
        img_filename = data_img[1]
        break

    return img_filename



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
