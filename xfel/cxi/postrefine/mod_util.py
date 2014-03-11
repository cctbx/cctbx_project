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

def get_overall_correlation (data_a, data_b) :
  """
  Correlate the averaged intensities to the intensities from the
  reference data set.
  """

  assert len(data_a) == len(data_b)
  corr = 0
  slope = 0
  try:
    sum_xx = 0
    sum_xy = 0
    sum_yy = 0
    sum_x  = 0
    sum_y  = 0
    N      = 0
    for i in xrange(len(data_a)):
      I_r       = data_a[i]
      I_o       = data_b[i]
      N      += 1
      sum_xx += I_r**2
      sum_yy += I_o**2
      sum_xy += I_r * I_o
      sum_x  += I_r
      sum_y  += I_o
    slope = (N * sum_xy - sum_x * sum_y) / (N * sum_xx - sum_x**2)
    corr  = (N * sum_xy - sum_x * sum_y) / (math.sqrt(N * sum_xx - sum_x**2) *
               math.sqrt(N * sum_yy - sum_y**2))
  except ZeroDivisionError:
    pass

  return corr, slope

class intensities_scaler(object):
  '''
  classdocs
  '''

  def __init__(self):
    '''
    Constructor
    '''

  def count_unique_frame_multiplicities(self, frame_no_list):
    from collections import Counter
    c = Counter(frame_no_list)

    return len(c.items())

  def output_mtz_files(self, observations_set,
        observations_weight_set,
        iph,
        output_mtz_file_prefix):


    #from observations_set, generate combine all
    miller_indices_all = flex.miller_index()
    I_all = flex.double()
    sigI_all = flex.double()
    for observations in observations_set:
      if observations is not None:
        for miller_index, i_obs, sigi_obs in zip(
            observations.indices(), observations.data(),
            observations.sigmas()):

          miller_indices_all.append(miller_index)
          I_all.append(i_obs)
          sigI_all.append(sigi_obs)
      else:
        print 'Something wrong'

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

    perm = miller_array_all.sort_permutation(by_value="packed_indices")
    miller_indices_all_sort = miller_array_all.indices().select(perm)
    I_obs_all_sort = miller_array_all.data().select(perm)
    sigI_obs_all_sort = miller_array_all.sigmas().select(perm)

    refl_now = 0
    miller_indices_merge = flex.miller_index()
    I_obs_merge_mean = flex.double()
    sigI_obs_merge_mean = flex.double()
    multiplicities = flex.double()
    r_merge_dif_all = flex.double()
    r_merge_w_all = flex.double()
    cn_miller_indices = 0
    I_even_merge_mean = flex.double()
    I_odd_merge_mean = flex.double()
    while refl_now < len(I_obs_all_sort)-1:
      miller_index_group = miller_indices_all_sort[refl_now]
      I_obs_group = flex.double()
      sigI_obs_group = flex.double()
      for i in range(refl_now, len(I_obs_all_sort)):
        if miller_indices_all_sort[i][0] == miller_index_group[0] and \
            miller_indices_all_sort[i][1] == miller_index_group[1] and \
            miller_indices_all_sort[i][2] == miller_index_group[2]:
          I_obs_group.append(I_obs_all_sort[i])
          sigI_obs_group.append(sigI_obs_all_sort[i])

        else:
          refl_now = i
          cn_miller_indices +=1
          break

      if len(I_obs_group) == 1:
        I_obs_merge = I_obs_group[0]
        sigI_obs_merge = sigI_obs_group[0]
        no_observations = 1
        r_dif = ((sigI_obs_merge**2)/sigI_obs_merge)*1.7321
        r_w = (I_obs_merge**2)/sigI_obs_merge
        I_even_merge = I_obs_merge
        I_odd_merge = I_obs_merge
      else:

        I_as_sigma = (I_obs_group - np.mean(I_obs_group))/np.std(I_obs_group)
        i_sel = (flex.abs(I_as_sigma) <= iph.sigma_max_merge)
        I_obs_group_sel = I_obs_group.select(i_sel)
        sigI_obs_group_sel = sigI_obs_group.select(i_sel)

        if len(I_obs_group_sel) == 1:
          I_obs_merge = I_obs_group_sel[0]
          sigI_obs_merge = sigI_obs_group_sel[0]
          no_observations = 1
          r_dif = ((sigI_obs_merge**2)/sigI_obs_merge)*1.7321
          r_w = (I_obs_merge**2)/sigI_obs_merge
          I_even_merge = I_obs_merge
          I_odd_merge = I_obs_merge
        else:
          I_obs_merge = np.mean(I_obs_group_sel)
          sigI_obs_merge = np.mean(sigI_obs_group_sel)
          no_observations = len(I_obs_group_sel)
          r_dif = flex.sum(((I_obs_group_sel - I_obs_merge)**2)/sigI_obs_group_sel)*math.sqrt(len(I_obs_group_sel)/(len(I_obs_group_sel)-1))
          r_w = flex.sum((I_obs_group_sel**2)/sigI_obs_group_sel)

          #sepearte the observations into two groups
          i_even = range(0,len(I_obs_group_sel),2)
          i_odd = range(1,len(I_obs_group_sel),2)
          I_even_sel = I_obs_group_sel.select(i_even)
          I_odd_sel = I_obs_group_sel.select(i_odd)
          if len(i_even) > len(i_odd):
            I_odd_sel.append(I_even_sel[len(I_even_sel)-1])
          I_even_merge = np.mean(I_even_sel)
          I_odd_merge = np.mean(I_odd_sel)


      if math.isnan(I_obs_merge) or math.isnan(sigI_obs_merge) or I_obs_merge==0 or sigI_obs_merge==0:
        print miller_index_group, I_obs_merge, sigI_obs_merge, ' - not merged'
      else:
        multiplicities.append(no_observations)
        miller_indices_merge.append(miller_index_group)
        I_obs_merge_mean.append(I_obs_merge)
        sigI_obs_merge_mean.append(sigI_obs_merge)
        r_merge_dif_all.append(r_dif)
        r_merge_w_all.append(r_w)
        I_even_merge_mean.append(I_even_merge)
        I_odd_merge_mean.append(I_odd_merge)

      if i==len(I_obs_all_sort)-1:
        break


    miller_set_merge=miller.set(
              crystal_symmetry=crystal_symmetry,
              indices=miller_indices_merge,
              anomalous_flag=iph.target_anomalous_flag)
    miller_array_merge_mean = miller_set_merge.array(data=I_obs_merge_mean,
              sigmas=sigI_obs_merge_mean).set_observation_type_xray_intensity()

    if output_mtz_file_prefix != '':
      #write out mtz file
      mtz_dataset_merge_mean = miller_array_merge_mean.as_mtz_dataset(column_root_label="IOBS")
      mtz_dataset_merge_mean.mtz_object().write(file_name=output_mtz_file_prefix+'_merge.mtz')

    cc_merge_mean = 0
    slope_merge_mean = 0
    r_all = 0
    if iph.miller_array_iso is not None:
      #calculate final cc
      matches_merge = miller.match_multi_indices(
                    miller_indices_unique=iph.miller_array_iso.indices(),
                    miller_indices=miller_indices_merge)

      I_ref_merge = flex.double([iph.miller_array_iso.data()[pair[0]] for pair in matches_merge.pairs()])
      I_obs_merge_mean = flex.double([miller_array_merge_mean.data()[pair[1]] for pair in matches_merge.pairs()])
      sigI_obs_merge_mean = flex.double([miller_array_merge_mean.sigmas()[pair[1]] for pair in matches_merge.pairs()])

      cc_merge_mean, slope_merge_mean = get_overall_correlation(I_obs_merge_mean, I_ref_merge)

      r_all = sum(((I_ref_merge - (slope_merge_mean*I_obs_merge_mean))**2)/sigI_obs_merge_mean)/sum(((slope_merge_mean*I_obs_merge_mean)**2)/sigI_obs_merge_mean)

      """
      plt.scatter(I_ref_merge, I_obs_merge_mean,s=10, marker='x', c='r')
      plt.title('CC=%6.5f Slope=%6.5f R=%6.5f <I_ref>=%6.5f <I_obs>=%6.5f'%(cc_merge_mean, slope_merge_mean, r_all, np.mean(I_ref_merge), np.mean(I_obs_merge_mean)))
      plt.xlabel('Reference intensity')
      plt.ylabel('Observed intensity')
      plt.show()
      """

    #report stat
    txt_out = ''
    if iph.flag_on_screen_output:
      binner_merge = miller_array_merge_mean.setup_binner(n_bins=iph.n_bins)
      binner_merge_indices = binner_merge.bin_indices()
      completeness_merge = miller_array_merge_mean.completeness(use_binning=True)


      txt_out = '------------------------------------------------------\n'
      txt_out += 'Summary for'+output_mtz_file_prefix+'_merge.mtz\n'
      txt_out += 'Bin Resolutions %Complete (obs/total) Obs# CCiso CC1/2 Qiso Qint  <I>  <sigma> <I/sigI> (median)\n'
      bin_multiplicities = []
      sum_bin_completeness = 0
      sum_bin_multiplicities = 0
      one_over_dsqr = flex.double()
      cc_bin = flex.double()
      n_ref_complete = flex.double()
      n_ref_observed = flex.double()
      for i in range(1,iph.n_bins+1):
        i_binner = (binner_merge_indices == i)

        corr_bin_mean = 0
        slope_bin_mean = 0
        r_bin = 0
        if iph.miller_array_iso is not None:
          #calculate C.C.iso per bin
          matches_bin = miller.match_multi_indices(
                    miller_indices_unique=iph.miller_array_iso.indices(),
                    miller_indices=miller_array_merge_mean.indices().select(i_binner))

          if len(matches_bin.pairs()) > 10:
            I_ref_bin = flex.double([iph.miller_array_iso.data()[pair[0]] for pair in matches_bin.pairs()])
            I_obs_bin_mean = flex.double([miller_array_merge_mean.data().select(i_binner)[pair[1]] for pair in matches_bin.pairs()])
            sigI_obs_bin_mean = flex.double([miller_array_merge_mean.sigmas().select(i_binner)[pair[1]] for pair in matches_bin.pairs()])
            corr_bin_mean, slope_bin_mean = get_overall_correlation(I_obs_bin_mean, I_ref_bin)

            r_bin = sum(((I_ref_bin - (slope_bin_mean*I_obs_bin_mean))**2)/sigI_obs_bin_mean)/sum(((slope_bin_mean*I_obs_bin_mean)**2)/sigI_obs_bin_mean)

        cc_bin.append(corr_bin_mean)

        if len(miller_array_merge_mean.data().select(i_binner)) == 0:
          mean_i_bin = 0
          mean_sigi_bin = 0
          mean_i_over_sigi_bin = 0
          median_i_over_sigi_bin = 0
          multiplicities_sel = 0
          r_merge_dif_sel = flex.double()
          r_merge_w_sel = flex.double()
          cc_onehalf_bin = 0
        else:
          mean_i_bin = np.mean(miller_array_merge_mean.data().select(i_binner))
          mean_sigi_bin = np.mean(miller_array_merge_mean.sigmas().select(i_binner))
          mean_i_over_sigi_bin = np.mean(miller_array_merge_mean.data().select(i_binner)/miller_array_merge_mean.sigmas().select(i_binner))
          median_i_over_sigi_bin = np.median(miller_array_merge_mean.data().select(i_binner)/miller_array_merge_mean.sigmas().select(i_binner))
          multiplicities_sel = multiplicities.select(i_binner)
          r_merge_dif_sel = r_merge_dif_all.select(i_binner)
          r_merge_w_sel = r_merge_w_all.select(i_binner)
          cc_onehalf_bin, slope_onehalf_bin = get_overall_correlation(I_even_merge_mean.select(i_binner), I_odd_merge_mean.select(i_binner))


        completeness = completeness_merge.data[i]

        sum_bin_completeness += completeness
        sum_bin_multiplicities += np.mean(multiplicities_sel)

        r_merge_bin = 0
        if len(r_merge_w_sel) > 0:
          if sum(r_merge_w_sel) > 0:
            r_merge_bin = (sum(r_merge_dif_sel)/sum(r_merge_w_sel))

        txt_out += '%02d %5.2f - %5.2f %6.2f (%05d/%05d) %5.2f %4.2f %4.2f %6.2f %6.2f %7.2f %7.2f %7.2f %7.2f' \
          %(i, binner_merge.bin_d_range(i)[0], binner_merge.bin_d_range(i)[1], completeness*100, \
          completeness_merge.binner.counts_given()[i], completeness_merge.binner.counts_complete()[i],\
          np.mean(multiplicities_sel), corr_bin_mean, cc_onehalf_bin, r_bin, r_merge_bin, \
          mean_i_bin, mean_sigi_bin, mean_i_over_sigi_bin, median_i_over_sigi_bin)
        txt_out += '\n'

        one_over_dsqr.append(1/(binner_merge.bin_d_range(i)[1]**2))
        n_ref_complete.append(completeness_merge.binner.counts_complete()[i])
        n_ref_observed.append(completeness_merge.binner.counts_given()[i])


      cc_onehalf_all, slope_onehalf_all = get_overall_correlation(I_even_merge_mean, I_odd_merge_mean)
      txt_out += 'Overall completeness: %3.4f no. of observations: %3.4f' %((sum_bin_completeness/iph.n_bins)*100, sum_bin_multiplicities/iph.n_bins)
      txt_out += '\n'
      txt_out += 'C.C.iso=%6.5f Slope=%6.5f' %(cc_merge_mean, slope_merge_mean)
      txt_out += '\n'
      txt_out += 'C.C.1/2=%6.5f Slope=%6.5f' %(cc_onehalf_all, slope_onehalf_all)
      txt_out += '\n'
      txt_out += 'Qiso=%6.5f Qint=%6.5f' %(r_all, (sum(r_merge_dif_all)/sum(r_merge_w_all)))
      txt_out += '\n'
      txt_out += 'No. of reflections after merge: %9.2f' %(len(miller_array_merge_mean.data()))
      txt_out += '\n'
      txt_out += '------------------------------------------------------'
      txt_out += '\n'
      print txt_out


    return miller_array_merge_mean, cc_merge_mean, slope_merge_mean, txt_out


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

    #fetch isomorphous structure
    if self.file_name_iso_mtz != '':
      reflection_file_iso = reflection_file_reader.any_reflection_file(self.file_name_iso_mtz)
      miller_arrays_iso=reflection_file_iso.as_miller_arrays()
      is_found_iso_as_intensity_array = False
      for miller_array_iso in miller_arrays_iso:
        if miller_array_iso.is_xray_intensity_array():
          self.miller_array_iso = miller_array_iso.generate_bijvoet_mates()
          is_found_iso_as_intensity_array = True
      if is_found_iso_as_intensity_array == False:
        print 'Cannot find intensity array in the isomorphous-reference mtz'
        exit()


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
