from __future__ import division, print_function

class misc_handler(object):
  def calculate_SE(self, results, iparams):
    """
    Take all post-refinement results, calculate covariances and new SE
    """
    if results[0].grad_set is None:
      return results
    else:
      #get miller array iso for comparison, if given
      mxh = mx_handler()
      flag_hklisoin_found, miller_array_iso = mxh.get_miller_array_from_reflection_file(iparams.hklisoin)
      #get reference set
      import os
      fileseq_list = flex.int()
      for file_in in os.listdir(iparams.run_no):
        if file_in.endswith('.mtz'):
          file_split = file_in.split('_')
          if len(file_split) > 3:
            fileseq_list.append(int(file_split[2]))
      if len(fileseq_list) == 0:
        hklref = 'mean_scaled_merge.mtz'
      else:
        hklref = 'postref_cycle_'+str(flex.max(fileseq_list))+'_merge.mtz'
      flag_hklref_found, miller_array_ref = mxh.get_miller_array_from_reflection_file(iparams.run_no+'/'+hklref)
      #calculate covariance
      X = np.array([[pres.G, pres.B, pres.tau, pres.rotx, pres.roty, pres.ry, pres.rz, pres.r0, pres.re, \
        pres.uc_params[0], pres.uc_params[1],pres.uc_params[2],pres.uc_params[3], \
        pres.uc_params[4],pres.uc_params[5]] for pres in results]).T
      COV = np.cov(X)
      COV_diag = flex.double([i for i in COV.diagonal()])
      #calculate standard errors
      output_results = []
      for pres in results:
        sigI = pres.observations.sigmas()
        observations_old = pres.observations.deep_copy()
        var_set = (pres.grad_set**2) * COV_diag
        sin_theta_over_lambda_sq = pres.observations.two_theta(wavelength=pres.wavelength).\
          sin_theta_over_lambda_sq().data()
        d_spacings = pres.observations.d_spacings().data()
        scale_factors_by_indices = pres.G * flex.exp(-2*pres.B*sin_theta_over_lambda_sq)
        var_scale_factors = var_set[0] + var_set[1]
        var_partiality = flex.sum(var_set[3:])
        err_scale_factors = var_scale_factors/(scale_factors_by_indices**2)
        err_partiality = var_partiality/(pres.partiality**2)
        err_tau = var_set[2]/((pres.tau**2)+1E-9)
        #determine weight
        I_o_full = ((4 * pres.rs_set * pres.observations.data())/(3 * pres.e_width_set * scale_factors_by_indices * pres.partiality)) + pres.tau
        observations_full = pres.observations.customized_copy(data=I_o_full)
        observations_err_scale = pres.observations.customized_copy(data=err_scale_factors)
        observations_err_p = pres.observations.customized_copy(data=err_partiality)
        observations_err_tau = pres.observations.customized_copy(data=flex.double([err_tau]*len(pres.observations.data())))
        flag_reset_weight = False
        if flag_hklref_found:
          ma_ref, ma_obs_full = miller_array_ref.common_sets(observations_full, assert_is_similar_symmetry=False)
          dummy, ma_err_scale = miller_array_ref.common_sets(observations_err_scale, assert_is_similar_symmetry=False)
          dummy, ma_err_p = miller_array_ref.common_sets(observations_err_p, assert_is_similar_symmetry=False)
          dummy, ma_err_tau = miller_array_ref.common_sets(observations_err_tau, assert_is_similar_symmetry=False)
          r1_factor = ma_ref.data()-ma_obs_full.data()
          SE_I_WEIGHT = max([flex.linear_correlation((ma_obs_full.sigmas()**2), flex.log(r1_factor**2)).coefficient(), 0])
          SE_SCALE_WEIGHT = max([flex.linear_correlation(ma_err_scale.data(), flex.log(r1_factor**2)).coefficient(), 0])
          SE_EOC_WEIGHT = max([flex.linear_correlation(ma_err_p.data(), flex.log(r1_factor**2)).coefficient(), 0])
          SE_TAU_WEIGHT = max([flex.linear_correlation(ma_err_tau.data(), flex.log(r1_factor**2)).coefficient(), 0])
        else:
          flag_reset_weight = True
        if True:
          SE_I_WEIGHT = self.CONST_SE_I_WEIGHT
          SE_SCALE_WEIGHT = self.CONST_SE_SCALE_WEIGHT
          SE_EOC_WEIGHT = self.CONST_SE_EOC_WEIGHT
          SE_TAU_WEIGHT = self.CONST_SE_TAU_WEIGHT
        #calculate new sigma
        new_sigI = flex.sqrt((SE_I_WEIGHT*(sigI**2)) + \
          (SE_SCALE_WEIGHT*(err_scale_factors * (flex.sum(sigI**2)/(flex.sum(err_scale_factors)+1E-9)))) + \
          (SE_TAU_WEIGHT*(err_tau * (flex.sum(sigI**2)/((err_tau*len(sigI))+1E-9)))) + \
          (SE_EOC_WEIGHT*(err_partiality * (flex.sum(sigI**2)/(flex.sum(err_partiality)+1E-9)))))
        #for i in new_sigI:
        #  if math.isnan(i) or i == float('inf') or i<0.1:
        #    print i, SE_I_WEIGHT, SE_SCALE_WEIGHT, SE_EOC_WEIGHT, SE_TAU_WEIGHT
        pres.set_params(observations=pres.observations.customized_copy(sigmas=new_sigI),
          observations_original=pres.observations_original.customized_copy(sigmas=new_sigI))
        output_results.append(pres)
        #for plotting
        if iparams.flag_plot_expert:
          plt.subplot(521)
          plt.scatter(1/(d_spacings**2), scale_factors_by_indices, s=10, marker='x', c='r')
          plt.title('Scale factors')
          plt.subplot(522)
          plt.scatter(1/(d_spacings**2), pres.partiality, s=10, marker='x', c='r')
          plt.title('Partiality')
          plt.subplot(523)
          plt.scatter(1/(d_spacings**2), err_scale_factors, s=10, marker='x', c='r')
          plt.title('Error in scale factors')
          plt.subplot(524)
          plt.scatter(1/(d_spacings**2), err_partiality, s=10, marker='x', c='r')
          plt.title('Error in partiality')
          plt.subplot(525)
          plt.scatter(1/(d_spacings**2), sigI, s=10, marker='x', c='r')
          plt.title('Original sigmas')
          plt.subplot(526)
          plt.scatter(1/(d_spacings**2), new_sigI, s=10, marker='x', c='r')
          plt.title('New sigmas')
          plt.subplot(527)
          plt.scatter(1/(d_spacings**2), flex.log(pres.observations.data()), s=10, marker='x', c='r')
          plt.title('Original I')
          plt.subplot(528)
          plt.scatter(1/(d_spacings**2), flex.log(observations_full.data()), s=10, marker='x', c='r')
          plt.title('New I')
          if miller_array_iso is not None:
            ma_iso, ma_obs_old = miller_array_iso.common_sets(observations_old, assert_is_similar_symmetry=False)
            ma_iso, ma_obs_full = miller_array_iso.common_sets(observations_full, assert_is_similar_symmetry=False)
            plt.subplot(529)
            plt.scatter(ma_obs_old.sigmas(), flex.log(flex.abs(ma_iso.data()-ma_obs_old.data())), s=10, marker='x', c='r')
            plt.title('Original SE vs Log Residual (R=%6.1f CC=%6.2f)'%(flex.sum(flex.abs(ma_iso.data()-ma_obs_old.data()))/flex.sum(flex.abs(ma_obs_old.data())), \
              flex.linear_correlation(ma_obs_old.sigmas(), flex.log(flex.abs(ma_iso.data()-ma_obs_old.data()))).coefficient()))
            plt.subplot(5,2,10)
            plt.scatter(ma_obs_full.sigmas(), flex.log(flex.abs(ma_iso.data()-ma_obs_full.data())), s=10, marker='x', c='r')
            plt.title('New SE vs Log Residual (R=%6.1f CC=%6.2f)'%(flex.sum(flex.abs(ma_iso.data()-ma_obs_full.data()))/flex.sum(flex.abs(ma_obs_full.data())), \
              flex.linear_correlation(ma_obs_full.sigmas(), flex.log(flex.abs(ma_iso.data()-ma_obs_full.data()))).coefficient()))
          plt.show()
      return output_results

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
    from scitbx.matrix import sqr
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
    """
    Constructor
    """

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

class wilson_plot_handler(object):
  """
  Take miller array and show Wilson Plot
  """
  def __init__(self):
    """
    Constructor
    """

  def show_plot(self, miller_array_in, n_bins=None):
    try:
      import matplotlib.pyplot as plt
    except Exception as e:
      print("Warning: error importing matplotlib.pyplot")
      print(e)
      return
    if n_bins is None:
      binner = miller_array_in.setup_binner(auto_binning=True)
    else:
      binner = miller_array_in.setup_binner(n_bins=n_bins)
    binner_indices = binner.bin_indices()
    avg_I_bin = flex.double()
    one_dsqr_bin = flex.double()
    for i in range(1, n_bins+1):
      i_binner = (binner_indices == i)
      I_sel = observations_original.data().select(i_binner)
      avg_I_bin.append(np.mean(I_sel))
      one_dsqr_bin.append(1/binner.bin_d_range(i)[1]**2)
    x_axis = one_dsqr_bin
    fig, ax1 = plt.subplots()
    ln1 = ax1.plot(x_axis, avg_I_bin, linestyle='-', linewidth=2.0, c='b')
    ax1.set_xlabel('1/d^2')
    ax1.set_ylabel('<I>', color='b')
    plt.grid()
    plt.show()
