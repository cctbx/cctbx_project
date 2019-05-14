from __future__ import absolute_import, division, print_function
from prime.postrefine import postref_handler
from prime.index_ambiguity.mod_lbfgs import lbfgs_handler
from cctbx.array_family import flex
import numpy as np
import random, ntpath, os

class isoform_cluster_handler(object):
  """
  handle indexing ambiguity main
  """
  def get_observations(self, pickle_filename, iparams):
    prh = postref_handler()
    #modify reflection filters
    iparams.merge.d_min = iparams.isoform_cluster.d_min
    iparams.merge.d_max = iparams.isoform_cluster.d_max
    iparams.merge.sigma_min = iparams.isoform_cluster.sigma_min
    pres, txt_out = prh.scale_frame_by_mean_I(-1, pickle_filename, iparams, 0, 'average')
    return pres.get_full_observations().merge_equivalents().array() if pres else None

  def calc_cc(self, frame_no, pickle_filename, main_obs, iparams, miller_array_ref_set):
    img_filename_only = ntpath.basename(pickle_filename)
    cc_set = []
    for miller_array_ref in miller_array_ref_set:
      corr = main_obs.correlation(miller_array_ref, assert_is_similar_symmetry=False)
      cc_set.append(corr.coefficient())
    i_best = np.argmax(cc_set)
    txt_out = ' {0:40} {1:2d} ==> CC:{2:6.2f} '.format(img_filename_only, i_best, cc_set[i_best])
    return i_best, txt_out

  def calc_r(self, frame_no, pickle_filename, main_obs, obs_results):
    #get main observation
    img_filename_only = ntpath.basename(pickle_filename)
    txt_out = ' {0:40} ==> '.format(img_filename_only)
    #calculate cc with other observations
    n_frames = len(obs_results)
    r_set = np.zeros((1, n_frames))
    n_refl_common_set = np.zeros((1, n_frames))
    n_frame_common = 0
    for j in range(frame_no+1, n_frames):
      obs_asu = obs_results[j]
      corr = main_obs.correlation(obs_asu, assert_is_similar_symmetry=False)
      c_main, c_other = main_obs.common_sets(obs_asu, assert_is_similar_symmetry=False)
      if len(c_main.indices()) > 5:
        r_set[0,j] = corr.coefficient()
        n_refl_common_set[0,j] = len(c_main.indices())
        n_frame_common += 1
    txt_out += ' <CC>:%6.2f <N_refl_common>: %6.1f N_frame_common: %6.0f'%(np.mean(r_set), np.mean(n_refl_common_set), n_frame_common)
    return r_set, main_obs, txt_out

  def optimize(self, r_matrix, flag_plot=False):
    xinp = flex.double([random.random() for i in range(len(r_matrix)*2)])
    xinp_copy = xinp[:]
    args = r_matrix
    lh = lbfgs_handler(current_x=xinp, args=args)
    x_set = np.array(lh.x).reshape((len(lh.x)/2,2))
    if flag_plot:
      import matplotlib.pyplot as plt
      xinp_set = np.array(xinp_copy).reshape((len(xinp_copy)/2,2))
      plt.subplot(2,1,1)
      plt.scatter(xinp_set[:,0], xinp_set[:,1], s=10, marker='x', c='r')
      plt.subplot(2,1,2)
      plt.scatter(x_set[:,0], x_set[:,1], s=10, marker='x', c='r')
      plt.show()
    return x_set

  def assign_cluster(self, frame_files, labels, k, output_path):
    sol_dict = {}
    cluster_files = [os.path.join(output_path,'cluster_'+str(i)+'.lst') for i in range(k)]
    for i in range(k):
      file_list_txt = ''
      frame_cluster = np.extract(labels==i, frame_files)
      for frame in frame_cluster:
        if frame not in sol_dict:
          sol_dict[frame] = i
          file_list_txt += frame+'\n'
      with open(cluster_files[i], 'w') as f: f.write(file_list_txt)
    return sol_dict, cluster_files
