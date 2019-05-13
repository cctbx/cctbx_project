from __future__ import division
from __future__ import print_function
from six.moves import cPickle as pickle
from prime.postrefine import postref_handler
from .mod_lbfgs import lbfgs_handler
import numpy as np
from cctbx import sgtbx
import random
from cctbx.array_family import flex
from prime.postrefine.mod_input import read_frame

class indamb_handler(object):
  """
  handle indexing ambiguity main
  """
  def __init__(self):
    """
    Constructor
    """

  def generate_twin_operators(self, obs_in, flag_all=False):
    #generate only true merohedral twin operators
    from mmtbx.scaling.twin_analyses import twin_laws
    TL = twin_laws(miller_array=obs_in)
    operators = []
    if flag_all:
      operators = TL.operators
    else:
      if TL.m > 0:
        operators = TL.operators
    return operators

  def generate_reindex_sets(self, obs_in):
    ops = self.generate_twin_operators(obs_in)
    alternates = {'h,k,l': obs_in}
    for op in ops:
      hkl = obs_in.indices()
      cb_op = sgtbx.change_of_basis_op(op.operator.r().as_hkl())
      hklrev = cb_op.apply(hkl)
      alternates[op.operator.r().as_hkl()] = obs_in.customized_copy(indices = hklrev).map_to_asu()
    return alternates

  def generate_forced_reindex_sets(self, obs_in, assigned_basis):
    alternates = {'h,k,l': obs_in}
    for ab in assigned_basis:
      cb_op = sgtbx.change_of_basis_op(ab)
      alternates[ab] = obs_in.change_basis(cb_op).map_to_asu()
    return alternates

  def get_observations(self, pickle_filename, iparams):
    main_obs_pickle = read_frame(pickle_filename)
    prh = postref_handler()
    avg_mode = 'average'
    try:
      inputs, txt_org = prh.organize_input(main_obs_pickle, iparams, avg_mode, pickle_filename=pickle_filename)
      main_obs = inputs[0]
    except Exception:
      print('Error reading input pickle.')
      return None
    main_asu = main_obs.map_to_asu().merge_equivalents().array()
    #get other indexing alternatives
    if iparams.indexing_ambiguity.mode == 'Forced':
      #generate other basis format according to the list
      alternates = self.generate_forced_reindex_sets(main_asu, iparams.indexing_ambiguity.assigned_basis)
    else:
      alternates = self.generate_reindex_sets(main_asu)
    return alternates

  def calc_cc(self, frame_no, pickle_filename, iparams, miller_array_ref):
    img_filename_only = ''
    pickle_filepaths = pickle_filename.split('/')
    img_filename_only = pickle_filepaths[len(pickle_filepaths)-1]
    try:
      alternates = self.get_observations(pickle_filename, iparams)
      cc_set = []
      for key in alternates.keys():
        corr = miller_array_ref.correlation(alternates[key], assert_is_similar_symmetry=False)
        cc_set.append(corr.coefficient())
      i_best = np.argmax(cc_set)
      txt_out = ' {0:40} ==> CC:{1:6.2f}'.format(img_filename_only+' '+alternates.keys()[i_best], cc_set[i_best])
      return alternates.keys()[i_best], txt_out
    except Exception:
      txt_out = ' {0:40} ==> CC:{1:6.2f}'.format(img_filename_only+' (h,k,l)', 0)
      return "h,k,l", txt_out

  def calc_r(self, frame_no, pickle_filename, index_basis, main_obs, obs_results):
    #get main observation
    img_filename_only = ''
    pickle_filepaths = pickle_filename.split('/')
    img_filename_only = pickle_filepaths[len(pickle_filepaths)-1]
    txt_out = ' {0:40} ==> '.format(img_filename_only+' '+index_basis)
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
    return r_set, txt_out

  def optimize(self, r_matrix, flag_plot=False):
    xinp = flex.double([random.random() for i in range(len(r_matrix)*2)])
    xinp_copy = xinp[:]
    args = r_matrix
    lh = lbfgs_handler(current_x=xinp, args=args)
    x_set = np.array(lh.x).reshape((len(lh.x)//2,2))
    if flag_plot:
      import matplotlib.pyplot as plt
      xinp_set = np.array(xinp_copy).reshape((len(xinp_copy)/2,2))
      ax1 = plt.subplot(2,1,1)
      ax1.scatter(xinp_set[:,0], xinp_set[:,1], s=10, marker='o', c='b')
      ax1.set_xlim([-0.2, 1.2])
      ax1.set_ylim([-0.2, 1.2])
      ax2 = plt.subplot(2,1,2)
      ax2.scatter(x_set[:,0], x_set[:,1], s=10, marker='o', c='b')
      ax2.set_xlim([-0.2, 1.2])
      ax2.set_ylim([-0.2, 1.2])
      plt.show()
    return x_set

  def assign_basis(self, frame_dup_files, basis_choices, labels, k, sample_fname):
    sol_dict = {}
    file_list_txt = ''
    for i in range(k):
      frame_cluster = np.extract(labels==i, frame_dup_files)
      basis_choice_cluster = np.extract(labels==i, basis_choices)
      for frame, basis_choice in zip(frame_cluster, basis_choice_cluster):
        if frame not in sol_dict:
          sol_dict[frame] = basis_choice
          file_list_txt += frame+'\n'
    f = open(sample_fname, 'w')
    f.write(file_list_txt)
    f.close()
    return sol_dict
