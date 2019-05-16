# LIBTBX_SET_DISPATCHER_NAME prime.postrefine
'''
Author      : Uervirojnangkoorn, M.
Created     : 7/13/2014
Description : Main commandline for prime.
'''
from __future__ import absolute_import, division, print_function

from libtbx.easy_mp import parallel_map
from cctbx.array_family import flex
from prime.command_line.solve_indexing_ambiguity import indexing_ambiguity_handler
from prime.postrefine.mod_mx import mx_handler
from prime.postrefine.mod_input import read_pickles
from prime.postrefine.mod_util import intensities_scaler
import os, sys, math
import numpy as np
from datetime import datetime, time
from libtbx.utils import Usage
from six.moves import range
from six.moves import zip

def determine_mean_I_mproc(args):
  from prime.postrefine import postref_handler
  prh = postref_handler()
  frame_file, iparams, avg_mode = args
  mean_I = prh.calc_mean_intensity(frame_file, iparams, avg_mode)
  return mean_I

def scale_frame_by_mean_I_mproc(args):
  from prime.postrefine import postref_handler
  prh = postref_handler()
  frame_no, frame_file, iparams, mean_of_mean_I, avg_mode = args
  pres = prh.scale_frame_by_mean_I(frame_no, frame_file, iparams, mean_of_mean_I, avg_mode)
  return pres

def postrefine_by_frame_mproc(args):
  from prime.postrefine import postref_handler
  prh = postref_handler()
  frame_no, frame_file, iparams, miller_array_ref, pres_in, avg_mode = args
  pres = prh.postrefine_by_frame(frame_no, frame_file, iparams, miller_array_ref, pres_in, avg_mode)
  return pres

def scale_frames(frames, frame_files, iparams):
  """scale frames"""
  avg_mode = 'average'
  if iparams.flag_apply_b_by_frame:
    mean_of_mean_I = 0
  else:
    #Calculate <I> for each frame
    frame_args = [(frame_file, iparams, avg_mode) for frame_file in frame_files]
    determine_mean_I_result = parallel_map(
            iterable=frame_args,
            func=determine_mean_I_mproc,
            processes=iparams.n_processors)
    frames_mean_I = flex.double()
    for result in determine_mean_I_result:
      if result is not None:
        mean_I, txt_out_result = result
        if mean_I is not None:
          frames_mean_I.append(mean_I)
    mean_of_mean_I = np.median(frames_mean_I)
  #use the calculate <mean_I> to scale each frame
  frame_args = [(frame_no, frame_file, iparams, mean_of_mean_I, avg_mode) for frame_no, frame_file in zip(frames, frame_files)]
  scale_frame_by_mean_I_result = parallel_map(
          iterable=frame_args,
          func=scale_frame_by_mean_I_mproc,
          processes=iparams.n_processors)
  observations_merge_mean_set = []
  for result in scale_frame_by_mean_I_result:
    if result is not None:
      pres, txt_out_result = result
      if pres is not None:
        observations_merge_mean_set.append(pres)
  return observations_merge_mean_set

def merge_frames(pres_set, iparams, avg_mode='average', mtz_out_prefix='mean_scaled'):
  """merge frames using average as the default"""
  miller_array_ref, txt_out = (None,'')
  intscal = intensities_scaler()
  if pres_set:
    prep_output = intscal.prepare_output(pres_set, iparams, avg_mode)
    if prep_output:
      mdh, _, txt_out_rejection  = intscal.calc_avg_I_cpp(prep_output, iparams, avg_mode)
      #select only indices with non-Inf non-Nan stats
      selections = flex.bool([False if (math.isnan(r0) or math.isinf(r0) or math.isnan(r1) or math.isinf(r1)) else True for r0, r1  in zip(mdh.r_meas_div, mdh.r_meas_divisor)])
      mdh.reduce_by_selection(selections)

      #handle rejected reflections
      rejections = {}
      for reject in txt_out_rejection.split('\n'):
        data = reject.split()
        if data:
          if not data[0] in rejections:
            rejections[data[0]] = flex.miller_index()
          rejections[data[0]].append(tuple([int(_d) for _d in data[1:4]]))

      if len(rejections) > 0:
        if not iparams.rejections:
          iparams.rejections = {}
        iparams.rejections.update(rejections)

      #merge all good indices
      mdh, txt_merge_mean_table = intscal.write_output(mdh,
          iparams, iparams.run_no+'/'+mtz_out_prefix, avg_mode)
      print(txt_merge_mean_table)
      print(prep_output[-1])
      txt_out = txt_merge_mean_table + prep_output[-1]
  return mdh, txt_out

def postrefine_frames(i_iter, frames, frame_files, iparams, pres_set, miller_array_ref, avg_mode):
  """postrefine given frames and previous postrefinement results"""
  miller_array_ref = miller_array_ref.generate_bijvoet_mates()
  txt_merge_postref = 'Post-refinement cycle '+str(i_iter+1)+' ('+avg_mode+')\n'
  txt_merge_postref += ' * R and CC show percent change.\n'
  print(txt_merge_postref)
  frame_args = [(frame_no, frame_file, iparams, miller_array_ref, pres_in, avg_mode) for frame_no, frame_file, pres_in in zip(frames, frame_files, pres_set)]
  postrefine_by_frame_result = parallel_map(
      iterable=frame_args,
      func=postrefine_by_frame_mproc,
      processes=iparams.n_processors)
  postrefine_by_frame_good = []
  postrefine_by_frame_pres_list = []
  for results in postrefine_by_frame_result:
    if results is not None:
      pres, txt_out_result = results
      postrefine_by_frame_pres_list.append(pres)
      if pres is not None:
        postrefine_by_frame_good.append(pres)
    else:
      postrefine_by_frame_pres_list.append(None)
  return postrefine_by_frame_good, postrefine_by_frame_pres_list, txt_merge_postref

def run(argv):
  #capture starting time
  time_global_start=datetime.now()
  import logging
  logging.captureWarnings(True)
  formatter = logging.Formatter('%(asctime)s\t%(levelname)s\t%(message)s')
  console_handler = logging.StreamHandler()
  console_handler.setLevel(logging.ERROR)
  console_handler.setFormatter(formatter)
  logging.getLogger().addHandler(console_handler)
  logging.getLogger('py.warnings').addHandler(console_handler)
  logging.basicConfig(format='%(asctime)s\t%(levelname)s\t%(message)s', level=logging.DEBUG)
  #0.1 determine indexing ambiguity and setup iparams
  txt_indexing_ambiguity = "Determine if there is an indexing ambiguity on the dataset"
  print(txt_indexing_ambiguity)
  idah = indexing_ambiguity_handler()
  sol_pickle, iparams = idah.run(argv)
  if sol_pickle is None:
    print("No ambiguity.")
    txt_indexing_ambiguity += "\nNo ambiguity."
  else:
    print("Ambiguity is solved. Solution file was saved to result folder.")
    txt_indexing_ambiguity += "Ambiguity is solved. Solution file was saved to result folder."
    iparams.indexing_ambiguity.index_basis_in = sol_pickle
  #0.2 setup parameters
  iparams.flag_volume_correction = False
  if iparams.partiality_model == "Lognormal":
    iparams.voigt_nu = 0.008 #use voigt_nu as lognpdf zero parameter
  #0.3 read frames
  frame_files = read_pickles(iparams.data)
  frames = range(len(frame_files))
  #1. prepare reference miller array
  txt_merge_mean = 'Generating a reference set (will not be used if hklrefin is set)'
  print(txt_merge_mean)
  #Always generate the mean-intensity scaled set.
  scaled_pres_set = scale_frames(frames, frame_files, iparams)
  mdh, _txt_merge_mean = merge_frames(scaled_pres_set, iparams)
  miller_array_ref = mdh.miller_array_merge
  txt_merge_mean += '\n' + _txt_merge_mean
  if not iparams.n_postref_cycle:
    with open(iparams.run_no+'/log.txt', 'a') as f:
      f.write(txt_indexing_ambiguity + txt_merge_mean)
    raise Usage("No. of post-refinement cycle was set to 0. Exit without post-refinement.")
  if iparams.hklrefin is not None:
    mxh = mx_handler()
    _, miller_array_ref = mxh.get_miller_array_from_reflection_file(iparams.hklrefin)
  if miller_array_ref is None:
    raise Usage("Problem with the assigned reference set. Try setting hklrefin=None and rerun the program.")
  #2. Post-refinement
  txt_merge_postref = ''
  postref_pres_set = [None]*len(frames)
  avg_mode = 'weighted'
  for i_iter in range(iparams.n_postref_cycle):
    if i_iter == (iparams.n_postref_cycle-1): avg_mode = 'final'
    postref_good_pres_set, postref_pres_set, _txt_merge_postref = postrefine_frames(i_iter, frames, frame_files, iparams, postref_pres_set, miller_array_ref, avg_mode)
    if postref_good_pres_set:
      mdh, _txt_merge_postref = merge_frames(postref_good_pres_set, iparams,
          avg_mode=avg_mode, mtz_out_prefix='postref_cycle_'+str(i_iter+1))
      miller_array_ref = mdh.miller_array_merge
      txt_merge_postref += _txt_merge_postref
    else:
      raise Usage("Problem with post-refinement. No images refined. Please check your input file.")
  #3. collect caculating time
  time_global_end=datetime.now()
  time_global_spent=time_global_end-time_global_start
  txt_out_time_spent = 'Total calculation time: '+'{0:.2f}'.format(time_global_spent.seconds)+ \
      ' seconds\nFinished: '+time_global_end.strftime("%A %d. %B %Y %H:%M:%S")+'\n'
  print(txt_out_time_spent)
  txt_out = txt_indexing_ambiguity + txt_merge_mean + txt_merge_postref + txt_out_time_spent
  with open(os.path.join(iparams.run_no,'log.txt'), 'a') as f:
    f.write(txt_out)
  with open(os.path.join(iparams.run_no,'.done'), 'w') as f:
    f.write('Done')
  return mdh

if (__name__ == "__main__"):
  run(sys.argv[1:] if len(sys.argv) > 1 else None)
