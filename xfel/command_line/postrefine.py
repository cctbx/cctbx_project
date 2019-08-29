# LIBTBX_SET_DISPATCHER_NAME cxi.postrefine
from __future__ import absolute_import, division, print_function
import os, sys
from libtbx.easy_mp import pool_map
import numpy as np
from cctbx.array_family import flex
from datetime import datetime, time
import logging
from six.moves import range

FORMAT = '%(levelname)s %(module)s.%(funcName)s: %(message)s'
logging.basicConfig(level=logging.WARNING, format=FORMAT)

def determine_mean_I_mproc(frame_no, frame_files, iph):
  from xfel.cxi.postrefine import postref_handler
  prh = postref_handler()
  mean_I = prh.calc_mean_intensity(frame_files[frame_no], iph)
  return mean_I

def scale_frame_by_mean_I_mproc(frame_no, frame_files, iph, mean_of_mean_I):
  from xfel.cxi.postrefine import postref_handler
  prh = postref_handler()
  pres = prh.scale_frame_by_mean_I(frame_no,frame_files[frame_no], iph, mean_of_mean_I)
  return pres

def postrefine_by_frame_mproc(frame_no, frame_files, iph, miller_array_ref):
  from xfel.cxi.postrefine import postref_handler
  prh = postref_handler()
  pres = prh.postrefine_by_frame(frame_no, frame_files[frame_no], iph, miller_array_ref)
  return pres

def read_input(args):
  file_name_input = ''
  frame_files = []
  for i in range(len(args)):
    pair=args[i].split('=')
    if pair[0]=='input':
      file_name_input = pair[1]

    #other args are considered as the pickle directory
    if len(pair) == 1:
      if pair[0].endswith('.pickle'):
        frame_files.append(pair[0])

  if file_name_input == '':
    print("Please provide input-parameters file (usage: input=yourinput.inp)")
    exit()

  from xfel.cxi.postrefine import postref_handler
  prh = postref_handler()
  iph = prh.read_input_parameters(file_name_input)

  #make run_no folder
  if not os.path.exists(iph.run_no):
    os.makedirs(iph.run_no)

  #check if pickle_dir is given in input file instead of from cmd arguments.
  if len(frame_files)==0:
    print('Path to pickle files is missing, please specify it at command line.')
    print('Usage: cxi.postrefine input=myinp.inp /path/to/pickles/*')
    exit()

  return iph, frame_files


if (__name__ == "__main__"):
  logging.info("Starting process.")
  #capture starting time
  time_global_start=datetime.now()

  #0 .read input parameters and frames (pickle files)
  iph, frame_files = read_input(args = sys.argv[1:])
  frames = list(range(len(frame_files)))

  #1. prepare reference miller array
  if iph.file_name_ref_mtz == '':
    #if iso. ref. is not given, use the <I> to scale each frame.

    #calculate mean intensity for each frame and determine median of the
    #distribution
    def determine_mean_I_mproc_wrapper(arg):
      return determine_mean_I_mproc(arg, frame_files, iph)

    determine_mean_I_result = pool_map(
            args=frames,
            func=determine_mean_I_mproc_wrapper,
            processes=None)

    frames_mean_I = flex.double()
    for result in determine_mean_I_result:
      if result is not None:
        frames_mean_I.append(result)

    mean_of_mean_I = np.median(frames_mean_I)

    #use the calculate <mean_I> to scale each frame
    def scale_frame_by_mean_I_mproc_wrapper(arg):
      return scale_frame_by_mean_I_mproc(arg, frame_files, iph, mean_of_mean_I)

    scale_frame_by_mean_I_result = pool_map(
            args=frames,
            func=scale_frame_by_mean_I_mproc_wrapper,
            processes=None)

    observations_merge_mean_set = []
    for pres in scale_frame_by_mean_I_result:
      if pres is not None:
        observations_merge_mean_set.append(pres)

    if len(observations_merge_mean_set) > 0:
      from xfel.cxi.postrefine import merge_observations
      miller_array_merge_mean, txt_merge_mean, csv_merge_mean = merge_observations(observations_merge_mean_set, iph, iph.run_no+'/mean_scaled','average')
      miller_array_ref = miller_array_merge_mean.expand_to_p1().generate_bijvoet_mates()
      if iph.flag_force_no_postrefine:
        txt_out = iph.txt_out + txt_merge_mean
        f = open(iph.run_no+'/log.txt', 'w')
        f.write(txt_out)
        f.close()
        exit()
    else:
      print("No frames merged as a reference set - exit without post-refinement")
      exit()
  else:
    miller_array_ref = iph.miller_array_iso
    txt_merge_mean = ''

  #2. Post-refinement
  n_iters = iph.n_postref_cycle
  txt_merge_postref = ''
  for i in range(n_iters):
    txt_merge_postref += 'Start post-refinement cycle '+str(i+1)+'\n'
    print(txt_merge_postref)
    def postrefine_by_frame_mproc_wrapper(arg):
      return postrefine_by_frame_mproc(arg, frame_files, iph, miller_array_ref)

    postrefine_by_frame_result = pool_map(
            args=frames,
            func=postrefine_by_frame_mproc_wrapper,
            processes=None)

    postrefine_by_frame_good = []
    for pres in postrefine_by_frame_result:
      if pres is not None:
        postrefine_by_frame_good.append(pres)

    if len(postrefine_by_frame_good) > 0:
      from xfel.cxi.postrefine import merge_observations
      miller_array_merge_postref, txt_merge_out, csv_merge_PR = merge_observations(postrefine_by_frame_good, iph, iph.run_no+'/postref_cycle_'+str(i+1),'weighted')
      miller_array_ref = miller_array_merge_postref.generate_bijvoet_mates()
      txt_merge_postref += txt_merge_out
    else:
      print("No frames merged as a reference set - exit without post-refinement")
      exit()

  #collect caculating time
  time_global_end=datetime.now()
  time_global_spent=time_global_end-time_global_start
  txt_out_time_spent = 'Total calculation time: '+'{0:.2f}'.format(time_global_spent.seconds)+' seconds\n'
  print(txt_out_time_spent)

  txt_out = iph.txt_out + txt_merge_mean + txt_merge_postref + txt_out_time_spent
  f = open(iph.run_no+'/log.txt', 'w')
  f.write(txt_out)
  f.close()

  with open("{}/mean_stats.csv".format(iph.run_no), 'wb') as mean_stats:
    mean_stats.write(csv_merge_mean)

  with open("{}/PostRef_stats.csv".format(iph.run_no), 'wb') as PR_stats:
    PR_stats.write(csv_merge_PR)

