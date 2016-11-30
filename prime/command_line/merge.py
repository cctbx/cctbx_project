from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME prime.merge
"""
Author      : Uervirojnangkoorn, M.
Created     : 7/14/2016
Description : prime.merge checks if run_no/pickles are found. If so, merge all
reflections and write out mtz file.
"""
import os, sys, math
from cctbx.array_family import flex
import cPickle as pickle
from prime.postrefine import postref_handler
from prime.postrefine.mod_partiality import partiality_handler

class merge_handler(object):
  """
  A wrapper class for the handler
  """
  def __init__(self):
    """
    Intialitze parameters
    """

  def run(self, args, avg_mode='average'):
    #read inputs
    from prime.postrefine.mod_input import process_input
    iparams, txt_out_input = process_input(argv=args, flag_check_exist=False)
    self.run_by_params(iparams, avg_mode=avg_mode)

  def run_by_params(self, iparams, avg_mode='average'):
    #read all result pickles
    try:
      DIR = iparams.run_no+'/pickles/'
      pickle_results = [pickle.load(open(DIR+fname, "rb")) for fname in os.listdir(DIR)]
      n_results = len(pickle_results)
      #for the last cycle merging with weak anomalous flag on
      #fetch the original observations from the integration pickles
      if iparams.flag_weak_anomalous and iparams.target_anomalous_flag and avg_mode=='final':
        prh = postref_handler()
        for pres in pickle_results:
          if pres is not None:
            observations_pickle = pickle.load(open(pres.pickle_filename,'rb'))
            inputs, _ = prh.organize_input(observations_pickle, iparams, avg_mode, pickle_filename=pres.pickle_filename)
            observations_original, alpha_angle_obs, spot_pred_x_mm, spot_pred_y_mm, detector_distance_mm, identified_isoform, mapped_predictions, xbeam, ybeam = inputs
            two_theta = observations_original.two_theta(wavelength=pres.wavelength).data()
            ph = partiality_handler()
            partiality_set, delta_xy_set, rs_set, rh_set = ph.calc_partiality_anisotropy_set(
              pres.unit_cell, pres.rotx, pres.roty, observations_original.indices(),
              pres.ry, pres.rz, pres.r0, pres.re, pres.voigt_nu,
              two_theta, alpha_angle_obs, pres.wavelength, pres.crystal_orientation,
              spot_pred_x_mm, spot_pred_y_mm, detector_distance_mm,
              iparams.partiality_model, iparams.flag_beam_divergence)
            observations, alt_hkl = prh.get_observations_non_polar(observations_original, pres.pickle_filename, iparams)
            pres.observations = observations
            pres.observations_original = observations_original
            pres.partiality = partiality_set
            pres.rs_set = rs_set
            pres.rh_set = rh_set
            pres.mapped_predictions = mapped_predictions
    except Exception:
      print "Error reading input pickles."
      print "Check if prime.run, prime.genref, or prime.postrefine was called prior to merge."
      exit()
    from prime.postrefine import prepare_output
    prep_output = prepare_output(pickle_results, iparams, avg_mode)
    txt_merge_mean = 'prime.merge\n'
    if prep_output is not None:
      #grab linear lists
      cn_group, group_id_list, miller_indices_all_sort, miller_indices_ori_all_sort,  I_obs_all_sort, \
      sigI_obs_all_sort, G_all_sort, B_all_sort, p_all_sort, rs_all_sort, wavelength_all_sort, \
      sin_sq_all_sort, SE_all_sort, uc_mean, wavelength_mean, pickle_filename_all_sort, txt_prep_out = prep_output
      #start averaging
      from prime.postrefine import calc_avg_I_cpp
      calc_average_I_result = calc_avg_I_cpp(cn_group, group_id_list, miller_indices_all_sort,
                                                 miller_indices_ori_all_sort, I_obs_all_sort,
                                                 sigI_obs_all_sort, G_all_sort, B_all_sort, p_all_sort,
                                                 rs_all_sort, wavelength_all_sort, sin_sq_all_sort,
                                                 SE_all_sort, avg_mode, iparams, pickle_filename_all_sort)
      #grab average results
      miller_index, I_merge, sigI_merge, stats, I_two_halves_tuple, txt_obs_out, txt_reject_out = calc_average_I_result
      I_even, I_odd, I_even_h, I_odd_h, I_even_k, I_odd_k, I_even_l, I_odd_l = I_two_halves_tuple
      #remove stat items with nan
      r_meas_w_top, r_meas_w_btm, r_meas_top, r_meas_btm, multiplicity = stats
      sel = flex.bool([math.isnan(r_top) or math.isinf(r_top)\
        or math.isnan(r_btm) or math.isinf(r_btm) for r_top, r_btm in zip(r_meas_w_top, r_meas_w_btm)])
      inverse_sel = ~sel
      inverse_isel = inverse_sel.iselection()
      miller_index = miller_index.select(inverse_isel)
      I_merge = I_merge.select(inverse_isel)
      sigI_merge = sigI_merge.select(inverse_isel)
      I_even = I_even.select(inverse_isel)
      I_odd = I_odd.select(inverse_isel)
      I_even_h = I_even_h.select(inverse_isel)
      I_odd_h = I_odd_h.select(inverse_isel)
      I_even_k = I_even_k.select(inverse_isel)
      I_odd_k = I_odd_k.select(inverse_isel)
      I_even_l = I_even_l.select(inverse_isel)
      I_odd_l = I_odd_l.select(inverse_isel)
      stat_all = (r_meas_w_top.select(inverse_isel), r_meas_w_btm.select(inverse_isel), \
        r_meas_top.select(inverse_isel), r_meas_btm.select(inverse_isel), multiplicity.select(inverse_isel))
      #write out rejected reflections
      f = open(iparams.run_no+'/rejections.txt', 'w')
      f.write(txt_reject_out)
      f.close()
      #get the latest no. of mtz file and write out the reflections
      DIR = iparams.run_no+'/mtz/'
      file_no = len([int(fname.split('.')[0]) for fname in os.listdir(DIR) if os.path.join(DIR, fname).endswith('mtz')])
      from prime.postrefine import write_output
      miller_array_ref, txt_merge_mean_table = write_output(miller_index,
                                                            I_merge, sigI_merge,
                                                            stat_all, (I_even, I_odd, I_even_h, I_odd_h,
                                                            I_even_k, I_odd_k, I_even_l, I_odd_l),
                                                            iparams, uc_mean,
                                                            wavelength_mean,
                                                            file_no,
                                                            avg_mode)
      txt_merge_mean +=  txt_merge_mean_table + txt_prep_out
      print txt_merge_mean
      f = open(iparams.run_no+'/log.txt', 'a')
      f.write(txt_merge_mean)
      f.close()

if __name__=="__main__":
  mh = merge_handler()
  mh.run(sys.argv[1:] if len(sys.argv)>1 else None)
