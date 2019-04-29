'''
Author      : Uervirojnangkoorn, M.
Created     : 7/13/2014
Detecting outliers with Wilson statistics (Read, 1999)
use the normalized intensity (E) limits (acentric 3.72, centric 4.89)
'''
from __future__ import division, print_function
from cctbx.array_family import flex
from cctbx import statistics

class outlier_handler(object):
  '''
  A wrapper class for finding outliers
  '''

  def __init__(self):
    '''
    Constructure
    '''

  def good_i_flags(self, miller_array_i, iparams, flag_show_summary=False):
    miller_array_f = miller_array_i.f_sq_as_f()
    good_i_flags = flex.bool([True]*len(miller_array_f.data()))

    try:
      miller_array_f.setup_binner(auto_binning=True)
      miller_array_e = miller_array_f.quasi_normalize_structure_factors()
      i_seq = flex.int([i for i in range(len(miller_array_e.data()))])

      e_lim_acentric = 3.72
      e_lim_centric = 4.89
      centric_flags_data = miller_array_e.centric_flags().data()
      data_e_acentric = miller_array_e.data().select(centric_flags_data==False)
      i_seq_acentric = i_seq.select(centric_flags_data==False)
      data_e_centric = miller_array_e.data().select(centric_flags_data==True)
      i_seq_centric = i_seq.select(centric_flags_data==True)

      i_seq_acentric_outliers = i_seq_acentric.select(data_e_acentric>=e_lim_acentric)
      i_seq_centric_outliers = i_seq_centric.select(data_e_centric>=e_lim_centric)

      i_seq_outliers = sorted(i_seq_acentric_outliers.concatenate(i_seq_centric_outliers))

      for i in i_seq_outliers:
        good_i_flags[i] = False

      if flag_show_summary:
        print('Acentric outliers:')
        for i in i_seq_acentric_outliers:
          print(miller_array_e.indices()[i], miller_array_e.data()[i])

        print('Centric outliers:')
        for i in i_seq_centric_outliers:
          print(miller_array_e.indices()[i], miller_array_e.data()[i])
    except Exception:
      dummy = 0

    return good_i_flags
