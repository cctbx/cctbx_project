# LIBTBX_SET_DISPATCHER_NAME prime.r_mtz_mtz
'''
Author      : Uervirojnangkoorn, M.
Created     : 6/29/2016
Description : Scale second mtz to the first (linear scale only) and calculate r-factors.
Note that all intensity array will be converted to amplitude.
'''
from __future__ import absolute_import, division, print_function

from iotbx import reflection_file_reader
import sys
from six.moves import range

def read_input(args):
  hkla = None
  hklb = ''
  d_min = 0
  d_max = 99
  n_bins = 20
  flag_b = False
  for i in range(len(args)):
    pair=args[i].split('=')
    if pair[0]=='hkla':
      hkla = pair[1]
    elif pair[0]=='hklb':
      hklb = pair[1]
    elif pair[0]=='d_min':
      d_min = float(pair[1])
    elif pair[0]=='d_max':
      d_max = float(pair[1])
    elif pair[0]=='n_bins':
      n_bins = int(pair[1])
    elif pair[0]=='flag_b':
      flag_b = bool(pair[1])
  if hklb == '':
    print("Please provide input hkl files.")
    exit()
  return hkla, hklb, d_min, d_max, n_bins, flag_b

def get_miller_f_from_reflection_file(hklin):
  flag_hklin_found = False
  miller_array_f = None
  if hklin is not None:
    flag_hklin_found = True
    reflection_file = reflection_file_reader.any_reflection_file(hklin)
    miller_arrays = reflection_file.as_miller_arrays()
    is_found_as_intensity_array = False
    is_found_as_amplitude_array = False
    is_found_as_complex_array = False
    for miller_array in miller_arrays:
      if miller_array.is_xray_amplitude_array():
        miller_array_f = miller_array.deep_copy()
        is_found_as_amplitude_array = True
        break
      elif miller_array.is_xray_intensity_array():
        is_found_as_intensity_array = True
        miller_array_converted_to_amplitude = miller_array.enforce_positive_amplitudes()
      elif miller_array.is_complex_array():
        is_found_as_complex_array = True
        miller_array_converted_to_amplitude = miller_array.amplitudes()
    if is_found_as_amplitude_array == False:
      if is_found_as_intensity_array or is_found_as_complex_array:
        miller_array_f = miller_array_converted_to_amplitude.deep_copy()
      else:
        flag_hklin_found = False
  return flag_hklin_found, miller_array_f


if (__name__ == "__main__"):
  #check input
  if len(sys.argv)==1:
    print('Use prime.r_mtz_mtz to calculate R-factor between two reflection files.')
    print('Only linear scale will be performed. Set flag_b=True to include B-factor scaling.')
    print('Usage: prime.r_mtz_mtz hkla=reflection1.mtz hklb=reflection2.mtz d_min=min_resolution d_max=max_resolution n_bins=no_of_bins flag_b=True_or_False.')
    exit()
  #read input parameters and frames (pickle files)
  hkla, hklb, d_min, d_max, n_bins, flag_b = read_input(args = sys.argv[1:])
  #first reflection file (reference)
  flag_hkla_found, miller_array_a = get_miller_f_from_reflection_file(hkla)
  #second reflection file
  flag_hklb_found, miller_array_b = get_miller_f_from_reflection_file(hklb)
  if flag_hkla_found and flag_hklb_found:
    #filter resolution
    miller_array_a = miller_array_a.resolution_filter(d_min=d_min, d_max=d_max)
    miller_array_b = miller_array_b.resolution_filter(d_min=d_min, d_max=d_max)
    #report
    print('First reflection file:', hkla)
    miller_array_a.show_summary()
    print('Second reflection file:', hklb)
    miller_array_b.show_summary()
    #scale b to a
    miller_array_b_scaled = miller_array_a.scale(miller_array_b, resolution_dependent=flag_b)
    #get common sets
    ma_common_a, ma_common_b = miller_array_a.common_sets(miller_array_b_scaled)
    #set up binning
    ma_common_a.setup_binner(n_bins=n_bins)
    #calculate R-factor
    r1_factor_bin = ma_common_a.r1_factor(ma_common_b, use_binning=True)
    r1_factor_bin.show()
    r1_factor = ma_common_a.r1_factor(ma_common_b, use_binning=False)
    print('Overall R-factor:', r1_factor)
