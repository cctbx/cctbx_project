from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from cctbx import maptbx

def split_into_hybrid_resolution_bins(miller_array, reflections_per_bin=100):
    """
    Takes a miller.array and splits it into resolution bins:
      - 4A to infinity: log binning
      - d_min to 4A: standard binning with a specified number of reflections per bin

    Returns a list of flex.bool arrays, where each flex.bool is a
    selection array for a single bin, matching the size and order of the original array.
    """
    d_spacings = miller_array.d_spacings().data()
    total_size = miller_array.size()
    # Define boolean selections for the two regimes
    sel_low_res = (d_spacings >= 4.0)
    sel_high_res = (d_spacings < 4.0)
    # Keep track of all combined bin selections
    all_bin_selections = []
    #
    # 1. Log Binning for 4A to Infinity
    #
    ma_low = miller_array.select(sel_low_res)
    if ma_low.size() > 0:
      # log_binning() returns a list of flex.bool arrays relative to ma_low
      low_bins_rel = ma_low.log_binning(
        n_reflections_in_lowest_resolution_bin = 100,
        max_number_of_bins = 99999,
        min_reflections_in_bin = 100)
      # Map the relative selections back to the original full-sized array
      indices_low = sel_low_res.iselection()
      for sel_rel in low_bins_rel:
        sel_full = flex.bool(total_size, False)
        selected_indices = indices_low.select(sel_rel)
        sel_full.set_selected(selected_indices, True)
        all_bin_selections.append(sel_full)
    #
    # 2. Standard Binning for d_min to 4A (Reflections Per Bin)
    #
    ma_high = miller_array.select(sel_high_res)
    if ma_high.size() > 0:
      # setup_binner creates a binner using the target reflections per bin
      ma_high.setup_binner(reflections_per_bin=reflections_per_bin)
      binner = ma_high.binner()
      # Map the binner selections back to the original full-sized array
      indices_high = sel_high_res.iselection()
      for i_bin in binner.range_used():
        sel_rel = binner.selection(i_bin)
        sel_full = flex.bool(total_size, False)
        selected_indices = indices_high.select(sel_rel)
        sel_full.set_selected(selected_indices, True)
        all_bin_selections.append(sel_full)
    return all_bin_selections

"""
Implementation of difference map calculation following the exact description in:

Yamashita, K., Palmer, C. M., Burnley, T. & Murshudov, G. N. (2021).
Acta Crystallographica Section D, 77, 1282–1291.

Hugely inefficient in favor of maximal transparency code <> paper. Though still
rather fast.

Actual implementation seems to use mask internally. They might be doing some
sharpening too. Both hypotheses are inspired by looking at the log file from
running Servalcat.

reflections_per_bin needs to be automated.
"""

def formula_3(fo1, fo2):
  x = fo1 - fo2
  diff = x - flex.mean(x)
  return flex.sum_sq(diff).real / diff.size() / 4

def formula_8(fo, fc):
  num = flex.sum( fo * flex.conj(fc) )
  den = flex.sum_sq(fc).real  # = sum |Fc|^2
  assert abs(den)>1.e-9
  return (num.real) / den   # should be real; use real() for safety

def formula_9(fo, fc, D, sig_n_sq):
  diff = fo - D * fc          # complex residual Fo - D Fc
  # Natively sum the squared moduli of the complex differences
  var_U_T = (flex.sum_sq(diff).real / diff.size()) - sig_n_sq
  return max(0.0, var_U_T)

def formula_4(fo1, fo2):
  return maptbx.cc_complex_complex(f_1 = fo1, f_2 = fo2)

def formula_5(FSChalf):
  return 2 * FSChalf / (FSChalf + 1)

def formula_13(sig_UT_sq, sig_n_sq):
  return sig_UT_sq / (sig_UT_sq + sig_n_sq)

def formula_17(FSCfull, fo, fc, w, D):
  num = w * (fo - D * fc)
  if FSCfull < 0 or abs(FSCfull) < 0.01: return num*0.0 # XXXX NEWNEWNEW
  den = ( FSCfull * flex.sum_sq(fo).real / fo.size() )**0.5
  return num/den

def one_bin(fo, fo1, fo2, fc):
  sig_n_sq  = formula_3(fo1 = fo1, fo2 = fo2)
  D         = formula_8(fo = fo, fc = fc)
  sig_UT_sq = formula_9(fo = fo, fc = fc, D = D, sig_n_sq = sig_n_sq)
  FSChalf   = formula_4(fo1 = fo1, fo2 = fo2)
  FSCfull   = formula_5(FSChalf = FSChalf)
  w         = formula_13(sig_UT_sq = sig_UT_sq, sig_n_sq = sig_n_sq)
  Fdiff     = formula_17(FSCfull = FSCfull, fo = fo, fc = fc, w = w, D = D)
  return Fdiff

def compute(fo1, fo2, fc, reflections_per_bin = 5000):
  fo1d = fo1.data()
  fo2d = fo2.data()
  fod  = (fo1d + fo2d) / 2.
  fcd  = fc.data()
  result = flex.complex_double(fod.size(), 0)

  ds = fo1.d_spacings()
  all_bins = split_into_hybrid_resolution_bins(miller_array=ds,
    reflections_per_bin=reflections_per_bin)

  #fo1.setup_binner(reflections_per_bin = reflections_per_bin)
  #ds = fo1.d_spacings().data()
  #for i_bin in fo1.binner().range_used():
  #  sel = fo1.binner().selection(i_bin)
  for sel in all_bins:
    Fdiff = one_bin(
      fo  = fod.select(sel),
      fo1 = fo1d.select(sel),
      fo2 = fo2d.select(sel),
      fc  = fcd.select(sel))
    result.set_selected(sel, Fdiff)
  return fo1.array(data = result)
