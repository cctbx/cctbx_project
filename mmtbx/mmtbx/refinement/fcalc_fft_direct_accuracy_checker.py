from cctbx.array_family import flex
import math, time
from cctbx import miller
from cctbx import crystal
from cctbx import adptbx
from libtbx.test_utils import approx_equal
from mmtbx import bulk_solvent
from cctbx import xray
from mmtbx.max_lik import maxlik
from cctbx import crystal

def fcalc_fft_direct_accuracy_checker(xray_structure,
                                      f_calc,
                                      nref_in_bin,
                                      r_lim,
                                      r_lim_high,
                                      scat_table):
  flag = "OK"
  xray_structure_copy =  xray_structure.deep_copy_scatterers()
  #xray_structure_copy.scattering_dict(table = scat_table)

  f_calc_fft = f_calc.structure_factors_from_scatterers(
                                 xray_structure = xray_structure_copy,
                                 algorithm      = "fft").f_calc()

  assert approx_equal(bulk_solvent.r_factor(flex.abs(f_calc_fft.data()),
                                            f_calc.data()), 0.0)

  f_calc_dir = f_calc.structure_factors_from_scatterers(
                                 xray_structure = xray_structure_copy,
                                 algorithm      = "direct").f_calc()

  f_calc_fft = abs(f_calc_fft)
  f_calc_dir =     f_calc_dir

  assert f_calc_fft.size() == f_calc_dir.size()

  f_calc_fft.setup_binner(reflections_per_bin = nref_in_bin)
  f_calc_dir.use_binning_of(f_calc_fft)
  for i_bin in f_calc_fft.binner().range_used():
    sel = f_calc_fft.binner().selection(i_bin)
    sel_f_calc_fft = f_calc_fft.select(sel)
    sel_f_calc_dir = f_calc_dir.select(sel)
    r = bulk_solvent.r_factor(sel_f_calc_fft.data(),sel_f_calc_dir.data())*100.
    if(r > r_lim and f_calc_dir.binner().bin_d_range(i_bin)[1] > 1.0):
      flag = "STOP"
      print "%s nref= %7d r-work= %8.4f r-lim= %5.2f" % \
             (f_calc_fft.binner().bin_legend(i_bin),sel_f_calc_fft.data().size(),r,r_lim)
    if(r > r_lim_high and f_calc_dir.binner().bin_d_range(i_bin)[1] <= 1.0):
      flag = "STOP"
      print "%s nref= %7d r-work= %8.4f r-lim= %5.2f" % \
             (f_calc_fft.binner().bin_legend(i_bin),sel_f_calc_fft.data().size(),r,r_lim_high)
  if(flag == "STOP"):
    print
    print "Problem: "
    print "    Insufficient accuracy of structure factors calculation for"
    print "    resolution zones listed above. Please change the parameters"
    print "    of structure factors generation in order to assert requested"
    print "    accuracy or make the accuracy criterion less strict. "
    print
    print "Requested accuracy:"
    print "    r-factor < %5.2f percent "%(1.0),"in any resolution zone > 1.0A"
    print "    r-factor < %5.2f percent "%(0.5),"in any resolution zone < 1.0A"
    print
    print "Approximate number of reflections requested for each zone: %5d"%(100)
    print
  return flag
