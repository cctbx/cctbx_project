from __future__ import absolute_import, division, print_function
import time, random
import numpy as np
from cctbx.development import random_structure
from cctbx.sgtbx import space_group_info
from cctbx.array_family import flex
#from cctbx.eltbx.discamb import tst_IAM_discamb_vs_cctbx as tst_disc
import pydiscamb
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
from matplotlib import pyplot as plt
import numpy as np

# ===========================================================================

tolerance_values = {
  "electron": {'score_max': 0.0003,
               'max_diff_max': 0.004},
  "it1992":   {'score_max': 0.0004,
               'max_diff_max': 0.008},
  "wk1995":   {'score_max': 0.0004,
               'max_diff_max': 0.008},
}

# ===========================================================================

def compare_structure_factors(x,y):
  x = flex.abs(x)
  y = flex.abs(y)
  scale = flex.sum(x*y)/flex.sum(y*y)
  num = flex.sum(flex.abs(x-scale*y))
  den = flex.sum(flex.abs(x+scale*y))
  diff = flex.abs(x-y)
  return num/den*2*100., flex.mean(diff), flex.max(diff)

# ===========================================================================

def compute_vals(table):
  # Get randomly 10 chiral spacegroup numbers
  # (computations are longer with higher dmin, so only run on subset of sgs)
  chiral_sgs = []
  for sg_number in range(1,231):
    if space_group_info(sg_number).group().is_chiral():
      chiral_sgs.append(sg_number)
  random_subset = random.sample(chiral_sgs, 10)
  plot_vals = []
  for sg_number in random_subset:
    sgi = space_group_info(sg_number)
    elements = ["C", "O", "N", "H", "S", "K", "Se", "Ca"] * 3
    xrs = random_structure.xray_structure(
      space_group_info       = sgi,
      elements               = elements,
      general_positions_only = False,
      use_u_aniso            = True,
      random_u_iso           = True,
      random_occupancy       = random.choice([True, False]),
    )
    xrs.scattering_type_registry(table = table)
    # Set high-resolution limit to a value between 0.5 and 2 Angstrom
    d_min = random.uniform(0.5, 2)
    fcalc_cctbx = xrs.structure_factors(
      d_min=d_min, algorithm="direct").f_calc().data()
    fcalc_discamb = flex.complex_double(
      pydiscamb.calculate_structure_factors_IAM(xrs, d_min))
    assert fcalc_discamb.size() == fcalc_cctbx.size()
    score, mean_diff, max_diff = compare_structure_factors(
      x=fcalc_cctbx, y=fcalc_discamb)
    #print (sgi, score, mean_diff, max_diff)
    plot_vals.append((score,mean_diff,max_diff))
    # Check if the differences are within the specified tolerance
    assert(score < tolerance_values[table]['score_max'])
    assert(max_diff < tolerance_values[table]['max_diff_max'])
  return plot_vals

# ===========================================================================

def run(make_plots):
  if make_plots: n_runs = 100
  else: n_runs = 1
  plot_data_electron = []
  plot_data_it1992 = []
  plot_data_wk1995 = []
  for i in range(0,n_runs):
    #print('run', i)
    #print('electron')
    plot_vals = compute_vals(table = "electron")
    plot_data_electron.extend(plot_vals)
    #print('it1992')
    plot_vals = compute_vals(table = "it1992")
    plot_data_it1992.extend(plot_vals)
    #print('wk1995')
    plot_vals = compute_vals(table = "wk1995")
    plot_data_wk1995.extend(plot_vals)
  #
  if make_plots:
    dtype1 = np.dtype([
                      ('score','float'),
                      ('mean_diff','float'),
                      ('max_diff','float')
                      ])
    v_electron = np.array(plot_data_electron, dtype = dtype1)
    v_it1992 = np.array(plot_data_it1992, dtype = dtype1)
    v_wk1995 = np.array(plot_data_wk1995, dtype = dtype1)
    print('Making plots...')
    for _table, _data in zip(['electron','it1992','wk1995'],
                             [v_electron, v_it1992, v_wk1995]):
      make_histogram(data = _data['score'],
                     table = _table,
                     ptype='score')
      make_histogram(data = _data['max_diff'],
                     table = _table,
                     ptype='max_diff')
      print(_table)
      #print('max score', round(max(_data['score']),5))
      #print('max max_diff', round(max(_data['max_diff']),5))

# ===========================================================================

def make_histogram(data, table, ptype):
  if ptype=='score':
    binwidth = 0.000005
    label = 'Score'
    maxlim = tolerance_values[table]['score_max']
  if ptype=='max_diff':
    binwidth = 0.0002
    label = 'Maximum difference'
    maxlim = tolerance_values[table]['max_diff_max']
  fig,ax = plt.subplots()
  bins = np.arange(min(data), max(data) + binwidth, binwidth)
  n, bins, patches = ax.hist(data, bins=bins)
  plt.axvline(x=maxlim, color='firebrick')
  ax.set_xlabel(label, fontsize=12, weight='bold')
  ax.set_ylabel('Counts', fontsize=12, weight='bold')
  ax.tick_params(colors='grey', labelsize=8)
  for spine in ax.spines.values():
    spine.set_color('grey')
  ax.xaxis.label.set_color('grey')
  ax.yaxis.label.set_color('grey')
  plt.title('table = %s, n_runs = 100' % table)
  filename = ptype+'_'+ table +'_high-res.png'
  plt.savefig(filename, dpi=600)
  plt.close()

# ===========================================================================

if (__name__ == "__main__"):
  t0 = time.time()
  run(make_plots=False)
  print("OK. Time: %8.3f"%(time.time()-t0))
