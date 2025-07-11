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

assertion_values = {
  "electron": {'score_max': 0.00015,
               'max_diff_max': 0.003},
  "it1992":   {'score_max': 0.00017,
               'max_diff_max': 0.011},
  "wk1995":   {'score_max': 0.0002,
               'max_diff_max': 0.011},
}

def compare_structure_factors(x,y):
  x = flex.abs(x)
  y = flex.abs(y)
  scale = flex.sum(x*y)/flex.sum(y*y)
  num = flex.sum(flex.abs(x-scale*y))
  den = flex.sum(flex.abs(x+scale*y))
  diff = flex.abs(x-y)
  return num/den*2*100., flex.mean(diff), flex.max(diff)

def test_aniso(table):
  """
  Test structure factor calculations for all chiral space groups using
  anisotropic atomic displacement.

  Parameters
  ----------
  table : str
    The scattering factor table to use ('electron', 'it1992', or 'wk1995').

  Returns
  -------
  plot_vals : list of tuples
    Each tuple contains (score, mean_diff, max_diff) for each space group.
    These values are used for plotting.

  Raises
  ------
  AssertionError
    If the relative difference score or maximum difference exceeds predefined
    thresholds for the specified table.

  Notes
  -----
  For each chiral space group, this function:
  - Generates random crystal structures with specified elements.
  - Computes structure factors using cctbx and DiSCaMB and compares the results.
  - Checks that the differences are within acceptable limits defined by
    `assertion_values`.
  """
  plot_vals = []
  for sg_number in range(1,231):
    sgi = space_group_info(sg_number)
    # Only iterate over 65 chiral space groups (so it is faster)
    if not sgi.group().is_chiral(): continue
    elements = ["C", "O", "N", "H", "S", "Na", "Fe", "Ca"] * 5
    random_occ = random.choice([True, False])
    xrs = random_structure.xray_structure(
      space_group_info       = sgi,
      elements               = elements,
      general_positions_only = False,
      use_u_aniso            = True,
      random_u_iso           = True,
      random_occupancy       = random_occ,
    )
    xrs.scattering_type_registry(table = table)
    # Set high-resolution limit to a value between 2 and 5 Angstrom
    d_min = random.uniform(2, 5)
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
    assert(score < assertion_values[table]['score_max'])
    assert(max_diff < assertion_values[table]['max_diff_max'])
  return plot_vals

def run(make_plots = False):
  """
  Run structure factor tests and optionally generate histograms for statistical
  analysis.

  Parameters
  ----------
  make_plots : bool, optional
    If True, generate and save histograms for score and max_diff for each table
    (default is False).

  Returns
  -------
  None

  Notes
  -----
  Repeatedly calls `test_aniso` for each table ('electron', 'it1992', 'wk1995')
  to collect score and difference data. If `make_plots` is True, histograms are
  created and saved for each table.
  """
  if make_plots: n_runs = 100
  else: n_runs = 1
  plot_data_electron = []
  plot_data_it1992 = []
  plot_data_wk1995 = []
  for i in range(0,n_runs):
    print('run', i)
    plot_vals = test_aniso(table = "electron")
    plot_data_electron.extend(plot_vals)
    plot_vals = test_aniso(table = "it1992")
    plot_data_it1992.extend(plot_vals)
    plot_vals = test_aniso(table = "wk1995")
    plot_data_wk1995.extend(plot_vals)
  #
  if make_plots:
    # create numpy array
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


def make_histogram(data, table, ptype):
  """
  Generate and save a histogram for the specified data.

  Parameters
  ----------
  data : array-like
      Data for histogram creation, typically score or max_diff values.
  table : str
      The scattering table type ('electron', 'it1992', or 'wk1995').
  ptype : str
      The type of histogram, either 'score' or 'max_diff'.

  Returns
  -------
  None

  Notes
  -----
  The function automatically adjusts bin width and label based on `ptype`, and
  adds a vertical line indicating the threshold from `assertion_values`.
  Histograms are saved as PNG files.
  """
  if ptype=='score':
    binwidth = 0.000003
    label = 'Score'
    maxlim = assertion_values[table]['score_max']
  if ptype=='max_diff':
    binwidth = 0.0002
    label = 'Maximum difference'
    maxlim = assertion_values[table]['max_diff_max']
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
  filename = ptype+'_'+ table +'.png'
  plt.savefig(filename, dpi=600)
  plt.close()

if (__name__ == "__main__"):
  """
  Run structure factor tests for all chiral space groups.
  Can optionally generate plots for lots of runs to gather statistics.
  """
  t0 = time.time()
  run(make_plots=False)
  print("OK. Time: %8.3f"%(time.time()-t0))
