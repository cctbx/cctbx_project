
"""
Rotamer distribution analysis tool for validation of models generated from
cryoEM data.  Written for use with EMRinger pkl output.

Author: Benjamin Barad
Reference:
  Barad BA, Echols N, Wang RYR, Cheng YC, DiMaio F, Adams PD, Fraser JS. (2015)
  Side-chain-directed model and map validation for 3D Electron Cryomicroscopy.
  Nature Methods, in press.
"""

from __future__ import absolute_import, division, print_function
from mmtbx.ringer import Peak, Peaklist
from libtbx import easy_pickle
from libtbx.utils import Sorry
from libtbx.math_utils import iceil
from collections import OrderedDict
import math
import os
import sys
from six.moves import range

# Residue_codes = ["PHE","TYR","TRP"]
Residue_codes = ["ARG","ASN","ASP","CYS","GLU","GLN","HIS",
"LEU","LYS","MET","PHE","SER","TRP","TYR","SEC","PYL"]

Ignored_codes = ["ALA","GLY","PRO","THR","ILE","VAL"]

def statistic(binned_peaks, n_angles=72):
  """
  This is the main pair of statistics used for the plots.  Normal approximation
  to the binomial theorem.
  """
  rotamer_count = sum(binned_peaks[0::2])
  total_count = sum(binned_peaks)
  stdev = math.sqrt(39.0/n_angles*33.0/n_angles*total_count)
  mean= total_count*39.0/n_angles
  # Hacky way to avoid zero division
  rotamer_ratio=rotamer_count/(total_count+0.000000000000000000001)
  zscore=(rotamer_count-mean)/(stdev+0.000000000000000000001)
  # print "\t Rotamer ratio: %.3f" % rotamer_ratio
  # print "\t Z-score = %.3f" % zscore
  # if (zscore>0):
  #   pscore_approx1=0.5-0.5*(math.erf(zscore/math.sqrt(2)))
  #   pscore_approx2=1.0/12*math.exp(-zscore*zscore/2)+1.0/4*math.exp(-zscore*zscore*2/3)
  #   # print "\t One approximation of the p-value is %g" % pscore_approx1
    # print "\t Another approximation of the p-value is %g" % pscore_approx2
  # else:
    # print "\t pscore greater than 0.5"

  return zscore, rotamer_ratio

def RMSD_statistic(peak_list):
  """
  Still not clear how useful RMSD is but angular deviations tend to be heavily
  dependent on sample size (as outliers are overweighted).
  """
  squared_deviations=[]
  for peak in peak_list:
    squared_deviations.append(min((i-peak.chi_value)**2 for i in [60,180,300]))
  RMSD = (sum(squared_deviations)/len(squared_deviations))**0.5
  return RMSD

def calculate_peaks(ringer,threshold):
  """
  Checks if something is greater than either of its neighbors (including
  wrapping) and returns if true and if above a threshold)
  """
  new_peaks=Peaklist()
  list = ringer._angles[1].densities
  for i in range(len(list)):
    if (list[i]==max(list) and list[i]>threshold):
      new_peaks.add_new(ringer.resname, ringer.resid, ringer.chain_id, 1, i, list[i])
  return new_peaks


def parse_pickle(filename, out=sys.stdout):
  print("===== Loading Pickle: %s =====" % filename)
  ringer_things = easy_pickle.load(filename)
  return process_raw_results(ringer_things, out=out)

def process_raw_results(ringer_result, out=sys.stdout):
  """
  All processes that require reading the pickle. Involves reading out the
  angles and calculating the thresholds.
  """
  import numpy as np
  chi = 1
  waves=[]
  averages=[]
  maxima=[]
  for residue in ringer_result :
    # TODO verify residue._angles is a dict
    if chi in residue._angles and residue.resname in Residue_codes:
      waves.append(residue)
      maxima.append(max(residue._angles[chi].densities))
      averages.append(np.average(residue._angles[chi].densities))
  max_max = max(maxima)
  avg_avg = np.average(averages)
  thresholds = [avg_avg+i*(max_max-avg_avg)/20 for i in range(20)]
  if min(thresholds) == max(thresholds):
    raise Sorry("""No features could be detected in the density around the model, so EMRinger can't
       proceed. This can be confirmed in pymol or coot. If features are present, raise
       the scale of the map values and try EMRinger again.""")
  print("Threshold list: %s" % thresholds, file=out)
  print("===== Pickle Parsed =====", file=out)
  return waves, thresholds

def calculate_binned_counts(peak_count, first=60, binsize=12,n_angles=72):
  """Bin peaks by rotamer regions for statistics."""
  first_loc = int(first/5)
  bins = int(n_angles/binsize)
  binned_output=[0]*bins
  for i in range(bins):
    for j in range(binsize):
      binned_output[i] += peak_count[int(first_loc+i*binsize-binsize/2+j)%n_angles]
  return binned_output

def calc_ratio(count_list, sampling_angle=5):
  """
  Calculate the same statistics as the "statistic" call, but do it without
  first binning the peaks.
  """
  # Calculate the same statistics as the "statistic" call, but do it without ifrst binning the peaks.
  total_angles=iceil(360.0/sampling_angle)
  binsize=int(total_angles/6)
  first_loc=60/sampling_angle
  binned_list=[0]*6
  for i in range(6):
    for j in range(binsize):
      binned_list[i] += count_list[int(first_loc+i*binsize-binsize/2+j)%total_angles]
  rotamer_count = sum(binned_list[0::2])
  total_count = sum(binned_list)
  stdev = math.sqrt((total_angles/2+3)*(total_angles/2-3)/(total_angles**2)*total_count)
  mean= total_count*(total_angles/2+3)/total_angles
  rotamer_ratio=rotamer_count/(total_count+0.000000000000000000001)
  zscore=(rotamer_count-mean)/(stdev+0.000000000000000000001)
  return rotamer_ratio, zscore

class main(object):
  def __init__(self,
      file_name,
      ringer_result=None,
      sampling_angle=5,
      out_dir=None,
      out=sys.stdout,
      quiet=False):
    self.threshold = waves = None
    if (ringer_result is not None):
      waves, self.thresholds = process_raw_results(ringer_result, out=out)
    else :
      assert (file_name is not None)
      waves, self.thresholds = parse_pickle(file_name, out=out)
    if not quiet:
      assert (out_dir is None) or os.path.isdir(out_dir)
      if (out_dir is None) and (not quiet):
        out_dir = file_name + ".output"
        if (not os.path.isdir(out_dir)):
          os.makedirs(file_name+'.output')
    Weird_residues=OrderedDict()
    self.peak_count={}
    residue_peak_count={}
    rotamer_ratios_residues={}
    zscores_residues={}
    for i in Residue_codes:
      residue_peak_count[i]={}
      rotamer_ratios_residues[i]=[]
      zscores_residues[i]=[]
    binned_peaks={}
    n_angles = iceil(360.0 / sampling_angle)
    self.zscores=[]
    self.rotamer_ratios=[]
    self.non_zero_thresholds=[]
    self.length = len(waves)
    self.peaks=OrderedDict()
        # calculate peaks and histogram
    for threshold in self.thresholds:
      if (not quiet):
        print("", file=out)
        print("===== Calculating Statistics for Threshold %.3f =====" %\
          threshold, file=out)
      self.peaks[threshold]=Peaklist()
      Weird_residues[threshold]=Peaklist()
      self.peak_count[threshold] = [0]*n_angles
      for i in Residue_codes:
        residue_peak_count[i][threshold]=[0]*n_angles
      for i in waves:
        self.peaks[threshold].append_lists(calculate_peaks(i, threshold))
      for peak in self.peaks[threshold].get_peaks():
        self.peak_count[threshold][peak.chi_value] += 1
        residue_peak_count[peak.resname][threshold][peak.chi_value]+=1
        if ((peak.chi_value<6) or (peak.chi_value>18 and peak.chi_value<30) or (peak.chi_value>42 and peak.chi_value<54) or (peak.chi_value>66)):
          Weird_residues[threshold].peaks.append(peak)
      # Calculate the binned peaks and ratios
      binned_peaks[threshold] = calculate_binned_counts(self.peak_count[threshold], 60)
      # print "For threshold %.3f" % threshold
      # print "Sample size = %d" % sum(binned_peaks[threshold])
      zscore_n, rotamer_ratio_n = statistic(binned_peaks[threshold], n_angles)
      if rotamer_ratio_n==0:
        break
      for i in Residue_codes:
        rotamer_ratios_residues_n, zscores_n = calc_ratio(residue_peak_count[i][threshold], sampling_angle)
        rotamer_ratios_residues[i].append(rotamer_ratios_residues_n)
        zscores_residues[i].append(zscores_n)
      self.non_zero_thresholds.append(threshold)
      self.zscores.append(zscore_n)
      self.rotamer_ratios.append(rotamer_ratio_n)
      if (not quiet):
        print("===== Plotting Histogram for Threshold %.3f =====" % \
          threshold, file=out)
        out_file = os.path.join(out_dir, "%.3f.histogram.png" % threshold)
        plot_peaks(
          peak_count=self.peak_count[threshold],
          file_name=out_file,
          threshold=threshold,
          first=60,
          title=RMSD_statistic(self.peaks[threshold].peaks),
          n_angles=n_angles)
        print("Saved plot to %s" % out_file, file=out)
      # plot_rotamers(binned_peaks[threshold], file, threshold, args.first_rotamer)
    #   print "Outliers at threshold %.2f: %s" % (threshold, str(Weird_residues[threshold]))
    if len(self.zscores) == 0:
      raise Sorry("""No scores could be calculated at any threshold for this map. This could be because the
       map is not sufficiently featured, or because of data corruption in the map.""")
    if (not quiet):
      print("", file=out)
      print("===== Plotting Statistics Across Thresholds =====", file=out)
      out_file = os.path.join(out_dir, "Total.threshold_scan.png")
      plot_progression(
        non_zero_thresholds=self.non_zero_thresholds,
        rotamer_ratios=self.rotamer_ratios,
        file_name=out_file,
        zscores=self.all_scores)
      print("Saved plot to %s" % out_file, file=out)
    # for i in Residue_codes:
    #   plot_progression(non_zero_thresholds, rotamer_ratios_residues[i], file, zscores_residues[i], i)
      print("", file=out)
      print("===== Writing Pickles Out =====", file=out)
      easy_pickle.dump(out_dir + '/Outliers.pkl',Weird_residues)
      print('Wrote ' + out_dir + '/Outliers.pkl', file=out)
      easy_pickle.dump(out_dir + '/rotamer_ratios.pkl', self.rotamer_ratios)
      print('Wrote ' + out_dir + '/rotamer_ratios.pkl', file=out)
      easy_pickle.dump(out_dir + '/zscores.pkl', self.zscores)
      print('Wrote ' + out_dir + '/zscores.pkl', file=out)
      easy_pickle.dump(out_dir + '/emringer_scores.pkl', self.all_scores)
      print('Wrote ' + out_dir + '/emringer_scores.pkl', file=out)
      easy_pickle.dump(out_dir + '/thresholds.pkl', self.thresholds)
      print('Wrote ' + out_dir + '/thresholds.pkl', file=out)
      easy_pickle.dump(out_dir + '/peak_counts.pkl', self.peak_count)
      print('Wrote ' + out_dir + '/peak_counts.pkl', file=out)
    self.zscore_max = max(self.zscores)
    self._zscore_max_index = self.zscores.index(self.zscore_max)

  @property
  def all_scores(self):
    return [ 10*z/math.sqrt(self.length) for z in self.zscores ]

  @property
  def optimal_threshold(self):
    return self.non_zero_thresholds[self._zscore_max_index]

  @property
  def rotamer_ratio(self):
    return self.rotamer_ratios[self._zscore_max_index]

  @property
  def score(self):
    """Overall EM-Ringer score"""
    return (10 * self.zscore_max / math.sqrt(self.length))

  def show_summary(self, out=sys.stdout):
    print("", file=out)
    print("=====Final Statistics for Model/Map Pair=====", file=out)
    print("Optimal Threshold: %.3f" % self.optimal_threshold, file=out)
    print("Rotamer-Ratio: %.3f" % self.rotamer_ratio, file=out)
    print("Max Zscore: %.3f" % self.zscore_max, file=out)
    print("Model Length: %d" % self.length, file=out)
    # print "Z-score/(length): %.8f" % (value/length)
    print("EMRinger Score: %8f" % self.score, file=out)
    return self

  def draw_wx_peaks_plot(self, plot, i_threshold):
    threshold = self.thresholds[i_threshold]
    plot.figure.clear()
    _plot_peaks(
      fig=plot.figure,
      peak_count=self.peak_count[threshold],
      threshold=threshold,
      first=60,
      title=RMSD_statistic(self.peaks[threshold].peaks))

  def draw_wx_progression_plot(self, plot):
    plot.figure.clear()
    _plot_progression(
      fig=plot.figure,
      non_zero_thresholds=self.non_zero_thresholds,
      rotamer_ratios=self.rotamer_ratios,
      zscores=self.all_scores)

########################################################################
# GUI and Output

def _plot_rotamers(fig, binned_output, threshold, first):
  """Binned histogram"""
  colors=['blue','red']*3
  angles = range(6)
  bin_angles = [(i*60+first)%360 for i in angles]
  ax = fig.add_subplot(111)
  ax.bar(bin_angles, binned_output, align='center', color=colors, width=60)

def plot_rotamers(binned_output, threshold, first):
  import matplotlib.pyplot as plt
  fig = plt.figure(1)
  _plot_rotamers(fig, binned_output, filename, threshold, first)
  plt.savefig('%s.output/%.3f.Phenixed_Histogram.png' % (filename,threshold))
  plt.close()

def _plot_peaks(fig, peak_count, threshold, first, title=0, n_angles=72):
  # rcParams.update({'figure.autolayout': True})
  colors = ['#F15854']*6+['#5DA5DA']*13+['#F15854']*11+['#5DA5DA']*13+['#F15854']*11+['#5DA5DA']*13+['#F15854']*5
  ax = fig.add_subplot(111)
  ax.axvspan((first-30), first+30, color='0.5', alpha=0.5)
  ax.axvspan(first+90, first+150, color='0.5', alpha=0.5)
  ax.axvspan(first+210, (first+270), color='0.5', alpha=0.5)
  ax.tick_params(axis='x',which='both',top='off')
  ax.tick_params(axis='y',which='both',right='off')
  peak_count_extra = peak_count+[peak_count[0]]
  angles = [i*5 for i in range(0,n_angles+1)]
  ax.bar(angles,peak_count_extra, width=5, align='center', color=colors)
  ax.set_title('Peak Counts', y=1.05) #  - Threshold %.3f' % (threshold)
  ax.set_xticks([i*60 for i in range(7)])
  ax.set_xlim(0,360)
  ax.set_xlabel(r'Chi1 Angle ($\degree$)', labelpad=10)
  ax.set_ylabel("Peak Count", labelpad=10)

def plot_peaks(peak_count, file_name, threshold, first, title=0, n_angles=72):
  import matplotlib.pyplot as plt
  fig = plt.figure(2, figsize=(5,4))
  _plot_peaks(fig, peak_count, threshold, first, title, n_angles=n_angles)
  plt.savefig(file_name)
  plt.close()

def _plot_progression(fig, non_zero_thresholds, rotamer_ratios, zscores,
    i="Total"):
  for j in range(len(zscores)):
    if zscores[j]>=0:
      non_zero_thresholds = non_zero_thresholds[j:]
      rotamer_ratios= rotamer_ratios[j:]
      zscores=zscores[j:]
      break
  ax1 = fig.add_subplot(111)
  ax1.plot(non_zero_thresholds, zscores, 'b-', linewidth=3.0, alpha=0.7)
  ax1.set_xlabel('Electron Potential Threshold')
  # Make the y-axis label and tick labels match the line color.
  ax1.set_ylabel('EMRinger Score', color='b')
  for tl in ax1.get_yticklabels():
    tl.set_color('b')
  ax1.set_ylim([-1,4])
  ax1.axhspan(-1,1,color='0.5',alpha=0.1)
  ax2 = ax1.twinx()
  ax2.grid()
  ax2.plot(non_zero_thresholds, rotamer_ratios, 'r-', label = i, linewidth=3.0, alpha=0.7)
  ax2.set_ylim([0.4,1])
  ax2.set_ylabel(r'% Rotameric Residues', color='r', labelpad=10)
  ax1.xaxis.set_ticks_position('bottom')
  # ax2.set_xlim([0.005,0.03])
  for tl in ax2.get_yticklabels():
    tl.set_color('r')
  if i != "Total":
    ax1.set_title("Threshold Scan - %s" % i, y=1.05)
  else:
    ax1.set_title("Threshold Scan", y=1.05)

def plot_progression(non_zero_thresholds, rotamer_ratios, file_name, zscores,
    i="Total"):
  import matplotlib.pyplot as plt
  fig = plt.figure(3, figsize=(5.5,4))
  _plot_progression(fig, non_zero_thresholds, rotamer_ratios, zscores, i)
  plt.savefig(file_name)
  plt.close()
