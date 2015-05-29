
"""
Rotamer distribution analysis tool for validation of models generated from
cryoEM data.  Written for use with EMRinger pkl output.

Author: Benjamin Barad
Reference:
  Barad BA, Echols N, Wang RYR, Cheng YC, DiMaio F, Adams PD, Fraser JS.
  Side-chain-directed model and map validation for 3D Electron Cryomicroscopy.
  Manuscript in preparation.
"""

from __future__ import division
from mmtbx.ringer import Peak, Peaklist
from libtbx import easy_pickle
import numpy as np
from collections import OrderedDict
import argparse
import math
import os
import sys

# Residue_codes = ["PHE","TYR","TRP"]
Residue_codes = ["ARG","ASN","ASP","CYS","GLU","GLN","HIS",
"LEU","LYS","MET","PHE","SER","TRP","TYR","SEC","PYL"]

Ignored_codes = ["ALA","GLY","PRO","THR","ILE","VAL"]

def statistic(binned_peaks):
  """
  This is the main pair of statistics used for the plots.  Normal approximation
  to the binomial theorem.
  """
  rotamer_count = sum(binned_peaks[0::2])
  total_count = sum(binned_peaks)
  stdev = math.sqrt(39.0/72*33.0/72*total_count)
  mean= total_count*39.0/72
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


def parse_pickle(filename):
  """
  All processes that require reading the pickle. Involves reading out the
  angles and calculating the thresholds.
  """
  print "===== Loading Pickle: %s =====" % filename
  chi = 1
  waves=[]
  averages=[]
  maxima=[]
  ringer_things = easy_pickle.load(filename)
  for i in ringer_things:
    if chi in i._angles.keys() and i.resname in Residue_codes:
      waves.append(i)
      maxima.append(max(i._angles[chi].densities))
      averages.append(np.average(i._angles[chi].densities))
  max_max = max(maxima)
  avg_avg = np.average(averages)
  thresholds = [avg_avg+i*(max_max-avg_avg)/20 for i in range(20)]
  print "Threshold list: %s" % thresholds
  print "===== Pickle Parsed ====="
  return waves, thresholds

def calculate_binned_counts(peak_count, first=60, binsize=12,n_angles=72):
  """Bin peaks by rotamer regions for statistics."""
  first_loc = int(first/5)
  bins = int(n_angles/binsize)
  binned_output=[0]*bins
  for i in range(bins):
    for j in range(binsize):
      binned_output[i] += peak_count[int(first_loc+i*binsize-binsize/2+j)%72]
  return binned_output

def calc_ratio(count_list, sampling_angle=5):
  """
  Calculate the same statistics as the "statistic" call, but do it without
  first binning the peaks.
  """
  # Calculate the same statistics as the "statistic" call, but do it without ifrst binning the peaks.
  total_angles=360/sampling_angle
  binsize=int(total_angles/6)
  first_loc=60/sampling_angle

  binned_list=[0]*6
  for i in range(6):
    for j in range(binsize):
      binned_list[i] += count_list[int(first_loc+i*binsize-binsize/2+j)%72]
  rotamer_count = sum(binned_list[0::2])
  total_count = sum(binned_list)
  stdev = math.sqrt((total_angles/2+3)*(total_angles/2-3)/(total_angles**2)*total_count)
  mean= total_count*(total_angles/2+3)/total_angles
  rotamer_ratio=rotamer_count/(total_count+0.000000000000000000001)
  zscore=(rotamer_count-mean)/(stdev+0.000000000000000000001)
  return rotamer_ratio, zscore

class main (object) :
  def __init__ (self, file_name, sampling_angle=5, out=sys.stdout,
      quiet=False) :
    if (not quiet) :
      out_dir = file_name + ".output"
      if (not os.path.isdir(out_dir)) :
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
    self.zscores=[]
    self.rotamer_ratios=[]
    self.non_zero_thresholds=[]
    waves, self.thresholds = parse_pickle(file_name)
    self.length = len(waves)
    self.peaks=OrderedDict()
        # calculate peaks and histogram
    for threshold in self.thresholds:
      if (not quiet) :
        print >>out, ""
        print >>out, "===== Calculating Statistics for Threshold %.3f =====" %\
          threshold
      self.peaks[threshold]=Peaklist()
      Weird_residues[threshold]=Peaklist()
      self.peak_count[threshold] = [0]*72
      for i in Residue_codes:
        residue_peak_count[i][threshold]=[0]*72
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
      zscore_n, rotamer_ratio_n = statistic(binned_peaks[threshold])
      if rotamer_ratio_n==0:
        break
      for i in Residue_codes:
        rotamer_ratios_residues_n, zscores_n = calc_ratio(residue_peak_count[i][threshold], sampling_angle)
        rotamer_ratios_residues[i].append(rotamer_ratios_residues_n)
        zscores_residues[i].append(zscores_n)
      self.non_zero_thresholds.append(threshold)
      self.zscores.append(zscore_n)
      self.rotamer_ratios.append(rotamer_ratio_n)
      if (not quiet) :
        print "===== Plotting Histogram for Threshold %.3f =====" % threshold
        plot_peaks(
          peak_count=self.peak_count[threshold],
          filename=file_name,
          threshold=threshold,
          first=60,
          title=RMSD_statistic(self.peaks[threshold].peaks))
      # plot_rotamers(binned_peaks[threshold], file, threshold, args.first_rotamer)
    #   print "Outliers at threshold %.2f: %s" % (threshold, str(Weird_residues[threshold]))
    if (not quiet) :
      print ""
      print "===== Plotting Statistics Across Thresholds ====="
      plot_progression(
        non_zero_thresholds=self.non_zero_thresholds,
        rotamer_ratios=self.rotamer_ratios,
        file_name=file_name,
        zscores=self.all_scores)
    # for i in Residue_codes:
    #   plot_progression(non_zero_thresholds, rotamer_ratios_residues[i], file, zscores_residues[i], i)
      print ""
      print "===== Writing Pickles Out ====="
      easy_pickle.dump(file_name+'.output/Outliers.pkl',Weird_residues)
      print 'Wrote ' + file_name+'.output/Outliers.pkl'
      easy_pickle.dump(file_name+'.output/rotamer_ratios.pkl', self.rotamer_ratios)
      print 'Wrote ' + file_name+'.output/rotamer_ratios.pkl'
      easy_pickle.dump(file_name+'.output/zscores.pkl', self.zscores)
      print 'Wrote ' + file_name+'.output/zscores.pkl'
      easy_pickle.dump(file_name+'.output/emringer_scores.pkl', self.all_scores)
      print 'Wrote ' + file_name+'.output/emringer_scores.pkl'
      easy_pickle.dump(file_name+'.output/thresholds.pkl', self.thresholds)
      print 'Wrote ' + file_name+'.output/thresholds.pkl'
      easy_pickle.dump(file_name+'.output/peak_counts.pkl', self.peak_count)
      print 'Wrote ' + file_name+'.output/peak_counts.pkl'
    self.zscore_max = max(self.zscores)
    self._zscore_max_index = self.zscores.index(self.zscore_max)

  @property
  def all_scores (self) :
    return [ 10*z/math.sqrt(self.length) for z in self.zscores ]

  @property
  def optimal_threshold (self) :
    return self.non_zero_thresholds[self._zscore_max_index]

  @property
  def rotamer_ratio (self) :
    return self.rotamer_ratios[self._zscore_max_index]

  @property
  def score (self) :
    """Overall EM-Ringer score"""
    return (10 * self.zscore_max / math.sqrt(self.length))

  def show_summary (self, out=sys.stdout) :
    print >> out, ""
    print >> out, "=====Final Statistics for Model/Map Pair====="
    print >> out, "Optimal Threshold: %.3f" % self.optimal_threshold
    print >> out, "Rotamer-Ratio: %.3f" % self.rotamer_ratio
    print >> out, "Max Zscore: %.3f" % self.zscore_max
    print >> out, "Model Length: %d" % self.length
    # print "Z-score/(length): %.8f" % (value/length)
    print >> out, "EMRinger Score: %8f" % self.score
    return self

  def draw_wx_peaks_plot (self, plot, i_threshold) :
    threshold = self.thresholds[i_threshold]
    plot.figure.clear()
    _plot_peaks(
      fig=plot.figure,
      peak_count=self.peak_count[threshold],
      threshold=threshold,
      first=60,
      title=RMSD_statistic(self.peaks[threshold].peaks))

########################################################################
# GUI and Output

def _plot_rotamers (fig, binned_output, threshold, first) :
  """Binned histogram"""
  colors=['blue','red']*3
  angles = range(6)
  bin_angles = [(i*60+first)%360 for i in angles]
  ax = fig.add_subplot(111)
  ax.bar(bin_angles, binned_output, align='center', color=colors, width=60)

def plot_rotamers (binned_output, threshold, first):
  import matplotlib.pyplot as plt
  fig = plt.figure(1)
  _plot_rotamers(fig, binned_output, filename, threshold, first)
  plt.savefig('%s.output/%.3f.Phenixed_Histogram.png' % (filename,threshold))
  plt.close()

def _plot_peaks (fig, peak_count, threshold, first, title=0):
  # rcParams.update({'figure.autolayout': True})
  colors = ['#F15854']*6+['#5DA5DA']*13+['#F15854']*11+['#5DA5DA']*13+['#F15854']*11+['#5DA5DA']*13+['#F15854']*5
  ax = fig.add_subplot(111)
  ax.axvspan((first-30), first+30, color='0.5', alpha=0.5)
  ax.axvspan(first+90, first+150, color='0.5', alpha=0.5)
  ax.axvspan(first+210, (first+270), color='0.5', alpha=0.5)
  ax.tick_params(axis='x',which='both',top='off')
  ax.tick_params(axis='y',which='both',right='off')
  peak_count_extra = peak_count+[peak_count[0]]
  angles = [i*5 for i in range(0,73)]
  ax.bar(angles,peak_count_extra, width=5, align='center', color=colors)
  ax.set_title('Peak Counts', y=1.05) #  - Threshold %.3f' % (threshold)
  ax.set_xticks([i*60 for i in range(7)])
  ax.set_xlim(0,360)
  ax.set_xlabel(r'Chi1 Angle ($\degree$)', labelpad=10)
  ax.set_ylabel("Peak Count", labelpad=10)

def plot_peaks (peak_count, filename, threshold, first, title=0) :
  import matplotlib.pyplot as plt
  fig = plt.figure(2, figsize=(5,4))
  _plot_peaks(fig, peak_count, threshold, first, title)
  plt.savefig('%s.output/%.3f.histogram.png' % (filename,threshold))
  print 'Saved plot to %s.output/%.3f.histogram.png' % (filename,threshold)
  # print 'RMSD at threshold %.3f is %.1f' % (threshold,title)
  # print 'Wrote '+filename+'/%.3f.Phenix_allpeaks.png' % threshold
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
  plt.savefig('%s.output/%s.threshold_scan.png' % (file_name, i))
  print 'Saved plot to %s.output/%s.threshold_scan.png' % (file_name, i)
  # print 'Wrote '+file+'/threshold_scan.png'
  plt.close()

def run (args, out=sys.stdout) :
    parser = argparse.ArgumentParser()
    parser.add_argument("files",nargs="*")
    parser.add_argument("-s", "--Sampling_Angle", dest="sampling_angle", help="Don't mess with this unless you've also made the corresponding change in ringer. By default it is 5, which is identical to the default in ringer.", nargs='?', default=5)
    parser.add_argument("-r", "--Residues", dest="residues")
    parser.add_argument("--gui", dest="show_gui", action="store_true",
      default=False)
    args = parser.parse_args(args)
    if (not args.show_gui) :
      import matplotlib
      matplotlib.use("Agg")
    app = None
    for file_name in args.files :
      result = main(
        file_name=file_name,
        sampling_angle=args.sampling_angle,
        out=out,
        quiet=False).show_summary(out=out)
      if (args.show_gui) :
        import wxtbx.plots.emringer
        import wxtbx.app
        if (app is None) :
          app = wxtbx.app.CCTBXApp(0)
        f1 = wxtbx.plots.emringer.peaks_plot_frame(
          parent=None,
          title="Histograms for %s" % file_name)
        f1.SetResult(result)
        f1.Show()
    if (args.show_gui) :
      app.MainLoop()

if __name__ == "__main__":
  run(sys.argv[1:])
