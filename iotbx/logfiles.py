
from __future__ import division
import cStringIO
import sys, os, re

def check_bin_format (bin) :
  try :
    d_max, d_min = float(bin[0]), float(bin[1])
  except ValueError, e :
    raise RuntimeError("%s\nOffending values: %s, %s"%(str(e),bin[0],bin[1]))

def float_or_none (n) :
  if n is None : return None
  else :         return float(n)

def percent_to_float (value) :
  assert value.endswith("%")
  return float(re.sub("\%$", "", value))

class experiment_info (object) :
  def extract_all_stats (self) :
    return self

class integration_info (object) :
  def __init__ (self, program_name="NULL") :
    self.program_name = program_name
    self.wavelength = None
    self.distance = None
    self.twotheta = None

  def set_wavelength (self, wavelength) :
    self.wavelength = float(wavelength)

  def set_distance (self, distance) :
    self.distance = float(distance)

  def set_2theta (self, twotheta) :
    self.twotheta = twotheta

  def extract_all_stats (self) :
    return self

class scaling_info (object) :
  def __init__ (self, program_name="NULL") :
    self.program_name = program_name
    self.stats_overall = {}
    self.binned_stats = {}
    self.bins = None
    self.d_max = None
    self.d_min = None
    self.n_refl = None
    self.n_refl_all = None

  def set_bins (self, bins) :
    for bin in bins :
      check_bin_format(bin)
    self.bins = bins
    if self.d_max is None :
      d_max = float(self.bins[0][0])
      d_min = float(self.bins[-1][1])
      self.set_d_max_min(d_max, d_min)

  def set_n_refl (self, n_refl, n_refl_all) :
    self.n_refl = n_refl
    self.n_refl_all = n_refl_all

  def add_bin_stat (self, bin, stat_name, value) :
    check_bin_format(bin)
    if stat_name in self.binned_stats :
      self.binned_stats[stat_name].append(value)
    else :
      self.binned_stats[stat_name] = [value]

  def add_overall_stat (self, stat_name, value) :
    self.stats_overall[stat_name] = value

  def set_d_max_min (self, d_max, d_min) :
    self.d_max = d_max
    self.d_min = d_min

  def extract_all_stats (self) :
    from libtbx import group_args
    d_min = float(self.bins[-1][1])
    d_max = float(self.bins[0][0])
    comp_overall = self.stats_overall.get("completeness", None)
    mult_overall = self.stats_overall.get("multiplicity", None)
    rmerg_overall = self.stats_overall.get("r_merge", None)
    s2n_overall = self.stats_overall.get("i/sigma", None)
    return group_args(d_max_min=(d_max, d_min),
                      n_refl=self.n_refl,
                      n_refl_all=self.n_refl_all,
                      completeness=float_or_none(comp_overall),
                      multiplicity=float_or_none(mult_overall),
                      r_sym=float_or_none(rmerg_overall),
                      i_over_sigma=float_or_none(s2n_overall))

  def extract_outer_shell_stats (self) :
    from libtbx import group_args
    d_min = float(self.bins[-1][1])
    d_max = float(self.bins[-1][0])
    comp_bin = self.binned_stats.get("completeness", [None])[-1]
    mult_bin = self.binned_stats.get("multiplicity", [None])[-1]
    rmerg_bin = self.binned_stats.get("r_merge", [None])[-1]
    s2n_bin = self.binned_stats.get("i/sigma", [None])[-1]
    return group_args(d_max_min=(d_max, d_min),
                      n_refl=None, # FIXME
                      n_refl_all=None,
                      completeness=float_or_none(comp_bin),
                      multiplicity=float_or_none(mult_bin),
                      r_sym=float_or_none(rmerg_bin),
                      i_over_sigma=float_or_none(s2n_bin))

class all_none (object) :
  def __getattr__ (self, name) :
    return None

class empty_info (object) :
  def extract_all_stats (self) :
    return all_none()
  def extract_outer_shell_stats (self) :
    return all_none()

class processing_info (object) :
  def __init__ (self, experiment, integration, scaling) :
    self.experiment = experiment
    self.integration = integration
    self.scaling = scaling

  def get_experiment_info (self) :
    if (self.experiment is not None) :
      return self.experiment
    return all_none() #empty_info()

  def get_integration_info (self) :
    if (self.integration is not None) :
      return self.integration
    return all_none() #empty_info()

  def get_scaling_info (self) :
    if (self.scaling is not None) :
      return self.scaling
    return empty_info()

  def format_remark_200 (self) :
    from libtbx.str_utils import format_value
    from libtbx.test_utils import approx_equal
    def format (obj, attr, fs="%.4f") :
      value = getattr(obj, attr, None)
      return format_value(fs, value, replace_none_with="NULL").strip()
    e = None
    if self.experiment is not None :
      e = self.experiment.extract_all_stats()
    i = None
    if self.integration is not None :
      i = self.integration.extract_all_stats()
    s = None
    if self.scaling is not None :
      s = self.scaling.extract_all_stats()
    lines = []
    lines.append("")
    lines.append("EXPERIMENTAL DETAILS")
    lines.append(" EXPERIMENT TYPE                : X-RAY DIFFRACTION")
    lines.append(" DATE OF DATA COLLECTION        : NULL")
    lines.append(" TEMPERATURE           (KELVIN) : NULL")
    lines.append(" PH                             : NULL")
    lines.append(" NUMBER OF CRYSTALS USED        : NULL")
    lines.append("")
    # TODO
    wavelength = getattr(e, "wavelength", None)
    if (wavelength is None) :
      wavelength = getattr(i, "wavelength", None)
    synchrotron = "NULL"
    if (wavelength is not None) :
      out = cStringIO.StringIO()
      if (not approx_equal(wavelength, 1.5418, eps=0.01, out=out) and
          not approx_equal(wavelength, 0.7107, eps=0.01, out=out)) :
        synchrotron = "Y"
      else :
        synchrotron = "N"
      wl = "%.4f" % wavelength
    else :
      wl = "NULL"
    lines.append(" SYNCHROTRON              (Y/N) : %s" % synchrotron)
    lines.append(" RADIATION SOURCE               : NULL")
    lines.append(" BEAMLINE                       : NULL")
    lines.append(" X-RAY GENERATOR MODEL          : NULL")
    lines.append(" MONOCHROMATIC OR LAUE    (M/L) : M")
    lines.append(" WAVELENGTH OR RANGE        (A) : %s" % wl)
    lines.append(" MONOCHROMATOR                  : NULL")
    lines.append(" OPTICS                         : NULL")
    lines.append("")
    int_software = getattr(self.integration, "program_name", "NULL")
    lines.append(" DETECTOR TYPE                  : NULL")
    lines.append(" DETECTOR MANUFACTURER          : NULL")
    lines.append(" INTENSITY-INTEGRATION SOFTWARE : %s" % int_software)
    scale_software = getattr(self.scaling, "program_name", "NULL")
    lines.append(" DATA SCALING SOFTWARE          : %s" % scale_software)
    lines.append("")
    lines.append("OVERALL.")
    comp_overall = format(s, "completeness", "%.1f")
    mult_overall = format(s, "multiplicity", "%.1f")
    rmerg_overall = format(s, "r_sym", "%.5f")
    s2n_overall = format(s, "i_over_sigma", "%.4f")
    lines.append(" COMPLETENESS FOR RANGE     (%%) : %s" % comp_overall)
    lines.append(" DATA REDUNDANCY                : %s" % mult_overall)
    lines.append(" R MERGE                    (I) : %s" % rmerg_overall)
    lines.append(" R SYM                      (I) : NULL")
    lines.append(" <I/SIGMA(I)> FOR THE DATA SET  : %s" % s2n_overall)
    lines.append("")
    lines.append("IN THE HIGHEST RESOLUTION SHELL.")
    shell = None
    if self.scaling is not None :
      shell = self.scaling.extract_outer_shell_stats()
    (_d_max, _d_min) = getattr(shell, "d_max_min", (None, None))
    d_max = format_value("%.2f", _d_max, replace_none_with="NULL").strip()
    d_min = format_value("%.2f", _d_min, replace_none_with="NULL").strip()
    comp_lastbin = format(shell, "completeness", "%.1f")
    mult_lastbin = format(shell, "multiplicity", "%.1f")
    rmerg_lastbin = format(shell, "r_sym", "%.5f")
    s2n_lastbin = format(shell, "i_over_sigma", "%.4f")
    lines.append(" HIGHEST RESOLUTION SHELL, RANGE HIGH (A) : %s" % d_min)
    lines.append(" HIGHEST RESOLUTION SHELL, RANGE LOW  (A) : %s" % d_max)
    lines.append(" COMPLETENESS FOR SHELL     (%%) : %s" % comp_lastbin)
    lines.append(" DATA REDUNDANCY IN SHELL       : %s" % mult_lastbin)
    lines.append(" R MERGE FOR SHELL          (I) : %s" % rmerg_lastbin)
    lines.append(" R SYM FOR SHELL            (I) : NULL")
    lines.append(" <I/SIGMA(I)> FOR SHELL         : %s" % s2n_lastbin)
    lines.append("")
    remark_lines = [ "REMARK 200 %s" % line for line in lines ]
    return "\n".join(remark_lines)

#-----------------------------------------------------------------------
# PARSERS
#
def parse_denzo (lines) :
  info = integration_info("HKL-2000")
  for i, line in enumerate(lines) :
    if line.strip().startswith("Wavelength ") :
      fields = line.strip().split()
      for field in fields :
        try :
          wavelength = float(field)
        except ValueError :
          pass
        else :
          info.set_wavelength(wavelength)
          break
    elif line.strip().startswith("Detector to crystal distance") :
      fields = line.strip().split()
      info.set_distance(float(fields[4]))
  return info

def parse_mosflm (lines) :
  info = integration_info("MOSFLM")
  for i, line in enumerate(lines) :
    line = line.strip()
    if line.startswith("Beam Parameters") :
      j = i
      while (j < len(lines)) :
        line = lines[j].strip()
        if line.startswith("Wavelength") :
          wavelength = float(line.split()[1])
          info.set_wavelength(wavelength)
          break
        j += 1
    elif line.startswith("Detector Parameters") :
      j = i
      while (j < (i + 100)) :
        line = lines[j].strip()
        if line.startswith("Crystal to detector distance") :
          fields = line.split()
          distance = float(fields[-2])
          info.set_distance(distance)
        elif line.startswith("Detector swing angle") :
          fields = line.split()
          twotheta = float(fields[-2])
          info.set_twotheta(twotheta)
        j += 1
      break
  return info

def parse_scalepack (lines) :
  n_lines = len(lines)
  mode = 0
  info = scaling_info("SCALA")
  def is_table_end (fields) :
    return (fields[0] == "All" and fields[1] == "hkl")
  n_refl_all = None
  n_refl = None
  for i, line in enumerate(lines) :
    if ("intensities and R-factors by batch number" in line) :
      j = i + 3
      while j < n_lines :
        line2 = lines[j].strip()
        if line2.startswith("All films") :
          n_refl_all = int(line2.split()[2])
          break
        j+= 1
    elif "Summary of observation redundancies by shells" in line :
      bins = []
      j = i + 3
      while (j < (i+100)) :
        line2 = lines[j]
        fields = line2.strip().split()
        if is_table_end(fields) :
          n_refl = int(fields[-1])
          info.set_n_refl(n_refl, n_refl_all)
          break
        else :
          bin_d_max_min = (fields[0], fields[1])
          bins.append(bin_d_max_min)
        j += 1
      assert (len(bins) > 0)
      info.set_bins(bins)
    elif "Average Redundancy Per Shell" in line :
      j = i + 3
      while (j < (i+100)) :
        line2 = lines[j]
        fields = line2.strip().split()
        if is_table_end(fields) :
          info.add_overall_stat("multiplicity", fields[-1])
          break
        else :
          bin = (fields[0], fields[1])
          info.add_bin_stat(bin, "multiplicity", fields[-1])
        j += 1
    elif "I/Sigma in resolution shells:" in line :
      j = i + 3
      while (j < (i+100)) :
        line2 = lines[j]
        fields = line2.strip().split()
        if is_table_end(fields) :
          info.add_overall_stat("completeness", fields[-1])
          break
        else :
          bin = (fields[0], fields[1])
          info.add_bin_stat(bin, "completeness", fields[-1])
        j += 1
    elif "Summary of reflections intensities and R-factors by shells" in line :
      j = i
      while (j < (i+100)) :
        line2 = lines[j]
        fields = line2.strip().split()
        j += 1
        if (len(fields) > 0 and fields[0] == "limit") :
          break
      while (j < (i+200)) :
        line2 = lines[j]
        fields = line2.strip().split()
        i_mean = float(fields[2])
        sig_i_mean = float(fields[3])
        r_merge = fields[-2] # XXX -1 (linear) or -2 (square) ???
        if (fields[0] == "All" and fields[1] == "reflections") :
          info.add_overall_stat("i/sigma", "%.2f" % (i_mean / sig_i_mean))
          info.add_overall_stat("r_merge", r_merge)
          break
        else :
          bin = (fields[0], fields[1])
          info.add_bin_stat(bin, "i/sigma", "%.2f" % (i_mean / sig_i_mean))
          info.add_bin_stat(bin, "r_merge", r_merge)
        j += 1
  return info

def parse_scala (lines) :
  from iotbx import data_plots
  info = scaling_info("SCALA")
  tables = data_plots.import_ccp4i_logfile(log_lines=lines)
  d_max = None
  for i, line in enumerate(lines) :
    if ("Summary data for " in line) :
      if (lines[i+1].startswith("</p>")) or ("<br" in line) :
        continue
      j = i
      n_refl = None
      n_refl_all = None
      while (j < len(lines)) :
        line = lines[j].strip()
        if line.startswith("Low resolution limit") :
          d_max = float(line.split()[3])
        elif line.startswith("Rmerge") and (not "bin" in line) :
          info.add_overall_stat("r_merge", float(line.split()[1]))
        elif line.startswith("Total number of observations") :
          n_refl_all = float(line.split()[4])
        elif line.startswith("Total number unique") :
          n_refl = float(line.split()[3])
          info.set_n_refl(n_refl, n_refl_all)
        elif (line.startswith("Mean((I)/sd(I))") or
              line.startswith("Mean(I)/sd(I)")) :
          info.add_overall_stat("i/sigma", float(line.split()[1]))
        elif line.startswith("Completeness") :
          info.add_overall_stat("completeness", float(line.split()[1]))
        elif line.startswith("Multiplicity") :
          info.add_overall_stat("multiplicity", float(line.split()[1]))
        elif ("Outlier rejection" in line) or ("$$" in line) :
          break
        j += 1
  assert (d_max is not None)
  for table in tables :
    if table.title.startswith("Analysis against resolution") :
      d_min_by_bin = table.get_column_by_label("Dmin(A)")
      bin_d_max = d_max
      bins = []
      for bin_d_min in d_min_by_bin :
        bins.append((bin_d_max, bin_d_min))
        bin_d_max = bin_d_min
      info.set_bins(bins)
      rmerge = table.get_column_by_label("Rmrg")
      for (rmerge_bin, bin) in zip(rmerge, bins) :
        info.add_bin_stat(bin, "r_merge", rmerge_bin)
      try :
        s2n = table.get_column_by_label("Mn(I/sd)")
      except Exception :
        s2n = table.get_column_by_label("Mn(I)/sd")
      for (s2n_bin, bin) in zip(s2n, bins) :
        info.add_bin_stat(bin, "i/sigma", s2n_bin)
    elif table.title.startswith("Completeness, multiplicity, Rmeas") :
      completeness = table.get_column_by_label("%poss")
      for (comp_bin, bin) in zip(completeness, bins) :
        info.add_bin_stat(bin, "completeness", comp_bin)
      multiplicity = table.get_column_by_label("Mlplct")
      for (mult_bin, bin) in zip(multiplicity, bins) :
        info.add_bin_stat(bin, "multiplicity", mult_bin)
      break
  return info

def parse_xds (lines) :
  info = integration_info("XDS")
  for i, line in enumerate(lines) :
    line = line.strip()
    if line.startswith("X-RAY_WAVELENGTH") :
      fields = line.split("=")[1].strip().split()
      info.set_wavelength(float(fields[0]))
      break
  return info

def parse_xscale (lines) :
  info = scaling_info("XSCALE")
  d_max = 0.0
  d_min = 999.99
  for i, line in enumerate(lines) :
    line = line.strip()
    if (line.startswith("INPUT_FILE=")) :
      fields = line.split()
      if (len(fields) < 4) :
        continue
      try :
        (d_max_file, d_min_file) = (float(fields[-2]), float(fields[-1]))
      except ValueError :
        pass
      except Exception, e :
        print line
        raise
      else :
        if (d_max_file > d_max) :
          d_max = d_max_file
        if (d_min_file < d_min) :
          d_min = d_min_file
        info.set_d_max_min(d_max, d_min)
    elif line.startswith("SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0") :
      j = i+3
      bins = []
      overall = []
      while (j < len(lines)) :
        line = lines[j].strip()
        fields = line.split()
        if (len(fields) < 12) :
          pass
        elif (fields[0] == "total") :
          overall = fields
          break
        elif re.match("^\d", line) :
          bins.append(fields)
        j += 1
      assert (len(bins) > 0) and (len(overall) == len(bins[0]))
      n_refl_all = int(overall[1])
      n_refl = int(overall[2])
      info.set_n_refl(n_refl, n_refl_all)
      info.add_overall_stat("completeness", percent_to_float(overall[4]))
      info.add_overall_stat("multiplicity", n_refl_all / n_refl)
      info.add_overall_stat("r_merge", percent_to_float(overall[5]) * 0.01)
      info.add_overall_stat("i/sigma", float(overall[8]))
      bins_d_max_min = []
      for bin in bins :
        d_min = float(bin[0])
        bins_d_max_min.append((d_max, d_min))
        d_max = d_min
      info.set_bins(bins_d_max_min)
      for fields, bin  in zip(bins, bins_d_max_min) :
        bin_n_refl_all = int(fields[1])
        bin_n_refl = int(fields[2])
        info.add_bin_stat(bin, "completeness", percent_to_float(fields[4]))
        info.add_bin_stat(bin, "multiplicity", bin_n_refl_all / bin_n_refl)
        info.add_bin_stat(bin, "r_merge", percent_to_float(fields[5]) * 0.01)
        info.add_bin_stat(bin, "i/sigma", float(fields[8]))
      break
  return info

def parse_all_files (args) :
  experiment = None
  integration = None
  scaling = None
  for arg in args :
    if os.path.isfile(arg) :
      lines = open(arg, "r").readlines()
      for line in lines :
        if "reading from a file" in line :
          scaling = parse_scalepack(lines)
          break
        elif "Oscillation starts at" in line :
          integration = parse_denzo(lines)
          break
        elif "SCALA - continuous scaling program" in line :
          scaling = parse_scala(lines)
          break
        elif "A.G.W. Leslie" in line :
          integration = parse_mosflm(lines)
          break
        elif "XSCALE (VERSION" in line :
          scaling = parse_xscale(lines)
          break
        elif "***** INTEGRATE *****" in line :
          integration = parse_xds(lines)
          break
  info = processing_info(experiment=experiment,
    integration=integration,
    scaling=scaling)
  return info

def exercise () :
  import libtbx.load_env
  denzo_log = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/tracking/denzo.log",
    test=os.path.isfile)
  scalepack_log = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/tracking/scalepack.log",
    test=os.path.isfile)
  if (denzo_log is None) :
    print "DENZO log not found, skipping test."
    return False
  info = parse_all_files([denzo_log, scalepack_log])
  output = info.format_remark_200().splitlines()
  assert ("\n".join([ line.strip() for line in output[-16:]]) == """\
REMARK 200 OVERALL.
REMARK 200  COMPLETENESS FOR RANGE     (%) : 88.4
REMARK 200  DATA REDUNDANCY                : 3.5
REMARK 200  R MERGE                    (I) : 0.09600
REMARK 200  R SYM                      (I) : NULL
REMARK 200  <I/SIGMA(I)> FOR THE DATA SET  : 11.5700
REMARK 200
REMARK 200 IN THE HIGHEST RESOLUTION SHELL.
REMARK 200  HIGHEST RESOLUTION SHELL, RANGE HIGH (A) : 1.98
REMARK 200  HIGHEST RESOLUTION SHELL, RANGE LOW  (A) : 2.05
REMARK 200  COMPLETENESS FOR SHELL     (%) : 34.1
REMARK 200  DATA REDUNDANCY IN SHELL       : 1.5
REMARK 200  R MERGE FOR SHELL          (I) : 0.93400
REMARK 200  R SYM FOR SHELL            (I) : NULL
REMARK 200  <I/SIGMA(I)> FOR SHELL         : 0.6000
REMARK 200""")

if __name__ == "__main__" :
  #info = parse_all_files(sys.argv[1:])
  #print info.format_remark_200()
  exercise()
  print "OK"
