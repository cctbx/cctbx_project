
def check_bin_format (bin) :
  try :
    d_max, d_min = float(bin[0]), float(bin[1])
  except ValueError, e :
    raise RuntimeError("%s\nOffending values: %s, %s"%(str(e),bin[0],bin[1]))

class processing_info (object) :
  def __init__ (self, bins) :
    self.bins = bins
    for bin in bins :
      check_bin_format(bin)
    self.stats_overall = {}
    self.binned_stats = {}

  def add_bin_stat (self, bin, stat_name, value) :
    check_bin_format(bin)
    if stat_name in self.binned_stats :
      self.binned_stats[stat_name].append(value)
    else :
      self.binned_stats[stat_name] = [value]

  def add_overall_stat (self, stat_name, value) :
    self.stats_overall[stat_name] = value

  def format_remark_200 (self) :
    lines = []
    lines.append("OVERALL.")
    comp_overall = self.stats_overall.get("completeness", "NULL")
    redu_overall = self.stats_overall.get("redundancy", "NULL")
    rmerg_overall = self.stats_overall.get("r_merge", "NULL")
    s2n_overall = self.stats_overall.get("i/sigma", "NULL")
    lines.append(" COMPLETENESS FOR RANGE     (%%) : %s" % comp_overall)
    lines.append(" DATA REDUNDANCY                : %s" % redu_overall)
    lines.append(" R MERGE                    (I) : %s" % rmerg_overall)
    lines.append(" R SYM                      (I) : NULL")
    lines.append(" <I/SIGMA(I)> FOR THE DATA SET  : %s" % s2n_overall)
    lines.append("")
    lines.append("IN THE HIGHEST RESOLUTION SHELL.")
    d_min = self.bins[-1][1]
    d_max = self.bins[-1][0]
    comp_lastbin = self.binned_stats.get("completeness", ["NULL"])[-1]
    redu_lastbin = self.binned_stats.get("redundancy", ["NULL"])[-1]
    rmerg_lastbin = self.binned_stats.get("r_merge", ["NULL"])[-1]
    s2n_lastbin = self.binned_stats.get("i/sigma", ["NULL"])[-1]
    lines.append(" HIGHEST RESOLUTION SHELL, RANGE HIGH (A) : %s" % d_min)
    lines.append(" HIGHEST RESOLUTION SHELL, RANGE LOW  (A) : %s" % d_max)
    lines.append(" COMPLETENESS FOR SHELL     (%%) : %s" % comp_lastbin)
    lines.append(" DATA REDUNDANCY IN SHELL       : %s" % redu_lastbin)
    lines.append(" R MERGE FOR SHELL          (I) : %s" % rmerg_lastbin)
    lines.append(" R SYM FOR SHELL            (I) : NULL")
    lines.append(" <I/SIGMA(I)> FOR SHELL         : %s" % s2n_lastbin)
    remark_lines = [ "REMARK 200 %s" % line for line in lines ]
    return "\n".join(remark_lines)

def parse_scalepack (lines) :
  mode = 0
  info = None
  def is_table_end (fields) :
    return (fields[0] == "All" and fields[1] == "hkl")
  for i, line in enumerate(lines) :
    if "Summary of observation redundancies by shells" in line :
      bins = []
      j = i + 3
      while (j < (i+100)) :
        line2 = lines[j]
        fields = line2.strip().split()
        if is_table_end(fields) :
          break
        else :
          bin_d_max_min = (fields[0], fields[1])
          bins.append(bin_d_max_min)
        j += 1
      assert (len(bins) > 0)
      info = processing_info(bins)
    elif "Average Redundancy Per Shell" in line :
      j = i + 3
      while (j < (i+100)) :
        line2 = lines[j]
        fields = line2.strip().split()
        if is_table_end(fields) :
          info.add_overall_stat("redundancy", fields[-1])
          break
        else :
          bin = (fields[0], fields[1])
          info.add_bin_stat(bin, "redundancy", fields[-1])
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
        r_merge = fields[-1] # XXX -1 (linear) or -2 (square) ???
        if (fields[0] == "All" and fields[1] == "reflections") :
          info.add_overall_stat("i/sigma", "%.2f" % (i_mean / sig_i_mean))
          info.add_overall_stat("r_merge", r_merge)
          break
        else :
          bin = (fields[0], fields[1])
          info.add_bin_stat(bin, "i/sigma", "%.2f" % (i_mean / sig_i_mean))
          info.add_bin_stat(bin, "r_merge", r_merge)
        j += 1
      break
  return info

def run (args) :
  from libtbx.utils import Sorry
  import os
  if not os.path.isfile(args[0]) :
    raise Sorry("First argument must be a valid file name.")
  lines = open(args[0], "r").readlines()
  for line in lines :
    if "reading from a file" in line :
      info = parse_scalepack(lines)
      break
  if info is not None :
    print info.format_remark_200()

if __name__ == "__main__" :
  import sys
  run(sys.argv[1:])
