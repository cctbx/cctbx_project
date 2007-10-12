import sys, os
from libtbx.str_utils import format_value
from libtbx import smart_open

class get_r_rfree_sigma(object):
  def __init__(self, remark_2_and_3_records, file_name):
    self.file_name = file_name
    self.r_work,self.r_free,self.sigma,self.high,self.low,self.resolution \
      = [None]*6
    start_DataUsedInRefinement = False
    start_FitToDataUsedInRefinement = False
    for line in remark_2_and_3_records:
      line = line.strip()
      if(not start_DataUsedInRefinement):
        start_DataUsedInRefinement = self.is_DataUsedInRefinement(line)
      if(start_DataUsedInRefinement and self.is_ResolutionRangeHigh(line)):
        try: self.high = float(line.split()[7])
        except ValueError: pass
      if(start_DataUsedInRefinement and self.is_ResolutionRangeLow(line)):
         try: self.low = float(line.split()[7])
         except ValueError: pass
      if(start_DataUsedInRefinement and self.is_DataCutoffSigma(line)):
         try: self.sigma = float(line.split()[6])
         except ValueError: pass
      if(not start_FitToDataUsedInRefinement):
        start_FitToDataUsedInRefinement = \
          self.is_FitToDataUsedInRefinement(line)
      if(start_FitToDataUsedInRefinement and self.is_RValueWorkingSet(line)):
        self.r_work = self.get_value(line)
      if(start_FitToDataUsedInRefinement and self.is_FreeRValue(line)):
        self.r_free = self.get_value(line)
      if(self.is_Resolution(line)):
        try: self.resolution = float(line.split()[3])
        except ValueError:
          try: self.resolution = float(line[22:28])
          except ValueError: pass

  def get_value(self, line):
    line = line.split()
    last = line[len(line)-1]
    value = None
    try: value = float(last)
    except ValueError:
      try: value = float(last[1:])
      except ValueError: pass
    return value

  def is_Resolution(self, line):
    r1 = line.startswith("REMARK   2 RESOLUTION")
    r2 = line.endswith("ANGSTROMS") or line.endswith("ANGSTROMS.")
    return r1 and r2

  def is_DataUsedInRefinement(self, line):
    r1 = line.startswith("REMARK   3  DATA USED IN REFINEMENT")
    return r1

  def is_ResolutionRangeHigh(self, line):
    r1 = line.startswith("REMARK   3   RESOLUTION RANGE HIGH")
    return r1

  def is_ResolutionRangeLow(self, line):
    r1 = line.startswith("REMARK   3   RESOLUTION RANGE LOW")
    return r1

  def is_DataCutoffSigma(self, line):
    r1 = line.startswith("REMARK   3   DATA CUTOFF            (SIGMA(F)) :")
    return r1

  def is_RValueWorkingSet(self, line):
    r1 = line.startswith("REMARK   3   R VALUE            (WORKING SET) :")
    r2 = line.startswith(
      "REMARK   3   R VALUE          (WORKING SET, NO CUTOFF) :")
    return r1 or r2

  def is_FreeRValue(self, line):
    r1 = line.startswith("REMARK   3   FREE R VALUE                     :")
    r2 = line.startswith(
      "REMARK   3   FREE R VALUE                  (NO CUTOFF) :")
    return r1 or r2

  def is_FitToDataUsedInRefinement(self, line):
    r1 = line.startswith("REMARK   3  FIT TO DATA USED IN REFINEMENT")
    r2 = line.startswith("REMARK   3  USING DATA ABOVE SIGMA CUTOFF.")
    r3 = line.startswith(
      "REMARK   3  FIT TO DATA USED IN REFINEMENT (NO CUTOFF).")
    result = r1 or r2 or r3
    return result

  def formatted_string(self):
    result = "%s %s %s %s %s %s %s" % (
      format_value("%6s",self.file_name),
      format_value("%6s",str(self.r_work)),
      format_value("%6s",str(self.r_free)),
      format_value("%6s",str(self.sigma)),
      format_value("%6s",str(self.high)),
      format_value("%6s",str(self.low)),
      format_value("%6s",str(self.resolution)))
    return result

  def show(self, log = None):
    if(log is None): log = sys.stdout
    print >> self.log, self.formatted_string()

def extract_remark_2_and_3_records(file_name):
  result = []
  file_lines = smart_open.for_reading(file_name = file_name).read().splitlines()
  for rec in file_lines:
    if(rec.startswith("REMARK   3 ") or rec.startswith("REMARK   2 ")):
      start = True
      result.append(rec)
    else:
      if(rec.startswith("ATOM ") or rec.startswith("HETATM ")):
        break
  return result

def extract(file_name):
  remarks_2_and_3 = extract_remark_2_and_3_records(file_name)
  if(len(remarks_2_and_3) == 0):
    print "No remarks 2 and 3:", file_name
    return None
  result = get_r_rfree_sigma(
    remark_2_and_3_records = remarks_2_and_3,
    file_name              = file_name)
  return result
