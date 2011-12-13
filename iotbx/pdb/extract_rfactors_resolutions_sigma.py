import sys
from libtbx.str_utils import format_value
from libtbx import smart_open

class get_r_rfree_sigma(object):
  def __init__(self, remark_2_and_3_records, file_name):
    self.file_name = file_name
    self.r_work,self.r_free,self.sigma,self.high,self.low,self.resolution \
      = [],[],[],[],[],[]
    start_DataUsedInRefinement = False
    start_FitToDataUsedInRefinement = False
    for line in remark_2_and_3_records:
      line = line.strip()
      flds = line.split()
      if(len(flds) == 0): continue
      if(not start_DataUsedInRefinement):
        start_DataUsedInRefinement = self.is_DataUsedInRefinement(line)
      def get_value_at(i):
        if (len(flds) <= i): return None
        try: return float(flds[i])
        except ValueError: return None
      if(start_DataUsedInRefinement and self.is_ResolutionRangeHigh(line)):
        self.high.append(get_value_at(i=7))
      if(start_DataUsedInRefinement and self.is_ResolutionRangeLow(line)):
        self.low.append(get_value_at(i=7))
      if(start_DataUsedInRefinement and self.is_DataCutoffSigma(line)):
        self.sigma.append(get_value_at(i=6))
      if(not start_FitToDataUsedInRefinement):
        start_FitToDataUsedInRefinement = \
          self.is_FitToDataUsedInRefinement(line)
      if(start_FitToDataUsedInRefinement and self.is_RValueWorkingSet(line)):
        self.r_work.append(self.get_value(flds=flds))
      if(start_FitToDataUsedInRefinement and self.is_FreeRValue(line)):
        self.r_free.append(self.get_value(flds=flds))
      if(self.is_Resolution(line)):
        tmp = get_value_at(i=3)
        if (self.resolution is None):
          try: tmp = float(line[22:28])
          except ValueError: pass
        self.resolution.append(tmp)
    if(len(self.r_work)>1 or len(self.r_work)==0): self.r_work = None
    else: self.r_work = self.r_work[0]
    if(len(self.r_free)>1 or len(self.r_free)==0): self.r_free = None
    else: self.r_free = self.r_free[0]
    if(len(self.sigma)>1 or len(self.sigma)==0): self.sigma = None
    else: self.sigma = self.sigma[0]
    if(len(self.high)>1 or len(self.high)==0): self.high = None
    else: self.high = self.high[0]
    if(len(self.low)>1 or len(self.low)==0): self.low = None
    else: self.low = self.low[0]
    if(len(self.resolution)>1 or len(self.resolution)==0): self.resolution = None
    else: self.resolution = self.resolution[0]

  def get_value(self, flds):
    last = flds[-1]
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
    r1 = line.startswith("REMARK   3   R VALUE            (WORKING SET) ")
    r2 = line.startswith(
      "REMARK   3   R VALUE          (WORKING SET, NO CUTOFF) ")
    return r1 or r2

  def is_FreeRValue(self, line):
    r1 = line.startswith("REMARK   3   FREE R VALUE                     ")
    r2 = line.startswith(
      "REMARK   3   FREE R VALUE                  (NO CUTOFF) ")
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
    return None
  result = get_r_rfree_sigma(
    remark_2_and_3_records = remarks_2_and_3,
    file_name              = file_name)
  return result
