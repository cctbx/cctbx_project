from iotbx.mtz import extract_from_symop_lib
from cctbx.array_family import flex

def setSpaceGroup(self, space_group_info, symbol=None):
  if (symbol == None):
    symbol = extract_from_symop_lib.ccp4_symbol(space_group_info)
    if (symbol == None):
      symbol = "No.%d" % space_group_info.type().number()
  assert not " " in symbol
  self.raw_setSpaceGroup(space_group_info.group(), symbol)

def label_phases(self, root_label):
  return self.phases_prefix() + root_label + self.phases_suffix()

def label_sigmas(self, root_label):
  return self.sigmas_prefix() + root_label + self.sigmas_suffix()

def label_plus(self, root_label):
  return root_label + self.plus_suffix()

def label_minus(self, root_label):
  return root_label + self.minus_suffix()

def set_phases_prefix(self, prefix):
  self._phases_prefix = prefix

def phases_prefix(self):
  if (not hasattr(self, "_phases_prefix")):
    self._phases_prefix = "PHI"
  return self._phases_prefix

def set_phases_suffix(self, suffix):
  self._phases_suffix = suffix

def phases_suffix(self):
  if (not hasattr(self, "_phases_suffix")):
    self._phases_suffix = ""
  return self._phases_suffix

def set_sigmas_prefix(self, prefix):
  self._sigmas_prefix = prefix

def sigmas_prefix(self):
  if (not hasattr(self, "_sigmas_prefix")):
    self._sigmas_prefix = "SIG"
  return self._sigmas_prefix

def set_sigmas_suffix(self, suffix):
  self._sigmas_suffix = suffix

def sigmas_suffix(self):
  if (not hasattr(self, "_sigmas_suffix")):
    self._sigmas_suffix = ""
  return self._sigmas_suffix

def set_plus_suffix(self, suffix):
  self._plus_suffix = suffix

def plus_suffix(self):
  if (not hasattr(self, "_plus_suffix")):
    self._plus_suffix = "(+)"
  return self._plus_suffix

def set_minus_suffix(self, suffix):
  self._minus_suffix = suffix

def minus_suffix(self):
  if (not hasattr(self, "_minus_suffix")):
    self._minus_suffix = "(-)"
  return self._minus_suffix

def _DoubleOrComplex(self,label,datatype,item_miller,item_data):
  if isinstance(item_data,flex.double):
    self.addColumn(label,datatype,item_miller,item_data)
  elif isinstance(item_data,flex.complex_double):
    self.addColumn(label,datatype,item_miller,flex.abs(item_data))
    self.addColumn(
      self.label_phases(label),"P",item_miller,flex.arg(item_data, True))
  else:
    raise RuntimeError, "Not implemented." # XXX

def _columnCombinations(self,label,datatype,carry_miller,carry_data):
  if len(carry_miller)==1:  # anomalous flag is False
    self._DoubleOrComplex(label,datatype,carry_miller[0],carry_data[0])
  elif len(carry_miller)==2:  # anomalous data
    self._DoubleOrComplex(
      self.label_plus(label),datatype,carry_miller[0],carry_data[0])
    self._DoubleOrComplex(
      self.label_minus(label),datatype,carry_miller[1],carry_data[1])
  else:
    raise RuntimeError, "Internal error."

def add_miller_array(self, miller_array, mtz_label):
  if (not miller_array.anomalous_flag()):
    self._columnCombinations(mtz_label,"F",
      [miller_array.indices()],[miller_array.data()])
    if (miller_array.sigmas()):
      self._columnCombinations(self.label_sigmas(mtz_label),"Q",
        [miller_array.indices()],[miller_array.sigmas()])
  else:
    carry_miller = []
    carry_data = []
    carry_sigma = []
    asu, matches = miller_array.match_bijvoet_mates()
    for i,s in ((0,"+"),(1,"-")):
      sel = matches.pairs_hemisphere_selection(s) \
          | matches.singles_hemisphere_selection(s)
      if (i == 0):
        carry_miller.append(asu.indices().select(sel))
      else:
        carry_miller.append(-asu.indices().select(sel))
      carry_data.append(asu.data().select(sel))
      if (asu.sigmas() != None):
        carry_sigma.append(asu.sigmas().select(sel))
    self._columnCombinations(mtz_label,"G",carry_miller,carry_data)
    if (asu.sigmas()):
      self._columnCombinations(
        self.label_sigmas(mtz_label),"L",carry_miller,carry_sigma)
