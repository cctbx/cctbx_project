from iotbx.mtz import extract_from_symop_lib
from cctbx import miller
from cctbx.array_family import flex

def setSpaceGroup(self, space_group_info, symbol=None):
  if (symbol is None):
    symbol = extract_from_symop_lib.ccp4_symbol(space_group_info)
    if (symbol is None):
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
    assert datatype in ("F", "G")
    self.addColumn(label,datatype,item_miller,flex.abs(item_data))
    self.addColumn(
      self.label_phases(label),"P",item_miller,flex.arg(item_data, True))
  else:
    raise RuntimeError, "Not implemented." # XXX

def add_miller_array(self, miller_array, mtz_label):
  if (not miller_array.anomalous_flag()):
    if (miller_array.is_xray_intensity_array()):
      column_types = "JQ"
    else:
      column_types = "FQ"
    self._DoubleOrComplex(
      mtz_label, column_types[0],
      miller_array.indices(), miller_array.data())
    if (miller_array.sigmas()):
      self._DoubleOrComplex(
        self.label_sigmas(mtz_label), column_types[1],
        miller_array.indices(), miller_array.sigmas())
  else:
    if (miller_array.is_xray_intensity_array()):
      column_types = "KM"
    else:
      column_types = "GL"
    asu, matches = miller_array.match_bijvoet_mates()
    for sign in ("+","-"):
      sel = matches.pairs_hemisphere_selection(sign)
      sel.extend(matches.singles_hemisphere_selection(sign))
      if (sign == "+"):
        label_sign = self.label_plus
        indices = asu.indices().select(sel)
      else:
        label_sign = self.label_minus
        indices = -asu.indices().select(sel)
      self._DoubleOrComplex(
        label_sign(mtz_label), column_types[0],
        indices, asu.data().select(sel))
      if (asu.sigmas() is not None):
        self._DoubleOrComplex(
          label_sign(self.label_sigmas(mtz_label)), column_types[1],
          indices, asu.sigmas().select(sel))

def miller_array_export_as_mtz(self, file_name, column_label,
                                     title=None,
                                     crystal_label="crystal",
                                     project_label="project",
                                     dataset_label="dataset",
                                     wavelength=1.0):
  from iotbx import mtz
  from cctbx import uctbx
  from cctbx import sgtbx
  w = mtz.MtzWriter()
  if (title is None):
    title = self.info()
  if (title is None):
    title = "cctbx.miller.array"
  unit_cell = self.unit_cell()
  if (unit_cell is None):
    unit_cell = uctbx.unit_cell((1,1,1,90,90,90))
  space_group_info = self.space_group_info()
  if (space_group_info is None):
    space_group_info = sgtbx.space_group_info(symbol="P 1")
  w.setTitle(title)
  w.setSpaceGroup(space_group_info)
  w.oneCrystal(crystal_label, project_label, unit_cell)
  w.oneDataset(dataset_label, wavelength)
  w.add_miller_array(self, column_label)
  w.write(file_name)

miller.array.export_as_mtz = miller_array_export_as_mtz
