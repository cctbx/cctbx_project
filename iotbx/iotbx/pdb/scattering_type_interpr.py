from iotbx.pdb import residue_info
import cctbx.eltbx.xray_scattering
from cctbx import eltbx
import cctbx.eltbx.chemical_elements
import string

class scan_atom_element_columns(object):

  def __init__(self, pdb_records):
    self.n_uninterpretable = 0
    self.n_interpretable = 0
    self.n_q = 0
    chemcial_elements = eltbx.chemical_elements.proper_and_isotopes_upper_set()
    for record in pdb_records:
      if (record.record_name in ("ATOM", "HETATM")):
        if (record.element == " Q"):
          self.n_q += 1
        elif (record.element.lstrip() in chemcial_elements):
          self.n_interpretable += 1
        else:
          self.n_uninterpretable += 1

def from_pdb_atom_record(record, have_useful_atom_element_columns=None):
  try:
    scattering_label = residue_info.get(
      residue_name=record.resName,
      atom_name=record.name).scattering_label
  except KeyError:
    pass
  else:
    return eltbx.xray_scattering.get_standard_label(
      label=scattering_label, exact=True)
  if (have_useful_atom_element_columns and len(record.element.strip()) > 0):
    result = eltbx.xray_scattering.get_standard_label(
      label=record.element+record.charge, exact=False, optional=True)
    if (result is not None): return result
    result = eltbx.xray_scattering.get_standard_label(
      label=record.element, exact=True, optional=True)
    if (result is not None): return result
  label = record.name[:2]
  if (label[0] in string.digits): label = label[1]
  result = eltbx.xray_scattering.get_standard_label(
    label=label.lstrip(), exact=True, optional=True)
  if (result is not None): return result
  if (record.name in ["PEAK", "SITE"]): return "const"
  raise RuntimeError(
    '%sUnknown x-ray scattering coefficients for "%s" "%s"' % (
      record.error_prefix(), record.name, record.resName))
