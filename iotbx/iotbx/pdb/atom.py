from iotbx import wildcard
from scitbx import stl
import scitbx.stl.map
from scitbx.python_utils.misc import adopt_init_args
import sys

class labels:

  def __init__(self,
        name=None,
        altLoc=None,
        resName=None,
        chainID=None,
        resSeq=None,
        iCode=None,
        segID=None,
        MODELserial=None):
    adopt_init_args(self, locals())

  def __str__(self):
    return ",".join([str(id) for id in [
      self.name, self.altLoc, self.resName, self.chainID,
      self.resSeq, self.iCode, self.segID, self.MODELserial]])

def labels_from_string(s):
  fields = []
  for field in s.split(","):
    if (field == "None"): field = None
    fields.append(field)
  assert len(fields) == 8
  return labels(*fields)

class attributes(labels):

  def __init__(self,
        name=None,
        altLoc=None,
        resName=None,
        chainID=None,
        resSeq=None,
        iCode=None,
        segID=None,
        MODELserial=None,
        coordinates=None,
        occupancy=None,
        tempFactor=None,
        element=None,
        charge=None,
        sigCoor=None,
        sigOcc=None,
        sigTemp=None,
        Ucart=None,
        sigUcart=None):
    adopt_init_args(self, locals())

  def set_from_ATOM_record(self, pdb_record):
    if (self.name is None):        self.name = pdb_record.name
    if (self.altLoc is None):      self.altLoc = pdb_record.altLoc
    if (self.resName is None):     self.resName = pdb_record.resName
    if (self.chainID is None):     self.chainID = pdb_record.chainID
    if (self.resSeq is None):      self.resSeq = pdb_record.resSeq
    if (self.iCode is None):       self.iCode = pdb_record.iCode
    if (self.segID is None):       self.segID = pdb_record.segID
    if (self.coordinates is None): self.coordinates = pdb_record.coordinates
    if (self.occupancy is None):   self.occupancy = pdb_record.occupancy
    if (self.tempFactor is None):  self.tempFactor = pdb_record.tempFactor
    if (self.element is None):     self.element = pdb_record.element
    if (self.charge is None):      self.charge = pdb_record.charge
    return self

  def set_from_SIGATM_record(self, pdb_record):
    if (self.sigCoor is None):     self.sigCoor = pdb_record.sigCoor
    if (self.sigOcc is None):      self.sigOcc = pdb_record.sigOcc
    if (self.sigTemp is None):     self.sigTemp = pdb_record.sigTemp
    return self

  def set_from_ANISOU_record(self, pdb_record):
    if (self.Ucart is None):       self.Ucart = pdb_record.Ucart
    return self

  def set_from_SIGUIJ_record(self, pdb_record):
    if (self.sigUcart is None):    self.sigUcart = pdb_record.sigUcart
    return self

  def show(self, f=None, prefix=""):
    if (f is None): f = sys.stdout
    for attr_name in [
        "name", "altLoc", "resName", "chainID", "resSeq", "iCode",
        "segID", "element", "charge",
        "coordinates", "sigCoor",
        "occupancy", "sigOcc",
        "tempFactor", "sigTemp",
        "Ucart", "sigUcart"]:
      value = getattr(self, attr_name)
      if (isinstance(value, list)):
        value = "(" + ", ".join(["%.6g" % v for v in value]) + ")"
      elif (isinstance(value, float)):
        value = "%.6g" % value
      elif (isinstance(value, str)):
        value = '"' + value.replace('"','\\"') + '"'
      print >> f, "%s%-12s" % (prefix, attr_name+":"), value

def _get_map_string(map, pattern):
  result = []
  for key,value in map.items():
    if (wildcard.is_match(string=key, pattern=pattern)):
      result.append(value)
  return result

class selection_cache:

  def __init__(self, atom_attributes_list):
    self.name = stl.map.stl_string_stl_vector_unsigned()
    self.altLoc = stl.map.stl_string_stl_vector_unsigned()
    self.resName = stl.map.stl_string_stl_vector_unsigned()
    self.chainID = stl.map.stl_string_stl_vector_unsigned()
    self.resSeq = stl.map.int_stl_vector_unsigned()
    self.iCode = stl.map.stl_string_stl_vector_unsigned()
    self.segID = stl.map.stl_string_stl_vector_unsigned()
    self.MODELserial = stl.map.int_stl_vector_unsigned()
    self.element = stl.map.stl_string_stl_vector_unsigned()
    self.charge = stl.map.stl_string_stl_vector_unsigned()
    for i_seq,atom_attributes in enumerate(atom_attributes_list):
      self.name.setdefault(atom_attributes.name).append(i_seq)
      self.altLoc.setdefault(atom_attributes.altLoc).append(i_seq)
      self.resName.setdefault(atom_attributes.resName).append(i_seq)
      self.chainID.setdefault(atom_attributes.chainID).append(i_seq)
      self.resSeq.setdefault(atom_attributes.resSeq).append(i_seq)
      self.iCode.setdefault(atom_attributes.iCode).append(i_seq)
      self.segID.setdefault(atom_attributes.segID).append(i_seq)
      self.MODELserial.setdefault(atom_attributes.MODELserial).append(i_seq)
      self.element.setdefault(atom_attributes.element).append(i_seq)
      self.charge.setdefault(atom_attributes.charge).append(i_seq)

  def get_name(self, pattern):
    return _get_map_string(map=self.name, pattern=pattern)

  def get_altLoc(self, pattern):
    return _get_map_string(map=self.altLoc, pattern=pattern)

  def get_resName(self, pattern):
    return _get_map_string(map=self.resName, pattern=pattern)

  def get_chainID(self, pattern):
    return _get_map_string(map=self.chainID, pattern=pattern)

  def get_resSeq(self, i):
    result = self.resSeq.get(i, None)
    if (result is None): return []
    return [result]

  def get_iCode(self, pattern):
    return _get_map_string(map=self.iCode, pattern=pattern)

  def get_segID(self, pattern):
    return _get_map_string(map=self.segID, pattern=pattern)

  def get_MODELserial(self, i):
    result = self.MODELserial.get(i, None)
    if (result is None): return []
    return [result]

  def get_element(self, pattern):
    return _get_map_string(map=self.element, pattern=pattern)

  def get_charge(self, pattern):
    return _get_map_string(map=self.charge, pattern=pattern)
