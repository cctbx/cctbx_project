from iotbx import wildcard
from scitbx.array_family import flex
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

  def residue_labels(self):
    return ",".join([str(id) for id in [
      self.resName, self.chainID,
      self.resSeq, self.iCode, self.segID]])

  def pdb_format(self):
    result = ""
    if (self.MODELserial is not None and self.MODELserial != 0):
      result += "MODEL     %4d: " % self.MODELserial
    result = '"%-4.4s%1.1s%-3.3s %1.1s%4d%1.1s"' % (
      self.name,self.altLoc,self.resName,self.chainID,self.resSeq,self.iCode)
    if (self.segID is not None and len(self.segID.strip()) != 0):
      result += ' ... segID="%s"' % self.segID
    return result

  def is_in_same_chain(self, other):
    if (self.chainID != other.chainID): return False
    if (self.chainID == " " and self.segID != other.segID): return False
    return True

  def is_label_equivalent(self, other):
    if (self.name != other.name): return False
    if (self.altLoc != other.altLoc): return False
    if (self.resName != other.resName): return False
    if (not self.is_in_same_chain(other)): return False
    if (self.resSeq != other.resSeq): return False
    if (self.iCode != other.iCode): return False
    return True

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
        sigUcart=None,
        line_number=None):
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
    self.n_seq = len(atom_attributes_list)
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

  def get_all_altLocs_sorted(self):
    result = self.altLoc.keys()
    if (" " in self.altLoc):
      result.remove(" ")
    result.sort()
    result.insert(0, " ")
    return result

  def get_model_indices(self):
    if (self.MODELserial.size() < 2):
      return None
    result = flex.size_t(self.n_seq, self.n_seq)
    for model_index,selection in self.MODELserial.items():
      assert model_index >= 0
      result.set_selected(selection, model_index)
    assert result.all_ne(self.n_seq)
    return result

  def get_conformer_indices(self):
    result = flex.size_t(self.n_seq, self.n_seq)
    for i,altLoc in enumerate(self.get_all_altLocs_sorted()):
      selection = self.altLoc.get(altLoc, None)
      if (selection is None):
        assert altLoc == " "
      else:
        result.set_selected(selection, i)
    assert result.all_ne(self.n_seq)
    return result

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

  def union(self, iselections):
    return flex.union(
      size=self.n_seq,
      iselections=iselections)

  def intersection(self, iselections):
    return flex.intersection(
      size=self.n_seq,
      iselections=iselections)

  def sel_name(self, pattern):
    return self.union(iselections=self.get_name(pattern=pattern))

  def sel_altLoc(self, pattern):
    return self.union(iselections=self.get_altLoc(pattern=pattern))

  def sel_resName(self, pattern):
    return self.union(iselections=self.get_resName(pattern=pattern))

  def sel_chainID(self, pattern):
    return self.union(iselections=self.get_chainID(pattern=pattern))

  def sel_resSeq(self, i):
    return self.union(iselections=self.get_resSeq(i=i))

  def sel_iCode(self, pattern):
    return self.union(iselections=self.get_iCode(pattern=pattern))

  def sel_segID(self, pattern):
    return self.union(iselections=self.get_segID(pattern=pattern))

  def sel_MODELserial(self, i):
    return self.union(iselections=self.get_MODELserial(i=i))

  def sel_element(self, pattern):
    return self.union(iselections=self.get_element(pattern=pattern))

  def sel_charge(self, pattern):
    return self.union(iselections=self.get_charge(pattern=pattern))

  def get_labels(self,
        name=None,
        altLoc=None,
        resName=None,
        chainID=None,
        resSeq=None,
        iCode=None,
        segID=None,
        MODELserial=None):
    result = []
    for arg,attr in [(name, self.name),
                     (altLoc, self.altLoc),
                     (resName, self.resName),
                     (chainID, self.chainID),
                     (resSeq, self.resSeq),
                     (iCode, self.iCode),
                     (segID, self.segID),
                     (MODELserial, self.MODELserial)]:
      if (arg is not None):
        isel = attr.get(arg, None)
        if (isel is not None): result.append(isel)
    return result
