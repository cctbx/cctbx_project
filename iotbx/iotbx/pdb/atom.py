from iotbx import simple_parser
from iotbx import simple_tokenizer
from iotbx import wildcard
from scitbx.array_family import flex
from scitbx import stl
import scitbx.stl.map
import sys

class labels:

  __slots__ = ["name", "altLoc", "resName", "chainID", "resSeq", "iCode",
               "segID", "MODELserial"]

  def __init__(self,
        name=None,
        altLoc=None,
        resName=None,
        chainID=None,
        resSeq=None,
        iCode=None,
        segID=None,
        MODELserial=None):
    self.name = name
    self.altLoc = altLoc
    self.resName = resName
    self.chainID = chainID
    self.resSeq = resSeq
    self.iCode = iCode
    self.segID = segID
    self.MODELserial = MODELserial

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

  __slots__ = ["is_hetatm",
               "coordinates", "occupancy", "tempFactor", "element", "charge",
               "sigCoor", "sigOcc", "sigTemp", "Ucart", "sigUcart",
               "line_number"]

  def __init__(self,
        is_hetatm=False,
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
    self.is_hetatm = is_hetatm
    self.name = name
    self.altLoc = altLoc
    self.resName = resName
    self.chainID = chainID
    self.resSeq = resSeq
    self.iCode = iCode
    self.segID = segID
    self.MODELserial = MODELserial
    self.coordinates = coordinates
    self.occupancy = occupancy
    self.tempFactor = tempFactor
    self.element = element
    self.charge = charge
    self.sigCoor = sigCoor
    self.sigOcc = sigOcc
    self.sigTemp = sigTemp
    self.Ucart = Ucart
    self.sigUcart = sigUcart
    self.line_number = line_number

  def set_from_ATOM_record(self, pdb_record):
    if (pdb_record.record_name == "HETATM"):
      self.is_hetatm = True
    else:
      self.is_hetatm = False
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

  def record_name(self):
    if (self.is_hetatm): return "HETATM"
    return "ATOM"

  def show(self, f=None, prefix=""):
    if (f is None): f = sys.stdout
    print >> f, "%s%-12s" % (prefix, "record name:"), self.record_name()
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

def _get_map_string(map, pattern, wildcard_escape_char='\\'):
  do_strip = False
  do_upper = False
  if (not isinstance(pattern, str)):
    if (pattern.quote_char is None):
      do_strip = True
      do_upper = True
      for c in pattern.value:
        if (c in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"):
          do_upper = False
          break
    pattern = pattern.value
    if (do_strip): pattern = pattern.strip()
    if (do_upper): pattern = pattern.upper()
  result = []
  for key,value in map.items():
    if (do_strip): key = key.strip()
    if (do_upper): key = key.upper()
    if (wildcard.is_match(
          string=key,
          pattern=pattern,
          escape_char=wildcard_escape_char)):
      result.append(value)
  return result

class selection_cache:

  def __init__(self, atom_attributes_list, wildcard_escape_char='\\'):
    self.n_seq = len(atom_attributes_list)
    self.wildcard_escape_char = wildcard_escape_char
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
    return _get_map_string(
      map=self.name,
      pattern=pattern,
      wildcard_escape_char=self.wildcard_escape_char)

  def get_altLoc(self, pattern):
    return _get_map_string(
      map=self.altLoc,
      pattern=pattern,
      wildcard_escape_char=self.wildcard_escape_char)

  def get_resName(self, pattern):
    return _get_map_string(
      map=self.resName,
      pattern=pattern,
      wildcard_escape_char=self.wildcard_escape_char)

  def get_chainID(self, pattern):
    return _get_map_string(
      map=self.chainID,
      pattern=pattern,
      wildcard_escape_char=self.wildcard_escape_char)

  def get_resSeq(self, i):
    result = self.resSeq.get(i, None)
    if (result is None): return []
    return [result]

  def get_resSeq_range(self, i, j):
    if (i is None): i = min(self.resSeq.keys())
    if (j is None): j = max(self.resSeq.keys())
    if (i > j):
      raise RuntimeError("resSeq range with first index > last index.")
    result = []
    for i in xrange(i,j+1):
      iselection = self.resSeq.get(i, None)
      if (iselection is not None): result.append(iselection)
    return result

  def get_iCode(self, pattern):
    return _get_map_string(
      map=self.iCode,
      pattern=pattern,
      wildcard_escape_char=self.wildcard_escape_char)

  def get_segID(self, pattern):
    return _get_map_string(
      map=self.segID,
      pattern=pattern,
      wildcard_escape_char=self.wildcard_escape_char)

  def get_MODELserial(self, i):
    result = self.MODELserial.get(i, None)
    if (result is None): return []
    return [result]

  def get_MODELserial_range(self, i, j):
    if (i is None): i = min(self.MODELserial.keys())
    if (j is None): j = max(self.MODELserial.keys())
    if (i > j):
      raise RuntimeError("MODELserial range with first index > last index.")
    result = []
    for i in xrange(i,j+1):
      iselection = self.MODELserial.get(i, None)
      if (iselection is not None): result.append(iselection)
    return result

  def get_element(self, pattern):
    return _get_map_string(
      map=self.element,
      pattern=pattern,
      wildcard_escape_char=self.wildcard_escape_char)

  def get_charge(self, pattern):
    return _get_map_string(
      map=self.charge,
      pattern=pattern,
      wildcard_escape_char=self.wildcard_escape_char)

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

  def sel_resSeq_range(self, i, j):
    return self.union(iselections=self.get_resSeq_range(i=i, j=j))

  def sel_iCode(self, pattern):
    return self.union(iselections=self.get_iCode(pattern=pattern))

  def sel_segID(self, pattern):
    return self.union(iselections=self.get_segID(pattern=pattern))

  def sel_MODELserial(self, i):
    return self.union(iselections=self.get_MODELserial(i=i))

  def sel_MODELserial_range(self, i, j):
    return self.union(iselections=self.get_MODELserial_range(i=i, j=j))

  def sel_element(self, pattern):
    return self.union(iselections=self.get_element(pattern=pattern))

  def sel_charge(self, pattern):
    return self.union(iselections=self.get_charge(pattern=pattern))

  def selection_tokenizer(self, string, contiguous_word_characters=None):
    if (contiguous_word_characters is None):
      contiguous_word_characters \
        = simple_tokenizer.default_contiguous_word_characters \
        + r"\*?[]^+-.:"
    word_stack = simple_tokenizer.split_into_words(
      input_string=string,
      contiguous_word_characters=contiguous_word_characters)
    word_stack.reverse()
    return word_stack

  def selection_parser(self,
        word_stack,
        callback=None,
        stop_word=None,
        expect_nonmatching_closing_parenthesis=False):
    result_stack = []
    for word,word_stack in simple_parser.infix_as_postfix(
          word_stack=word_stack,
          stop_word=stop_word,
          expect_nonmatching_closing_parenthesis
            =expect_nonmatching_closing_parenthesis):
      if (word.value == "not"):
        assert len(result_stack) >= 1
        arg = result_stack.pop()
        result_stack.append(~arg)
      elif (word.value in ["and", "or"]):
        assert len(result_stack) >= 2
        rhs = result_stack.pop()
        lhs = result_stack.pop()
        if (word.value == "and"):
          result_stack.append(lhs & rhs)
        else:
          result_stack.append(lhs | rhs)
      else:
        lword = word.value.lower()
        if (lword == "all"):
          result_stack.append(flex.bool(self.n_seq, True))
        elif (lword == "none"):
          result_stack.append(flex.bool(self.n_seq, False))
        elif (lword == "name"):
          if (len(word_stack) == 0): raise RuntimeError("Missing argument.")
          result_stack.append(self.sel_name(pattern=word_stack.pop()))
        elif (lword in ["altloc", "altid"]):
          if (len(word_stack) == 0): raise RuntimeError("Missing argument.")
          result_stack.append(self.sel_altLoc(pattern=word_stack.pop()))
        elif (lword == "resname"):
          if (len(word_stack) == 0): raise RuntimeError("Missing argument.")
          result_stack.append(self.sel_resName(pattern=word_stack.pop()))
        elif (lword == "chain"):
          if (len(word_stack) == 0): raise RuntimeError("Missing argument.")
          result_stack.append(self.sel_chainID(pattern=word_stack.pop()))
        elif (lword in ["resseq", "resid", "model"]):
          if (len(word_stack) == 0): raise RuntimeError("Missing argument.")
          arg = word_stack.pop()
          i_colon_or_dash = arg.value.find(":")
          if (i_colon_or_dash < 0): i_colon_or_dash = arg.value.find("-")
          if (i_colon_or_dash < 0):
            try: i = int(arg.value)
            except ValueError: raise RuntimeError("Value error.")
            if (lword != "model"):
              result_stack.append(self.sel_resSeq(i=i))
            else:
              result_stack.append(self.sel_MODELserial(i=i))
          else:
            i,j = arg.value[:i_colon_or_dash], arg.value[i_colon_or_dash+1:]
            if (len(i) == 0):
              i = None
            else:
              try: i = int(i)
              except ValueError: raise RuntimeError("Value error.")
            if (len(j) == 0):
              j = None
            else:
              try: j = int(j)
              except ValueError: raise RuntimeError("Value error.")
            if (lword != "model"):
              result_stack.append(self.sel_resSeq_range(i=i, j=j))
            else:
              result_stack.append(self.sel_MODELserial_range(i=i, j=j))
        elif (lword == "icode"):
          if (len(word_stack) == 0): raise RuntimeError("Missing argument.")
          result_stack.append(self.sel_iCode(pattern=word_stack.pop()))
        elif (lword == "segid"):
          if (len(word_stack) == 0): raise RuntimeError("Missing argument.")
          result_stack.append(self.sel_segID(pattern=word_stack.pop()))
        elif (lword == "element"):
          if (len(word_stack) == 0): raise RuntimeError("Missing argument.")
          result_stack.append(self.sel_element(pattern=word_stack.pop()))
        elif (lword == "charge"):
          if (len(word_stack) == 0): raise RuntimeError("Missing argument.")
          result_stack.append(self.sel_charge(pattern=word_stack.pop()))
        elif (callback is not None):
          if (not callback(
                    word=word,
                    word_stack=word_stack,
                    result_stack=result_stack)):
            raise RuntimeError("Syntax error.")
        else:
          raise RuntimeError("Syntax error.")
    if (len(result_stack) == 0):
      return flex.bool(self.n_seq, False)
    selection = result_stack[0]
    for result in result_stack[1:]:
      selection &= result
    return selection

  def selection(self, string, contiguous_word_characters=None, callback=None):
    return self.selection_parser(
      word_stack=self.selection_tokenizer(
        string=string,
        contiguous_word_characters=contiguous_word_characters),
      callback=callback)

  def iselection(self, string, contiguous_word_characters=None):
    return self.selection(
      string=string,
      contiguous_word_characters=contiguous_word_characters).iselection()

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
