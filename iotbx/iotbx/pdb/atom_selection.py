from iotbx import simple_parser
from iotbx import wildcard
from scitbx.array_family import flex
from scitbx import stl
import scitbx.stl.map
from libtbx.phil import tokenizer
from libtbx.utils import Sorry
from cStringIO import StringIO
import traceback

def _character_case_id(strings):
  have_upper = False
  have_lower = False
  for s in strings:
    for c in s:
      if   (c in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"):
        if (have_lower): return 0
        have_upper = True
      elif (c in "abcdefghijklmnopqrstuvwxyz"):
        if (have_upper): return 0
        have_lower = True
  if (have_upper): return 1
  if (have_lower): return -1
  return 0

def _get_map_string(map, pattern, wildcard_escape_char='\\'):
  pattern_was_quoted = True
  if (not isinstance(pattern, str)):
    if (pattern.quote_token is None):
      pattern_was_quoted = False
    pattern = pattern.value
    if (not pattern_was_quoted): pattern = pattern.strip()
  result = []
  def match():
    for key,value in map.items():
      if (not pattern_was_quoted): key = key.strip()
      if (wildcard.is_match(
            string=key,
            pattern=pattern,
            escape_char=wildcard_escape_char)):
        result.append(value)
  match()
  if (    len(result) == 0
      and not pattern_was_quoted
      and _character_case_id(strings=[pattern]) != 0):
    keys_case_id = _character_case_id(strings=map.keys())
    if (keys_case_id != 0):
      if (keys_case_id > 0):
        pattern = pattern.upper()
      else:
        pattern = pattern.lower()
      match()
  return result

class selection_tokenizer(tokenizer.word_iterator):

  def __init__(self, string, contiguous_word_characters=None):
    if (contiguous_word_characters is None):
      contiguous_word_characters \
        = tokenizer.default_contiguous_word_characters \
        + r"\*?[]^+-.:"
    tokenizer.word_iterator.__init__(self,
      input_string=string,
      list_of_settings=[tokenizer.settings(
        contiguous_word_characters=contiguous_word_characters)])

  def pop_argument(self, keyword):
    word = self.try_pop()
    if (word is None): raise RuntimeError("Missing argument for %s." % keyword)
    return word

class AtomSelectionError(Sorry):
  __module__ = "exceptions"

class cache(object):

  def __init__(self, atoms, wildcard_escape_char='\\'):
    self.atoms = atoms
    self.n_seq = len(atoms)
    self.wildcard_escape_char = wildcard_escape_char
    self.name = stl.map.stl_string_stl_vector_unsigned()
    self.altLoc = stl.map.stl_string_stl_vector_unsigned()
    self.resName = stl.map.stl_string_stl_vector_unsigned()
    self.chainID = stl.map.stl_string_stl_vector_unsigned()
    self.resSeq = stl.map.stl_string_stl_vector_unsigned()
    self.iCode = stl.map.stl_string_stl_vector_unsigned()
    self.resid = stl.map.stl_string_stl_vector_unsigned()
    self.segID = stl.map.stl_string_stl_vector_unsigned()
    self.MODELserial = stl.map.int_stl_vector_unsigned()
    self.element = stl.map.stl_string_stl_vector_unsigned()
    self.charge = stl.map.stl_string_stl_vector_unsigned()
    self.anisou = stl.vector.unsigned()
    for i_seq,atom in enumerate(atoms):
      resname = "   "
      chainid = " "
      resseq = "    "
      icode = " "
      model_serial = "    "
      ag = atom.parent()
      if (ag is not None):
        altloc = ag.altloc
        if (altloc == ""): altloc = " "
        resname = ag.resname
        rg = ag.parent()
        if (rg is not None):
          resseq = rg.resseq
          icode = rg.icode
          ch = rg.parent()
          if (ch is not None):
            chain_id = ch.id
            mo = ch.parent()
            if (mo is not None):
              model_id = mo.id
      self.name.setdefault(atom.name).append(i_seq)
      self.altLoc.setdefault(altloc).append(i_seq)
      self.resName.setdefault(resname).append(i_seq)
      self.chainID.setdefault(chain_id).append(i_seq)
      self.resSeq.setdefault(resseq).append(i_seq)
      self.iCode.setdefault(icode).append(i_seq)
      self.resid.setdefault(resseq+icode).append(i_seq)
      self.segID.setdefault(atom.segid).append(i_seq)
      self.MODELserial.setdefault(int(model_id)).append(i_seq)
      self.element.setdefault(atom.element).append(i_seq)
      self.charge.setdefault(atom.charge).append(i_seq)
      if (atom.uij_is_defined()): self.anisou.append(i_seq)

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

  def get_resSeq(self, pattern):
    return _get_map_string(
      map=self.resSeq,
      pattern=pattern,
      wildcard_escape_char=self.wildcard_escape_char)

  def get_resSeq_range(self, start, stop):
    from iotbx.pdb import utils_base_256_ordinal as o
    o_start = None
    o_stop = None
    if (start is not None and start.count(" ") != len(start)):
      o_start = o(start)
    if (stop is not None and stop.count(" ") != len(stop)):
      o_stop = o(stop)
    if (    o_start is not None
        and o_stop is not None
        and o_start > o_stop):
      raise RuntimeError(
        "range with first index > last index: resseq %s:%s" % (start, stop))
    result = []
    for s,iselection in self.resSeq.items():
      os = o(s)
      if (o_start is not None and os < o_start): continue
      if (o_stop  is not None and os > o_stop): continue
      result.append(iselection)
    return result

  def get_iCode(self, pattern):
    return _get_map_string(
      map=self.iCode,
      pattern=pattern,
      wildcard_escape_char=self.wildcard_escape_char)

  def get_resid(self, pattern):
    return _get_map_string(
      map=self.resid,
      pattern=pattern,
      wildcard_escape_char=self.wildcard_escape_char)

  def get_resid_range(self, start, stop):
    from iotbx.pdb import utils_base_256_ordinal as o
    def shift(s):
      if (len(s) < 5 and s[-1] in "0123456789"): return s + " "
      return s
    o_start = None
    o_stop = None
    if (start is not None and start.count(" ") != len(start)):
      o_start = o(shift(start))
    if (stop is not None and stop.count(" ") != len(stop)):
      o_stop = o(shift(stop))
    if (    o_start is not None
        and o_stop is not None
        and o_start > o_stop):
      raise RuntimeError(
        "range with first index > last index: resid %s:%s" % (start, stop))
    result = []
    for s,iselection in self.resid.items():
      os = o(s)
      if (o_start is not None and os < o_start): continue
      if (o_stop  is not None and os > o_stop): continue
      result.append(iselection)
    return result

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

  def get_anisou(self):
    return [self.anisou]

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

  def sel_resSeq(self, pattern):
    return self.union(iselections=self.get_resSeq(pattern=pattern))

  def sel_resSeq_range(self, start, stop):
    return self.union(iselections=self.get_resSeq_range(start=start,stop=stop))

  def sel_iCode(self, pattern):
    return self.union(iselections=self.get_iCode(pattern=pattern))

  def sel_resid(self, pattern):
    return self.union(iselections=self.get_resid(pattern=pattern))

  def sel_resid_range(self, start, stop):
    return self.union(iselections=self.get_resid_range(start=start,stop=stop))

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

  def sel_anisou(self):
    return self.union(iselections=self.get_anisou())

  def selection_tokenizer(self, string, contiguous_word_characters=None):
    return selection_tokenizer(string, contiguous_word_characters)

  def selection_parser(self,
        word_iterator,
        callback=None,
        stop_word=None,
        expect_nonmatching_closing_parenthesis=False):
    result_stack = []
    for word,word_iterator in simple_parser.infix_as_postfix(
          word_iterator=word_iterator,
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
          result_stack.append(
            self.sel_name(pattern=word_iterator.pop_argument(word.value)))
        elif (lword in ["altloc", "altid"]):
          result_stack.append(
            self.sel_altLoc(pattern=word_iterator.pop_argument(word.value)))
        elif (lword == "resname"):
          result_stack.append(
            self.sel_resName(pattern=word_iterator.pop_argument(word.value)))
        elif (lword == "chain"):
          result_stack.append(
            self.sel_chainID(pattern=word_iterator.pop_argument(word.value)))
        elif (lword in ["resseq", "resid", "resi", "model"]):
          arg = word_iterator.pop_argument(word.value)
          i_colon_or_dash = arg.value.find(":")
          if (i_colon_or_dash < 0): i_colon_or_dash = arg.value.find("-")
          if (i_colon_or_dash < 0):
            if (lword == "resseq"):
              result_stack.append(self.sel_resSeq(pattern=arg))
            elif (lword in ["resid", "resi"]):
              result_stack.append(self.sel_resid(pattern=arg))
            else:
              try: i = int(arg.value)
              except ValueError: raise RuntimeError("Value error.")
              result_stack.append(self.sel_MODELserial(i=i))
          else:
            start = arg.value[:i_colon_or_dash]
            stop = arg.value[i_colon_or_dash+1:]
            if (lword == "resseq"):
              result_stack.append(
                self.sel_resSeq_range(start=start, stop=stop))
            elif (lword in ["resid", "resi"]):
              result_stack.append(
                self.sel_resid_range(start=start, stop=stop))
            else:
              if (len(start) == 0):
                i = None
              else:
                try: i = int(start)
                except ValueError: raise RuntimeError("Value error.")
              if (len(stop) == 0):
                j = None
              else:
                try: j = int(stop)
                except ValueError: raise RuntimeError("Value error.")
              result_stack.append(self.sel_MODELserial_range(i=i, j=j))
        elif (lword == "icode"):
          result_stack.append(
            self.sel_iCode(pattern=word_iterator.pop_argument(word.value)))
        elif (lword == "segid"):
          result_stack.append(
            self.sel_segID(pattern=word_iterator.pop_argument(word.value)))
        elif (lword == "element"):
          result_stack.append(
            self.sel_element(pattern=word_iterator.pop_argument(word.value)))
        elif (lword == "charge"):
          result_stack.append(
            self.sel_charge(pattern=word_iterator.pop_argument(word.value)))
        elif (lword == "anisou"):
          result_stack.append(self.sel_anisou())
        elif (callback is not None):
          if (not callback(
                    word=word,
                    word_iterator=word_iterator,
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
    try:
      return self.selection_parser(
        word_iterator=self.selection_tokenizer(
          string=string,
          contiguous_word_characters=contiguous_word_characters),
        callback=callback)
    except (AtomSelectionError, KeyboardInterrupt): raise
    except:
      s = StringIO()
      traceback.print_exc(file=s)
      msg = s.getvalue().splitlines()
      msg.extend([
        "Atom selection string leading to error:",
        "  " + string])
      raise AtomSelectionError("\n".join(msg))

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

  def link_iselections(self, link_record):
    sel_null = stl.vector.unsigned()
    fs = flex.size_t
    return [
                    fs(self.name.get(link_record.name1, sel_null))
      .intersection(fs(self.altLoc.get(link_record.altLoc1, sel_null)))
      .intersection(fs(self.resName.get(link_record.resName1, sel_null)))
      .intersection(fs(self.chainID.get(link_record.chainID1, sel_null)))
      .intersection(fs(self.resSeq.get(link_record.resSeq1, sel_null)))
      .intersection(fs(self.iCode.get(link_record.iCode1, sel_null))),
                    fs(self.name.get(link_record.name2, sel_null))
      .intersection(fs(self.altLoc.get(link_record.altLoc2, sel_null)))
      .intersection(fs(self.resName.get(link_record.resName2, sel_null)))
      .intersection(fs(self.chainID.get(link_record.chainID2, sel_null)))
      .intersection(fs(self.resSeq.get(link_record.resSeq2, sel_null)))
      .intersection(fs(self.iCode.get(link_record.iCode2, sel_null)))]
