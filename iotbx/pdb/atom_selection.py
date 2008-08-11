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

def _get_serial_range(sel_keyword, map, start, stop):
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
      "range with first index > last index: %s %s:%s" % (
        sel_keyword, start, stop))
  result = []
  for s,iselection in map.items():
    os = o(s)
    if (o_start is not None and os < o_start): continue
    if (o_stop  is not None and os > o_stop): continue
    result.append(iselection)
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

  def __init__(self, root, wildcard_escape_char='\\'):
    self.wildcard_escape_char = wildcard_escape_char
    root.get_atom_selection_cache(self)

  def get_name(self, pattern):
    return _get_map_string(
      map=self.name,
      pattern=pattern,
      wildcard_escape_char=self.wildcard_escape_char)

  def get_altloc(self, pattern):
    return _get_map_string(
      map=self.altloc,
      pattern=pattern,
      wildcard_escape_char=self.wildcard_escape_char)

  def get_resname(self, pattern):
    return _get_map_string(
      map=self.resname,
      pattern=pattern,
      wildcard_escape_char=self.wildcard_escape_char)

  def get_chain_id(self, pattern):
    return _get_map_string(
      map=self.chain_id,
      pattern=pattern,
      wildcard_escape_char=self.wildcard_escape_char)

  def get_resseq(self, pattern):
    return _get_map_string(
      map=self.resseq,
      pattern=pattern,
      wildcard_escape_char=self.wildcard_escape_char)

  def get_resseq_range(self, start, stop):
    return _get_serial_range(
      sel_keyword="resseq", map=self.resseq, start=start, stop=stop)

  def get_icode(self, pattern):
    return _get_map_string(
      map=self.icode,
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

  def get_segid(self, pattern):
    return _get_map_string(
      map=self.segid,
      pattern=pattern,
      wildcard_escape_char=self.wildcard_escape_char)

  def get_model_id(self, pattern):
    return _get_map_string(
      map=self.model_id,
      pattern=pattern,
      wildcard_escape_char=self.wildcard_escape_char)

  def get_model_id_range(self, start, stop):
    return _get_serial_range(
      sel_keyword="model", map=self.model_id, start=start, stop=stop)

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

  def sel_altloc(self, pattern):
    return self.union(iselections=self.get_altloc(pattern=pattern))

  def sel_resname(self, pattern):
    return self.union(iselections=self.get_resname(pattern=pattern))

  def sel_chain_id(self, pattern):
    return self.union(iselections=self.get_chain_id(pattern=pattern))

  def sel_resseq(self, pattern):
    return self.union(iselections=self.get_resseq(pattern=pattern))

  def sel_resseq_range(self, start, stop):
    return self.union(iselections=self.get_resseq_range(start=start,stop=stop))

  def sel_icode(self, pattern):
    return self.union(iselections=self.get_icode(pattern=pattern))

  def sel_resid(self, pattern):
    return self.union(iselections=self.get_resid(pattern=pattern))

  def sel_resid_range(self, start, stop):
    return self.union(iselections=self.get_resid_range(start=start,stop=stop))

  def sel_segid(self, pattern):
    return self.union(iselections=self.get_segid(pattern=pattern))

  def sel_model_id(self, pattern):
    return self.union(iselections=self.get_model_id(pattern=pattern))

  def sel_model_id_range(self, start, stop):
    return self.union(iselections=self.get_model_id_range(
      start=start,stop=stop))

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
            self.sel_altloc(pattern=word_iterator.pop_argument(word.value)))
        elif (lword == "resname"):
          result_stack.append(
            self.sel_resname(pattern=word_iterator.pop_argument(word.value)))
        elif (lword == "chain"):
          result_stack.append(
            self.sel_chain_id(pattern=word_iterator.pop_argument(word.value)))
        elif (lword in ["resseq", "resid", "resi", "model"]):
          arg = word_iterator.pop_argument(word.value)
          i_colon_or_dash = arg.value.find(":")
          if (i_colon_or_dash < 0): i_colon_or_dash = arg.value.find("-")
          if (i_colon_or_dash < 0):
            if (lword == "resseq"):
              result_stack.append(self.sel_resseq(pattern=arg))
            elif (lword in ["resid", "resi"]):
              result_stack.append(self.sel_resid(pattern=arg))
            else:
              result_stack.append(self.sel_model_id(pattern=arg))
          else:
            start = arg.value[:i_colon_or_dash]
            stop = arg.value[i_colon_or_dash+1:]
            if (lword == "resseq"):
              result_stack.append(
                self.sel_resseq_range(start=start, stop=stop))
            elif (lword in ["resid", "resi"]):
              result_stack.append(
                self.sel_resid_range(start=start, stop=stop))
            else:
              result_stack.append(
                self.sel_model_id_range(start=start, stop=stop))
        elif (lword == "icode"):
          result_stack.append(
            self.sel_icode(pattern=word_iterator.pop_argument(word.value)))
        elif (lword == "segid"):
          result_stack.append(
            self.sel_segid(pattern=word_iterator.pop_argument(word.value)))
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
        altloc=None,
        resname=None,
        chain_id=None,
        resseq=None,
        icode=None,
        segid=None,
        model_id=None):
    result = []
    for arg,attr in [(name, self.name),
                     (altloc, self.altloc),
                     (resname, self.resname),
                     (chain_id, self.chain_id),
                     (resseq, self.resseq),
                     (icode, self.icode),
                     (segid, self.segid),
                     (model_id, self.model_id)]:
      if (arg is not None):
        isel = attr.get(arg, None)
        if (isel is not None): result.append(isel)
    return result

  def link_iselections(self, link_record):
    sel_null = stl.vector.unsigned()
    fs = flex.size_t
    return [
                    fs(self.name.get(link_record.name1, sel_null))
      .intersection(fs(self.altloc.get(link_record.altloc1, sel_null)))
      .intersection(fs(self.resname.get(link_record.resname1, sel_null)))
      .intersection(fs(self.chain_id.get(link_record.chain_id1, sel_null)))
      .intersection(fs(self.resseq.get(link_record.resseq1, sel_null)))
      .intersection(fs(self.icode.get(link_record.icode1, sel_null))),
                    fs(self.name.get(link_record.name2, sel_null))
      .intersection(fs(self.altloc.get(link_record.altloc2, sel_null)))
      .intersection(fs(self.resname.get(link_record.resname2, sel_null)))
      .intersection(fs(self.chain_id.get(link_record.chain_id2, sel_null)))
      .intersection(fs(self.resseq.get(link_record.resseq2, sel_null)))
      .intersection(fs(self.icode.get(link_record.icode2, sel_null)))]
