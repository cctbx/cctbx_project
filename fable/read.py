from fable \
  import unsigned_integer_scan, \
  identifier_scan, \
  find_closing_parenthesis, \
  SemanticError
from fable import tokenization
from fable import intrinsics
from fable import equivalence
from fable import utils
import sys

class Error(Exception): pass

class raise_errors_mixin(object):

  __slots__ = []

  def text_location(O, i):
    sl, i = O.stmt_location(i)
    if (i is not None and sl.stmt_offs is not None):
      i += sl.stmt_offs
    return sl, i

  def format_error(O, i, msg, prefix=""):
    sl, i = O.stmt_location(i)
    from libtbx.str_utils import expandtabs_track_columns
    t, js = expandtabs_track_columns(s=sl.text)
    if (i is None):
      ptr = ""
    else:
      if (i < 0):
        j = -i - 1
        t += " " * (j - len(t) + 1)
      else:
        j = js[sl.stmt_offs + i]
      ptr = "\n" + "-"*(3+j) + "^"
    if (msg is None): intro = ""
    else:             intro = "%s:\n  at " % msg
    result = "%s%s:\n  |%s|%s" % (
      intro, sl.format_file_name_and_line_number(), t, ptr)
    if (prefix is None or prefix == ""):
      return result
    return "\n".join([(prefix + line).rstrip()
      for line in result.splitlines()])

  def raise_error(O, msg, i=None, ErrorType=None):
    if (ErrorType is None): ErrorType = Error
    raise ErrorType(O.format_error(i=i, msg=msg))

  def raise_syntax_error(O, i=None):
    O.raise_error(msg="Syntax error", i=i)

  def raise_syntax_error_or_not_implemented(O, i=None):
    O.raise_error(msg="Syntax error or not implemented", i=i)

  def raise_semantic_error(O, msg=None, i=None):
    O.raise_error(msg=msg, i=i, ErrorType=SemanticError)

  def raise_internal_error(O, i=None):
    O.raise_error(
      msg="Sorry: fable internal error", i=i, ErrorType=AssertionError)

class source_line(raise_errors_mixin):

  __slots__ = [
    "global_line_index",
    "file_name",
    "line_number",
    "text",
    "label",
    "stmt",
    "stmt_offs",
    "is_cont",
    "index_of_exclamation_mark"]

  def format_file_name_and_line_number(O):
    return "%s(%d)" % (O.file_name, O.line_number)

  def stmt_location(O, i):
    return O, i

  def __init__(O, global_line_index_generator, file_name, line_number, text):
    O.global_line_index = global_line_index_generator.next()
    O.file_name = file_name
    O.line_number = line_number
    O.text = text
    O.label = None
    O.stmt = ""
    O.stmt_offs = None
    i = text.find("\t", 0, 6)
    if (i >= 0):
      soff = i + 1
      O.is_cont = False
      l = text[:i].strip()
      s = text[soff:72]
    else:
      soff = 6
      c = text[5:6]
      O.is_cont = (c != " " and c != "\t" and c != "")
      l = text[:5].strip()
      s = text[6:72]
    if (len(l) == 0):
      if (len(s) != 0):
        O.stmt = s
        O.stmt_offs = soff
    else:
      i = unsigned_integer_scan(code=l)
      if (i < 0 or i != len(l)):
        O.is_cont = False
      else:
        if (O.is_cont):
          O.raise_error(
            msg="A continuation character is illegal on a line with"
                " a statement label",
            i=-6)
        O.label = l
        if (len(s) == 0):
          O.raise_error(msg="Labelled statement is empty", i=-7)
        O.stmt = s
        O.stmt_offs = soff
    if (not O.is_cont and len(O.stmt.rstrip()) == 0):
      O.stmt_offs = None
    O.index_of_exclamation_mark = None

class stripped_source_line(raise_errors_mixin):

  __slots__ = [
    "source_line_cluster",
    "label",
    "code0_locations",
    "start",
    "code",
    "strings",
    "strings_locs",
    "string_indices"]

  def __init__(O,
        source_line_cluster,
        code0_locations,
        code,
        start,
        strings,
        strings_locs,
        string_indices):
    if (source_line_cluster is None):
      assert code0_locations is None
    else:
      assert len(source_line_cluster) != 0
      assert len(code0_locations) >= start + len(code)
    assert len(strings) == len(string_indices)
    O.source_line_cluster = source_line_cluster
    O.label = None
    for sl in O.source_line_cluster:
      if (sl.label is not None):
        O.label = sl.label
        break
    O.code0_locations = code0_locations
    O.start = start
    O.code = code
    O.strings = strings
    O.strings_locs = strings_locs
    O.string_indices = string_indices

  def code_with_strings(O):
    result = []
    j = 0
    for c in O.code:
      if (c == "'"):
        result.append("'" + O.strings[j].replace("'","''") + "'")
        j += 1
      else:
        result.append(c)
    assert j == len(O.strings)
    return "".join(result)

  def stmt_location(O, i):
    if (i is None):
      return O.source_line_cluster[0], None
    if (i < 0):
      return O.source_line_cluster[0], i
    return O.code0_locations[O.start + i]

  def is_comment(O):
    return (O.source_line_cluster[0].stmt_offs is None)

  def __getitem__(O, key):
    if (isinstance(key, slice)):
      start, stop, step = key.indices(len(O.code))
      assert step == 1
      del step
    else:
      start = key
      stop = key + 1
    slice_strings = []
    slice_strings_locs = []
    slice_string_indices = []
    for s,locs,si in zip(O.strings, O.strings_locs, O.string_indices):
      if (si < start): continue
      if (si >= stop): break
      slice_strings.append(s)
      slice_strings_locs.append(locs)
      slice_string_indices.append(si-start)
    return stripped_source_line_slice(
      source_line_cluster=O.source_line_cluster,
      code0_locations=O.code0_locations,
      start=O.start+start,
      code=O.code[key],
      strings=slice_strings,
      strings_locs=slice_strings_locs,
      string_indices=slice_string_indices)

  def raise_if_not_identifier(O):
    i = identifier_scan(code=O.code)
    if (i < 0 or i != len(O.code)):
      O.raise_error("Not an identifier: %s" % repr(O.code), i=0)

  def extract_identifier(O):
    O.raise_if_not_identifier()
    return O.code

  def index_of_closing_parenthesis(O, start=0):
    i = find_closing_parenthesis(code=O.code, start=start)
    if (i < 0):
      O.raise_error(msg='Missing a closing ")"', i=max(0, start-1))
    return i

  def comma_scan(O, start=0):
    code = O.code
    n = len(code)
    i = start
    while (i < n):
      c = code[i]
      if (c == ","):
        return i
      i += 1
      if (c == "("):
        i = O.index_of_closing_parenthesis(start=i) + 1
    return -1

def get_hollerith_count_index(code):
  i = len(code)
  while (i != 0):
    i -= 1
    c = code[i]
    digit = "0123456789".find(c)
    if (digit < 0):
      if (i+1 == len(code)):
        return None
      if (",(/$".find(c) >= 0):
        return i+1
      return None
  return None

class stripped_source_line_slice(stripped_source_line):

  __slots__ = stripped_source_line.__slots__

def strip_spaces_separate_strings(source_line_cluster):
  code = []
  locs = []
  strings = []
  strings_locs = []
  string_indices = []
  ca = code.append
  la = locs.append
  n_sl = len(source_line_cluster)
  i_sl = 0
  while (i_sl < n_sl):
    sl = source_line_cluster[i_sl]
    s = sl.stmt
    n = len(s)
    i = 0
    while (i < n):
      c = s[i]
      if (c == "!"):
        sl.index_of_exclamation_mark = i
        break
      if (c == "'" or c == '"'):
        opening_quote = c
        string_indices.append(len(code))
        ca("'")
        la((sl,i))
        i += 1
        string_chars = []
        string_chars_locs = []
        in_string = True
        while in_string:
          while (i < n):
            c = s[i]
            ci = i
            i += 1
            if (c == opening_quote):
              if (not s.startswith(opening_quote, i)):
                in_string = False
                break
              i += 1
            string_chars.append(c)
            string_chars_locs.append((sl,ci))
          else:
            i_sl += 1
            if (i_sl == n_sl):
              locs[-1][0].raise_error(
                msg="Missing terminating %s character" % opening_quote,
                i=locs[-1][1])
            sl = source_line_cluster[i_sl]
            s = sl.stmt
            n = len(s)
            i = 0
        strings.append("".join(string_chars))
        strings_locs.append(string_chars_locs)
      elif (" \t".find(c) < 0):
        c = c.lower()
        if (c == 'h'):
          j = get_hollerith_count_index(code)
        else:
          j = None
        if (j is None):
          ca(c.lower())
          la((sl,i))
          i += 1
        else:
          hollerith_count = int("".join(code[j:]))
          del code[j:]
          del locs[j:]
          string_indices.append(len(code))
          ca("'")
          la((sl,i))
          i += 1
          string_chars = []
          string_chars_locs = []
          while True:
            if (i < n):
              string_chars.append(s[i])
              string_chars_locs.append((sl,i))
              i += 1
              if (len(string_chars) == hollerith_count):
                break
            else:
              i_sl += 1
              if (i_sl == n_sl):
                break
              sl = source_line_cluster[i_sl]
              s = sl.stmt
              n = len(s)
              i = 0
          if (len(string_chars) != hollerith_count):
            locs[-1][0].raise_error(
              msg="Missing characters for Hollerith constant",
              i=locs[-1][1])
          strings.append("".join(string_chars))
          strings_locs.append(string_chars_locs)
      else:
        i += 1
    i_sl += 1
  return stripped_source_line(
    source_line_cluster=source_line_cluster,
    code0_locations=locs,
    code="".join(code),
    start=0,
    strings=strings,
    strings_locs=strings_locs,
    string_indices=string_indices)

class fmt_string_stripped(raise_errors_mixin):

  __slots__ = ["code", "locs", "strings", "strings_locs", "string_indices"]

  def __init__(O, fmt_tok):
    ssl = fmt_tok.ssl
    i = ssl.string_indices.index(fmt_tok.i_code)
    fmt_string = ssl.strings[i]
    fmt_string_locs = ssl.strings_locs[i]
    assert len(fmt_string) == len(fmt_string_locs)
    code = []
    O.locs = []
    O.strings = []
    O.strings_locs = []
    O.string_indices = []
    ca = code.append
    la = O.locs.append
    n = len(fmt_string)
    have_leading_parenthesis = False
    i = 0
    while (i < n):
      c = fmt_string[i]
      loc = fmt_string_locs[i]
      if (c == "'" or c == '"'):
        if (not have_leading_parenthesis):
          raise_must_start()
        opening_quote = c
        O.string_indices.append(len(code))
        ca("'")
        la(loc)
        i += 1
        string_chars = []
        string_chars_locs = []
        in_string = True
        while in_string:
          while (i < n):
            c = fmt_string[i]
            loc = fmt_string_locs[i]
            i += 1
            if (c == opening_quote):
              if (not fmt_string.startswith(opening_quote, i)):
                in_string = False
                break
              i += 1
            string_chars.append(c)
            string_chars_locs.append(loc)
          else:
            loc = O.locs[-1]
            loc[0].raise_error(
              msg='Missing terminating %s within character format'
                  ' specifier "%s"' % (opening_quote, fmt_string),
              i=loc[1])
        O.strings.append("".join(string_chars))
        O.strings_locs.append(string_chars_locs)
      else:
        if (" \t".find(c) < 0):
          if (have_leading_parenthesis):
            ca(c.lower())
            la(loc)
          else:
            if (c != "("):
              raise_must_start()
            have_leading_parenthesis = True
        i += 1
    def raise_must_start():
      fmt_tok.raise_error(msg='Format string must start with "("')
    def raise_must_end():
      fmt_tok.raise_error(msg='Format string must end with ")"')
    if (len(code) == 0):
      if (have_leading_parenthesis):
        raise_must_end()
      raise_must_start()
    elif (code[-1] != ")"):
      raise_must_end()
    code.pop()
    O.locs.pop()
    O.code = "".join(code)

  def stmt_location(O, i):
    if (i is None): i = 0
    return O.locs[i]

def combine_continuation_lines_and_strip_spaces(source_lines):
  result = []
  rapp = result.append
  n_sl = len(source_lines)
  i_sl = 0
  while (i_sl < n_sl):
    sl = source_lines[i_sl]
    if (sl.stmt_offs is None):
      rapp(strip_spaces_separate_strings(source_line_cluster=[sl]))
      i_sl += 1
    else:
      assert not sl.is_cont
      code_sls = [sl]
      k_sl = i_sl
      for j_sl in xrange(i_sl+1, n_sl):
        sl = source_lines[j_sl]
        if (sl.is_cont):
          code_sls.append(sl)
          k_sl = j_sl
        elif (sl.stmt_offs is not None):
          break
      for j_sl in xrange(i_sl+1, k_sl):
        sl = source_lines[j_sl]
        if (not sl.is_cont):
          rapp(strip_spaces_separate_strings(source_line_cluster=[sl]))
      rapp(strip_spaces_separate_strings(source_line_cluster=code_sls))
      i_sl = k_sl + 1
  return result

def load_includes(global_line_index_generator, stripped_source_lines):
  import os.path as op
  result = []
  for ssl in stripped_source_lines:
    if (ssl.code == "include'"):
      assert len(ssl.strings) == 1
      file_name = ssl.strings[0]
      if (op.isabs(file_name)):
        file_path = file_name
      else:
        sl = ssl.code0_locations[-1][0]
        file_path = op.join(op.dirname(sl.file_name), file_name)
      if (not op.isfile(file_path)):
        ssl.raise_semantic_error(msg="Missing include file", i=7)
      # TODO potential performance problem if deeply nested includes
      result.extend(load(
        global_line_index_generator=global_line_index_generator,
        file_name=file_path))
    else:
      result.append(ssl)
  return result

def load(global_line_index_generator, file_name, skip_load_includes=False):
  source_lines = []
  for i_line,line in enumerate(open(file_name).read().splitlines()):
    source_lines.append(source_line(
      global_line_index_generator=global_line_index_generator,
      file_name=file_name,
      line_number=i_line+1,
      text=line))
  stripped_source_lines = combine_continuation_lines_and_strip_spaces(
    source_lines=source_lines)
  if (skip_load_includes):
    return stripped_source_lines
  return load_includes(
    global_line_index_generator=global_line_index_generator,
    stripped_source_lines=stripped_source_lines)

def tokenize_expression(
      ssl,
      start=0,
      stop=None,
      allow_commas=False,
      allow_equal_signs=False):
  result = []
  if (stop is None): stop = len(ssl.code)
  tokenize_expression_impl(
    tokens=result,
    tokenizer=tokenization.ssl_iterator(ssl=ssl, start=start, stop=stop),
    allow_commas=allow_commas,
    allow_equal_signs=allow_equal_signs,
    tok_opening_parenthesis=None)
  return result

def tokenize_expression_impl(
      tokens,
      tokenizer,
      allow_commas,
      allow_equal_signs,
      tok_opening_parenthesis):
  from tokenization import tk_seq, tk_parentheses
  if (allow_commas):
    tlist = []
    tokens.append(tk_seq(ssl=tokenizer.ssl, i_code=tokenizer.i, value=tlist))
  else:
    tlist = tokens
  tapp = tlist.append
  for tok in tokenizer:
    if (tok.is_op()):
      tv = tok.value
      if (tv == "("):
        nested_tokens = []
        tokenize_expression_impl(
          tokens=nested_tokens,
          tokenizer=tokenizer,
          allow_commas=True,
          allow_equal_signs=allow_equal_signs,
          tok_opening_parenthesis=tok)
        tapp(tk_parentheses(
          ssl=tok.ssl, i_code=tok.i_code, value=nested_tokens))
        continue
      if (tv == ")"):
        if (tok_opening_parenthesis is None):
          tok.raise_missing_opening()
        return
      if (tv == ","):
        if (not allow_commas):
          tok.ssl.raise_syntax_error(i=tok.i_code)
        tlist = []
        tokens.append(tk_seq(
          ssl=tokenizer.ssl, i_code=tokenizer.i, value=tlist))
        tapp = tlist.append
        continue
      if (tv == "="):
        if (not allow_equal_signs):
          tok.ssl.raise_syntax_error(i=tok.i_code)
        tapp(tok)
        continue
    tapp(tok)
    continue
  if (tok_opening_parenthesis is not None):
    tok_opening_parenthesis.raise_missing_closing()

def indices_of_tokenized_equal_signs(tokens):
  result = []
  for i,tok in enumerate(tokens):
    if (tok.is_op() and tok.value == "="):
      result.append(i)
  return result

# variable types
class vt_used(object): pass
class vt_scalar(object): pass
class vt_string(object): pass
class vt_array(object): pass
class vt_intrinsic(object): pass
class vt_external(object): pass
class vt_function(object): pass
class vt_subroutine(object): pass

# variable storage
class vs_fproc_name(object): pass
class vs_argument(object): pass
class vs_common(object): pass
class vs_save(object): pass
class vs_local(object): pass
class vs_parameter(object): pass

class fdecl_info(object):

  __slots__ = [
    "id_tok",
    "var_type",
    "var_storage",
    "data_type",
    "size_tokens",
    "dim_tokens",
    "parameter_assignment_tokens",
    "f90_decl",
    "is_modified",
    "use_count",
    "passed_as_arg",
    "passed_as_arg_plain"]

  def __init__(O,
        id_tok,
        var_type,
        var_storage,
        data_type,
        size_tokens,
        dim_tokens,
        f90_decl=None):
    assert size_tokens is None or isinstance(size_tokens, list)
    assert dim_tokens is None or isinstance(dim_tokens, list)
    O.id_tok = id_tok
    O.var_type = var_type
    O.var_storage = var_storage
    O.data_type = data_type
    O.size_tokens = size_tokens
    O.dim_tokens = dim_tokens
    O.parameter_assignment_tokens = None
    O.f90_decl = f90_decl
    O.is_modified = False
    O.use_count = 0
    O.passed_as_arg = {}
    O.passed_as_arg_plain = {}

  def is_used(O): return O.var_type is vt_used
  def is_scalar(O): return O.var_type is vt_scalar
  def is_string(O): return O.var_type is vt_string
  def is_array(O): return O.var_type is vt_array
  def is_intrinsic(O): return O.var_type is vt_intrinsic
  def is_external(O): return O.var_type is vt_external
  def is_function(O): return O.var_type is vt_function
  def is_subroutine(O): return O.var_type is vt_subroutine
  def is_user_defined_callable(O):
    vt = O.var_type
    return (vt is vt_external or vt is vt_function or vt is vt_subroutine)

  def is_fproc_name(O): return O.var_storage is vs_fproc_name
  def is_argument(O): return O.var_storage is vs_argument
  def is_common(O): return O.var_storage is vs_common
  def is_save(O): return O.var_storage is vs_save
  def is_local(O): return O.var_storage is vs_local
  def is_parameter(O): return O.var_storage is vs_parameter

  def required_parameter_assignment_tokens(O):
    result = O.parameter_assignment_tokens
    if (result is None):
      O.id_tok.raise_internal_error()
    return result

def extract_size_tokens(ssl, start):
  code = ssl.code
  c = code[start]
  if (c == "("):
    i_clp = ssl.index_of_closing_parenthesis(start=start+1)
    return i_clp+1, tokenize_expression(
      ssl=ssl,
      start=start+1,
      stop=i_clp)
  i_size = unsigned_integer_scan(code=code, start=start)
  if (i_size < 0):
    ssl.raise_syntax_error(i=start)
  return \
    i_size, \
    [tokenization.tk_integer(ssl=ssl, i_code=start, value=code[start:i_size])]

data_types = """\
byte
character
complex
doublecomplex
doubleprecision
integer
logical
real
""".splitlines()

def extract_data_type(ssl, start=0, optional=False):
  sw = ssl.code.startswith
  for data_type in data_types:
    if (sw(data_type, start)):
      return (
        start + len(data_type),
        tokenization.tk_identifier(ssl=ssl, i_code=start, value=data_type))
  if (not optional):
    ssl.raise_syntax_error()
  return None, None

def extract_data_type_and_size(ssl, start=0, optional=False):
  i_code, data_type = extract_data_type(
    ssl=ssl, start=start, optional=optional)
  if (optional and i_code is None):
    return None, None, None
  if (not ssl.code.startswith("*", i_code)):
    return i_code, data_type, None
  i_code, size_tokens = extract_size_tokens(ssl=ssl, start=i_code+1)
  return i_code, data_type, size_tokens

def extract_f90_decl(ssl, start):
  code = ssl.code
  if (start == len(code)):
    ssl.raise_syntax_error(i=start)
  i_cc = code.find("::", start)
  if (i_cc >= 0):
    return i_cc+2, ssl[start:i_cc]
  if (code.startswith("(", start)):
    i_clp = ssl.index_of_closing_parenthesis(start=start+1)
    return i_clp+1, ssl[start+1:i_clp]
  return start, None

def extract_fdecl(
      result,
      ssl,
      start,
      data_type,
      size_tokens,
      allow_size,
      f90_decl=None):
  code = ssl.code
  stop = len(code)
  def parse_decl(start):
    item_size_tokens = None
    dim_tokens = None
    i_id = identifier_scan(code=code, start=start)
    if (i_id < 0):
      ssl.raise_syntax_error(i=start)
    i_code = i_id
    while True:
      if (i_code == stop):
        break
      c = code[i_code]
      if (c == ","):
        i_code += 1
        break
      if (c == "("):
        if (dim_tokens is not None):
          ssl.raise_syntax_error(i=i_code)
        i_clp = ssl.index_of_closing_parenthesis(start=i_code+1)
        dim_tokens = tokenize_expression(
          ssl=ssl,
          start=i_code+1,
          stop=i_clp,
          allow_commas=True)
        i_code = i_clp + 1
      elif (c == "*"):
        if (not allow_size or item_size_tokens is not None):
          ssl.raise_syntax_error(i=i_code)
        i_code, item_size_tokens = extract_size_tokens(ssl=ssl, start=i_code+1)
      else:
        ssl.raise_syntax_error(i=i_code)
    if (item_size_tokens is None):
      item_size_tokens = size_tokens
    result.append(fdecl_info(
      id_tok=tokenization.tk_identifier(
        ssl=ssl, i_code=start, value=code[start:i_id]),
      var_type=None,
      var_storage=None,
      data_type=data_type,
      size_tokens=item_size_tokens,
      dim_tokens=dim_tokens,
      f90_decl=f90_decl))
    return i_code
  if (ssl.code.startswith(",", start)):
    start += 1
  while (start < stop):
    start = parse_decl(start=start)

def dimensions_are_simple(dim_tokens):
  is_star = tokenization.tok_seq_is_star
  for tok_seq in dim_tokens:
    if (is_star(tok_seq=tok_seq)):
      return False
    for tok in tok_seq.value:
      if (tok.is_op_with(value=":")):
        return False
  return True

def process_labels_list(ssl, start, stop, len_min, len_max):
  if (start == stop):
    ssl.raise_syntax_error(i=start)
  result = []
  code = ssl.code
  i = start
  while True:
    if (len_max is not None and len(result) >= len_max):
      ssl.raise_syntax_error(i=i)
    j = unsigned_integer_scan(code=code, start=i, stop=stop)
    if (j < 0):
      ssl.raise_syntax_error(i=i)
    result.append(
      tokenization.tk_integer(ssl=ssl, i_code=i, value=code[i:j]))
    if (j == stop):
      break
    if (code[j] != ","):
      ssl.raise_syntax_error(i=j)
    i = j + 1
  if (len(result) < len_min):
    ssl.raise_syntax_error(i=i)
  return result

class executable_info(object):

  __slots__ = []

  def __init__(O, **ks):
    O.key = O.__class__.__name__[3:]
    O.ssl = ks["ssl"]
    O.start = ks["start"]
    for k,v in ks.items():
      setattr(O, k, v)

  def s4it(O, callback, tokens):
    tokenization.search_for_id_tokens(
      callback=callback, tokens=tokens, with_next_token=True)

  def s4it_slots(O, callback, obj_with_slots):
    for s in obj_with_slots.__slots__:
      attr = getattr(obj_with_slots, s)
      if (attr is not None):
        O.s4it(callback, attr)

  def set_is_modified(O, fdecl_by_identifier):
    pass

def mksl(*names): return ("key", "ssl", "start") + names

class ei_allocate(executable_info):
  __slots__ = mksl()

  def search_for_id_tokens(O, callback):
    pass # TODO

class ei_assign(executable_info):
  __slots__ = mksl()

  def search_for_id_tokens(O, callback):
    pass # TODO

class ei_assignment(executable_info):
  __slots__ = mksl("lhs_tokens", "rhs_tokens")

  def search_for_id_tokens(O, callback):
    O.s4it(callback, O.lhs_tokens)
    O.s4it(callback, O.rhs_tokens)

  def set_is_modified(O, fdecl_by_identifier):
    assert len(O.lhs_tokens) != 0
    id_tok = O.lhs_tokens[0]
    if (not id_tok.is_identifier()):
      id_tok.raise_syntax_error()
    tf = fdecl_by_identifier.get(id_tok.value)
    assert tf is not None
    tf.is_modified = True

class ei_file_positioning(executable_info):
  __slots__ = mksl("io_function", "alist")

  def search_for_id_tokens(O, callback):
    if (O.alist is not None):
      O.s4it_slots(callback, O.alist)

class ei_call(executable_info):
  __slots__ = mksl("subroutine_name", "arg_token")

  def search_for_id_tokens(O, callback):
    callback(O.subroutine_name, O.arg_token)
    if (O.arg_token is not None):
      O.s4it(callback, O.arg_token.value)

class ei_close(executable_info):
  __slots__ = mksl("cllist")

  def search_for_id_tokens(O, callback):
    O.s4it_slots(callback, O.cllist)

class ei_continue(executable_info):
  __slots__ = mksl()

  def search_for_id_tokens(O, callback):
    pass

class ei_cycle(executable_info):
  __slots__ = mksl()

  def search_for_id_tokens(O, callback):
    pass

class ei_deallocate(executable_info):
  __slots__ = mksl()

  def search_for_id_tokens(O, callback):
    pass # TODO

class ei_do(executable_info):
  __slots__ = mksl("label", "id_tok", "tokens")

  def search_for_id_tokens(O, callback):
    callback(O.id_tok, None)
    O.s4it(callback, O.tokens)

  def set_is_modified(O, fdecl_by_identifier):
    fdecl = fdecl_by_identifier[O.id_tok.value]
    fdecl.is_modified = True

class ei_dowhile(executable_info):
  __slots__ = mksl("label", "cond_tokens")

  def search_for_id_tokens(O, callback):
    O.s4it(callback, O.cond_tokens)

class ei_else(executable_info):
  __slots__ = mksl()

  def search_for_id_tokens(O, callback):
    pass

class ei_elseif_then(executable_info):
  __slots__ = mksl("cond_tokens")

  def search_for_id_tokens(O, callback):
    O.s4it(callback, O.cond_tokens)

class ei_enddo(executable_info):
  __slots__ = mksl()

  def search_for_id_tokens(O, callback):
    pass

class ei_endif(executable_info):
  __slots__ = mksl()

  def search_for_id_tokens(O, callback):
    pass

class ei_entry(executable_info):
  __slots__ = mksl()

  def search_for_id_tokens(O, callback):
    pass # TODO

class ei_exit(executable_info):
  __slots__ = mksl()

  def search_for_id_tokens(O, callback):
    pass

class ei_goto(executable_info):
  __slots__ = mksl("label")

  def search_for_id_tokens(O, callback):
    pass

class ei_goto_computed(executable_info):
  __slots__ = mksl("labels", "tokens")

  def search_for_id_tokens(O, callback):
    O.s4it(callback, O.tokens)

class ei_if(executable_info):
  __slots__ = mksl("cond_tokens")

  def search_for_id_tokens(O, callback):
    O.s4it(callback, O.cond_tokens)

class ei_if_then(executable_info):
  __slots__ = mksl("cond_tokens")

  def search_for_id_tokens(O, callback):
    O.s4it(callback, O.cond_tokens)

class ei_if_arithmetic(executable_info):
  __slots__ = mksl("cond_tokens", "labels")

  def search_for_id_tokens(O, callback):
    O.s4it(callback, O.cond_tokens)

class ei_inquire(executable_info):
  __slots__ = mksl("iuflist")

  def search_for_id_tokens(O, callback):
    O.s4it_slots(callback, O.iuflist)

class ei_open(executable_info):
  __slots__ = mksl("olist")

  def search_for_id_tokens(O, callback):
    O.s4it_slots(callback, O.olist)

class ei_print(executable_info):
  __slots__ = mksl("cilist", "iolist", "fmt_tokens")

  def search_for_id_tokens(O, callback):
    O.s4it_slots(callback, O.cilist)
    O.s4it(callback, O.iolist)

class ei_read(executable_info):
  __slots__ = mksl("cilist", "iolist", "fmt_tokens")

  def search_for_id_tokens(O, callback):
    if (O.cilist is not None):
      O.s4it_slots(callback, O.cilist)
    if (O.iolist is not None):
      O.s4it(callback, O.iolist)

  def set_is_modified(O, fdecl_by_identifier):
    if (O.iolist is not None):
      def callback(tok):
        fdecl = fdecl_by_identifier[tok.value]
        fdecl.is_modified = True
      tokenization.search_for_data_or_read_target_tokens(
        callback=callback, tokens=O.iolist)

class ei_return(executable_info):
  __slots__ = mksl("return_label")

  def search_for_id_tokens(O, callback):
    pass

class ei_stop(executable_info):
  __slots__ = mksl("arg_token")

  def search_for_id_tokens(O, callback):
    pass

class ei_write(executable_info):
  __slots__ = mksl("cilist", "iolist", "fmt_tokens")

  def search_for_id_tokens(O, callback):
    O.s4it_slots(callback, O.cilist)
    O.s4it(callback, O.iolist)

  def set_is_modified(O, fdecl_by_identifier):
    if (    O.cilist is not None
        and O.cilist.unit is not None
        and len(O.cilist.unit) != 0):
      first_tok = O.cilist.unit[0]
      if (first_tok.is_identifier()):
        fdecl = fdecl_by_identifier[first_tok.value]
        if (    fdecl.data_type is not None
            and fdecl.data_type.value == "character"):
          fdecl.is_modified = True

del mksl

class fproc_p_methods(object):
  "Separated from class fproc for clarity and a minor getattr speed gain."

  __slots__ = []

  def p_allocate(O, ssl, start):
    O.executable.append(ei_allocate(ssl=ssl, start=start)) # TODO

  def p_assign(O, ssl, start):
    O.executable.append(ei_assign(ssl=ssl, start=start)) # TODO

  def p_file_positioning(O, ssl, start, io_function):
    liof = len(io_function)
    if (ssl.code.startswith("(", start+liof)):
      tz = tokenization.ssl_iterator(ssl=ssl, start=start+liof)
      alist = collect_io_alist(tz=tz, unit=None)
      if (alist.unit is None):
        ssl.raise_semantic_error(
          msg="Required UNIT information is not defined", i=start)
    else:
      alist = collect_io_alist(
        tz=None,
        unit=tokenize_expression(ssl=ssl, start=start+liof))
    O.executable.append(ei_file_positioning(
      ssl=ssl, start=start, io_function=io_function, alist=alist))

  def p_backspace(O, ssl, start):
    O.p_file_positioning(ssl=ssl, start=start, io_function="backspace")

  def p_call(O, ssl, start):
    tokens = tokenize_expression(ssl=ssl, start=start+4)
    if (   len(tokens) == 0
        or len(tokens) > 2
        or not tokens[0].is_identifier()):
      ssl.raise_syntax_error()
    subroutine_name = tokens[0]
    if (len(tokens) == 1):
      arg_token = None
    else:
      if (not tokens[1].is_parentheses()):
        ssl.raise_syntax_error()
      arg_token = tokens[1]
    O.executable.append(ei_call(ssl=ssl, start=start,
      subroutine_name=subroutine_name,
      arg_token=arg_token))

  def p_close(O, ssl, start):
    tz = tokenization.ssl_iterator(ssl=ssl, start=start+5)
    cllist = collect_io_cllist(tz=tz)
    tok = tz.look_ahead(optional=True)
    if (tok is not None):
      tok.raise_syntax_error()
    O.executable.append(ei_close(ssl=ssl, start=start, cllist=cllist))
    O.uses_io = True

  def p_common(O, ssl, start):
    assert start == 0
    code = ssl.code
    if (len(code) == 6):
      ssl.raise_syntax_error()
    c = code[6]
    if (c == "/"):
      i = code.find("/", 7)
      if (i < 0):
        ssl.raise_syntax_error_or_not_implemented()
      if (i == 7):
        common_name = "commonymous"
        i_code = 8
      else:
        common_name = ssl[7:i].extract_identifier()
        i_code = i + 1
    else:
      common_name = "commonymous"
      i_code = 6
    extract_fdecl(
      result=O.common.get(key=common_name),
      ssl=ssl,
      start=i_code,
      data_type=None,
      size_tokens=None,
      allow_size=False)

  def p_continue(O, ssl, start):
    O.executable.append(ei_continue(ssl=ssl, start=start))

  def p_cycle(O, ssl, start):
    if (len(ssl.code) != start+5):
      ssl.raise_syntax_error(i=start+5)
    O.executable.append(ei_cycle(ssl=ssl, start=start))

  def p_data(O, ssl, start):
    assert start == 0
    tz = tokenization.ssl_iterator(ssl=ssl, start=4)
    tok = None
    while True: # loop over nlist, clist pairs
      nlist = []
      while True:
        if (tok is None):
          tok = tz.get()
        if (tok.is_identifier()):
          ntoks = [tok]
          nlist.append(
            tokenization.tk_seq(ssl=tz.ssl, i_code=tz.i, value=ntoks))
          tok = tz.get()
          if (tok.is_op_with(value="(")):
            tz.collect_to_matching_parenthesis(
              callback=ntoks.append, opening_token=tok)
            tok = tz.get()
            if (tok.is_op_with(value="(")):
              tz.collect_to_matching_parenthesis(
                callback=ntoks.append, opening_token=tok)
              tok = tz.get()
        elif (tok.is_op_with(value="(")):
          ntoks = []
          nlist.append(
            tokenization.tk_seq(ssl=tz.ssl, i_code=tz.i, value=ntoks))
          ntoks.append(tz.get_implied_do(opening_token=tok))
          tok = tz.get()
        else:
          tok.raise_syntax_error()
        if (not tok.is_op_with(value=",")):
          break
        tok = None
      if (not tok.is_op_with(value="/")):
        tok.raise_syntax_error()
      clist = []
      repetition_tok = None
      sign_count = 0
      ctoks = []
      while True:
        tok = tz.get()
        if (    len(ctoks) == 0
            and repetition_tok is None
            and (tok.is_integer() or tok.is_identifier())
            and tz.look_ahead().is_op_with(value="*")):
          repetition_tok = tok
          tz.get()
          tok = tz.get()
        if (tok.is_op()):
          if (tok.value in ["+", "-"]):
            if (sign_count != 0 or len(ctoks) != 0):
              tok.raise_syntax_error()
            sign_count = 1
            ctoks.append(tok)
          elif (tok.value == "("):
            if (len(ctoks) != sign_count):
              tok.raise_syntax_error()
            ctoks.append(tz.get_complex_literal(opening_token=tok))
          elif (tok.value == "/"):
            if (len(ctoks) == sign_count):
              tok.raise_syntax_error()
            clist.append((repetition_tok, ctoks))
            break
          elif (tok.value == ","):
            if (len(ctoks) == sign_count):
              tok.raise_syntax_error()
            clist.append((repetition_tok, ctoks))
            repetition_tok = None
            sign_count = 0
            ctoks = []
          else:
            tok.raise_syntax_error()
        else:
          if (len(ctoks) != sign_count):
            tok.raise_syntax_error()
          ctoks.append(tok)
      O.data.append((nlist, clist))
      tok = tz.get(optional=True)
      if (tok is None):
        break
      if (tok.is_op_with(value=",")):
        tok = None

  def p_deallocate(O, ssl, start):
    O.executable.append(ei_deallocate(ssl=ssl, start=start)) # TODO

  def p_dimension(O, ssl, start):
    assert start == 0
    extract_fdecl(
      result=O.dimension,
      ssl=ssl,
      start=9,
      data_type=None,
      size_tokens=None,
      allow_size=False)

  def p_do(O, ssl, start):
    assert start == 0
    code = ssl.code
    i = unsigned_integer_scan(code=code, start=2)
    if (i < 3):
      i = 2
      label = None
    else:
      label = code[2:i]
      if (code[i] == ","):
        i += 1
    j = identifier_scan(code=code, start=i)
    assert j >= 3
    assert code[j] == "="
    tokens = tokenize_expression(ssl=ssl, start=j+1, allow_commas=True)
    if (not (2 <= len(tokens) <= 3)):
      ssl.raise_syntax_error(i=j+1)
    O.executable.append(ei_do(ssl=ssl, start=start,
      label=label,
      id_tok=tokenization.tk_identifier(ssl=ssl, i_code=i, value=code[i:j]),
      tokens=tokens))

  def p_dowhile(O, ssl, start, label_end=None):
    assert start == 0
    if (label_end is None):
      i = 7
      label = None
    else:
      i = label_end + 5
      label = ssl.code[2:label_end]
    cond_tokens = tokenize_expression(ssl=ssl, start=i)
    if (len(cond_tokens) != 1):
      ssl.raise_syntax_error(i=i)
    O.executable.append(ei_dowhile(ssl=ssl, start=start,
      label=label,
      cond_tokens=cond_tokens))

  def p_else(O, ssl, start):
    assert start == 0
    O.executable.append(ei_else(ssl=ssl, start=start))

  def p_elseif(O, ssl, start):
    assert start == 0
    O.p_if_elseif(ssl=ssl, keyword="elseif", start=0)

  def p_enddo(O, ssl, start):
    assert start == 0
    O.executable.append(ei_enddo(ssl=ssl, start=start))

  def p_endfile(O, ssl, start):
    O.p_file_positioning(ssl=ssl, start=start, io_function="endfile")

  def p_endif(O, ssl, start):
    assert start == 0
    O.executable.append(ei_endif(ssl=ssl, start=start))

  def p_entry(O, ssl, start):
    assert start == 0
    O.executable.append(ei_entry(ssl=ssl, start=start)) # TODO

  def p_equivalence(O, ssl, start):
    assert start == 0
    buffer = []
    tz = tokenization.ssl_iterator(ssl=ssl, start=11)
    while True:
      tok = tz.get()
      if (not tok.is_op_with(value="(")):
        tok.raise_syntax_error()
      def callback(tok):
        if (len(tok.value) == 0):
          tok.raise_syntax_error()
        for tok_seq in tok.value:
          if (len(tok_seq.value) == 0):
            tok.raise_syntax_error()
          id_tok = tok_seq.value[0]
          if (not id_tok.is_identifier()):
            id_tok.raise_syntax_error()
        O.equivalence.append(tok)
      tz.collect_to_matching_parenthesis(
        callback=callback,
        opening_token=tok)
      tok = tz.get(optional=True)
      if (tok is None):
        break
      if (not tok.is_op_with(value=",")):
        tok.raise_syntax_error()

  def p_exit(O, ssl, start):
    if (len(ssl.code) != start+4):
      ssl.raise_syntax_error(i=start+4)
    O.executable.append(ei_exit(ssl=ssl, start=start))

  def p_external(O, ssl, start):
    assert start == 0
    tokenization.ssl_iterator(
      ssl=ssl, start=8).collect_comma_separated_identifiers(
        callback=O.external.append, one_required=True)

  def p_format(O, ssl, start):
    assert start == 0
    code = ssl.code
    assert code.startswith("format(")
    assert code.endswith(")")
    if (ssl.label is None):
      ssl.raise_error(
        msg="FORMAT without a statement label in columns 1-5", i=0)
    if (ssl.label in O.format):
      ssl.raise_error(
        msg="Duplicate statement label in columns 1-5", i=-1)
    O.format[ssl.label] = list(tokenization.fss_iterator(fss=ssl[7:-1]))

  def p_goto(O, ssl, start):
    code = ssl.code
    i = start + 4
    if (i == len(code)):
      ssl.raise_syntax_error(i=i)
    j = unsigned_integer_scan(code=code, start=i)
    if (j == len(code)):
      O.executable.append(ei_goto(ssl=ssl, start=start,
        label=tokenization.tk_integer(ssl=ssl, i_code=i, value=code[i:])))
      return
    if (j > 0):
      ssl.raise_syntax_error(i=i)
    if (code[i] == "("):
      # GO TO (s [,s]...)[,] i
      j = ssl.index_of_closing_parenthesis(start=i+1)
      labels = process_labels_list(
        ssl=ssl, start=i+1, stop=j, len_min=1, len_max=None)
      j += 1
      if (j == len(code)):
        ssl.raise_syntax_error(i=j)
      if (code[j] == ","):
        j += 1
        if (j == len(code)):
          ssl.raise_syntax_error(i=j)
      tokens = tokenize_expression(ssl=ssl, start=j)
      O.executable.append(ei_goto_computed(ssl=ssl, start=start,
        labels=labels,
        tokens=tokens))
      return
    # GO TO i [[,](s [,s]...)]
    if (code[-1] != ")"):
      ssl.raise_syntax_error()
    k = code.rfind("(")
    if (k < 0):
      ssl.raise_syntax_error()
    j = k - 1
    if (code[j] == ","):
      j -= 1
    tokens = tokenize_expression(ssl=ssl, start=i+1, stop=j+1)
    labels = process_labels_list(
      ssl=ssl, start=k+1, stop=len(code)-1, len_min=1, len_max=None)
    O.executable.append(ei_goto_computed(ssl=ssl, start=start,
      labels=labels,
      tokens=tokens))

  def p_if(O, ssl, start):
    i = O.p_if_elseif(ssl=ssl, keyword="if", start=start)
    if (i is not None):
      O.process_body_line(ssl=ssl, start=i)

  def p_if_elseif(O, ssl, keyword, start):
    i_open = start + len(keyword) + 1
    i_clp = ssl.index_of_closing_parenthesis(start=i_open)
    cond_tokens = tokenize_expression(ssl=ssl, start=i_open, stop=i_clp)
    if (len(cond_tokens) == 0):
      ssl.raise_syntax_error(i=i_open)
    code = ssl.code
    if (code.startswith("then", i_clp+1) and len(code) == i_clp+5):
      if (start != 0):
        ssl.raise_syntax_error()
      if (keyword == "if"): ei = ei_if_then
      else:                 ei = ei_elseif_then
      O.executable.append(ei(ssl=ssl, start=start, cond_tokens=cond_tokens))
      return None
    if (keyword != "if"):
      ssl.raise_syntax_error()
    i = i_clp + 1
    if (i == len(code)):
      ssl.raise_syntax_error(i=i)
    j = unsigned_integer_scan(code=code, start=i, stop=i+1)
    if (j < 0):
      if (start != 0):
        ssl.raise_syntax_error()
      O.executable.append(ei_if(ssl=ssl, start=start, cond_tokens=cond_tokens))
      return i
    labels = process_labels_list(
      ssl=ssl, start=i, stop=len(code), len_min=3, len_max=3)
    O.executable.append(ei_if_arithmetic(ssl=ssl, start=start,
      cond_tokens=cond_tokens,
      labels=labels))
    return None

  def p_implicit(O, ssl, start):
    assert start == 0
    if (ssl.code == "implicitnone"):
      O.implicit = {}
      return
    i_code, data_type = extract_data_type(ssl=ssl, start=8)
    if (   not ssl.code.startswith("(", i_code)
        or not ssl.code.endswith(")")):
      ssl.raise_syntax_error_or_not_implemented()
    letters = "abcdefghijklmnopqrstuvwxyz"
    def get(c):
      i = letters.find(c)
      if (i < 0):
        ssl.raise_syntax_error_or_not_implemented()
      return i
    for part in ssl.code[i_code+1:-1].split(","):
      if (len(part) == 3 and part[1] == "-"):
        i = get(part[0])
        j = get(part[2])
        for c in letters[i:j+1]:
          O.implicit[c] = data_type
      else:
        for c in part:
          if (c not in "abcdefghijklmnopqrstuvwxyz"):
            ssl.raise_syntax_error_or_not_implemented()
          O.implicit[c] = data_type

  def p_inquire(O, ssl, start):
    tz = tokenization.ssl_iterator(ssl=ssl, start=start+7)
    iuflist = collect_io_iuflist(tz=tz)
    tok = tz.look_ahead(optional=True)
    if (tok is not None):
      tok.raise_syntax_error()
    O.executable.append(ei_inquire(ssl=ssl, start=start, iuflist=iuflist))
    O.uses_io = True

  def p_intrinsic(O, ssl, start):
    assert start == 0
    tokenization.ssl_iterator(
      ssl=ssl, start=9).collect_comma_separated_identifiers(
        callback=O.intrinsic.append, one_required=True)

  def p_open(O, ssl, start):
    tz = tokenization.ssl_iterator(ssl=ssl, start=start+4)
    olist = collect_io_olist(tz=tz)
    tok = tz.look_ahead(optional=True)
    if (tok is not None):
      tok.raise_syntax_error()
    O.executable.append(ei_open(ssl=ssl, start=start, olist=olist))
    O.uses_io = True

  def p_parameter(O, ssl, start):
    assert start == 0
    code = ssl.code
    if (   not code.startswith("(", 9)
        or not code.endswith(")")):
      ssl.raise_syntax_error()
    tokens_ll = tokenize_expression(
      ssl=ssl,
      start=10,
      stop=len(code)-1,
      allow_commas=True,
      allow_equal_signs=True)
    for tokens_l in tokens_ll:
      i_equal_signs = indices_of_tokenized_equal_signs(tokens=tokens_l.value)
      if (len(i_equal_signs) != 1 or i_equal_signs[0] != 1):
        ssl.raise_syntax_error()
      key_token = tokens_l.value[0]
      if (not key_token.is_identifier()):
        key_token.raise_syntax_error()
      O.parameter.append((key_token, tokens_l.value[2:]))

  def p_print(O, ssl, start):
    tz = tokenization.ssl_iterator(ssl=ssl, start=start+5)
    fmt_buffer = []
    tz.collect_comma_separated_expressions(
      callback=fmt_buffer.append,
      first_get_optional=False,
      stop_after_given_number_of_commas=1)
    assert len(fmt_buffer) == 1
    cilist = collect_io_cilist(tz=None, fmt=fmt_buffer[0].value)
    fmt_tokens = None
    if (len(cilist.fmt) == 1 and cilist.fmt[0].is_string()):
      fmt_tokens = list(tokenization.fss_iterator(
        fss=fmt_string_stripped(fmt_tok=cilist.fmt[0])))
    iolist = collect_iolist(tz=tz)
    O.executable.append(ei_print(ssl=ssl, start=start,
      cilist=cilist, fmt_tokens=fmt_tokens, iolist=iolist))
    O.uses_io = True
    O.uses_write = True

  def p_read_write(O, ssl, start, ei_type):
    tz = tokenization.ssl_iterator(ssl=ssl, start=start)
    cilist = collect_io_cilist(tz=tz)
    if (cilist.unit is None):
      ssl.raise_semantic_error(
        msg="Required UNIT information is not defined", i=start)
    fmt_tokens = None
    if (cilist.fmt is not None):
      if (len(cilist.fmt) == 1 and cilist.fmt[0].is_string()):
        fmt_tokens = list(tokenization.fss_iterator(
          fss=fmt_string_stripped(fmt_tok=cilist.fmt[0])))
    if (ei_type is ei_write and cilist.end is not None):
      cilist.end[0].raise_semantic_error(
        msg="END is invalid for WRITE statements")
    iolist = collect_iolist(tz=tz)
    O.executable.append(ei_type(ssl=ssl, start=start,
      cilist=cilist, fmt_tokens=fmt_tokens, iolist=iolist))

  def p_read(O, ssl, start):
    code = ssl.code
    if (code.startswith("*", start+4)):
      if (code.startswith(",", start+5)):
        tz = tokenization.ssl_iterator(ssl=ssl, start=start+6)
        iolist = []
        tz.collect_comma_separated_expressions(
          callback=iolist.append,
          enable_implied_do=1)
        iolist = tokenization.remove_redundant_parentheses(tokens=iolist)
        if (len(iolist) == 0):
          ssl.raise_syntax_error(i=start+5)
      elif (len(code) == start+5):
        iolist = None
      else:
        ssl.raise_syntax_error(i=start+5)
      O.executable.append(ei_read(ssl=ssl, start=start,
        cilist=None, fmt_tokens=None, iolist=iolist))
    elif (code.startswith("(", start+4)):
      O.p_read_write(ssl=ssl, start=start+4, ei_type=ei_read)
    else:
      ssl.raise_syntax_error(i=start+4)
    O.uses_io = True
    O.uses_read = True

  def p_return(O, ssl, start):
    O.executable.append(ei_return(ssl=ssl, start=start,
      return_label=ssl[start:]))

  def p_rewind(O, ssl, start):
    O.p_file_positioning(ssl=ssl, start=start, io_function="rewind")

  def p_save(O, ssl, start):
    assert start == 0
    if (O.save is not None):
      if (tokenization.ssl_iterator(
            ssl=ssl, start=4).collect_comma_separated_identifiers(
              callback=O.save.append, enable_common=True) == 0):
        O.save = None

  def p_stop(O, ssl, start):
    tz = tokenization.ssl_iterator(ssl=ssl, start=start+4)
    tok = tz.get(optional=True)
    if (tok is not None):
      if (not (tok.is_integer() or tok.is_string())):
        tok.raise_syntax_error()
      next_tok = tz.get(optional=True)
      if (next_tok is not None):
        next_tok.raise_syntax_error()
    O.executable.append(ei_stop(ssl=ssl, start=start, arg_token=tok))

  def p_write(O, ssl, start):
    O.p_read_write(ssl=ssl, start=start+5, ei_type=ei_write)
    O.uses_io = True
    O.uses_write = True

def collect_keyword_arguments(O, tz, n_implied):
  for known in O.__slots__:
    setattr(O, known, None)
  tok = tz.get()
  if (not tok.is_op()):
    tok.raise_syntax_error()
  if (tok.value == "*"):
    tok = tz.get()
    if (not tok.is_op_with(value=",")):
      tok.raise_syntax_error()
    O.fmt = [tok]
  else:
    if (tok.value != "("):
      tok.raise_syntax_error()
    while True: # loop over comma-separated arguments
      tok = tz.get()
      if (tok.is_op_with(value=")")):
        break
      value_tokens = []
      next_tok = tz.look_ahead()
      if (next_tok.is_op_with(value="=")):
        tz.get()
        if (not tok.is_identifier()):
          tok.raise_syntax_error()
        key = tok.value
        if (key not in O.__slots__):
          tok.raise_syntax_error()
      else:
        value_tokens.append(tok)
        for key in O.__slots__[:n_implied]:
          if (getattr(O, key) is None):
            break
        else:
          tok.raise_syntax_error()
      while True:
        tok = tz.look_ahead()
        if (tok.is_op_with(value=")")):
          break
        tz.get()
        if (tok.is_op_with(value=",")):
          break
        if (tok.is_op_with(value="(")):
          nested_tokens = []
          tz.collect_comma_separated_expressions(
            callback=nested_tokens.append,
            opening_token=tok)
          value_tokens.append(tokenization.tk_parentheses(
            ssl=tok.ssl, i_code=tok.i_code, value=nested_tokens))
        else:
          value_tokens.append(tok)
      setattr(O, key, value_tokens)

class collect_io_cilist(object):
  "Control Information List f77_std 12.8"

  __slots__ = ["unit", "fmt", "rec", "iostat", "err", "end"]

  def __init__(O, tz, fmt=None):
    if (fmt is None):
      collect_keyword_arguments(O=O, tz=tz, n_implied=2)
    else:
      for known in O.__slots__:
        setattr(O, known, None)
      O.fmt = fmt

class collect_io_olist(object):
  "Open List f77_std 12.10.1"

  chain = ["access", "form", "recl", "blank", "status", "iostat"]

  __slots__ = ["unit", "file", "err"] + chain

  def __init__(O, tz):
    collect_keyword_arguments(O=O, tz=tz, n_implied=1)

class collect_io_cllist(object):
  "Close List f77_std 12.10.2"

  chain = ["iostat", "status"]

  __slots__ = ["unit", "err"] + chain

  def __init__(O, tz):
    collect_keyword_arguments(O=O, tz=tz, n_implied=1)

class collect_io_iuflist(object):
  "iulist or iflist f77_std 12.10.3"

  chain = [
    "iostat", "exist", "opened", "number", "named", "name", "access",
    "sequential", "direct", "form", "formatted", "unformatted", "recl",
    "nextrec", "blank"]

  __slots__ = ["unit", "file", "err"] + chain

  def __init__(O, tz):
    collect_keyword_arguments(O=O, tz=tz, n_implied=1)

class collect_io_alist(object):
  "f77_std 12.10.4"

  chain = ["iostat"]

  __slots__ = ["unit", "err"] + chain

  def __init__(O, tz, unit):
    if (tz is not None):
      assert unit is None
      collect_keyword_arguments(O=O, tz=tz, n_implied=1)
    else:
      O.unit = unit
      O.iostat = None
      O.err = None

def collect_iolist(tz):
  result = []
  tok = tz.look_ahead(optional=True)
  if (tok is not None):
    if (tok.is_op_with(value=",")):
      tz.get()
    tz.collect_comma_separated_expressions(
      callback=result.append,
      enable_implied_do=1)
    result = tokenization.remove_redundant_parentheses(tokens=result)
  return result

class fproc(fproc_p_methods):

  __slots__ = [
    "leading_comments",
    "trailing_comments",
    "top_ssl", "fproc_type", "data_type", "size_tokens",
    "body_lines", "end_ssl",
    "name_plain", "name", "args",
    "body_lines_processed_already",
    "common",
    "data",
    "declarations",
    "dimension",
    "equivalence",
    "executable",
    "external",
    "format",
    "implicit",
    "intrinsic",
    "parameter",
    "save",
    "fdecl_by_identifier",
    "args_fdecl",
    "uses_common",
    "uses_save",
    "uses_io",
    "uses_read",
    "uses_write",
    "_fmt_counts_by_statement_label",
    "_common_name_by_identifier",
    "_equivalence_info",
    "_classified_equivalence_info",
    "_target_statement_labels",
    "dynamic_parameters",
    "needs_cmn",
    "needs_sve_dynamic_parameters",
    "ignore_common_and_save",
    "variant_common_names",
    "needs_is_called_first_time",
    "needs_variant_bind",
    "data_init_after_variant_bind",
    "is_passed_as_external",
    "externals_passed_by_arg_identifier"]

  def __init__(O,
        leading_comments,
        top_ssl,
        fproc_type,
        i_code,
        data_type,
        size_tokens,
        body_lines,
        end_ssl):
    assert fproc_type in ["program", "function", "subroutine", "blockdata"]
    O.leading_comments = leading_comments
    O.trailing_comments = []
    O.top_ssl = top_ssl
    O.fproc_type = fproc_type
    O.body_lines = body_lines
    O.end_ssl = end_ssl
    O.data_type = data_type
    O.size_tokens = size_tokens
    O.set_name_and_args(i_code=i_code)
    O.body_lines_processed_already = False
    O.common = utils.keyed_lists()
    O.data = []
    O.declarations = []
    O.dimension = []
    O.equivalence = []
    O.executable = []
    O.external = []
    O.format = {}
    O.init_implicit()
    O.intrinsic = []
    O.parameter = []
    O.save = []
    O.fdecl_by_identifier = None
    O.args_fdecl = None
    O.uses_common = None
    O.uses_save = None
    O.uses_io = False
    O.uses_read = False
    O.uses_write = False
    O._fmt_counts_by_statement_label = None
    O._common_name_by_identifier = None
    O._equivalence_info = None
    O._classified_equivalence_info = None
    O._target_statement_labels = None
    O.dynamic_parameters = set()
    O.needs_cmn = None
    O.needs_sve_dynamic_parameters = False
    O.ignore_common_and_save = False
    O.variant_common_names = None
    O.needs_is_called_first_time = None
    O.needs_variant_bind = None
    O.data_init_after_variant_bind = None
    O.is_passed_as_external = False
    O.externals_passed_by_arg_identifier = {}

  def is_program(O): return (O.fproc_type == "program")
  def is_function(O): return (O.fproc_type == "function")
  def is_subroutine(O): return (O.fproc_type == "subroutine")
  def is_blockdata(O): return (O.fproc_type == "blockdata")

  def first_body_source_line(O):
    assert len(O.body_lines) != 0
    assert len(O.body_lines[0].source_line_cluster) != 0
    return O.body_lines[0].source_line_cluster[0]

  def set_name_and_args(O, i_code):
    O.name_plain = None
    O.name = None
    O.args = []
    if (O.top_ssl is None):
      assert O.is_program()
      assert i_code == 0
      O.name = tokenization.tk_identifier(
        ssl=None, i_code=None, value=O.fproc_type+"_unnamed")
      return
    j_code = i_code + len(O.fproc_type)
    tz = tokenization.ssl_iterator(ssl=O.top_ssl, start=j_code)
    O.name = tz.get(optional=True)
    if (O.name is None):
      if (not O.is_program() and not O.is_blockdata()):
        O.top_ssl.raise_syntax_error(i=j_code-1)
      O.name = tokenization.tk_identifier(
        ssl=O.top_ssl, i_code=0, value=O.fproc_type+"_unnamed")
      return
    opening_token = tz.get(optional=True)
    if (opening_token is None):
      if (O.is_program() or O.is_blockdata()):
        O.name_plain = O.name
        O.name = tokenization.tk_identifier(
          ssl=O.name.ssl,
          i_code=O.name.i_code,
          value=O.fproc_type+"_"+O.name.value)
      return
    if (not opening_token.is_op_with(value="(") or O.is_blockdata()):
      opening_token.raise_syntax_error()
    need_arg = False
    while True:
      tok = tz.get_inside_parentheses(opening_token)
      if (tok.is_identifier()):
        O.args.append(tok)
      elif (tok.is_op_with(value="*")):
        if (O.fproc_type != "subroutine"):
          tok.raise_syntax_error()
        O.args.append(tok)
      elif (need_arg):
        tok.raise_syntax_error()
      elif (tok.is_op_with(value=")")):
        break
      else:
        tok.raise_syntax_error()
      tok = tz.get_inside_parentheses(opening_token)
      if (tok.is_op_with(value=")")):
        break
      if (not tok.is_op_with(value=",")):
        tok.raise_syntax_error()
      need_arg = True
    tok = tz.get(optional=True)
    if (tok is not None):
      tok.raise_syntax_error()

  def all_ssl(O):
    result = list(O.leading_comments)
    if (O.top_ssl is not None):
      result.append(O.top_ssl)
    result.extend(O.body_lines)
    result.append(O.end_ssl)
    result.extend(O.trailing_comments)
    return result

  def init_implicit(O):
    O.implicit = {}
    data_type = tokenization.tk_identifier(
      ssl=None, i_code=None, value="real")
    for c in "abcdefghopqrstuvwxyz":
      O.implicit[c] = data_type
    data_type = tokenization.tk_identifier(
      ssl=None, i_code=None, value="integer")
    for c in "ijklmn":
      O.implicit[c] = data_type

  def process_body_line(O, ssl, start):
    code = ssl.code
    if (len(code) == start): return
    i_lid = identifier_scan(code, start=start) # i_leading_identifier
    if (i_lid < 0):
      ssl.raise_syntax_error()
    if (i_lid == len(code)):
      if (code.endswith("continue", start)):
        O.p_continue(ssl=ssl, start=start)
        return
      for s in [
            "assign",
            "backspace",
            "call",
            "cycle",
            "endfile",
            "exit",
            "goto",
            "print",
            "return",
            "rewind",
            "stop"]:
        if (code.startswith(s, start)):
          p = getattr(fproc_p_methods, "p_"+s)
          p(O, ssl=ssl, start=start)
          return
      if (start != 0):
        ssl.raise_syntax_error(i=start)
      if (code in ["else", "enddo", "endif"]):
        p = getattr(fproc_p_methods, "p_"+code)
        p(O, ssl=ssl, start=start)
        return
      for s in [
            "common",
            "external",
            "entry",
            "implicit",
            "intrinsic",
            "save"]:
        if (code.startswith(s)):
          p = getattr(fproc_p_methods, "p_"+s)
          p(O, ssl=ssl, start=start)
          return
      O.process_declaration(ssl=ssl, start=start, enable_size=False)
      return
    c = code[i_lid]
    if (c == "="):
      i = ssl.comma_scan(start=i_lid+1)
      if (i < 0):
        O.process_assignment(ssl=ssl, start=start, i_equal_sign=i_lid)
        return
      if (start != 0):
        ssl.raise_syntax_error()
      if (code.startswith("do")):
        O.p_do(ssl=ssl, start=start)
        return
      ssl.raise_syntax_error()
    if (c == "("):
      i_clp = ssl.index_of_closing_parenthesis(start=i_lid+1)
      if (i_clp+1 == len(code)):
        cid = code[start:i_lid]
        if (cid in [
              "allocate",
              "backspace",
              "close",
              "deallocate",
              "endfile",
              "inquire",
              "open",
              "read",
              "rewind",
              "write"]):
          p = getattr(fproc_p_methods, "p_"+cid)
          p(O, ssl=ssl, start=start)
          return
        for s in ["call", "goto"]:
          if (code.startswith(s, start)):
            p = getattr(fproc_p_methods, "p_"+s)
            p(O, ssl=ssl, start=start)
            return
        if (start != 0):
          ssl.raise_syntax_error(i=start)
        if (cid.startswith("do") and cid.endswith("while")):
          label_end = unsigned_integer_scan(code=cid, start=2)
          if (label_end == len(cid)-5):
            O.p_dowhile(ssl=ssl, start=start, label_end=label_end)
            return
        for s in [
              "common",
              "dimension",
              "dowhile",
              "entry",
              "equivalence",
              "format",
              "implicit",
              "parameter"]:
          if (code.startswith(s)):
            p = getattr(fproc_p_methods, "p_"+s)
            p(O, ssl=ssl, start=start)
            return
        O.process_declaration(ssl=ssl, start=start, enable_size=True)
        return
      c = code[i_clp+1]
      if (c == "="):
        O.process_assignment(ssl=ssl, start=start, i_equal_sign=i_clp+1)
        return
      if (c == "("):
        i_clp2 = ssl.index_of_closing_parenthesis(start=i_clp+2)
        if (i_clp2+1 < len(code) and code[i_clp2+1] == "="):
          O.process_assignment(ssl=ssl, start=start, i_equal_sign=i_clp2+1)
          return
        for s in ["allocate(", "backspace(", "deallocate(", "read(", "write("]:
          if (code.startswith(s, start)):
            p = getattr(fproc_p_methods, "p_"+s[:-1])
            p(O, ssl=ssl, start=start)
            return
        if (code.startswith("data", start)):
          O.p_data(ssl=ssl, start=start)
          return
        ssl.raise_syntax_error_or_not_implemented()
      if (c == ","):
        cid = code[start:i_lid]
        if (cid == "goto"):
          O.p_goto(ssl=ssl, start=start)
          return
        if (cid in [
              "allocate",
              "backspace",
              "equivalence",
              "deallocate",
              "read",
              "write"]):
          p = getattr(fproc_p_methods, "p_"+cid)
          p(O, ssl=ssl, start=start)
          return
        if (start != 0):
          ssl.raise_syntax_error(i=start)
        for s in ["common", "data", "dimension", "print"]:
          if (code.startswith(s)):
            p = getattr(fproc_p_methods, "p_"+s)
            p(O, ssl=ssl, start=start)
            return
        O.process_declaration(ssl=ssl, start=start, enable_size=True)
        return
      for s in [
            "allocate(",
            "backspace(",
            "deallocate(",
            "goto(",
            "if(",
            "read(",
            "write("]:
        if (code.startswith(s, start)):
          p = getattr(fproc_p_methods, "p_"+s[:-1])
          p(O, ssl=ssl, start=start)
          return
      if (start != 0):
        ssl.raise_syntax_error(i=start)
      if (code.startswith("elseif(")):
        O.p_elseif(ssl=ssl, start=start)
        return
      if (code.startswith("data")):
        O.p_data(ssl=ssl, start=start)
        return
      O.process_declaration(ssl=ssl, start=start, enable_size=True)
      return
    if (c == "/"):
      if (start != 0):
        ssl.raise_syntax_error(i=start)
      for s in ["common", "data", "save"]:
        if (code.startswith(s)):
          p = getattr(fproc_p_methods, "p_"+s)
          p(O, ssl=ssl, start=start)
          return
      ssl.raise_syntax_error_or_not_implemented()
    if (c == ","):
      if (code.startswith("goto", start)):
        O.p_goto(ssl=ssl, start=start)
        return
      if (code.startswith("print", start)):
        O.p_print(ssl=ssl, start=start)
        return
      if (start != 0):
        ssl.raise_syntax_error(i=start)
      for s in ["common", "data", "external", "intrinsic", "save"]:
        if (code.startswith(s)):
          p = getattr(fproc_p_methods, "p_"+s)
          p(O, ssl=ssl, start=start)
          return
      if (    code.startswith("do")
          and unsigned_integer_scan(code=code, start=2) == i_lid):
        O.p_do(ssl=ssl, start=start)
        return
      O.process_declaration(ssl=ssl, start=start, enable_size=False)
      return
    if (code.endswith("stop'", start)):
      O.p_stop(ssl=ssl, start=start)
      return
    for s in ["backspace", "print", "read", "rewind"]:
      if (code.startswith(s, start)):
        p = getattr(fproc_p_methods, "p_"+s)
        p(O, ssl=ssl, start=start)
        return
    if (start != 0):
      ssl.raise_syntax_error(i=start)
    O.process_declaration(ssl=ssl, start=start, enable_size=True)

  def process_declaration(O, ssl, start, enable_size):
    assert start == 0
    if (enable_size):
      i_code, data_type, size_tokens = extract_data_type_and_size(ssl=ssl)
    else:
      i_code, data_type = extract_data_type(ssl=ssl)
      size_tokens = None
    code = ssl.code
    if (i_code == len(code)):
      ssl.raise_syntax_error(i=start)
    i_code, f90_decl = extract_f90_decl(ssl=ssl, start=i_code)
    extract_fdecl(
      result=O.declarations,
      ssl=ssl,
      start=i_code,
      data_type=data_type,
      size_tokens=size_tokens,
      allow_size=True,
      f90_decl=f90_decl)

  def process_assignment(O, ssl, start, i_equal_sign):
    if (i_equal_sign+1 == len(ssl.code)):
      ssl.raise_syntax_error()
    lhs_tokens = tokenize_expression(ssl=ssl, start=start, stop=i_equal_sign)
    rhs_tokens = tokenize_expression(ssl=ssl, start=i_equal_sign+1)
    O.executable.append(ei_assignment(ssl=ssl, start=start,
      lhs_tokens=lhs_tokens, rhs_tokens=rhs_tokens))

  def process_body_lines(O):
    assert not O.body_lines_processed_already
    O.body_lines_processed_already = True
    for ssl in O.body_lines:
      O.process_body_line(ssl=ssl, start=0)

  def show_fdecl(O):
    "for debugging; not exercised"
    assert O.fdecl_by_identifier is not None
    print O.name.value
    for key in sorted(O.fdecl_by_identifier.keys()):
      fdecl = O.fdecl_by_identifier[key]
      print " ", fdecl.id_tok.value
      print "   ", fdecl.var_type
      if (fdecl.var_storage is not None):
        print "   ", fdecl.var_storage
      if (fdecl.data_type is not None):
        print "   ", fdecl.data_type.value
      if (fdecl.parameter_assignment_tokens is not None):
        print "    parameter"
    print

  def build_fdecl_by_identifier(O):
    if (not O.body_lines_processed_already):
      O.process_body_lines()
    assert O.fdecl_by_identifier is None
    O.fdecl_by_identifier = {}
    def make_fdecl(
         id_tok,
         var_type=None,
         var_storage=None,
         data_type=None,
         size_tokens=None,
         dim_tokens=None):
      O.fdecl_by_identifier[id_tok.value] = result = fdecl_info(
        id_tok=id_tok,
        var_type=var_type,
        var_storage=var_storage,
        data_type=data_type,
        size_tokens=size_tokens,
        dim_tokens=dim_tokens)
      return result
    if (O.is_function()):
      vt = vt_function
    elif (O.is_subroutine() or O.is_blockdata()):
      vt = vt_subroutine
    else:
      vt = None
    if (vt is not None):
      make_fdecl(
        id_tok=O.name,
        var_type=vt,
        var_storage=vs_fproc_name,
        data_type=O.data_type,
        size_tokens=O.size_tokens)
    def raise_confl_decl(id_tok):
      id_tok.raise_semantic_error(
        msg="Conflicting declaration: %s" % id_tok.value)
    def raise_confl_or_repeated_decl(id_tok):
      id_tok.raise_semantic_error(
        msg="Conflicting or repeated declaration: %s" % id_tok.value)
    for fdecl in O.declarations:
      id_tok = fdecl.id_tok
      tf = O.fdecl_by_identifier.get(id_tok.value)
      if (tf is None):
        make_fdecl(
          id_tok=id_tok,
          data_type=fdecl.data_type,
          size_tokens=fdecl.size_tokens,
          dim_tokens=fdecl.dim_tokens)
      elif (tf.var_storage is vs_fproc_name):
        if (tf.data_type is not None):
          raise_confl_or_repeated_decl(id_tok=id_tok)
        if (fdecl.dim_tokens is not None):
          raise_confl_or_repeated_decl(id_tok=id_tok)
        tf.data_type = fdecl.data_type
        tf.size_tokens = fdecl.size_tokens
      elif (tf.data_type is not None):
        raise_confl_or_repeated_decl(id_tok=id_tok)
    for id_tok in O.args:
      if (id_tok.value == "*"): continue
      tf = O.fdecl_by_identifier.get(id_tok.value)
      if (tf is not None):
        if (tf.var_storage is not None):
          raise_confl_or_repeated_decl(id_tok=id_tok)
        tf.var_storage = vs_argument
      else:
        make_fdecl(id_tok=id_tok, var_storage=vs_argument)
    for id_tok,assignment_tokens in O.parameter:
      tf = O.fdecl_by_identifier.get(id_tok.value)
      if (tf is None):
        tf = make_fdecl(id_tok=id_tok, var_storage=vs_parameter)
      else:
        if (tf.var_storage is not None):
          raise_confl_or_repeated_decl(id_tok=id_tok)
        if (tf.dim_tokens is not None):
          raise_confl_or_repeated_decl(id_tok=fdecl.id_tok)
        tf.var_storage = vs_parameter
      tf.parameter_assignment_tokens = assignment_tokens
    def get_implicit_data_type(id_tok, optional=False):
      result = O.implicit.get(id_tok.value[0])
      if (result is None and not optional):
        id_tok.raise_semantic_error("Unknown data type: %s" % id_tok.value)
      return result
    def set_dim_tokens(fdecl):
      tf = O.fdecl_by_identifier.get(fdecl.id_tok.value)
      if (tf is None):
        tf = make_fdecl(
          id_tok=fdecl.id_tok,
          size_tokens=fdecl.size_tokens,
          dim_tokens=fdecl.dim_tokens)
      elif (fdecl.dim_tokens is not None):
        if (tf.dim_tokens is not None):
          fdecl.id_tok.raise_semantic_error(
            msg="Conflicting or repeated dimension: %s" % fdecl.id_tok.value)
        tf.dim_tokens = fdecl.dim_tokens
      return tf
    for fdecl in O.dimension:
      set_dim_tokens(fdecl=fdecl)
    for fdecl_list in O.common.lists:
      for fdecl in fdecl_list:
        tf = set_dim_tokens(fdecl=fdecl)
        if (tf.var_storage is not None):
          raise_confl_or_repeated_decl(id_tok=fdecl.id_tok)
        tf.var_storage = vs_common
        if (tf.data_type is None):
          tf.data_type = get_implicit_data_type(id_tok=fdecl.id_tok)
    if (O.save is not None):
      for id_tok in O.save:
        tf = O.fdecl_by_identifier.get(id_tok.value)
        if (tf is None):
          make_fdecl(id_tok=id_tok, var_storage=vs_save)
        else:
          vs = tf.var_storage
          if (vs is None):
            tf.var_storage = vs_save
          elif (vs is not vs_common):
            raise_confl_or_repeated_decl(id_tok=id_tok)
    for id_tok in O.external:
      tf = O.fdecl_by_identifier.get(id_tok.value)
      if (tf is None):
        make_fdecl(id_tok=id_tok, var_type=vt_external)
      elif (tf.dim_tokens is not None):
        raise_confl_or_repeated_decl(id_tok=id_tok)
      else:
        vs = tf.var_storage
        if (vs is not None and vs is not vs_argument):
          raise_confl_or_repeated_decl(id_tok=id_tok)
        vt = tf.var_type
        if (    vt is not None
            and vt is not vs_external
            and vt is not vs_function):
          raise_confl_or_repeated_decl(id_tok=id_tok)
        if (tf.data_type is None):
          tf.var_type = vt_external
        else:
          tf.var_type = vt_function
    for id_tok in O.intrinsic:
      tf = O.fdecl_by_identifier.get(id_tok.value)
      if (tf is None):
        make_fdecl(id_tok=id_tok, var_type=vt_intrinsic)
      elif (tf.dim_tokens is not None):
        raise_confl_or_repeated_decl(id_tok=id_tok)
      else:
        vs = tf.var_storage
        if (vs is not None):
          raise_confl_or_repeated_decl(id_tok=id_tok)
    #
    for nlist,clist in O.data:
      id_toks = tokenization.extract_identifiers(tokens=nlist)
      for id_tok in id_toks:
        tf = O.fdecl_by_identifier.get(id_tok.value)
        if (tf is None):
          make_fdecl(id_tok=id_tok, var_type=vt_scalar)
      def callback(tok):
        tf = O.fdecl_by_identifier.get(tok.value)
        if (tf is not None):
          if (tf.var_type is None):
            tf.var_type = vt_scalar
          if (tf.var_storage is None or tf.var_storage is vs_local):
            tf.var_storage = vs_save
          tf.is_modified = True
          tf.use_count += 1
      tokenization.search_for_data_or_read_target_tokens(
        callback=callback, tokens=nlist)
    #
    for id_tok in O.args:
      tf = O.fdecl_by_identifier.get(id_tok.value)
      if (tf is not None and tf.dim_tokens is not None):
        dim_id_toks = tokenization.extract_identifiers(tokens=tf.dim_tokens)
        for dim_id_tok in dim_id_toks:
          dim_tf = O.fdecl_by_identifier.get(dim_id_tok.value)
          if (dim_tf is not None and dim_tf.var_type is None):
            dim_tf.var_type = vt_used
            dim_tf.use_count += 1
    #
    for equiv_tok in O.equivalence:
      for tok_seq in equiv_tok.value:
        id_tok = tok_seq.value[0]
        tf = O.fdecl_by_identifier.get(id_tok.value)
        if (tf is None):
          make_fdecl(
            id_tok=id_tok,
            data_type=get_implicit_data_type(id_tok=id_tok))
    #
    for ei in O.executable:
      if (ei.key == "call"):
        id_tok = ei.subroutine_name
        tf = O.fdecl_by_identifier.get(id_tok.value)
        if (tf is None):
          make_fdecl(id_tok=id_tok, var_type=vt_subroutine)
        else:
          vt = tf.var_type
          if (vt is None or vt is vt_external):
            if (tf.data_type is not None):
              raise_confl_decl(id_tok=id_tok)
            if (tf.dim_tokens is not None):
              raise_confl_decl(id_tok=id_tok)
            tf.var_type = vt_subroutine
            tf.use_count += 1
          elif (vt is not vt_subroutine):
            raise_confl_decl(id_tok=id_tok)
      def search_for_id_tokens_callback(id_tok, next_tok):
        followed_by_parenthesis = (
          next_tok is not None and next_tok.is_parentheses())
        tf = O.fdecl_by_identifier.get(id_tok.value)
        if (tf is None):
          if (not followed_by_parenthesis):
            tf = make_fdecl(id_tok=id_tok, var_type=vt_used)
          elif (id_tok.value in intrinsics.set_lower):
            tf = make_fdecl(id_tok=id_tok, var_type=vt_intrinsic)
          else:
            tf = make_fdecl(id_tok=id_tok, var_type=vt_external)
          tf.use_count += 1
          return
        tf.use_count += 1
        if (tf.var_type is vt_intrinsic):
          if (id_tok.value not in intrinsics.set_lower):
            id_tok.raise_semantic_error(
              msg="Unknown intrinsic: %s" % id_tok.value)
          if (not followed_by_parenthesis):
            id_tok.raise_semantic_error(
              msg="Improper use of intrinsic: %s" % id_tok.value)
        elif (followed_by_parenthesis):
          vt = tf.var_type
          vs = tf.var_storage
          if (tf.dim_tokens is not None):
            if (tf.var_type is None):
              tf.var_type = vt_used
            if (tf.data_type is None):
              tf.data_type = get_implicit_data_type(id_tok=id_tok)
          elif (    tf.data_type is not None
                and tf.data_type.value == "character"):
            if (tf.var_type is None):
              tf.var_type = vt_used
          elif (vs is vs_argument):
            if (vt is None):
              tf.var_type = vt_external
            elif (    vt is not vt_external
                  and vt is not vt_function
                  and vt is not vt_subroutine):
              raise_confl_decl(id_tok=id_tok)
          elif (   vt is vt_external
                or vt is vt_function
                or vt is vt_subroutine):
            pass
          elif (vt is vt_intrinsic):
            if (id_tok.value not in intrinsics.set_lower):
              id_tok.raise_semantic_error(
                msg="Unknown intrinsic: %s" % id_tok.value)
          elif (vt is vt_used):
            raise_confl_decl(id_tok=id_tok)
          else:
            if (tf.var_storage is not None):
              pass # XXX should be error; ignored due to
                   #     lack of proper handling of f90 declarations
            tf.var_storage = None
            if (id_tok.value in intrinsics.set_lower):
              tf.var_type = vt_intrinsic
            else:
              tf.var_type = vt_function
        else:
          if (tf.var_type is None):
            tf.var_type = vt_used
          if (tf.data_type is None
                and tf.var_type is not vt_external
                and tf.var_type is not vt_subroutine):
            tf.data_type = get_implicit_data_type(
              id_tok=tf.id_tok, optional=True)
      ei.search_for_id_tokens(callback=search_for_id_tokens_callback)
      ei.set_is_modified(fdecl_by_identifier=O.fdecl_by_identifier)
    #
    O.uses_common = False
    for tf in O.fdecl_by_identifier.values():
      vt = tf.var_type
      vs = tf.var_storage
      if (   vt is vt_external
          or vt is vt_function
          or vt is vt_subroutine):
        if (not (vs is None or vs is vs_fproc_name or vs is vs_argument)):
          tf.id_tok.raise_internal_error()
      elif (vt is vt_intrinsic):
        if (vs is not None):
          tf.id_tok.raise_internal_error()
      else:
        if (vs is None):
          if (O.save is None):
            tf.var_storage = vs_save
          else:
            tf.var_storage = vs_local
        if (vt is vt_used):
          tf.var_type = vt_scalar
          if (tf.data_type is None):
            tf.data_type = get_implicit_data_type(id_tok=tf.id_tok)
        elif (vt is None and vs is vs_argument and tf.data_type is None):
          tf.data_type = get_implicit_data_type(id_tok=tf.id_tok)
      if (tf.is_common()):
        O.uses_common = True
    #
    for ei in O.executable:
      def search_for_id_tokens_callback(id_tok, next_tok):
        if (   next_tok is None
            or not next_tok.is_parentheses()):
          return
        tf = O.fdecl_by_identifier.get(id_tok.value)
        assert tf is not None
        if (not tf.is_user_defined_callable()):
          return
        called_identifier = id_tok.value
        for i_arg,tok_seq in enumerate(next_tok.value):
          assert tok_seq.is_seq()
          if (len(tok_seq.value) == 0):
            continue
          first_arg_tok = tok_seq.value[0]
          if (not first_arg_tok.is_identifier()):
            continue
          tf_arg = O.fdecl_by_identifier.get(first_arg_tok.value)
          assert tf_arg is not None
          if (tf_arg.is_fproc_name()):
            return
          tf_arg.passed_as_arg.setdefault(
            called_identifier, set()).add(i_arg)
          if (len(tok_seq.value) == 1):
            tf_arg.passed_as_arg_plain.setdefault(
              called_identifier, set()).add(i_arg)
      ei.search_for_id_tokens(callback=search_for_id_tokens_callback)
    #
    assert O.args_fdecl is None
    O.args_fdecl = []
    for id_tok in O.args:
      if (id_tok.value == "*"): continue
      tf = O.fdecl_by_identifier.get(id_tok.value)
      assert tf is not None
      O.args_fdecl.append(tf)
    #
    equiv_info = O.equivalence_info()
    for equiv_tok_cluster in equiv_info.equiv_tok_clusters:
      cluster_is_modified = False
      tf_cluster = []
      for equiv_tok in equiv_tok_cluster:
        for tok_seq in equiv_tok.value:
          id_tok = tok_seq.value[0]
          tf = O.fdecl_by_identifier[id_tok.value]
          tf_cluster.append(tf)
          if (tf.is_modified):
            cluster_is_modified = tf.is_modified
      if (cluster_is_modified):
        for tf in tf_cluster:
          tf.is_modified = True

  def get_fdecl(O, id_tok):
    return O.fdecl_by_identifier[id_tok.value]

  def fmt_counts_by_statement_label(O):
    assert O.body_lines_processed_already
    result = O._fmt_counts_by_statement_label
    if (result is None):
      from libtbx import dict_with_default_0
      result = dict_with_default_0()
      for ei in O.executable:
        if (ei.key in ["read", "write", "print"] and ei.fmt_tokens is None):
          tl = ei.cilist.fmt
          if (tl is not None and len(tl) == 1):
            tok = tl[0]
            if (tok.is_integer()):
              result[tok.value] += 1
    return result

  def common_name_by_identifier(O):
    result = O._common_name_by_identifier
    if (result is None):
      result = {}
      for common_name,fdecl_list in O.common.items():
        for fdecl in fdecl_list:
          identifier = fdecl.id_tok.value
          if (identifier in result):
            fdecl.id_tok.raise_semantic_error(
              msg="Identifier appears in multiple COMMON statements: %s" %
                identifier)
          result[identifier] = common_name
      O._common_name_by_identifier = result
    return result

  def equivalence_info(O):
    result = O._equivalence_info
    if (result is None):
      cu = equivalence.cluster_unions()
      for equiv_tok in O.equivalence:
        cu.add(
          key_cluster=[tok_seq.value[0].value for tok_seq in equiv_tok.value])
      cu.tidy()
      result = equivalence_info()
      for i in xrange(len(cu.unions)):
        result.equiv_tok_clusters.append([])
      for equiv_tok in O.equivalence:
        result.equiv_tok_clusters[
          cu.indices[equiv_tok.value[0].value[0].value]].append(
            equiv_tok)
      result.set_derived()
      O._equivalence_info = result
    return result

  def classified_equivalence_info(O):
    assert O.fdecl_by_identifier is not None
    result = O._classified_equivalence_info
    if (result is None):
      result = classified_equivalence_info()
      equiv_info = O.equivalence_info()
      for equiv_tok_cluster in equiv_info.equiv_tok_clusters:
        highest_priority = 0
        for equiv_tok in equiv_tok_cluster:
          for tok_seq in equiv_tok.value:
            identifier = tok_seq.value[0].value
            fdecl = O.fdecl_by_identifier[identifier]
            if (fdecl.is_common()):
              priority = 3
            elif (fdecl.is_save()):
              priority = 2
            elif (fdecl.is_local()):
              priority = 1
            else:
              tok_seq.raise_semantic_error(msg="Invalid EQUIVALENCE")
            highest_priority = max(highest_priority, priority)
        assert highest_priority != 0
        slot = getattr(result, ["local", "save", "common"][highest_priority-1])
        slot.equiv_tok_clusters.append(
          equiv_info.equiv_tok_cluster_by_identifier[identifier])
      result.set_derived()
      O._classified_equivalence_info = result
    return result

  def set_uses_save(O):
    cei = O.classified_equivalence_info()
    O.uses_save = (len(O.data) != 0)
    if (not O.uses_save):
      for fdecl in O.fdecl_by_identifier.values():
        if (fdecl.is_save()):
         equiv_tok_cluster = cei.common.equiv_tok_cluster_by_identifier.get(
           fdecl.id_tok.value)
         if (equiv_tok_cluster is None):
           O.uses_save = True
           break

  def target_statement_labels(O):
    result = O._target_statement_labels
    if (O._target_statement_labels is None):
      result = {}
      for ei in O.executable:
        if (ei.key == "goto"):
          result.setdefault(ei.label.value, []).append(ei.label)
        elif (ei.key in ["goto_computed", "if_arithmetic"]):
          for tok in ei.labels:
            result.setdefault(tok.value, []).append(tok)
        elif (ei.key == "open"):
          if (ei.olist.err is not None):
            tok = tokenization.get_statement_label_token(tokens=ei.olist.err)
            result.setdefault(tok.value, []).append(tok)
        elif (ei.key == "close"):
          if (ei.cllist.err is not None):
            tok = tokenization.get_statement_label_token(tokens=ei.cllist.err)
            result.setdefault(tok.value, []).append(tok)
        elif (ei.key == "inquire"):
          if (ei.iuflist.err is not None):
            tok = tokenization.get_statement_label_token(tokens=ei.iuflist.err)
            result.setdefault(tok.value, []).append(tok)
        elif (ei.key == "file_positioning"):
          if (ei.alist.err is not None):
            tok = tokenization.get_statement_label_token(tokens=ei.alist.err)
            result.setdefault(tok.value, []).append(tok)
        elif (ei.key in ["read", "write"]):
          cilist = ei.cilist
          if (cilist is not None):
            for slot in ["end", "err"]:
              tokens = getattr(cilist, slot)
              if (tokens is not None):
                tok = tokenization.get_statement_label_token(tokens=tokens)
                result.setdefault(tok.value, []).append(tok)
      O._target_statement_labels = result
    return result

  def _eval_const_expression_simple_identifier(O,
        identifier, buffer, allow_power):
    if (identifier in O.dynamic_parameters):
      return False
    fdecl = O.fdecl_by_identifier.get(identifier)
    if (fdecl is None):
      return False
    tokens = fdecl.parameter_assignment_tokens
    if (tokens is None):
      return False
    code = tokenization.tokens_as_python_code(
      tokens=tokens, allow_power=allow_power)
    if (code is None):
      return False
    expr = "%s = %s" % (identifier, code)
    buffer.append(expr)
    return O._eval_const_expression_simple_tokens(
      tokens=tokens, buffer=buffer, allow_power=allow_power)

  def _eval_const_expression_simple_tokens(O,
        tokens, buffer, allow_power):
    for id_tok in tokenization.extract_identifiers(tokens=tokens):
      if (not O._eval_const_expression_simple_identifier(
        identifier=id_tok.value, buffer=buffer, allow_power=allow_power)):
          return False
    return True

  def eval_const_expression_simple(O,
        identifier=None, tokens=None, allow_power=True):
    assert O.fdecl_by_identifier is not None
    assert "_" not in O.fdecl_by_identifier # not supported
    assert [identifier, tokens].count(None) == 1
    buffer = []
    if (identifier is None):
      code = tokenization.tokens_as_python_code(
        tokens=tokens, allow_power=allow_power)
      if (code is None):
        return None
      buffer.append("_ = %s" % code)
      if (not O._eval_const_expression_simple_tokens(
            tokens=tokens, buffer=buffer, allow_power=allow_power)):
        return None
    else:
      if (not O._eval_const_expression_simple_identifier(
            identifier=identifier, buffer=buffer, allow_power=allow_power)):
        return None
    buffer.reverse()
    code = "\n".join(buffer)
    exec_globals = {}
    exec_locals = {}
    exec(code, exec_globals, exec_locals)
    if (identifier is None):
      return exec_locals["_"]
    return exec_locals[identifier]

  def eval_dimensions_simple(O, dim_tokens, allow_power=True):
    vals = []
    for tok_seq in dim_tokens:
      if (tokenization.tok_seq_is_star(tok_seq=tok_seq)):
        vals.append(None)
      else:
        for i,tok in enumerate(tok_seq.value):
          if (tok.is_op_with(value=":")):
            fl = []
            for tokens in (tok_seq.value[:i], tok_seq.value[i+1:]):
              fl.append(O.eval_const_expression_simple(
                tokens=tokens, allow_power=allow_power))
            f,l = fl
            if (f is None or l is None):
              vals.append(None)
            else:
              vals.append(l-f+1)
            break
        else:
          vals.append(O.eval_const_expression_simple(
            tokens=tok_seq.value, allow_power=allow_power))
    return vals

class equivalence_info(object):

  __slots__ = [
    "equiv_tok_clusters",
    "equiv_tok_cluster_by_identifier",
    "identifier_clusters",
    "identifier_cluster_by_identifier"]

  def __init__(O):
    O.equiv_tok_clusters = []

  def set_derived(O):
    O.equiv_tok_cluster_by_identifier = {}
    O.identifier_clusters = []
    O.identifier_cluster_by_identifier = {}
    for equiv_tok_cluster in O.equiv_tok_clusters:
      identifier_cluster = []
      O.identifier_clusters.append(identifier_cluster)
      for equiv_tok in equiv_tok_cluster:
        for tok_seq in equiv_tok.value:
          identifier = tok_seq.value[0].value
          O.equiv_tok_cluster_by_identifier[identifier] = equiv_tok_cluster
          identifier_cluster.append(identifier)
          O.identifier_cluster_by_identifier[identifier] = identifier_cluster

class classified_equivalence_info(object):

  __slots__ = ["common", "save", "local"]

  def __init__(O):
    for slot in O.__slots__:
      setattr(O, slot, equivalence_info())

  def set_derived(O):
    for slot in O.__slots__:
      getattr(O, slot).set_derived()

  def has_save(O):
    return (len(O.save.equiv_tok_clusters) != 0)

class split_fprocs(object):

  __slots__ = [
    "program",
    "subroutine",
    "function",
    "blockdata",
    "all_in_input_order",
    "_fprocs_by_name",
    "_fprocs_by_name_plain"]

  def __init__(O):
    O.program = []
    O.subroutine = []
    O.function = []
    O.blockdata = []
    O.all_in_input_order = []
    O._fprocs_by_name = None
    O._fprocs_by_name_plain = None

  def by_type(O):
    return [O.program, O.blockdata, O.subroutine, O.function]

  def process(O, stripped_source_lines):
    ssls = iter(stripped_source_lines)
    leading_comments = []
    for curr_ssl in ssls:
      if (curr_ssl.is_comment()):
        leading_comments.append(curr_ssl)
        continue
      assert len(curr_ssl.code) != 0
      def collect_until_end(
            fproc_type, top_ssl, i_code, data_type, size_tokens,
            first_body_line=None):
        body_lines = []
        if (first_body_line is not None):
          body_lines.append(first_body_line)
        specific_end = "end"+fproc_type
        for ssl in ssls:
          if (ssl.code in ["end", specific_end]):
            result = fproc(
              leading_comments=leading_comments,
              top_ssl=top_ssl,
              fproc_type=fproc_type,
              i_code=i_code,
              data_type=data_type,
              size_tokens=size_tokens,
              body_lines=body_lines,
              end_ssl=ssl)
            O.all_in_input_order.append(result)
            return result
          body_lines.append(ssl)
        if (top_ssl is None):
          top_ssl = first_body_line
        top_ssl.raise_error(msg="Missing END for %s" % (fproc_type.upper()))
      for fproc_type in ["program", "blockdata", "subroutine", "function"]:
        if (curr_ssl.code.startswith(fproc_type)):
          getattr(O, fproc_type).append(collect_until_end(
            fproc_type=fproc_type,
            top_ssl=curr_ssl,
            i_code=0,
            data_type=None,
            size_tokens=None))
          break
      else:
        i_code, data_type, size_tokens = extract_data_type_and_size(
          ssl=curr_ssl, optional=True)
        if (i_code is None or
              not curr_ssl.code.startswith("function", i_code)):
          O.program.append(collect_until_end(
            fproc_type="program",
            top_ssl=None,
            i_code=0,
            data_type=None,
            size_tokens=None,
            first_body_line=curr_ssl))
        else:
          O.function.append(collect_until_end(
            fproc_type="function",
            top_ssl=curr_ssl,
            i_code=i_code,
            data_type=data_type,
            size_tokens=size_tokens))
      leading_comments = []
    if (len(leading_comments) != 0 and len(O.all_in_input_order) != 0):
      O.all_in_input_order[-1].trailing_comments = leading_comments

  def show_counts_by_type(O, out=None, prefix=""):
    if (out is None): out = sys.stdout
    print >> out, prefix + "Counts by Fortran procedure type:"
    for attr in O.__slots__[:4]:
      print >> out, prefix + "  %s: %s" % (attr, len(getattr(O, attr)))

  def process_body_lines(O):
    for fproc in O.all_in_input_order:
      fproc.process_body_lines()

  def build_fdecl_by_identifier(O):
    for fproc in O.all_in_input_order:
      fproc.build_fdecl_by_identifier()

  def fprocs_by_name(O, plain=False):
    if (O._fprocs_by_name is None):
      O._fprocs_by_name = {}
      O._fprocs_by_name_plain = {}
      for fprocs in O.by_type():
        for fproc in fprocs:
          other = O._fprocs_by_name.get(fproc.name.value)
          if (other is not None):
            msg = ["Fortran procedure name conflict:"]
            for name in [other.name, fproc.name]:
              if (name.ssl is None):
                msg.append(
                  "  %d. definition: %s (implied)\n"
                  "    before %s" % (
                    len(msg),
                    fproc.name.value,
                    fproc.first_body_source_line()
                      .format_file_name_and_line_number()))
              else:
                msg.append(name.format_error(
                  msg="%d. definition" % len(msg), prefix="  "))
            from libtbx.utils import Sorry
            raise Sorry("\n".join(msg))
          O._fprocs_by_name[fproc.name.value] = fproc
          if (fproc.name_plain is not None):
            O._fprocs_by_name_plain[fproc.name_plain.value] = fproc
    if (plain):
      return O._fprocs_by_name_plain
    return O._fprocs_by_name

  def build_bottom_up_fproc_list_following_calls(O, top_procedures=None):
    return build_bottom_up_fproc_list_following_calls(
      all_fprocs=O, top_procedures=top_procedures)

class build_bottom_up_fproc_list_following_calls(object):

  __slots__ = [
    "all_fprocs",
    "top_procedures",
    "deps_by_fproc_identifier",
    "bottom_up_list",
    "forward_uses_by_identifier",
    "dependency_cycles",
    "intrinsics_extra",
    "missing_external_fdecls_by_identifier"]

  def __init__(O, all_fprocs, top_procedures=None):
    O.all_fprocs = all_fprocs
    O.top_procedures = top_procedures
    fprocs_by_name = O.all_fprocs.fprocs_by_name()
    #
    for fproc in O.all_fprocs.all_in_input_order:
      for fdecl in fproc.fdecl_by_identifier.values():
        if (    fdecl.is_user_defined_callable()
            and not fdecl.var_storage is vs_argument
            and len(fdecl.passed_as_arg_plain) != 0):
          def recursively_update_externals_passed(
                primary_external_identifier,
                procs_visited_already,
                caller_fdecl):
            for called_name,i_arg_set in \
                  caller_fdecl.passed_as_arg_plain.items():
              called_fproc = fprocs_by_name.get(called_name)
              if (called_fproc is None):
                continue
              for i_arg in i_arg_set:
                arg_identifier = called_fproc.args[i_arg].value
                primaries = called_fproc.externals_passed_by_arg_identifier \
                  .setdefault(arg_identifier, set())
                if (not primary_external_identifier in primaries):
                  primaries.add(primary_external_identifier)
                  if (called_name not in procs_visited_already):
                    procs_visited_already.add(called_name)
                    recursively_update_externals_passed(
                      primary_external_identifier=primary_external_identifier,
                      procs_visited_already=procs_visited_already,
                      caller_fdecl=called_fproc.fdecl_by_identifier[
                        arg_identifier])
          primary_external_identifier = fdecl.id_tok.value
          primary_fproc = fprocs_by_name.get(primary_external_identifier)
          if (primary_fproc is not None):
            primary_fproc.is_passed_as_external = True
          recursively_update_externals_passed(
            primary_external_identifier=primary_external_identifier,
            procs_visited_already=set([fproc.name.value]),
            caller_fdecl=fdecl)
    #
    O.deps_by_fproc_identifier = {}
    external_fdecls = {}
    def get_dependencies(fproc):
      deps = set()
      for primaries in fproc.externals_passed_by_arg_identifier.values():
        deps.update(primaries)
      for identifier in sorted(fproc.fdecl_by_identifier.keys()):
        if (identifier == fproc.name.value): continue
        fdecl = fproc.fdecl_by_identifier[identifier]
        if (    fdecl.is_user_defined_callable()
            and fdecl.var_storage is not vs_argument):
          deps.add(fdecl.id_tok.value)
          external_fdecls.setdefault(identifier, []).append(fdecl)
      if (fproc.is_program()):
        for b in O.all_fprocs.blockdata:
          deps.add(b.name.value)
      result = sorted(deps)
      O.deps_by_fproc_identifier[fproc.name.value] = result
      return result
    if (O.top_procedures is None or len(O.top_procedures) == 0):
      connections_for_topological_sort = []
      for fproc in O.all_fprocs.all_in_input_order:
        connections_for_topological_sort.append(
          (fproc.name.value, get_dependencies(fproc=fproc)))
    else:
      top_procedures_tidy = []
      for top_procedure_or_procedures in O.top_procedures:
        for top_procedure in top_procedure_or_procedures.split(","):
          top_fproc = fprocs_by_name.get(top_procedure)
          if (top_fproc is None):
            top_fproc = fprocs_by_name.get("program_"+top_procedure)
          if (top_fproc is None):
            raise RuntimeError(
              "Unknown Fortran procedure name: %s" % top_procedure)
          top_procedures_tidy.append(top_procedure)
          def recurse(fproc):
            for identifier in get_dependencies(fproc=fproc):
              if (identifier in O.deps_by_fproc_identifier):
                continue
              next_fproc = fprocs_by_name.get(identifier)
              if (next_fproc is not None):
                recurse(fproc=next_fproc)
          recurse(fproc=top_fproc)
      O.top_procedures = top_procedures_tidy
      connections_for_topological_sort = []
      for fproc in O.all_fprocs.all_in_input_order:
        if (fproc.name.value in O.deps_by_fproc_identifier):
          connections_for_topological_sort.append(
            (fproc.name.value, get_dependencies(fproc=fproc)))
    #
    from libtbx import topological_sort
    successors_by_node = dict(connections_for_topological_sort)
    O.bottom_up_list = []
    bottom_up_set = set()
    O.forward_uses_by_identifier = {}
    forward_uses_set = set()
    O.intrinsics_extra = set()
    O.missing_external_fdecls_by_identifier = {}
    for identifier in topological_sort.stable(
                        connections=connections_for_topological_sort):
      fproc = fprocs_by_name.get(identifier)
      if (fproc is not None):
        O.bottom_up_list.append(fproc)
        bottom_up_set.add(identifier)
        for dep in successors_by_node[identifier]:
          if (    dep in successors_by_node
              and dep not in bottom_up_set
              and dep not in forward_uses_set):
            O.forward_uses_by_identifier.setdefault(identifier, []).append(dep)
            forward_uses_set.add(dep)
      elif (identifier in intrinsics.extra_set_lower):
        O.intrinsics_extra.add(identifier)
      elif (identifier not in all_fprocs.fprocs_by_name(plain=True)):
        O.missing_external_fdecls_by_identifier[identifier] = \
          external_fdecls[identifier]
    O.dependency_cycles = topological_sort.strongly_connected_components(
      successors_by_node=successors_by_node)

  def each_fproc_update_is_modified(O):
    fprocs_by_name = O.all_fprocs.fprocs_by_name()
    for caller_fproc in O.bottom_up_list:
      for caller_fdecl in caller_fproc.fdecl_by_identifier.values():
        for called_identifier, i_args in caller_fdecl.passed_as_arg.items():
          primaries = caller_fproc.externals_passed_by_arg_identifier.get(
            called_identifier)
          if (primaries is None):
            primaries = [called_identifier]
          for called_identifier in primaries:
            called_fproc = fprocs_by_name.get(called_identifier)
            if (called_fproc is not None):
              for i_arg in sorted(i_args):
                if (i_arg >= len(called_fproc.args_fdecl)):
                  continue
                arg_fdecl = called_fproc.args_fdecl[i_arg]
                if (arg_fdecl.is_modified):
                  caller_fdecl.is_modified = True
    return O

  def each_fproc_update_needs_cmn(O):
    fprocs_by_name = O.all_fprocs.fprocs_by_name()
    have_blockdata = (len(O.all_fprocs.blockdata) != 0)
    for caller_fproc in O.bottom_up_list:
      caller_fproc.needs_cmn = \
           caller_fproc.uses_common \
        or caller_fproc.uses_save \
        or caller_fproc.uses_io \
        or (have_blockdata and caller_fproc.is_program()) \
        or len(caller_fproc.dynamic_parameters) != 0
      if (not caller_fproc.needs_cmn):
        for dependency in O.deps_by_fproc_identifier.get(
                            caller_fproc.name.value, []):
          called_fproc = fprocs_by_name.get(dependency)
          if (called_fproc is not None and called_fproc.needs_cmn):
            caller_fproc.needs_cmn = True
            break
    return O

def process(file_names, basic_only=False, skip_load_includes=False):
  assert not skip_load_includes or basic_only
  all_fprocs = split_fprocs()
  import itertools
  global_line_index_generator = itertools.count()
  for file_name in file_names:
    all_fprocs.process(stripped_source_lines=load(
      global_line_index_generator=global_line_index_generator,
      file_name=file_name,
      skip_load_includes=skip_load_includes))
  if (not basic_only):
    all_fprocs.build_fdecl_by_identifier()
    for fproc in all_fprocs.all_in_input_order:
      fproc.common_name_by_identifier()
      fproc.set_uses_save()
      fproc.target_statement_labels()
  return all_fprocs
