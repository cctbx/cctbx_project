from fable \
  import unsigned_integer_scan, \
  identifier_scan, \
  floating_point_scan_after_exponent_char, \
  floating_point_scan_after_dot

class token(object):

  __slots__ = ["typeid", "ssl", "i_code", "value"]

  def __init__(O, typeid, ssl, i_code, value):
    O.typeid = typeid
    O.ssl = ssl
    O.i_code = i_code
    O.value = value

  def type(O):
    return [
      "identifier",
      "integer",
      "hexadecimal",
      "real",
      "double_precision",
      "logical",
      "string",
      "op",
      "complex",
      "seq",
      "parentheses",
      "implied_do",
      "power",
      "format"][O.typeid]

  def is_identifier(O): return O.typeid == 0
  def is_integer(O): return O.typeid == 1
  def is_hexadecimal(O): return O.typeid == 2
  def is_real(O): return O.typeid == 3
  def is_double_precision(O): return O.typeid == 4
  def is_logical(O): return O.typeid == 5
  def is_string(O): return O.typeid == 6
  def is_op(O): return O.typeid == 7
  def is_complex(O): return O.typeid == 8
  def is_seq(O): return O.typeid == 9
  def is_parentheses(O): return O.typeid == 10
  def is_implied_do(O): return O.typeid == 11
  def is_power(O): return O.typeid == 12
  def is_format(O): return O.typeid == 13

  def is_identifier_or_scalar_number(O):
    return O.typeid <= 4

  def is_op_with(O, value):
    return O.typeid == 7 and O.value == value

  def is_unary_plus_or_minus(O):
    return O.typeid == 7 and O.value in ["+", "-"]

  def is_seq_or_parentheses(O):
    return 9 <= O.typeid <= 10

  def is_seq_or_parentheses_or_implied_do(O):
    return 9 <= O.typeid <= 11

  def stmt_location(O):
    return O.ssl.stmt_location(i=O.i_code)

  def format_error(O, msg, prefix=""):
    return O.ssl.format_error(msg=msg, i=O.i_code, prefix=prefix)

  def raise_error(O, msg, ErrorType=None):
    O.ssl.raise_error(msg=msg, i=O.i_code, ErrorType=ErrorType)

  def raise_syntax_error(O):
    O.ssl.raise_syntax_error(i=O.i_code)

  def raise_semantic_error(O, msg=None):
    O.ssl.raise_semantic_error(msg=msg, i=O.i_code)

  def raise_not_supported(O):
    O.raise_error(
      msg="Sorry: not supported",
      ErrorType=RuntimeError)

  def raise_internal_error(O):
    O.ssl.raise_internal_error(i=O.i_code)

  def raise_if_not_identifier(O):
    i = identifier_scan(code=O.value)
    if (i < 0 or i != len(O.value)):
      O.raise_error(msg="Not an identifier: %s" % repr(O.value))

  def raise_missing_closing(O):
    O.raise_error(msg='Missing a closing ")"')

  def raise_missing_opening(O):
    O.raise_error(msg='Closing ")" without a matching opening parenthesis')

  def fquote(O):
    return "'" + O.value.replace("'", "''") + "'"

def tk_identifier(ssl, i_code, value): return token(0, ssl, i_code, value)
def tk_integer(ssl, i_code, value): return token(1, ssl, i_code, value)
def tk_hexadecimal(ssl, i_code, value): return token(2, ssl, i_code, value)
def tk_real(ssl, i_code, value): return token(3, ssl, i_code, value)
def tk_double_precision(ssl, i_code, value): return token(4, ssl, i_code, value)
def tk_logical(ssl, i_code, value): return token(5, ssl, i_code, value)
def tk_string(ssl, i_code, value): return token(6, ssl, i_code, value)
def tk_op(ssl, i_code, value): return token(7, ssl, i_code, value)
def tk_complex(ssl, i_code, value): return token(8, ssl, i_code, value)
def tk_seq(ssl, i_code, value): return token(9, ssl, i_code, value)
def tk_parentheses(ssl, i_code, value): return token(10, ssl, i_code, value)
def tk_implied_do(ssl, i_code, value): return token(11, ssl, i_code, value)
def tk_power(ssl, i_code, value): return token(12, ssl, i_code, value)
def tk_format(ssl, i_code, value): return token(13, ssl, i_code, value)

class ssl_iterator(object):
  "stripped source line iterator"

  __slots__ = ["ssl", "start", "stop", "i", "buffer"]

  def __init__(O, ssl, start, stop=None):
    O.ssl = ssl
    O.start = start
    if (stop is None): O.stop = len(ssl.code)
    else:              O.stop = stop
    O.i = start
    O.buffer = []

  def __iter__(O): return O

  def next(O):
    tok = O.get(optional=True)
    if (tok is None):
      raise StopIteration
    return tok

  def get(O, optional=False):
    ssl = O.ssl
    code = ssl.code
    stop = O.stop
    def return_none():
      if (not optional):
        if (stop == 0): ssl.raise_syntax_error()
        else:           ssl.raise_syntax_error(i=stop-1)
      O.i = None
      return None
    if (len(O.buffer) != 0):
      O.i, result = O.buffer.pop(0)
      if (result is None):
        return_none()
      return result
    if (O.i is None):
      if (optional):
        return None
      ssl.raise_internal_error(i=stop-1)
    while (O.i < stop):
      i_code = O.i
      c = code[i_code]
      if ("(),=:+-".find(c) >= 0):
        O.i += 1
        return tk_op(ssl=ssl, i_code=i_code, value=c)
      if (c == "'"):
        O.i += 1
        return tk_string(
          ssl=ssl,
          i_code=i_code,
          value=ssl.strings[ssl.string_indices.index(i_code)])
      if (c == "x" and code.startswith("'", i_code+1)):
        O.i += 2
        return tk_hexadecimal(
          ssl=ssl,
          i_code=i_code,
          value=ssl.strings[ssl.string_indices.index(i_code+1)])
      j = identifier_scan(code=code, start=i_code)
      if (j > 0):
        O.i = j
        return tk_identifier(ssl=ssl, i_code=i_code, value=code[i_code:j])
      if (c == "*"):
        if (code.startswith("*", i_code+1)):
          O.i += 2
          return tk_op(ssl=ssl, i_code=i_code, value="**")
        O.i += 1
        return tk_op(ssl=ssl, i_code=i_code, value="*")
      if (c == "/"):
        if (code.startswith("/", i_code+1)):
          O.i += 2
          return tk_op(ssl=ssl, i_code=i_code, value="//")
        if (code.startswith("=", i_code+1)):
          O.i += 2
          return tk_op(ssl=ssl, i_code=i_code, value="/=")
        O.i += 1
        return tk_op(ssl=ssl, i_code=i_code, value="/")
      if (c == "."):
        return O.__after_dot(i_fld=i_code, i_dot=i_code)
      j = unsigned_integer_scan(code=code, start=i_code, stop=stop)
      if (j > 0):
        if (j == stop):
          O.i = j
          return tk_integer(ssl=ssl, i_code=i_code, value=code[i_code:j])
        cj = code[j]
        if (cj == "."):
          if (j + 1 == stop):
            O.i = stop
            return tk_real(ssl=ssl, i_code=i_code, value=code[i_code:stop])
          return O.__after_dot(i_fld=i_code, i_dot=j)
        if (cj == "e" or cj == "d"):
          k = floating_point_scan_after_exponent_char(
            code=code, start=j+1, stop=stop)
          if (k < 0):
            ssl.raise_error(
              msg="Invalid floating-point literal", i=j+1)
          O.i = k
          if (cj == "d"):
            return tk_double_precision(
              ssl=ssl, i_code=i_code, value=code[i_code:k])
          return tk_real(ssl=ssl, i_code=i_code, value=code[i_code:k])
        O.i = j
        return tk_integer(ssl=ssl, i_code=i_code, value=code[i_code:j])
      ssl.raise_syntax_error(i=i_code)
    return_none()

  def __after_dot(O, i_fld, i_dot):
    ssl = O.ssl
    stop = O.stop
    if (i_dot + 1 == stop):
      ssl.raise_error(
        msg="Expression unexpectedly ends with a dot", i=i_dot)
    code = ssl.code
    csw = code.startswith
    c = code[i_dot+1]
    if (c == "a"):
      if (not csw("nd.", i_dot+2)):
        ssl.raise_syntax_error(i=i_dot+2)
      tok = tk_op(ssl=ssl, i_code=i_dot, value=".and.")
      if (i_dot == i_fld):
        O.i = i_dot+5
        return tok
      O.buffer.append((i_dot+5, tok))
      O.i = i_dot
      return tk_integer(ssl=ssl, i_code=i_fld, value=code[i_fld:i_dot])
    if (c == "o"):
      if (not csw("r.", i_dot+2)):
        ssl.raise_syntax_error(i=i_dot+2)
      tok = tk_op(ssl=ssl, i_code=i_dot, value=".or.")
      if (i_dot == i_fld):
        O.i = i_dot+4
        return tok
      O.buffer.append((i_dot+4, tok))
      O.i = i_dot
      return tk_integer(ssl=ssl, i_code=i_fld, value=code[i_fld:i_dot])
    if (c == "e" or c == "d"):
      if (c == "e"):
        if (csw("q.", i_dot+2)):
          tok = tk_op(ssl=ssl, i_code=i_dot, value=".eq.")
          if (i_dot == i_fld):
            O.i = i_dot+4
            return tok
          O.buffer.append((i_dot+4, tok))
          O.i = i_dot
          return tk_integer(ssl=ssl, i_code=i_fld, value=code[i_fld:i_dot])
        if (csw("qv.", i_dot+2)):
          tok = tk_op(ssl=ssl, i_code=i_dot, value=".eqv.")
          if (i_dot == i_fld):
            O.i = i_dot+5
            return tok
          O.buffer.append((i_dot+5, tok))
          return tk_integer(ssl=ssl, i_code=i_fld, value=code[i_fld:i_dot])
      if (i_dot == i_fld):
        ssl.raise_syntax_error(i=i_dot+1)
      j = floating_point_scan_after_exponent_char(
        code=code, start=i_dot+2, stop=stop)
      if (j < 0):
        ssl.raise_syntax_error(i=i_dot+1)
      O.i = j
      if (c == "d"):
        return tk_double_precision(
          ssl=ssl, i_code=i_fld, value=code[i_fld:j])
      return tk_real(ssl=ssl, i_code=i_fld, value=code[i_fld:j])
    if (c == "f"):
      if (not csw("alse.", i_dot+2) or i_dot != i_fld):
        ssl.raise_syntax_error(i=i_dot+1)
      O.i = i_dot+7
      return tk_logical(ssl=ssl, i_code=i_dot, value=".false.")
    if (c == "t"):
      if (not csw("rue.", i_dot+2) or i_dot != i_fld):
        ssl.raise_syntax_error(i=i_dot+1)
      O.i = i_dot+6
      return tk_logical(ssl=ssl, i_code=i_dot, value=".true.")
    if (c == "n"):
      if (csw("ot.", i_dot+2)):
        if (i_dot != i_fld):
          ssl.raise_syntax_error(i=i_dot+1)
        O.i = i_dot+5
        return tk_op(ssl=ssl, i_code=i_dot, value=".not.")
      if (csw("e.", i_dot+2)):
        tok = tk_op(ssl=ssl, i_code=i_dot, value=".ne.")
        if (i_dot == i_fld):
          O.i = i_dot+4
          return tok
        O.buffer.append((i_dot+4, tok))
        O.i = i_dot
        return tk_integer(ssl=ssl, i_code=i_fld, value=code[i_fld:i_dot])
      if (csw("eqv.", i_dot+2)):
        tok = tk_op(ssl=ssl, i_code=i_dot, value=".neqv.")
        if (i_dot == i_fld):
          O.i = i_dot+6
          return tok
        O.buffer.append((i_dot+6, tok))
        O.i = i_dot
        return tk_integer(ssl=ssl, i_code=i_fld, value=code[i_fld:i_dot])
      ssl.raise_syntax_error(i=i_dot+1)
    if (c == "g" or c == "l"):
      if (not csw("t.", i_dot+2) and not csw("e.", i_dot+2)):
        ssl.raise_syntax_error(i=i_dot+1)
      tok = tk_op(ssl=ssl, i_code=i_dot, value=code[i_dot:i_dot+4])
      if (i_dot == i_fld):
        O.i = i_dot+4
        return tok
      O.buffer.append((i_dot+4, tok))
      O.i = i_dot
      return tk_integer(ssl=ssl, i_code=i_fld, value=code[i_fld:i_dot])
    j = floating_point_scan_after_dot(code=code, start=i_dot+1, stop=stop)
    if (j < 0):
      ssl.raise_syntax_error(i=i_dot+1)
    O.i = j
    if (code.find("d", i_dot+1, O.i) < 0): tk_type = tk_real
    else:                                  tk_type = tk_double_precision
    return tk_type(ssl=ssl, i_code=i_fld, value=code[i_fld:j])

  def look_ahead(O, optional=False):
    if (len(O.buffer) != 0):
      return O.buffer[0][1]
    prev_i = O.i
    tok = O.get(optional=optional)
    O.buffer.append((O.i, tok))
    O.i = prev_i
    return tok

  def get_complex_literal(O, opening_token):
    assert opening_token.ssl is O.ssl
    invalid_message = "Invalid complex number literal"
    result = [opening_token]
    def get_part():
      tok = O.get()
      sign_tok = None
      if (tok.is_op()):
        if (tok.value not in ["+", "-"]):
          tok.raise_error(msg=invalid_message)
        sign_tok = tok
        tok = O.get()
      if (not tok.is_identifier_or_scalar_number()):
        tok.raise_error(msg=invalid_message)
      return (sign_tok, tok)
    result.append(get_part())
    tok = O.get()
    if (not tok.is_op_with(value=",")):
      tok.raise_error(msg=invalid_message)
    result.append(get_part())
    tok = O.get()
    if (not tok.is_op_with(value=")")):
      tok.raise_error(msg=invalid_message)
    result.append(tok)
    return tk_complex(ssl=O.ssl, i_code=opening_token.i_code, value=result)

  def get_inside_parentheses(O, opening_token):
    tok = O.get(optional=True)
    if (tok is None):
      opening_token.raise_missing_closing()
    return tok

  def collect_comma_separated_identifiers(O,
        callback,
        one_required=False,
        enable_common=False):
    n_list = 0
    while True:
      tok = O.get(optional=not one_required)
      if (tok is None):
        break
      if (tok.is_identifier()):
        callback(tok)
      elif (not enable_common):
        tok.raise_syntax_error()
      else:
        if (not tok.is_op_with("/")):
          tok.raise_syntax_error()
        tok = O.get()
        if (not tok.is_identifier()):
          tok.raise_syntax_error()
        tok = O.get()
        if (not tok.is_op_with("/")):
          tok.raise_syntax_error()
        # SAVE of COMMON are ignored (considered redundant)
      n_list += 1
      tok = O.get(optional=True)
      if (tok is None):
        break
      if (not tok.is_op_with(value=",")):
        tok.raise_syntax_error()
      one_required = True
    return n_list

  def collect_comma_separated_expressions(O,
        callback,
        opening_token=None,
        first_get_optional=True,
        stop_after_given_number_of_commas=None,
        enable_implied_do=0):
    tokens = tk_seq(ssl=O.ssl, i_code=O.i, value=[])
    tapp = tokens.value.append
    last_comma = None
    n_callbacks = 0
    i_assignment_op = None
    get_optional = first_get_optional
    while True:
      tok = O.get(optional=get_optional)
      get_optional = True
      if (tok is None):
        if (opening_token is not None):
          opening_token.raise_missing_closing()
        if (len(tokens.value) == 0):
          if (last_comma is None):
            return None
          last_comma.raise_syntax_error()
        callback(tokens)
        return None
      if (tok.is_op()):
        tv = tok.value
        if (tv ==")"):
          if (opening_token is None):
            tok.raise_missing_opening()
          if (i_assignment_op is None):
            if (enable_implied_do == 2):
              tok.raise_syntax_error()
          elif (i_assignment_op == n_callbacks):
            tok.raise_syntax_error()
          if (len(tokens.value) == 0):
            if (last_comma is None):
              return i_assignment_op
            last_comma.raise_syntax_error()
          callback(tokens)
          return i_assignment_op
        if (tv == ","):
          if (len(tokens.value) == 0):
            tok.raise_syntax_error()
          callback(tokens)
          n_callbacks += 1
          if (stop_after_given_number_of_commas is not None
                and n_callbacks == stop_after_given_number_of_commas):
            return
          tokens = tk_seq(ssl=O.ssl, i_code=O.i, value=[])
          tapp = tokens.value.append
          last_comma = tok
        elif (tv == "("):
          nested_tokens = []
          if (O.collect_comma_separated_expressions(
                callback=nested_tokens.append,
                opening_token=tok,
                enable_implied_do=int(enable_implied_do!=0)) is None):
            tk_type = tk_parentheses
          else:
            tk_type = tk_implied_do
          tapp(tk_type(
            ssl=tok.ssl, i_code=tok.i_code, value=nested_tokens))
        elif (tv == "="):
          if (   opening_token is None
              or enable_implied_do == 0
              or n_callbacks == 0
              or i_assignment_op is not None
              or len(tokens.value) != 1
              or not tokens.value[0].is_identifier()):
            tok.raise_syntax_error()
          i_assignment_op = n_callbacks
          tapp(tok)
        else:
          tapp(tok)
      else:
        tapp(tok)

  def collect_to_matching_parenthesis(O, callback, opening_token):
    nested_tokens = []
    O.collect_comma_separated_expressions(
      callback=nested_tokens.append,
      opening_token=opening_token)
    callback(tk_parentheses(
      ssl=opening_token.ssl, i_code=opening_token.i_code, value=nested_tokens))

  def get_implied_do(O, opening_token):
    result = []
    O.collect_comma_separated_expressions(
      callback=result.append,
      opening_token=opening_token,
      enable_implied_do=2)
    return tk_implied_do(ssl=O.ssl, i_code=opening_token.i_code, value=result)

def remove_redundant_parentheses(tokens):
  result = []
  for tok in tokens:
    def handle_seq_with_one_parentheses():
      if (tok.is_seq()):
        tv = tok.value
        if (len(tv) == 1):
          inner_tok = tv[0]
          if (inner_tok.is_parentheses()):
            result.extend(
              remove_redundant_parentheses(tokens=inner_tok.value))
            return
      return result.append(tok)
    handle_seq_with_one_parentheses()
  return result

def group_power(tokens):
  result = []
  i_tok = 0
  n_toks = len(tokens)
  while (i_tok < n_toks):
    tok = tokens[i_tok]; i_tok += 1
    if (tok.is_op_with(value="**")):
      if (len(result) == 0):
        tok.raise_syntax_error()
      if (i_tok == len(tokens)):
        tok.raise_syntax_error()
      if (    result[-1].is_parentheses()
          and len(result) != 1
          and result[-2].is_identifier()):
        base_tok = tk_seq(
          ssl=result[-2].ssl,
          i_code=result[-2].i_code,
          value=result[-2:])
        result.pop()
      else:
        base_tok = result[-1]
      result.pop()
      exponent_tok = tokens[i_tok]; i_tok += 1
      if (    exponent_tok.is_identifier()
          and i_tok < len(tokens)
          and tokens[i_tok].is_parentheses()):
        next_tok = tokens[i_tok]; i_tok += 1
        exponent_tok = tk_seq(
          ssl=exponent_tok.ssl,
          i_code=exponent_tok.i_code,
          value=[exponent_tok, next_tok])
      result.append(tk_power(
        ssl=exponent_tok.ssl,
        i_code=exponent_tok.i_code,
        value=[base_tok, exponent_tok]))
    else:
      result.append(tok)
  return result

class implied_do_info(object):

  __slots__ = ["dlist_size", "id_tok", "fls_tokens"]

  def __init__(O, tokens):
    assert len(tokens) >= 3
    def get_first_tokens(i):
      tok = tokens[i]
      if (not tok.is_seq()): return False
      tv = tok.value
      if (len(tv) < 3): return False
      if (not tv[0].is_identifier()): return False
      if (not tv[1].is_op_with(value="=")): return False
      O.dlist_size = len(tokens) + i
      O.id_tok = tv[0]
      O.fls_tokens = [tk_seq(ssl=tv[2].ssl, i_code=tv[2].i_code, value=tv[2:])]
      return True
    if (get_first_tokens(i=-2)):
      O.fls_tokens.append(tokens[-1])
    else:
      assert get_first_tokens(i=-3)
      O.fls_tokens.extend(tokens[-2:])

def tok_seq_is_star(tok_seq):
  return (len(tok_seq.value) == 1 and tok_seq.value[0].is_op_with(value="*"))

def tokens_as_strings(tokens, result=None):
  if (result is None):
    result = []
  for tok in tokens:
    if (tok.is_seq()):
      result.append("[")
      tokens_as_strings(tokens=tok.value, result=result)
      result.append("]")
    elif (tok.is_parentheses()):
      result.append("(")
      tokens_as_strings(tokens=tok.value, result=result)
      result.append(")")
    elif (tok.is_implied_do()):
      result.append("{")
      tokens_as_strings(tokens=tok.value, result=result)
      result.append("}")
    elif (tok.is_string()):
      result.append(tok.fquote())
    elif (tok.is_power()):
      result.append("power(")
      tokens_as_strings(tokens=tok.value, result=result)
      result.append(")")
    elif (tok.is_format()):
      result.append("format(" + tok.fquote() + ")")
    else:
      result.append(tok.value)
  return result

def tokens_as_string(tokens):
  return " ".join(tokens_as_strings(tokens=tokens))

def tokens_as_python_code(tokens, result=None, allow_power=True):
  join_result = (result is None)
  if (join_result):
    result = []
  for tok in tokens:
    if (tok.is_seq()):
      if (tokens_as_python_code(tokens=tok.value, result=result) is None):
        return None
    elif (tok.is_parentheses()):
      result.append("(")
      if (tokens_as_python_code(tokens=tok.value, result=result) is None):
        return None
      result.append(")")
    elif (tok.is_implied_do()):
      return None
    elif (tok.is_string()):
      return None
    elif (tok.is_power()):
      return None
    elif (tok.is_format()):
      return None
    elif (not allow_power and tok.is_op_with(value="**")):
      return None
    else:
      result.append(tok.value)
  if (join_result):
    return " ".join(result)
  return result

def search_for_id_tokens(callback, tokens, with_next_token):
  for i_tok,tok in enumerate(tokens):
    if (tok.is_seq_or_parentheses_or_implied_do()):
      search_for_id_tokens(
        callback=callback, tokens=tok.value, with_next_token=with_next_token)
    elif (tok.is_identifier()):
      if (not with_next_token):
        callback(tok)
      else:
        next_tok = None
        j_tok = i_tok + 1
        if (j_tok != len(tokens)):
          next_tok = tokens[j_tok]
        callback(tok, next_tok)

def extract_identifiers(tokens, result=None):
  if (result is None): result = []
  search_for_id_tokens(
    callback=result.append, tokens=tokens, with_next_token=False)
  return result

def search_for_data_or_read_target_tokens(callback, tokens):
  for tok in tokens:
    assert tok.is_seq()
    if (len(tok.value) == 0):
      continue
    first_target_tok = tok.value[0]
    if (first_target_tok.is_identifier()):
      callback(tok=first_target_tok)
    elif (first_target_tok.is_implied_do()):
      idi = implied_do_info(tokens=first_target_tok.value)
      search_for_data_or_read_target_tokens(
        callback=callback, tokens=first_target_tok.value[:idi.dlist_size])

class fss_iterator(object):
  "format string stripped iterator"

  __slots__ = ["fss", "i"]

  def __init__(O, fss):
    O.fss = fss
    O.i = 0

  def __iter__(O): return O

  def next(O):
    tok = O.get()
    if (tok is None):
      raise StopIteration
    return tok

  def get(O):
    fss = O.fss
    code = fss.code
    stop = len(code)
    assert O.i is not None
    def raise_invalid():
      fss.raise_error(msg="Invalid FORMAT specification", i=min(stop-1, O.i))
    while (O.i < stop):
      i_code = O.i
      c = code[i_code]
      if (c == ","):
        O.i += 1
        continue
      if (c == "x"):
        O.i += 1
        return tk_format(ssl=fss, i_code=i_code, value="1x")
      if ("():/$".find(c) >= 0):
        O.i += 1
        return tk_op(ssl=fss, i_code=i_code, value=c)
      if (c == "'"):
        O.i += 1
        return tk_string(
          ssl=fss,
          i_code=i_code,
          value=fss.strings[fss.string_indices.index(i_code)])
      if (c == "+" or c == "-"):
        j = unsigned_integer_scan(code=code, start=i_code+1, stop=stop)
        if (j < 0 or not code.startswith("p", j)):
          raise_invalid()
        O.i = j + 1
        return tk_format(
          ssl=fss,
          i_code=i_code,
          value=code[i_code:O.i])
      j = unsigned_integer_scan(code=code, start=i_code, stop=stop)
      if (j > 0):
        O.i = j
        if (code.startswith("h", O.i)):
          # extract Hollerith edit descriptor if it did not confuse the
          # previous ignorant parsing algorithms
          sl, i_text = fss.text_location(i=O.i)
          assert sl.text[i_text].lower() == "h"
          i_text += 1
          nh = int(code[i_code:j])
          j_text = i_text + nh
          tv = sl.text[i_text:j_text]
          if (len(tv) != nh):
            raise RuntimeError(fss.format_error(
              i=i_code,
              msg="FATAL: Not supported: FORMAT Hollerith edit descriptor"
                  " spanning continuation lines"))
          if (tv.find("'") >= 0 or tv.find('"') >= 0):
            raise RuntimeError(fss.format_error(
              i=i_code,
              msg="FATAL: Not supported:"
                  " FORMAT Hollerith edit descriptor with quotes"))
          while (O.i < stop):
            sl, i_text = fss.text_location(i=O.i)
            if (i_text >= j_text):
              break
            O.i += 1
          return tk_string(ssl=fss, i_code=i_code, value=tv)
        if (code.startswith("x", O.i) or code.startswith("p", O.i)):
          O.i += 1
          return tk_format(ssl=fss, i_code=i_code, value=code[i_code:O.i])
        return tk_integer(ssl=fss, i_code=i_code, value=code[i_code:j])
      if ("defgiz".find(c) >= 0):
        j = unsigned_integer_scan(code=code, start=i_code+1, stop=stop)
        if (j > 0):
          O.i = j
          if (code.startswith(".", j)):
            j = unsigned_integer_scan(code=code, start=j+1, stop=stop)
            if (j < 0):
              raise_invalid()
            O.i = j
        return tk_format(ssl=fss, i_code=i_code, value=code[i_code:O.i])
      if (c == "a" or c == "l"):
        O.i += 1
        j = unsigned_integer_scan(code=code, start=i_code+1, stop=stop)
        if (j > 0):
          O.i = j
        return tk_format(ssl=fss, i_code=i_code, value=code[i_code:O.i])
      if (code.startswith("bn", i_code) or code.startswith("bz", i_code)):
        O.i += 2
        return tk_format(ssl=fss, i_code=i_code, value=code[i_code:O.i])
      if (c == "s"):
        O.i += 1
        if (code.startswith("p", O.i) or code.startswith("s", O.i)):
          O.i += 1
        return tk_format(ssl=fss, i_code=i_code, value=code[i_code:O.i])
      if (c == "t"):
        O.i += 1
        if (code.startswith("l", O.i) or code.startswith("r", O.i)):
          O.i += 1
        j = unsigned_integer_scan(code=code, start=O.i, stop=stop)
        if (j < 0):
          raise_invalid()
        O.i = j
        return tk_format(ssl=fss, i_code=i_code, value=code[i_code:O.i])
      raise_invalid()
    O.i = None
    return None

def get_statement_label_token(tokens):
  assert len(tokens) != 0
  tok = tokens[0]
  if (len(tokens) != 1 or not tok.is_integer()):
    tok.raise_error(msg="Invalid statement label")
  return tok

def fmt_tokens_as_string(tokens, comma=","):
  result = []
  def result_ends_with_comma():
    return (len(result) != 0 and result[-1] == comma)
  for tok in tokens:
    if (tok.is_integer()):
      result.append(tok.value)
    elif (tok.is_string()):
      result.append("'" + tok.value.replace("'","''") + "'")
      result.append(comma)
    elif (tok.is_format()):
      result.append(tok.value)
      result.append(comma)
    elif (tok.is_op()):
      tv = tok.value
      if (tv == "(" or tv == ")"):
        if (result_ends_with_comma()):
          result.pop()
        result.append(tv)
        if (tv == ")"): result.append(comma)
      else:
        result.append(tv)
        result.append(comma)
    else:
      raise AssertionError
  if (result_ends_with_comma()):
    result.pop()
  return "".join(result)
