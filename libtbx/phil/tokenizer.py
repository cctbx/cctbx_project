from libtbx.str_utils import string_representation
from libtbx import Auto, slots_getstate_setstate

def escape_python_str(quote_char, string):
  return string.replace("\\", "\\\\").replace(quote_char, "\\"+quote_char)

def quote_python_str(quote_token, string):
  return quote_token \
       + escape_python_str(quote_char=quote_token[0], string=string) \
       + quote_token

class character_iterator(slots_getstate_setstate):

  __slots__ = [
    "input_string",
    "i_char",
    "line_number",
    "backup_i_char",
    "backup_line_number"]

  def __init__(O, input_string):
    O.input_string = input_string
    O.i_char = 0
    O.line_number = 1
    O.backup_i_char = None
    O.backup_line_number = None

  def mark_for_backup(O):
    O.backup_i_char = O.i_char
    O.backup_line_number = O.line_number

  def restore_backup(O):
    O.i_char = O.backup_i_char
    O.line_number = O.backup_line_number

  def look_ahead(O, n):
    end = O.i_char + n
    if (end > len(O.input_string)): return None
    return O.input_string[O.i_char:end]

  def look_ahead_1(O):
    if (O.i_char == len(O.input_string)): return None
    return O.input_string[O.i_char]

  def skip_ahead_1(O):
    if (O.input_string[O.i_char] == "\n"): O.line_number += 1
    O.i_char += 1

  def next(O):
    result = O.look_ahead_1()
    if (result is not None):
      O.i_char += 1
      if (result == "\n"): O.line_number += 1
    return result

  def scan_for_start(O, intro, followups):
    while True:
      if (O.i_char == len(O.input_string)): return 0
      if (O.input_string[O.i_char] != "\n"):
        O.i_char += 1
      else:
        while True:
          O.line_number += 1
          O.i_char += 1
          if (O.i_char == len(O.input_string)): return 0
          if (O.input_string[O.i_char] != "\n"): break
        if (O.input_string.find(
              intro, O.i_char, O.i_char+len(intro)) != O.i_char):
          O.i_char += 1
        else:
          O.i_char += len(intro)
          while True:
            if (O.i_char == len(O.input_string)): return 0
            if (O.input_string[O.i_char].isspace()): break
            O.i_char += 1
          while True:
            if (O.i_char == len(O.input_string)): return 0
            if (not O.input_string[O.i_char].isspace()): break
            O.i_char += 1
          for i_followup,followup in enumerate(followups):
            if (O.input_string.find(
                  followup,
                  O.i_char,
                  O.i_char+len(followup)) == O.i_char):
              O.i_char += len(followup)
              while True:
                if (O.i_char == len(O.input_string)): return i_followup
                c = O.input_string[O.i_char]
                O.i_char += 1
                if (c == "\n"):
                  O.line_number += 1
                  return i_followup
                if (not c.isspace()): break
              break
    return 0

def where(source_info, line_number):
  result = []
  if (source_info is not None):
    result.append(source_info)
    if (line_number is not None):
      result.append("line %d" % line_number)
  elif (line_number is not None):
    result.append("input line %d" % line_number)
  if (len(result) == 0): return None
  return ", ".join(result)

def where_str(source_info, line_number):
  s = where(source_info, line_number)
  if (s is None): return ""
  return " (%s)" % s

class word(slots_getstate_setstate):

  __slots__ = ["value", "quote_token", "line_number", "source_info"]

  def __init__(O,
        value,
        quote_token=None,
        line_number=None,
        source_info=None):
    O.value = value
    O.quote_token = quote_token
    O.line_number = line_number
    O.source_info = source_info

  def __str__(O):
    if (O.quote_token is None):
      return O.value
    if (O.quote_token is Auto):
      return string_representation(
        string=O.value, preferred_quote='"', alternative_quote='"')
    return quote_python_str(quote_token=O.quote_token, string=O.value)

  def where(O):
    return where(O.source_info, O.line_number)

  def where_str(O):
    return where_str(O.source_info, O.line_number)

  def raise_syntax_error(O, message):
    raise RuntimeError(
      'Syntax error: %s"%s"%s' % (message, O.value, O.where_str()))

  def assert_expected(O, value):
    if (O.value != value):
      O.raise_syntax_error('expected "%s", found ' % value)

  def raise_unquoted_word_expected(O):
    raise RuntimeError('Unquoted word expected, found %s%s' % (
      str(O), O.where_str()))

default_contiguous_word_characters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ" \
                                   + "abcdefghijklmnopqrstuvwxyz" \
                                   + "0123456789" \
                                   + "_"
class settings(slots_getstate_setstate):

  __slots__ = [
    "unquoted_single_character_words",
    "contiguous_word_characters",
    "enable_quoted_t_n_r_x_escapes",
    "enable_unquoted_embedded_quotes",
    "comment_characters",
    "meta_comment"]

  def __init__(O,
        unquoted_single_character_words="",
        contiguous_word_characters=None,
        enable_quoted_t_n_r_x_escapes=False,
        enable_unquoted_embedded_quotes=True,
        comment_characters="",
        meta_comment=None):
    O.unquoted_single_character_words = unquoted_single_character_words
    if (contiguous_word_characters is None):
      O.contiguous_word_characters = default_contiguous_word_characters
    else:
      O.contiguous_word_characters = contiguous_word_characters
    O.enable_quoted_t_n_r_x_escapes = enable_quoted_t_n_r_x_escapes
    O.enable_unquoted_embedded_quotes = enable_unquoted_embedded_quotes
    O.comment_characters = comment_characters
    O.meta_comment = meta_comment

class word_iterator(slots_getstate_setstate):

  __slots__ = [
    "char_iter",
    "source_info",
    "list_of_settings"]

  def __init__(O,
        input_string,
        source_info=None,
        file_name=None,
        list_of_settings=None):
    assert source_info is None or file_name is None
    O.char_iter = character_iterator(input_string)
    if (source_info is None and file_name is not None):
      O.source_info = 'file "%s"' % file_name
    else:
      O.source_info = source_info
    if (list_of_settings is None):
      O.list_of_settings = [settings()]
    else:
      O.list_of_settings = list_of_settings

  def __iter__(O):
    return O

  def next(O, settings_index=0):
    settings = O.list_of_settings[settings_index]
    char_iter = O.char_iter
    char_iter.mark_for_backup()
    while True:
      c = char_iter.next()
      if (c is None): break
      if (c.isspace()): continue
      if (c in settings.comment_characters
          and (settings.meta_comment is None
               or char_iter.look_ahead(n=len(settings.meta_comment))
                  != settings.meta_comment)):
        while True:
          c = char_iter.next()
          if (c is None or c == "\n"): break
      elif (c in ['"', "'"]):
        quote_char = c
        word_value = []
        word_line_number = char_iter.line_number
        if (char_iter.look_ahead(n=2) == quote_char+quote_char):
          char_iter.skip_ahead_1()
          char_iter.skip_ahead_1()
          quote_token = quote_char+quote_char+quote_char
        else:
          quote_token = quote_char
        while True:
          c = char_iter.next()
          if (c == quote_char):
            if (quote_token is quote_char): break
            if (char_iter.look_ahead(n=2) == quote_char+quote_char):
              char_iter.skip_ahead_1()
              char_iter.skip_ahead_1()
              break
          if (c == "\\"):
            _ = char_iter.look_ahead_1()
            if (_ == "\\"):
              char_iter.skip_ahead_1()
            elif (_ == quote_char):
              c = char_iter.next()
            elif (_ == "\n"):
              char_iter.skip_ahead_1()
              continue
            elif (settings.enable_quoted_t_n_r_x_escapes):
              if (_ in "tnr"):
                char_iter.skip_ahead_1()
                c = eval('"'+c+_+'"')
              elif (_ == "x"):
                hex = char_iter.look_ahead(n=3)
                hex_chars = "0123456789ABCDEFabcdef"
                if (hex is not None
                      and hex[1] in hex_chars
                      and hex[2] in hex_chars):
                  for _ in xrange(3):
                    char_iter.skip_ahead_1()
                  c = eval('"'+c+hex+'"')
          if (c is None):
            raise RuntimeError(
              "Syntax error: missing closing quote%s" % (
                where_str(O.source_info, char_iter.line_number)))
          word_value.append(c)
        return word(
          value="".join(word_value),
          quote_token=quote_token,
          line_number=word_line_number,
          source_info=O.source_info)
      else:
        word_value = [c]
        word_line_number = char_iter.line_number
        if (c not in settings.unquoted_single_character_words
            and (settings.contiguous_word_characters == ""
                 or c in settings.contiguous_word_characters)):
          while True:
            c = char_iter.look_ahead_1()
            if (c is None): break
            if (c.isspace()): break
            if (c in settings.unquoted_single_character_words): break
            if (settings.contiguous_word_characters != ""
                and c not in settings.contiguous_word_characters
                and (not settings.enable_unquoted_embedded_quotes
                     or c not in ['"', "'"])):
              break
            word_value.append(c)
            char_iter.skip_ahead_1()
        return word(
          value="".join(word_value),
          line_number=word_line_number,
          source_info=O.source_info)
    raise StopIteration

  def try_pop(O, settings_index=0):
    try: return O.next(settings_index)
    except StopIteration: return None

  def pop(O, settings_index=0):
    try: return O.next(settings_index)
    except StopIteration: raise RuntimeError("Unexpected end of input.")

  def try_pop_unquoted(O, settings_index=0):
    word = O.try_pop(settings_index)
    if (word is not None):
      if (word.quote_token is not None): word.raise_unquoted_word_expected()
    return word

  def pop_unquoted(O, settings_index=0):
    word = O.pop(settings_index)
    if (word.quote_token is not None): word.raise_unquoted_word_expected()
    return word

  def backup(O):
    O.char_iter.restore_backup()

  def scan_for_start(O, intro, followups):
    O.char_iter.scan_for_start(intro, followups)
