class character_iterator:

  def __init__(self, input_string):
    self.input_string = input_string
    self.i_char = 0
    self.line_number = 1

  def mark_for_backup(self):
    self.backup_i_char = self.i_char
    self.backup_line_number = self.line_number

  def restore_backup(self):
    self.i_char = self.backup_i_char
    self.line_number = self.backup_line_number

  def look_ahead(self, n):
    end = self.i_char + n
    if (end > len(self.input_string)): return None
    return self.input_string[self.i_char:end]

  def look_ahead_1(self):
    if (self.i_char == len(self.input_string)): return None
    return self.input_string[self.i_char]

  def skip_ahead_1(self):
    if (self.input_string[self.i_char] == "\n"): self.line_number += 1
    self.i_char += 1

  def next(self):
    result = self.look_ahead_1()
    if (result is not None):
      self.i_char += 1
      if (result == "\n"): self.line_number += 1
    return result

def where(file_name, line_number):
  result = []
  if (file_name is not None):
    result.append('file "%s"' % file_name)
    if (line_number is not None):
      result.append("line %d" % line_number)
  elif (line_number is not None):
    result.append("input line %d" % line_number)
  if (len(result) == 0): return None
  return ", ".join(result)

class word: # FUTURE word(object)

  __slots__ = ["value", "quote_token", "line_number", "file_name"]

  def __init__(self,
        value,
        quote_token=None,
        line_number=None,
        file_name=None):
    self.value = value
    self.quote_token = quote_token
    self.line_number = line_number
    self.file_name = file_name

  def __str__(self):
    if (self.quote_token is None):
      return self.value
    return self.quote_token \
         + self.value \
            .replace("\\", "\\\\") \
            .replace(self.quote_token, "\\"+self.quote_token) \
         + self.quote_token

  def where(self):
    return where(self.file_name, self.line_number)

  def where_str(self):
    where = self.where()
    if (self.where() is None): return ""
    return " (%s)" % where

  def raise_syntax_error(self, message):
    raise RuntimeError(
      'Syntax error: %s"%s"%s' % (message, self.value, self.where_str()))

  def assert_expected(self, value):
    if (self.value != value):
      self.raise_syntax_error('expected "%s", found ' % value)

  def raise_unquoted_word_expected(self):
    raise RuntimeError('Unquoted word expected, found %s%s' % (
      str(self), self.where_str()))

default_contiguous_word_characters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ" \
                                   + "abcdefghijklmnopqrstuvwxyz" \
                                   + "0123456789" \
                                   + "_"
class settings:

  def __init__(self,
        unquoted_single_character_words="",
        contiguous_word_characters=None,
        enable_unquoted_embedded_quotes=True,
        comment_characters=""):
    self.unquoted_single_character_words = unquoted_single_character_words
    if (contiguous_word_characters is None):
      self.contiguous_word_characters = default_contiguous_word_characters
    else:
      self.contiguous_word_characters = contiguous_word_characters
    self.enable_unquoted_embedded_quotes = enable_unquoted_embedded_quotes
    self.comment_characters = comment_characters

class word_iterator:

  def __init__(self, input_string, file_name=None, list_of_settings=None):
    self.char_iter = character_iterator(input_string)
    self.file_name = file_name
    if (list_of_settings is None):
      self.list_of_settings = [settings()]
    else:
      self.list_of_settings = list_of_settings

  def __iter__(self):
    return self

  def next(self, settings_index=0):
    settings = self.list_of_settings[settings_index]
    char_iter = self.char_iter
    char_iter.mark_for_backup()
    while True:
      c = char_iter.next()
      if (c is None): break
      if (c.isspace()): continue
      if (c in settings.comment_characters
          and char_iter.look_ahead_1().isspace()):
        while True:
          c = char_iter.next()
          if (c is None or c == "\n"): break
      elif (c in ['"', "'"]):
        quote_char = c
        word_value = ""
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
            if (char_iter.look_ahead_1() == "\\"):
              char_iter.skip_ahead_1()
            elif (char_iter.look_ahead_1() == quote_char):
              c = char_iter.next()
            elif (char_iter.look_ahead_1() == "\n"):
              char_iter.skip_ahead_1()
              c = char_iter.next()
          if (c is None):
            raise RuntimeError(
              "Syntax error: missing closing quote (%s)" % (
                where(self.file_name, char_iter.line_number)))
          word_value += c
        return word(
          value=word_value,
          quote_token=quote_token,
          line_number=word_line_number,
          file_name=self.file_name)
      else:
        word_value = c
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
            word_value += c
            char_iter.skip_ahead_1()
        return word(
          value=word_value,
          line_number=word_line_number,
          file_name=self.file_name)
    raise StopIteration

  def try_pop(self, settings_index=0):
    try: return self.next(settings_index)
    except StopIteration: return None

  def pop(self, settings_index=0):
    try: return self.next(settings_index)
    except StopIteration: raise RuntimeError("Unexpected end of input.")

  def try_pop_unquoted(self, settings_index=0):
    word = self.try_pop(settings_index)
    if (word is not None):
      if (word.quote_token is not None): word.raise_unquoted_word_expected()
    return word

  def pop_unquoted(self, settings_index=0):
    word = self.pop(settings_index)
    if (word.quote_token is not None): word.raise_unquoted_word_expected()
    return word

  def backup(self):
    self.char_iter.restore_backup()
