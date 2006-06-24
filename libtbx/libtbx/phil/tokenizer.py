class character_iterator(object):

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

  def scan_for_start(self, intro, followups):
    while True:
      if (self.i_char == len(self.input_string)): return 0
      if (self.input_string[self.i_char] != "\n"):
        self.i_char += 1
      else:
        while True:
          self.line_number += 1
          self.i_char += 1
          if (self.i_char == len(self.input_string)): return 0
          if (self.input_string[self.i_char] != "\n"): break
        if (self.input_string.find(
              intro, self.i_char, self.i_char+len(intro)) != self.i_char):
          self.i_char += 1
        else:
          self.i_char += len(intro)
          while True:
            if (self.i_char == len(self.input_string)): return 0
            if (self.input_string[self.i_char].isspace()): break
            self.i_char += 1
          while True:
            if (self.i_char == len(self.input_string)): return 0
            if (not self.input_string[self.i_char].isspace()): break
            self.i_char += 1
          for i_followup,followup in enumerate(followups):
            if (self.input_string.find(
                  followup,
                  self.i_char,
                  self.i_char+len(followup)) == self.i_char):
              self.i_char += len(followup)
              while True:
                if (self.i_char == len(self.input_string)): return i_followup
                c = self.input_string[self.i_char]
                self.i_char += 1
                if (c == "\n"):
                  self.line_number += 1
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

class word: # FUTURE word(object)

  __slots__ = ["value", "quote_token", "line_number", "source_info"]

  def __init__(self,
        value,
        quote_token=None,
        line_number=None,
        source_info=None):
    self.value = value
    self.quote_token = quote_token
    self.line_number = line_number
    self.source_info = source_info

  def __str__(self):
    if (self.quote_token is None):
      return self.value
    return self.quote_token \
         + self.value \
            .replace("\\", "\\\\") \
            .replace(self.quote_token, "\\"+self.quote_token) \
         + self.quote_token

  def where(self):
    return where(self.source_info, self.line_number)

  def where_str(self):
    return where_str(self.source_info, self.line_number)

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
class settings(object):

  def __init__(self,
        unquoted_single_character_words="",
        contiguous_word_characters=None,
        enable_unquoted_embedded_quotes=True,
        comment_characters="",
        meta_comment=None):
    self.unquoted_single_character_words = unquoted_single_character_words
    if (contiguous_word_characters is None):
      self.contiguous_word_characters = default_contiguous_word_characters
    else:
      self.contiguous_word_characters = contiguous_word_characters
    self.enable_unquoted_embedded_quotes = enable_unquoted_embedded_quotes
    self.comment_characters = comment_characters
    self.meta_comment = meta_comment

class word_iterator(object):

  def __init__(self,
        input_string,
        source_info=None,
        file_name=None,
        list_of_settings=None):
    assert source_info is None or file_name is None
    self.char_iter = character_iterator(input_string)
    if (source_info is None and file_name is not None):
      self.source_info = 'file "%s"' % file_name
    else:
      self.source_info = source_info
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
          and (settings.meta_comment is None
               or char_iter.look_ahead(n=len(settings.meta_comment))
                  != settings.meta_comment)):
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
              "Syntax error: missing closing quote%s" % (
                where_str(self.source_info, char_iter.line_number)))
          word_value += c
        return word(
          value=word_value,
          quote_token=quote_token,
          line_number=word_line_number,
          source_info=self.source_info)
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
          source_info=self.source_info)
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

  def scan_for_start(self, intro, followups):
    self.char_iter.scan_for_start(intro, followups)
