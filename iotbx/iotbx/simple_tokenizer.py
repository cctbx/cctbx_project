class character_iterator:

  def __init__(self, input_string):
    self.input_string = input_string
    self.i_char = 0
    self.line_number = 1

  def next(self):
    if (self.i_char == len(self.input_string)): return None
    result = self.input_string[self.i_char]
    self.i_char += 1
    if (result == "\n"): self.line_number += 1
    return result

  def look_ahead(self, n=1):
    end = self.i_char + n
    if (end > len(self.input_string)): return None
    return self.input_string[self.i_char:end]

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

class word(object):

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

default_contiguous_word_characters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ" \
                                   + "abcdefghijklmnopqrstuvwxyz" \
                                   + "0123456789" \
                                   + "_"
def split_into_words(
      input_string,
      file_name=None,
      unquoted_single_character_words="",
      contiguous_word_characters=None,
      enable_unquoted_embedded_quotes=True):
  if (contiguous_word_characters is None):
    contiguous_word_characters = default_contiguous_word_characters
  words = []
  char_iter = character_iterator(input_string)
  c = char_iter.next()
  while (c is not None):
    if (c.isspace()):
      c = char_iter.next()
      continue
    if (c in ['"', "'"]):
      quote_char = c
      word_value = ""
      word_line_number = char_iter.line_number
      if (char_iter.look_ahead(n=2) == quote_char+quote_char):
        char_iter.next()
        char_iter.next()
        quote_token = quote_char+quote_char+quote_char
      else:
        quote_token = quote_char
      while True:
        c = char_iter.next()
        if (c == quote_char):
          if (quote_token is quote_char): break
          if (char_iter.look_ahead(n=2) == quote_char+quote_char):
            char_iter.next()
            char_iter.next()
            break
        if (c == "\\"):
          if (char_iter.look_ahead() == "\\"):
            char_iter.next()
          elif (char_iter.look_ahead() == quote_char):
            c = char_iter.next()
          elif (char_iter.look_ahead() == "\n"):
            char_iter.next()
            c = char_iter.next()
        if (c is None):
          raise RuntimeError(
            "Syntax error: missing closing quote (%s)" % (
              where(file_name, char_iter.line_number)))
        word_value += c
      words.append(word(
        value=word_value,
        quote_token=quote_token,
        line_number=word_line_number,
        file_name=file_name))
      c = char_iter.next()
    else:
      word_value = c
      word_line_number = char_iter.line_number
      if (c in unquoted_single_character_words
          or (contiguous_word_characters != ""
              and c not in contiguous_word_characters)):
        c = char_iter.next()
      else:
        while True:
          c = char_iter.next()
          if (c is None): break
          if (c.isspace()): break
          if (c in unquoted_single_character_words): break
          if (contiguous_word_characters != ""
              and c not in contiguous_word_characters
              and (not enable_unquoted_embedded_quotes
                   or c not in ['"', "'"])):
            break
          word_value += c
      words.append(word(
        value=word_value,
        line_number=word_line_number,
        file_name=file_name))
  return words

class word_stack:

  def __init__(self, words):
    self.words = words
    self.i_word = 0

  def __len__(self):
    return len(self.words) - self.i_word

  def pop(self):
    if (self.i_word == len(self.words)):
      raise RuntimeError("Unexpected end of input.")
    word = self.words[self.i_word]
    self.i_word += 1
    return word

  def pop_unquoted(self):
    word = self.pop()
    if (word.quote_token is not None):
      raise RuntimeError('Unquoted word expected, found %s%s' % (
        str(word), word.where_str()))
    return word

  def push_back(self, word):
    self.i_word -= 1
    assert self.words[self.i_word] is word

def as_word_stack(
      input_string,
      file_name=None,
      unquoted_single_character_words="",
      contiguous_word_characters=None,
      enable_unquoted_embedded_quotes=True):
  return word_stack(split_into_words(
    input_string=input_string,
    file_name=file_name,
    unquoted_single_character_words=unquoted_single_character_words,
    contiguous_word_characters=contiguous_word_characters,
    enable_unquoted_embedded_quotes=enable_unquoted_embedded_quotes))
