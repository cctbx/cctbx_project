class character_iterator:

  def __init__(self, input_string):
    self.remaining = input_string
    self.i_char = 0
    self.line_number = 1

  def next(self):
    if (len(self.remaining) == 0): return None
    result = self.remaining[0]
    self.remaining = self.remaining[1:]
    self.i_char += 1
    if (result == "\n"): self.line_number += 1
    return result

  def look_ahead(self, n=1):
    if (len(self.remaining) < n): return None
    return self.remaining[:n]

class word:

  def __init__(self,
        value,
        quote_token=None,
        line_number=None):
    self.value = value
    self.quote_token = quote_token
    self.line_number = line_number

  def __str__(self):
    if (self.quote_token is None):
      return self.value
    return self.quote_token \
         + self.value \
            .replace("\\", "\\\\") \
            .replace(self.quote_token, "\\"+self.quote_token) \
         + self.quote_token

  def raise_syntax_error(self, message):
    raise RuntimeError(
      'Syntax error: %s"%s" (input line %d)' % (
        message, self.value, self.line_number))

  def assert_expected(self, value):
    if (self.value != value):
      self.raise_syntax_error('expected "{", found ')

default_contiguous_word_characters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ" \
                                   + "abcdefghijklmnopqrstuvwxyz" \
                                   + "0123456789" \
                                   + "_"
def split_into_words(
      input_string,
      contiguous_word_characters=None,
      enable_unquoted_embedded_quotes=True,
      auto_split_unquoted={}):
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
            "Syntax error: missing closing quote (input line %d)." % (
              char_iter.line_number))
        word_value += c
      words.append(word(
        value=word_value,
        quote_token=quote_token,
        line_number=word_line_number))
      c = char_iter.next()
    else:
      word_value = c
      word_line_number = char_iter.line_number
      if (contiguous_word_characters != ""
          and c not in contiguous_word_characters):
        c = char_iter.next()
      else:
        while True:
          c = char_iter.next()
          if (c is None): break
          if (c.isspace()): break
          if (contiguous_word_characters != ""
              and c not in contiguous_word_characters
              and (not enable_unquoted_embedded_quotes
                   or c not in ['"', "'"])):
            break
          word_value += c
      if (word_value not in auto_split_unquoted):
        words.append(word(
          value=word_value,
          line_number=word_line_number))
      else:
        for value in auto_split_unquoted[word_value]:
          words.append(word(
            value=value,
            line_number=word_line_number))
  return words

class word_stack:

  def __init__(self, words):
    self.stack = words[:]
    self.stack.reverse()

  def __len__(self):
    return len(self.stack)

  def pop(self):
    if (len(self.stack) == 0):
      raise RuntimeError("Unexpected end of input.")
    return self.stack.pop()

  def pop_unquoted(self):
    result = self.pop()
    if (result.quote_token is not None):
      raise RuntimeError('Unquoted word expected: line=%d, word="%s"' % (
        result.line_number, result.value))
    return result

  def push_back(self, word):
    self.stack.append(word)

def as_word_stack(
      input_string,
      contiguous_word_characters=None,
      enable_unquoted_embedded_quotes=True,
      auto_split_unquoted={}):
  return word_stack(split_into_words(
    input_string=input_string,
    contiguous_word_characters=contiguous_word_characters,
    enable_unquoted_embedded_quotes=enable_unquoted_embedded_quotes,
    auto_split_unquoted=auto_split_unquoted))
