class character_iterator:

  def __init__(self, input_string):
    self.remaining = input_string
    self.line_number = 1

  def next(self):
    if (len(self.remaining) == 0): return None
    result = self.remaining[0]
    self.remaining = self.remaining[1:]
    if (result == "\n"): self.line_number += 1
    return result

  def look_ahead(self, n=1):
    if (len(self.remaining) < n): return None
    return self.remaining[:n]

class word:

  def __init__(self,
        value,
        quote_char=None,
        line_number=None):
    self.value = value
    self.quote_char = quote_char
    self.line_number = line_number

default_contiguous_word_characters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ" \
                                   + "abcdefghijklmnopqrstuvwxyz" \
                                   + "0123456789" \
                                   + "_"
def split_into_words(
      input_string,
      contiguous_word_characters=None,
      enable_unquoted_embedded_quotes=True):
  if (contiguous_word_characters is None):
    contiguous_word_characters = default_contiguous_word_characters
  words = []
  char_iter = character_iterator(input_string)
  c = char_iter.next()
  while True:
    if (c is None): break
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
        triple_quote = True
      else:
        triple_quote = False
      while True:
        c = char_iter.next()
        if (c == quote_char):
          if (not triple_quote): break
          if (char_iter.look_ahead(n=2) == quote_char+quote_char):
            char_iter.next()
            char_iter.next()
            break
        if (c == "\\"):
          if (char_iter.look_ahead() == quote_char):
            c = char_iter.next()
          elif (char_iter.look_ahead() == "\n"):
            char_iter.next()
            c = char_iter.next()
        if (c is None): raise RuntimeError("No closing quote.")
        word_value += c
      words.append(word(
        value=word_value,
        quote_char=quote_char,
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
      words.append(word(
        value=word_value,
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
    if (result.quote_char is not None):
      raise RuntimeError('Unquoted word expected: line=%d, word="%s"' % (
        result.line_number, result.value))
    return result

  def push_back(self, word):
    self.stack.append(word)

def as_word_stack(
      input_string,
      contiguous_word_characters=None,
      enable_unquoted_embedded_quotes=True):
  return word_stack(split_into_words(
    input_string=input_string,
    contiguous_word_characters=contiguous_word_characters,
    enable_unquoted_embedded_quotes=enable_unquoted_embedded_quotes))
