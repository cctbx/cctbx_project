class character_iterator:

  def __init__(self, input_string):
    self.remaining = input_string

  def next(self):
    if (len(self.remaining) == 0): return None
    result = self.remaining[0]
    self.remaining = self.remaining[1:]
    return result

  def look_ahead(self):
    if (len(self.remaining) == 0): return None
    return self.remaining[0]

class word:

  def __init__(self, value, quote_char=None):
    self.value = value
    self.quote_char = quote_char

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
      while True:
        c = char_iter.next()
        if (c == quote_char): break
        if (c == "\\" and char_iter.look_ahead() == quote_char):
          c = char_iter.next()
        if (c is None): raise RuntimeError("No closing quote.")
        word_value += c
      words.append(word(value=word_value, quote_char=quote_char))
      c = char_iter.next()
    else:
      word_value = c
      if (c not in contiguous_word_characters):
        c = char_iter.next()
      else:
        while True:
          c = char_iter.next()
          if (c is None): break
          if (c not in contiguous_word_characters
              and (not enable_unquoted_embedded_quotes
                   or c not in ['"', "'"])):
            break
          word_value += c
      words.append(word(value=word_value))
  return words
