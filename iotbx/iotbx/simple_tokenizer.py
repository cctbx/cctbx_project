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

def split_into_words(
      input_string,
      contiguous_word_characters="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                 "abcdefghijklmnopqrstuvwxyz"
                                 "0123456789"
                                 "_"):
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
          if (c not in contiguous_word_characters): break
          word_value += c
      words.append(word(value=word_value))
  return words

def exercise():
  tests = [
  ["",
    []],
  ["resname=a and chain=b",
    ['resname', '=', 'a', 'and', 'chain', '=', 'b']],
  ["resname a and chain b",
    ['resname', 'a', 'and', 'chain', 'b']],
  ["resname resname and chain chain",
    ['resname', 'resname', 'and', 'chain', 'chain']],
  ["resname \"a b\"",
    ['resname', 'a b']],
  ["resname a",
    ['resname', 'a']],
  ["resname ala and backbone",
    ['resname', 'ala', 'and', 'backbone']],
  ["resname ala or backbone",
    ['resname', 'ala', 'or', 'backbone']],
  ["name x and x > 10",
    ['name', 'x', 'and', 'x', '>', '10']],
  ["((expr or expr) and expr)",
    ['(', '(', 'expr', 'or', 'expr', ')', 'and', 'expr', ')']],
  ["resname and and chain b",
    ['resname', 'and', 'and', 'chain', 'b']],
  ["resname ( and chain b",
    ['resname', '(', 'and', 'chain', 'b']],
  ["resname \"(\" and chain b",
    ['resname', '(', 'and', 'chain', 'b']],
  ["all_hydrophobic_within(5) and resname ALA",
    ['all_hydrophobic_within', '(', '5', ')', 'and', 'resname', 'ALA']],
  ["something(a, b)",
    ['something', '(', 'a', ',', 'b', ')']],
  ["something(a b)",
    ['something', '(', 'a', 'b', ')']],
  ["something(\"a\"\"b\")",
    ['something', '(', 'a', 'b', ')']],
  ["resname 'a'",
    ['resname', 'a']],
  ["resname '\"'",
    ['resname', '"']],
  ["resname '\"\\''",
    ['resname', '"\'']],
  ["resname \"'\\\"\"",
    ['resname', '\'"']],
  ]
  for input_string,expected_result in tests:
    assert [word.value for word in split_into_words(input_string=input_string)
           ] == expected_result
  print "OK"

if (__name__ == "__main__"):
  exercise()
