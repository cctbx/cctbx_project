from __future__ import generators

class operator_priority_evaluator:

  def __init__(self, operator_dict):
    self.operator_dict = operator_dict

  def __call__(self, word):
    if (word.quote_char is not None): return 0
    return self.operator_dict.get(word.value, 0)

def infix_as_postfix(
      words=None,
      word_stack=None,
      operator_dict={"not": 3, "and": 2, "or": 1},
      stop_word=None,
      expect_nonmatching_closing_parenthesis=False):
  """http://www.programmersheaven.com/2/Art_Expressions_p1"""
  assert (words is None) != (word_stack is None)
  operator_priority = operator_priority_evaluator(operator_dict=operator_dict)
  if (word_stack is None):
    word_stack = words[:]
    word_stack.reverse()
  parse_stack = []
  while (len(word_stack) > 0):
    word = word_stack.pop()
    if (stop_word is not None and word.value == stop_word):
      break
    if (word.value == "("):
      parse_stack.append(word)
    elif (word.value == ")"):
      while True:
        if (len(parse_stack) == 0):
          if (expect_nonmatching_closing_parenthesis): return
          raise RuntimeError(
            "Closing parenthesis without a matching opening parenthesis.")
        item = parse_stack.pop()
        if (item.value == "("): break
        yield item, word_stack
    else:
      word_priority = operator_priority(word)
      if (word_priority == 0):
        yield word, word_stack
      else:
        while True:
          if (len(parse_stack) == 0
              or parse_stack[-1].value == "("
              or word_priority > operator_priority(parse_stack[-1])):
            parse_stack.append(word)
            break
          else:
            yield parse_stack.pop(), None
  assert not expect_nonmatching_closing_parenthesis
  while (len(parse_stack) > 0):
    yield parse_stack.pop(), None
