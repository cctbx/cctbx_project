"""Tools for parsing
"""
from __future__ import absolute_import, division, print_function
class operator_priority_evaluator(object):

  def __init__(self, operator_dict):
    self.operator_dict = operator_dict

  def __call__(self, word):
    if (word.quote_token is not None): return 0
    return self.operator_dict.get(word.value.lower(), 0)

def infix_as_postfix(
      word_iterator,
      operator_dict={"not": 3, "and": 2, "or": 1},
      stop_if_parse_stack_is_empty=False,
      stop_word=None,
      expect_nonmatching_closing_parenthesis=False):
  """http://www.programmersheaven.com/2/Art_Expressions_p1"""
  operator_priority = operator_priority_evaluator(operator_dict=operator_dict)
  parse_stack = []
  while True:
    word = word_iterator.try_pop()
    if (word is None): break
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
        yield item, word_iterator
      if (len(parse_stack) == 0 and stop_if_parse_stack_is_empty): return
    else:
      word_priority = operator_priority(word)
      if (word_priority == 0):
        yield word, word_iterator
        if (len(parse_stack) == 0 and stop_if_parse_stack_is_empty): return
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
