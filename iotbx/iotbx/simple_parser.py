from __future__ import generators
from iotbx import simple_tokenizer
import sys

class operator_priority_evaluator:

  def __init__(self, operator_dict):
    self.operator_dict = operator_dict

  def __call__(self, word):
    if (word.quote_char is not None): return 0
    return self.operator_dict.get(word.value, 0)

def infix_as_postfix(words, operator_dict={"not": 3, "and": 2, "or": 1}):
  """http://www.programmersheaven.com/2/Art_Expressions_p1"""
  operator_priority = operator_priority_evaluator(operator_dict=operator_dict)
  word_stack = words[:]
  word_stack.reverse()
  parse_stack = []
  while (len(word_stack) > 0):
    word = word_stack.pop()
    if (word.value == "("):
      parse_stack.append(word)
    elif (word.value == ")"):
      while True:
        if (len(parse_stack) == 0):
          raise RuntimeError(
            "Closing parenthesis without a matching opening parenthesis.");
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
  while (len(parse_stack) > 0):
    yield parse_stack.pop(), None

def exercise():
  tests = [
  ["a",
    ['a']
  ],
  ["a and b",
    ['a', 'b', 'and']
  ],
  ["a or b",
    ['a', 'b', 'or']
  ],
  ["not a or b",
    ['a', 'not', 'b', 'or']
  ],
  ["not a or b and c",
    ['a', 'not', 'b', 'c', 'and', 'or']
  ],
  ["not (a or b) and c",
    ['a', 'b', 'or', 'not', 'c', 'and']
  ],
  ["(not (a or b) and c)",
    ['a', 'b', 'or', 'not', 'c', 'and']
  ],
  ["not ((a or b) and c)",
    ['a', 'b', 'or', 'c', 'and', 'not']
  ],
  ]
  verbose = "--verbose" in sys.argv[1:]
  for input_string,expected_result in tests:
    infix = simple_tokenizer.split_into_words(input_string=input_string)
    if (verbose): print input_string
    postfix = [word for word,word_stack in infix_as_postfix(infix)]
    if (verbose): print [word.value for word in postfix]
    assert [word.value for word in postfix] == expected_result
    if (verbose): print
  print "OK"

if (__name__ == "__main__"):
  exercise()
