from __future__ import absolute_import, division, print_function
from iotbx import simple_parser
from libtbx.phil import tokenizer
import sys

def exercise_basic():
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
    infix = tokenizer.word_iterator(input_string=input_string)
    if (verbose): print(input_string)
    postfix = [word
      for word,word_iterator in simple_parser.infix_as_postfix(infix)]
    if (verbose): print([word.value for word in postfix])
    assert [word.value for word in postfix] == expected_result
    if (verbose): print()

def rewrite_parser(
      word_iterator,
      stop_if_parse_stack_is_empty=False,
      stop_word=None,
      expect_nonmatching_closing_parenthesis=False):
  result_stack = []
  for word,word_iterator in simple_parser.infix_as_postfix(
         word_iterator=word_iterator,
         stop_if_parse_stack_is_empty=stop_if_parse_stack_is_empty,
         stop_word=stop_word,
         expect_nonmatching_closing_parenthesis
           =expect_nonmatching_closing_parenthesis):
    if (word.value == "not"):
      arg = result_stack.pop()
      result_stack.append("(!%s)" % arg)
    elif (word.value in ["and", "or"]):
      rhs = result_stack.pop()
      lhs = result_stack.pop()
      if (word.value == "and"):
        result_stack.append("(%s&%s)" % (lhs, rhs))
      else:
        result_stack.append("(%s|%s)" % (lhs, rhs))
    elif (word.value == "within"):
      assert word_iterator.pop().value == "("
      radius = float(word_iterator.pop().value)
      assert word_iterator.pop().value == ","
      nested_result = rewrite_parser(
        word_iterator=word_iterator,
        expect_nonmatching_closing_parenthesis=True)
      if (nested_result == ""): raise RuntimeError("Missing argument.")
      result_stack.append("@(%.2f,%s)" % (radius, nested_result))
    elif (word.value == "around"):
      assert word_iterator.pop().value == "("
      nested_result = rewrite_parser(
        word_iterator=word_iterator,
        stop_word=",")
      if (nested_result == ""): raise RuntimeError("Missing argument.")
      radius = float(word_iterator.pop().value)
      assert word_iterator.pop().value == ")"
      result_stack.append("@(%.2f,%s)" % (radius, nested_result))
    elif (word.value == "for"):
      var = word_iterator.pop().value
      assert word_iterator.pop().value == "in"
      nested_result = rewrite_parser(
        word_iterator=word_iterator,
        stop_if_parse_stack_is_empty=True)
      if (nested_result == ""): raise RuntimeError("Missing argument.")
      result_stack.append("(for %s in %s)" % (var, nested_result))
    else:
      result_stack.append(word.value)
  if (len(result_stack) == 0):
    return ""
  result = result_stack[0]
  for item in result_stack[1:]:
    result = "(%s&%s)" % (result, item)
  return result

def rewrite(input_string):
  word_iterator = tokenizer.word_iterator(input_string=input_string)
  return rewrite_parser(word_iterator=word_iterator)

def exercise_nested():
  tests = [
  ["a",
    "a",
  ],
  ["a and b",
    "(a&b)",
  ],
  ["a and b or c",
    "((a&b)|c)",
  ],
  ["a or b and c",
    "(a|(b&c))",
  ],
  ["within(5, a or b and c)",
    "@(5.00,(a|(b&c)))",
  ],
  ["within(5, (a or b) and c)",
    "@(5.00,((a|b)&c))",
  ],
  ["around(a or b and c, 5)",
    "@(5.00,(a|(b&c)))",
  ],
  ["around((a or b) and c, 5)",
    "@(5.00,((a|b)&c))",
  ],
  ["around((a or b) and within(3, c or d), 5)",
    "@(5.00,((a|b)&@(3.00,(c|d))))"
  ],
  ["for i in a",
    "(for i in a)",
  ],
  ["for i in not a",
    "(for i in (!a))",
  ],
  ["for i in a or b",
    "((for i in a)|b)",
  ],
  ["for i in (a or b)",
    "(for i in (a|b))",
  ],
  ["for i in (a or b) and c",
    "((for i in (a|b))&c)",
  ],
  ]
  verbose = "--verbose" in sys.argv[1:]
  for input_string,expected_result in tests:
    show = verbose or expected_result is None
    if (show): print(input_string)
    result = rewrite(input_string=input_string)
    if (show): print(result)
    if (expected_result is not None):
      assert result == expected_result
    if (show): print()

def exercise():
  exercise_basic()
  exercise_nested()
  print("OK")

if (__name__ == "__main__"):
  exercise()
