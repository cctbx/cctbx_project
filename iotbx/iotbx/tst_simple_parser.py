from iotbx import simple_parser
from iotbx import simple_tokenizer
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
    infix = simple_tokenizer.split_into_words(input_string=input_string)
    if (verbose): print input_string
    postfix = [word
      for word,word_stack in simple_parser.infix_as_postfix(infix)]
    if (verbose): print [word.value for word in postfix]
    assert [word.value for word in postfix] == expected_result
    if (verbose): print

def rewrite_parser(
      word_stack,
      stop_word=None,
      expect_nonmatching_closing_parenthesis=False):
  result_stack = []
  for word,word_stack in simple_parser.infix_as_postfix(
         word_stack=word_stack,
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
      assert word_stack.pop().value == "("
      radius = float(word_stack.pop().value)
      assert word_stack.pop().value == ","
      sel = rewrite_parser(
        word_stack=word_stack,
        expect_nonmatching_closing_parenthesis=True)
      if (sel == ""): raise RuntimeError("Missing argument.")
      result_stack.append("@(%.2f,%s)" % (radius, sel))
    elif (word.value == "around"):
      assert word_stack.pop().value == "("
      sel = rewrite_parser(
        word_stack=word_stack,
        stop_word=",")
      if (sel == ""): raise RuntimeError("Missing argument.")
      radius = float(word_stack.pop().value)
      assert word_stack.pop().value == ")"
      result_stack.append("@(%.2f,%s)" % (radius, sel))
    else:
      result_stack.append(word.value)
  if (len(result_stack) == 0):
    return ""
  result = result_stack[0]
  for item in result_stack[1:]:
    result = "(%s&%s)" % (result, item)
  return result

def rewrite(input_string):
  word_stack = simple_tokenizer.split_into_words(input_string=input_string)
  word_stack.reverse()
  return rewrite_parser(word_stack=word_stack)

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
  ]
  verbose = "--verbose" in sys.argv[1:]
  for input_string,expected_result in tests:
    if (verbose): print input_string
    result = rewrite(input_string=input_string)
    if (verbose): print result
    if (expected_result is not None):
      assert result == expected_result
    if (verbose): print

def exercise():
  exercise_basic()
  exercise_nested()
  print "OK"

if (__name__ == "__main__"):
  exercise()
