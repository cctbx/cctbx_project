from __future__ import absolute_import, division, print_function
from libtbx.phil import tokenizer

def exercise_basic(verbose):
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
  ["resname 'a \\\\'",
    ['resname', 'a \\']],
  ["resname 'a'",
    ['resname', 'a']],
  ["resname '\"'",
    ['resname', '"']],
  ["resname '\"\\''",
    ['resname', '"\'']],
  ["resname \"'\\\"\"",
    ['resname', '\'"']],
  ["name o1'",
    ['name', 'o1\'']],
  ['name """o1\'"""',
    ['name', 'o1\'']],
  ['name """o1\n  o2\'"""',
    ['name', "o1\n  o2'"]],
  ['name """o1\\\n  o2\'"""',
    ['name', "o1  o2'"]],
  ]
  for input_string,expected_result in tests:
    show = verbose or expected_result is None
    if (show): print(input_string)
    result = [word.value
      for word in tokenizer.word_iterator(input_string=input_string)]
    if (show): print(result)
    if (expected_result is not None):
      assert result == expected_result
    if (show): print()

def exercise_pickle():
  # TODO: verify this is intended change for py2/3 compat
  from six.moves import cPickle as pickle
  for p in [pickle]:
    o = tokenizer.word(value="hello")
    l = p.loads(p.dumps(o))
    assert l.value == "hello"
    o = tokenizer.settings(meta_comment="%")
    l = p.loads(p.dumps(o))
    assert l.meta_comment == "%"
    o = tokenizer.word_iterator(input_string="all")
    l = p.loads(p.dumps(o))
    assert l.char_iter.input_string == "all"

def run(args):
  assert args in [[], ["--verbose"]]
  verbose = len(args) != 0
  exercise_basic(verbose=verbose)
  exercise_pickle()
  print("OK")

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
