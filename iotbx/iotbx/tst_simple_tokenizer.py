from iotbx import simple_tokenizer

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
  ["name o1'",
    ['name', 'o1\'']],
  ]
  for input_string,expected_result in tests:
    assert [word.value
      for word in simple_tokenizer.split_into_words(input_string=input_string)
           ] == expected_result
  print "OK"

if (__name__ == "__main__"):
  exercise()
