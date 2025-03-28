from __future__ import absolute_import, division, print_function

def test_1():
  """
  Test generate_n_char_string method, iterator that generates strings of
  length n_char
  """

  from iotbx.pdb.utils import generate_n_char_string
  string_list = []
  g = generate_n_char_string(n_chars=2)
  for i in range(1000):
    string_list.append(g.next())
  assert [string_list[55], string_list[900], string_list[87]] == "A3 Og BZ".split(), "%s %s %s" %(string_list[55], string_list[900], string_list[87])

  string_list = []
  g = generate_n_char_string(n_chars=2, include_upper=False)
  for i in range(1000):
    string_list.append(g.next())
  assert [string_list[55], string_list[900], string_list[87]] == "bt za cp".split(), "%s %s %s" %(string_list[55], string_list[900], string_list[87])

  string_list = []
  g = generate_n_char_string(n_chars=2, include_numbers=False)
  for i in range(1000):
    string_list.append(g.next())
  assert [string_list[55], string_list[900], string_list[87]] == "BD RQ Bj".split(), "%s %s %s" %(string_list[55], string_list[900], string_list[87])

  string_list = []
  g = generate_n_char_string(n_chars=1, include_special_chars=True)
  for i in range(120):
    string_list.append(g.next())
    if string_list[-1] is None:
      string_list = string_list[:-1]
      break
  assert [string_list[55], string_list[89], string_list[86]] == " 3 + %".split(), "%s %s %s" %(string_list[55], string_list[89], string_list[86])

  string_list = []
  g = generate_n_char_string(n_chars=3, include_upper=False,
    include_lower=False,
   include_numbers=True,
   include_special_chars=False,
   end_with_tilde=True)
  for i in range(1000):
    string_list.append(g.next())
    if string_list[-1] is None: break
  assert [string_list[55], string_list[57], string_list[23]] == "55~ 57~ 23~".split(), "%s %s %s" %(string_list[55], string_list[57], string_list[23])

  print("OK")

if (__name__ == "__main__"):
  test_1()
