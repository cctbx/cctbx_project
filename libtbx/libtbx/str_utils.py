from __future__ import generators

def size_as_string_with_commas(sz):
  if (sz is None): return "unknown"
  if (sz < 0):
    sz = -sz
    sign = "-"
  else:
    sign = ""
  result = []
  while True:
    if (sz >= 1000):
      result.insert(0, "%03d" % (sz % 1000))
      sz //= 1000
      if (sz == 0): break
    else:
      result.insert(0, "%d" % sz)
      break
  return sign + ",".join(result)

def show_string(s):
  if (s is None): return None
  if (s.find('"') < 0): return '"'+s+'"'
  if (s.find("'") < 0): return "'"+s+"'"
  return '"'+s.replace('"','\\"')+'"'

def overwrite_at(s, offset, replacement):
  return s[:offset] + replacement + s[offset+len(replacement):]

def contains_one_of(label, patterns):
  for pattern in patterns:
    if (label.find(pattern) >= 0):
      return True
  return False

def line_breaker(string, width):
  if (width <= 0 or len(string) == 0):
    yield string
  else:
    i_block_start = 0
    i_last_space = None
    for i,c in enumerate(string):
      if (c == " "):
        i_last_space = i
      if (i-i_block_start >= width and i_last_space is not None):
        yield string[i_block_start:i_last_space]
        i_block_start = i_last_space + 1
        i_last_space = None
    if (i_block_start < len(string)):
      yield string[i_block_start:]

def exercise():
  import libtbx.forward_compatibility
  assert size_as_string_with_commas(0) == "0"
  assert size_as_string_with_commas(1) == "1"
  assert size_as_string_with_commas(-1) == "-1"
  assert size_as_string_with_commas(10) == "10"
  assert size_as_string_with_commas(100) == "100"
  assert size_as_string_with_commas(1000) == "1,000"
  assert size_as_string_with_commas(12345) == "12,345"
  assert size_as_string_with_commas(12345678) == "12,345,678"
  assert size_as_string_with_commas(-12345678) == "-12,345,678"
  assert show_string("abc") == '"abc"'
  assert show_string("a'c") == '"a\'c"'
  assert show_string('a"c') == "'a\"c'"
  assert show_string('\'"c') == '"\'\\"c"'
  for string, expected_result in [
    ("", [""]),
    ("this is", ["this is"]),
    ("this is a", ["this is", "a"]),
    ("this is a sentence", ["this is", "a", "sentence"]),
    ("this is a longer sentence", ["this is", "a", "longer", "sentence"]),
    ("this is a very long sentence indeed",
      ["this is", "a very", "long", "sentence", "indeed"])]:
    assert [block for block in line_breaker(string, width=7)]==expected_result
  print "OK"

if (__name__ == "__main__"):
  exercise()
