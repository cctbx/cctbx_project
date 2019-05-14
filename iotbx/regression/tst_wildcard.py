from __future__ import absolute_import, division, print_function
from iotbx import wildcard

def exercise():
  assert wildcard.is_match(string="", pattern="")
  assert wildcard.is_match("", "*")
  assert wildcard.is_match("a", "?")
  assert wildcard.is_match("a", "[a-z]")
  assert wildcard.is_match("wildcard", "wi*card")
  assert not wildcard.is_match("wildcard", "wi*cart")
  assert wildcard.is_match("wildcard", "wi*c*d")
  assert not wildcard.is_match("wildcard", "wi*z*d")
  assert wildcard.is_match("wildcard", "w*c??d")
  assert not wildcard.is_match("wildcard", "w*d??d")
  assert wildcard.is_match("wildcard", "w*c[a]?d")
  assert not wildcard.is_match("wildcard", "w*c[^a]?d")
  assert not wildcard.is_match("wild*ard", r"wild\*ard")
  assert wildcard.is_match("wild*ard", r"wild\*ard", escape_char='\\')
  assert wildcard.is_match("wild*ard", r"wild\**", escape_char='\\')
  assert wildcard.is_match(
    r"\*?[^a-z]", r"\\\*\?\[\^\a\-\z\]", escape_char='\\')
  assert not wildcard.is_match(
    r"\*?[^a-z]", r"\\\*\?\[\^\a\-\z\]")
  print("OK")

if (__name__ == "__main__"):
  exercise()
