from __future__ import absolute_import, division, print_function

import boost_adaptbx_char_array_ext as ext
import sys

def exercise_char_n():
  assert ext.char_3_holder().value == "bar"
  assert ext.char_5_holder().value == "barte"
  assert ext.use_char_n("rab") == "bar"
  assert ext.use_char_n("etrab") == "barte"

def run(args):
  iterations = 100
  if (len(args) > 0):
    iterations = int(args[0])
  i = 0
  while (iterations == 0 or i < iterations):
    exercise_char_n()
    i += 1
  print("OK")

if (__name__ == "__main__"):
  run(sys.argv[1:])
