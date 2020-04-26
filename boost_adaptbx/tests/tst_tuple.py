from __future__ import absolute_import, division, print_function

import boost_adaptbx.tuple

def exercise():
  doc = boost_adaptbx.tuple.exercise.__doc__
  assert doc.replace("\n","").startswith("exercise(")
  assert boost_adaptbx.tuple.exercise(1) == (2, 0.5)
  assert boost_adaptbx.tuple.exercise(2) == (4, 1)

def run():
  exercise()
  print('OK')

if __name__ == '__main__':
  run()
