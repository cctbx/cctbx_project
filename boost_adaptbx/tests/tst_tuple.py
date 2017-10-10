from __future__ import division
from __future__ import print_function
import boost.tuple

def exercise():
  doc = boost.tuple.exercise.__doc__
  assert doc.replace("\n","").startswith("exercise(")
  assert boost.tuple.exercise(1) == (2, 0.5)
  assert boost.tuple.exercise(2) == (4, 1)

def run():
  exercise()
  print('OK')

if __name__ == '__main__':
  run()
