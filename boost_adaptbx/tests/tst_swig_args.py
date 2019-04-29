from __future__ import absolute_import, division, print_function
import example
import boost_python_swig_args_ext
import sys

def exercise():
  c = example.Circle(10)
  c.x = 20
  c.y = 30

  s = example.Square(10)
  s.x = -10
  s.y = 5

  forever = "--forever" in sys.argv[1:]
  while True:
    boost_python_swig_args_ext.show(c.this)
    boost_python_swig_args_ext.show(s.this)
    if (not forever): break

if (__name__ == "__main__"):
  exercise()
