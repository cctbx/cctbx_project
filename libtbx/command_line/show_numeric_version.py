from __future__ import absolute_import, division, print_function
def run():
  try: import Numeric
  except ImportError: print("None")
  else: print(Numeric.__version__)

if (__name__ == "__main__"):
  run()
