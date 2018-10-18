from __future__ import division
from __future__ import print_function
def run():
  try: import Numeric
  except ImportError: print("None")
  else: print(Numeric.__version__)

if (__name__ == "__main__"):
  run()
