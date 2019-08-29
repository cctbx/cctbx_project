from __future__ import absolute_import, division, print_function
import sys
from time_cmd import run_string

def run():
  cmd = "fftpacktimer"
  factor = 1
  for arg in sys.argv[1:]:
    try:
      factor = int(arg)
    except Exception:
      cmd = arg
  for fft_type in ("cc", "rc"):
    for N, iter in (
      (2*3*4*5, 86000),
      (2*3*3*4*4*5*5, 620),
      (2*3*4*5*7, 6400),
      (2*3*3*4*4*5*5*7*7, 2),
                    ):
      t = run_string("%s %s %d %d %d" % (cmd, fft_type, N, iter, factor))
      print("u+s:", t, "t/s:", iter * factor / t)
      sys.stdout.flush()

if (__name__ == "__main__"):
  run()
