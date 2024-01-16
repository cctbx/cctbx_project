# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
from crys3d.regression import tests_HKLviewer
from libtbx.test_utils import contains_substring
import os, sys, subprocess

os.environ['PYTHONIOENCODING'] = 'UTF-8'

# Using the HKLviewer Qt GUI exerciseQtGUI() runs HKLviewer to enact the settings in philstr3 and
# eventually asserts that the visible reflections in the browser match the indices in
# reflections2match3. Due to occasional instability of websockets on virtual machines the test
# is run in a loop until it passes but no longer than maxruns times.

# This test is creating the F/SigF dataset on the fly from iotbx/regression/data/phaser_1.mtz,
# expands data to P1 with Friedel pairs, slices with a clip plane at l=9 and only shows
# reflections with F/SigF<=1. Then saves those reflections to a new file OSbrowserLowValueBinFSigF.mtz.
# Then checks that this file matches the info in expectedstr

expectedstr = """
Starting job
===============================================================================
1 Miller arrays in this dataset:
 Labels          |       Type      |   λ/Å   |  #HKLs  |               Span              |     min,max data       |     min,max sigmas     |  d_min,d_max/Å   |Anomalous|Sym.uniq.|Data compl.|
  LowValuesFSigF |       Amplitude |       1 |     363 |           (-9, 0, 0), (9, 9, 9) |    0.78407,       11.26|        nan,         nan|     2.5,    8.771|   False |    True |   0.33925 |

===============================================================================
Job complete
"""

def run():
  if 'linux' in sys.platform and os.environ.get("DISPLAY", None) != "0:0":
    return # Only run test if on a local linux machine with a connected monitor
  count = 0
  while True:
    print("running %d" %count)
    # exerciseQtGUI() is unstable and might yield a bogus failure. If so, repeat the test
    # at most maxruns times or until it passes
    if not tests_HKLviewer.runagain(tests_HKLviewer.exerciseQtGUI,
                                    tests_HKLviewer.philstr3 %"QtGuiLowValueBinFSigF.mtz",
                                    tests_HKLviewer.reflections2match3,
                                    "QtGuiBinFSigF"):
      break
    count +=1
    assert(count < tests_HKLviewer.maxruns)

  # Now check that the produced mtz file matches the info in expectedstr
  obj = subprocess.Popen("cctbx.HKLinfo QtGuiLowValueBinFSigF.mtz",
                          shell=True,
                          stdin = subprocess.PIPE,
                          stdout = subprocess.PIPE,
                          stderr = subprocess.STDOUT)
  souterr,err = obj.communicate()
  tests_HKLviewer.Append2LogFile("QtGuiBinFSigFHKLviewer.log", souterr)
  souterr = souterr.decode().replace("\r\n", "\n") # omit \r\n line endings on Windows
  assert (contains_substring( souterr, expectedstr ) )


if __name__ == '__main__':
  run()
