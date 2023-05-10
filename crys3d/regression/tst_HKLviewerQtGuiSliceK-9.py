from __future__ import absolute_import, division, print_function
from crys3d.regression import tests_HKLviewer
import os, sys

# Using the HKLviewer Qt GUI exerciseQtGUI() runs HKLviewer to enact the settings in philstr1 and
# eventually asserts that the visible reflections in the browser match the indices in
# reflections2match1. Due to occasional instability of websockets on virtual machines the test
# is run in a loop until it passes but no longer than maxruns times.

# With the HKLviewer Qt GUI test expanding amplitude data from iotbx/regression/data/phaser_1.mtz
# to P1 with Friedel pairs, slice with a clip plane at k= -9

def run():
  if 'linux' in sys.platform and os.environ.get("DISPLAY", None) != "0:0":
    return # Only run test if on a local linux machine with a connected monitor
  count = 0
  while True:
    print("running %d" %count)
    # exerciseQtGUI() is unstable and might yield a bogus failure. If so, repeat the test
    # at most maxruns times or until it passes
    if not tests_HKLviewer.runagain(tests_HKLviewer.exerciseQtGUI,
                                    tests_HKLviewer.philstr1,
                                    tests_HKLviewer.reflections2match1,
                                    "QtGuiSliceK-9" ):
      break
    count +=1
    assert(count < tests_HKLviewer.maxruns)

if __name__ == '__main__':
  run()
