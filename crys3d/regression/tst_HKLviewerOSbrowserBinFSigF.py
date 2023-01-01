from __future__ import absolute_import, division, print_function
from crys3d.regression import tests_HKLviewer

# Using a webbrowser
# Test creating the F/SigF dataset on the fly from iotbx/regression/data/phaser_1.mtz,
# expand data to P! with Friedel pairs, slice with a clip plane at l=9 and only show
# reflections with F/SigF<=1

def run():
  count = 0
  while True:
    print("running %d" %count)
    # exercise_OSbrowser() is unstable and might yield a bogus failure. If so, repeat the test
    # at most maxruns times or until it passes
    if not tests_HKLviewer.runagain(tests_HKLviewer.exercise_OSbrowser,
                                    tests_HKLviewer.philstr3,
                                    tests_HKLviewer.reflections2match3,
                                    "OSbrowserBinFSigF"):
      break
    count +=1
    assert(count < tests_HKLviewer.maxruns)


if __name__ == '__main__':
  run()
