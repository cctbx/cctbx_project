from __future__ import absolute_import, division, print_function
from crys3d.regression import tests_HKLviewer

# Using a webbrowser test expanding amplitude data from iotbx/regression/data/phaser_1.mtz
# to P! with Friedel pairs, slice with a clip plane at k= -9

def run():
  count = 0
  while True:
    print("running %d" %count)
    # exercise_OSbrowser() is unstable and might yield a bogus failure. If so, repeat the test
    # at most maxruns times or until it passes
    if not tests_HKLviewer.runagain(tests_HKLviewer.exercise_OSbrowser,
                                    tests_HKLviewer.philstr1,
                                    tests_HKLviewer.reflections2match1,
                                    "OSbrowserSliceK-9" ):
      break
    count +=1
    assert(count < tests_HKLviewer.maxruns)

if __name__ == '__main__':
  run()
