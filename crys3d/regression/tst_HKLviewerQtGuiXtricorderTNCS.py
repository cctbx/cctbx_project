from __future__ import absolute_import, division, print_function
from crys3d.regression import tests_PhenixHKLviewer as tsthkl

# With the HKLviewer Qt GUI
# test creating the F/SigF dataset on the fly from iotbx/regression/data/phaser_1.mtz,
# expand data to P! with Friedel pairs, slice with a clip plane at l=9 and only show
# reflections with F/SigF<=1

def run():
  count = 0
  while True:
    print("running %d" %count)
    # exerciseQtGUI() is unstable and might yield a bogus failure. If so, repeat the test
    # at most maxruns times or until it passes
    if not tsthkl.runagain(tsthkl.exerciseQtGUI,
                                    tsthkl.philstr1,
                                    tsthkl.reflections2match1,
                                    "QtGuiXtricorderTNCS"):
      break
    count +=1
    assert(count < tsthkl.maxruns)


if __name__ == '__main__':
  run()
