from __future__ import absolute_import, division, print_function
from crys3d.regression import tests_HKLviewer

# With the HKLviewer Qt GUI
# test creating the F/SigF dataset on the fly from iotbx/regression/data/phaser_1.mtz,
# expand data to P! with Friedel pairs, slice with a clip plane at l=9 and only show
# reflections with F/SigF<=1

def run():
  count = 0
  while True:
    print("running %d" %count)
    if not tests_HKLviewer.runagain(tests_HKLviewer.exerciseQtGUI,
                                    tests_HKLviewer.philstr3,
                                    tests_HKLviewer.reflections2match3,
                                    "QtGuiBinFSigF"):
      break
    count +=1
    assert(count < 3)


if __name__ == '__main__':
  run()
