from __future__ import absolute_import, division, print_function
from crys3d.regression import tests_HKLviewer


def run():
  tests_HKLviewer.exerciseQtGUI(tests_HKLviewer.philstr1, tests_HKLviewer.reflections2match1,
                               "QtGuiSliceK-9")
  print("OK")


if __name__ == '__main__':
  run()
