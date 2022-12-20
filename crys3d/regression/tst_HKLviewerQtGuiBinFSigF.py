from __future__ import absolute_import, division, print_function
from crys3d.regression import tests_HKLviewer


def run():
  tests_HKLviewer.exerciseQtGUI(tests_HKLviewer.philstr3, tests_HKLviewer.reflections2match3,
                               "QtGuiBinFSigF")
  print("OK")


if __name__ == '__main__':
  run()
