import UnitTest, regression
import os, sys
from sys import stderr

workDir = os.path.dirname(sys.argv[0]).strip()
if workDir:
  os.chdir(workDir)
  # else we must be there already
print("working directory=", os.getcwd())


def fullTest():
  print("----------------------------  unit test  -----------------------------", file=stderr)
  UnitTest.testAll()

  print("--------------------------  regression test  -------------------------", file=stderr)
  regression.test()


if __name__ == '__main__':
  fullTest()
