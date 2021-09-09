from __future__ import nested_scopes, generators, division, absolute_import
from __future__ import  with_statement, print_function, unicode_literals
import UnitTest, regression
import os, sys

workDir = os.path.dirname(sys.argv[0]).strip()
if workDir:
  os.chdir(workDir)
  # else we must be there already

def fullTest():
  UnitTest.testAll()
  regression.test()

if __name__ == '__main__':
  fullTest()
