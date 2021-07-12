import UnitTest, download
import os, sys, subprocess
from sys import stderr

workDir = os.path.dirname(sys.argv[0]).strip()
if workDir:
    os.chdir(workDir)
    # else we must be there already
print("working directory=", os.getcwd())


def execute(program, args=None):
    argList = [program] + args
    code = subprocess.call(argList)
    return code


def regressionTest(molecule):
    code = 9
    if download.download(molecule):
      if execute("convert.bat", [molecule]) == 0:
          code = execute("regression.bat", [molecule])
          os.remove(molecule + ".py-out.txt")
      else:
          print("unable to convert", molecule)
      os.remove(molecule + ".pdb")
    else:
      print("unable to download", molecule)
    return code


def fullTest():
    print("----------------------------  unit test  -----------------------------", file=stderr)
    UnitTest.test()

    print("--------------------------  regression test  -------------------------", file=stderr)
    code = regressionTest("1ehz")
    if code == 0:
        print("success", file=stderr)


if __name__ == '__main__':
    fullTest()