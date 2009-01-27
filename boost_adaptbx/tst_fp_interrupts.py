import boost.python
ext = boost.python.import_ext("tst_fp_interrupts_ext")

class should_have_crashed:
  """ Execution flow should never have reached that point"""

def run(case):
  arguments = {"division_by_zero": (1,0),
               "inexact": (0,0),
               "overflow": (1.e300, 1.e-300)}[case]
  a = ext.division(*arguments)
  raise should_have_crashed

if __name__ == '__main__':
  import sys
  run(sys.argv[1])
