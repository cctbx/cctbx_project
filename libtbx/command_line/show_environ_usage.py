import sys

def run(args):
  assert len(args) == 0
  print """
LIBTBX_DISABLE_TRACEBACKLIMIT
  If set, Sorry and Usage exceptions are shown with the full traceback.

LIBTBX_VALGRIND
  Run "libtbx.valgrind python" for more information.

LIBTBX_PRINT_TRACE
  If set, print trace of all Python code executed.
  This can lead to very large output.

LIBTBX_NATIVE_TAR
  Inspected by libtbx.bundle_as_selfx to find alternative tar command.
  Example: setenv LIBTBX_NATIVE_TAR $HOME/bin/tar

LIBTBX_FULL_TESTING
  If set, forces libtbx.env.full_testing = True.
"""

if (__name__ == "__main__"):
  run(sys.argv[1:])
