from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME boost_adaptbx.show_environ_usage

import sys

def run(args):
  assert len(args) == 0
  print("""
BOOST_ADAPTBX_TRAP_FPE
BOOST_ADAPTBX_TRAP_INVALID
BOOST_ADAPTBX_TRAP_OVERFLOW
  If set to a non-empty value (0 included), trap FE_DIVBYZERO,
  FE_INVALID or FE_OVERFLOW respectively; set to an empty value, switch
  the trap off
  (see boost_adaptbx.boost.python.floating_point_exceptions and man fenv).
BOOST_ADAPTBX_ENABLE_TRACE
  If set to a non-empty value, print a libc call stack trace on a
  crash, if possible.

BOOST_ADAPTBX_DOCSTRING_OPTIONS
  Example:
    setenv BOOST_ADAPTBX_DOCSTRING_OPTIONS "show_user_defined=True,show_signatures=False"
""")

if (__name__ == "__main__"):
  run(sys.argv[1:])
