# LIBTBX_SET_DISPATCHER_NAME boost_adaptbx.show_environ_usage

import sys

def run(args):
  assert len(args) == 0
  print """
BOOST_ADAPTBX_SIGNALS_DEFAULT
  If NOT set, enable Python and libc call stack traces if possible.
BOOST_ADAPTBX_FPE_DEFAULT
  If NOT set, trap floating-point exceptions if possible.
BOOST_ADAPTBX_FE_DIVBYZERO_DEFAULT
BOOST_ADAPTBX_FE_INVALID_DEFAULT
BOOST_ADAPTBX_FE_OVERFLOW_DEFAULT
  If NOT set, trap FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW if possible.
  See also: man fenv

BOOST_ADAPTBX_DOCSTRING_OPTIONS
  Example:
    setenv BOOST_ADAPTBX_DOCSTRING_OPTIONS "show_user_defined=True,show_signatures=False"
"""

if (__name__ == "__main__"):
  run(sys.argv[1:])
