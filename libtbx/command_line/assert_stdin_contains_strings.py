from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME libtbx.assert_stdin_contains_strings
# LIBTBX_SET_DISPATCHER_NAME libtbx.assert_stdin_does_not_contain_strings
import libtbx.load_env
import sys

def run(args):
  dn = libtbx.env.dispatcher_name
  does_not = (dn.find("does_not") >= 0)
  input = "\n".join(sys.stdin.read().splitlines()) + "\n"
  for arg in args:
    s = "\n".join(arg.splitlines())
    if ((input.find(s) < 0) == does_not): continue
    print("BEGIN OF INPUT (%s)" % dn)
    print("v" * 79)
    sys.stdout.write(input)
    print("^" * 79)
    print("END OF INPUT (%s)" % dn)
    if (does_not): s = "found"
    else:          s = "not"
    raise AssertionError("String %s in input: %s" % (s, arg))

if (__name__ == "__main__"):
  run(sys.argv[1:])
