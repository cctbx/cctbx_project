from __future__ import absolute_import, division, print_function
from types import FrameType
from six.moves import cStringIO as StringIO
import gc
import sys

# =============================================================================
# This method is written by Michael Droettboom and is found at
#
#   http://code.activestate.com/recipes/523004/
#
# -----------------------------------------------------------------------------
def print_cycles(objects, outstream=sys.stdout, show_progress=False):
    """
    objects:       A list of objects to find cycles in.  It is often useful
                   to pass in gc.garbage to find the cycles that are
                   preventing some objects from being garbage collected.
    outstream:     The stream for output.
    show_progress: If True, print the number of objects reached as they are
                   found.
    """
    def print_path(path):
        for i, step in enumerate(path):
            # next "wraps around"
            next = path[(i + 1) % len(path)]

            outstream.write("   %s -- " % str(type(step)))
            if isinstance(step, dict):
                for key, val in step.items():
                    if val is next:
                        outstream.write("[%s]" % repr(key))
                        break
                    if key is next:
                        outstream.write("[key] = %s" % repr(val))
                        break
            elif isinstance(step, list):
                outstream.write("[%d]" % step.index(next))
            elif isinstance(step, tuple):
                outstream.write("[%d]" % list(step).index(next))
            else:
                outstream.write(repr(step))
            outstream.write(" ->\n")
        outstream.write("\n")

    def recurse(obj, start, all, current_path):
        if show_progress:
            outstream.write("%d\r" % len(all))

        all[id(obj)] = None

        referents = gc.get_referents(obj)
        for referent in referents:
            # If we've found our way back to the start, this is
            # a cycle, so print it out
            if referent is start:
                print_path(current_path)

            # Don't go back through the original list of objects, or
            # through temporary references to the object, since those
            # are just an artifact of the cycle detector itself.
            elif referent is objects or isinstance(referent, FrameType):
                continue

            # We haven't seen this object before, so recurse
            elif id(referent) not in all:
                recurse(referent, start, all, current_path + [obj])

    for obj in objects:
        outstream.write("Examining: %r\n" % obj)
        recurse(obj, obj, { }, [])

# =============================================================================
# test for print_cycles
# -----------------------------------------------------------------------------
def exercise():

  class bad_class(object):

    def __init__(self):
      self.a = self.b

    def b():
      return 1+2

  t = bad_class()
  sio = StringIO()
  print_cycles([t],outstream=sio)
  lines = sio.getvalue().splitlines()

  # 'Examining: <__main__.bad_class object at 0xb7eb5d2c>'
  # "   <class '__main__.bad_class'> -- <__main__.bad_class object at 0xb7eb5d2c> ->"
  # "   <type 'dict'> --  ->"
  # ''
  # For Python 3.11, the <type 'dict'> line is missing
  if sys.version_info.major == 3 and sys.version_info.minor >= 11:
    assert len(lines) == 3
  else:
    assert len(lines) == 4
  assert lines[0].startswith("Examining: <")
  assert lines[0].count("bad_class") == 1
  assert lines[1].count("bad_class") == 2
  if sys.version_info.major == 3 and sys.version_info.minor < 11:
    assert lines[2].count("dict") == 1

  print("OK")

# =============================================================================
if (__name__ == "__main__"):
  exercise()
