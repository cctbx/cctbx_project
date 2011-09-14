import random
import sys

def run(args, command_name="libtbx.obfuscate"):
  if (len(args) != 2):
    print
    print "Interleaves string with random digits."
    print "Main purpose: maintenance of libtbx/windows_dispatcher.c"
    print
    print "usage:   %s # STRING" % command_name
    print "example: %s 20 ABC" % command_name
    print
    return
  l = int(args[0])
  s = ""
  for c in args[1]:
    i = random.randrange(10)
    s += str(i)
    s += c
  i = random.randrange(10)
  s += str(i)
  while (len(s) < l): s += "_"+s
  s = s[:l]
  n = 77
  i = 0
  while (i < len(s)):
    print '"%s"' % s[i:i+n]
    i += n

if (__name__ == "__main__"):
    run(sys.argv[1:])
