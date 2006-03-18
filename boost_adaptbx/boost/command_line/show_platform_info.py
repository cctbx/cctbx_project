import boost.python
import sys

def run():
  sys.stdout.write(boost.python.platform_info)
  print "sys.byteorder:", sys.byteorder
  try: import thread
  except ImportError: print "import thread: NO"
  else: print "import thread: OK"

if (__name__ == "__main__"):
  run()
