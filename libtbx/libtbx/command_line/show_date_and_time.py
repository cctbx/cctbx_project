import libtbx.utils
import time
import sys
s = libtbx.utils.date_and_time()
if ("-v" in sys.argv[1:]):
  s += " Seconds since the Epoch %.2f" % time.time()
print s
