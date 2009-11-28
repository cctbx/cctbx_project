def usage():
  from libtbx.utils import Usage
  import libtbx.load_env
  raise Usage("%s timeout file..." % libtbx.env.dispatcher_name)

def run(args):
  if (len(args) < 2): usage()
  timeout = float(args[0])
  if (timeout <= 0): usage()
  file_names = args[1:]
  #
  import libtbx.load_env
  import time
  import os
  op = os.path
  time_start = time.time()
  while True:
    for file_name in file_names:
      if (not op.exists(file_name)):
        if (time.time() - time_start > timeout):
          raise RuntimeError(
            "%s timeout exceeded." % libtbx.env.dispatcher_name)
        time.sleep(0.1)
        break
    else:
      break

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
