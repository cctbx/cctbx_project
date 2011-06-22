import sys
sys.stderr = sys.stdout

import cgitb
if (__name__ == "__main__"):
  cgitb.enable()

from cctbx.web import cgi_utils
import os, cgi

if (0):
  cgi.test()
  sys.exit(0)

def run():
  print "Content-type: text/html"
  print

  server_info = cgi_utils.server_info()

  form = cgi.FieldStorage()
  target_module = form["target_module"].value

  exec "import " + target_module + " as target"
  inp = target.interpret_form_data(form)

  # optionally capture input to facilitate debugging
  capture_input_dir = "/var/tmp/cctbx_web"
  if (capture_input_dir is not None and os.path.isdir(capture_input_dir)):
    import time, pickle
    time_stamp = "%d_%02d_%02d_%02d_%02d_%02d" % (
      time.localtime(time.time())[:6])
    f = open(capture_input_dir+"/"+target_module+"_"+time_stamp, "wb")
    pickle.dump([server_info, target_module, inp], f)
    f.close()

  print '[<a href="'+server_info.base()+'">Index of services</a>]'
  print '[<a href="'+server_info.file(target_module+".html")+'">New input</a>]'
  print "<hr>"

  import traceback
  class empty: pass
  status = empty()
  status.in_table = False
  try:
    target.run(server_info, inp, status)
  except Exception:
    if (status.in_table): print "</table><pre>"
    ei = sys.exc_info()
    print traceback.format_exception_only(ei[0], ei[1])[0]
    print
    print
    print "Details:"
    print
    traceback.print_exc()
  else:
    print "<hr>"
    print '[<a href="'+server_info.base()+'">Index of services</a>]'
    print \
      '[<a href="'+server_info.file(target_module+".html")+'">New input</a>]'

if (__name__ == "__main__"):
  run()
