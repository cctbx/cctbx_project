import sys
sys.stderr = sys.stdout

import cgitb
cgitb.enable()

from cctbx.web import utils
import os, cgi, urlparse

if (0):
  cgi.test()
  sys.exit(0)

print "Content-type: text/html"
print

server_info = utils.server_info()

form = cgi.FieldStorage()
target_module = form["target_module"].value

exec "import " + target_module + " as target"
inp = target.interpret_form_data(form)

if (1):
  # capture input to facilitate debugging
  import time, pickle
  time_stamp = "%d_%02d_%02d_%02d_%02d_%02d" % time.localtime(time.time())[:6]
  f = open("/var/tmp/cctbx_web/"+target_module+"_"+time_stamp, "wb")
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
except:
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
  print '[<a href="'+server_info.file(target_module+".html")+'">New input</a>]'
