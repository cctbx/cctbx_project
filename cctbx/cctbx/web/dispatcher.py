import sys
sys.stderr = sys.stdout

import cgitb
cgitb.enable()

import os, cgi, urlparse

if (0):
  cgi.test()

if (0):
  print "Content-type: text/html"
  print

cctbx_url = list(urlparse.urlsplit(os.environ["HTTP_REFERER"]))
cctbx_url[2] = "/".join(cctbx_url[2].split("/")[:-1])

form = cgi.FieldStorage()
target_module = form["target_module"].value

exec "import " + target_module + " as target"
inp = target.interpret_form_data(form)

if (1):
  # capture input to facilitate debugging
  import time, pickle
  time_stamp = "%d_%02d_%02d_%02d_%02d_%02d" % time.localtime(time.time())[:6]
  f = open("/var/tmp/cctbx_web/"+target_module+"_"+time_stamp, "wb")
  pickle.dump([cctbx_url, target_module, inp], f)
  f.close()

target.run(cctbx_url, inp)
