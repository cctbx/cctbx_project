import sys
sys.stderr = sys.stdout

import os, cgi, urlparse

if (0):
  print "Content-type: text/plain"
  print
  import os
  for k,v in os.environ.items():
    print k, v

cctbx_url = list(urlparse.urlsplit(os.environ["HTTP_REFERER"]))
cctbx_url[2] = "/".join(cctbx_url[2].split("/")[:-1])

form = cgi.FieldStorage()
target_module = form["target_module"].value
exec "import " + target_module + " as target"
target.run(cctbx_url, form)
