from __future__ import absolute_import, division, print_function
from cctbx.web import cgi_utils
import pydoc
import cgi
from six.moves import cStringIO as StringIO
import sys

def interpret_form_data(form):
  inp = cgi_utils.inp_from_form(form,
     (("query", ""),))
  return inp

def run(server_info, inp, status):
  print("<pre>")
  sys.argv = ["libtbx.help"] + inp.query.split()
  s = StringIO()
  sys.stdout = s
  pydoc.cli()
  sys.stdout = sys.__stdout__
  s = s.getvalue()
  if sys.version_info.major >= 3:
    import html
    sys.stdout.write(html.escape(s))
  else:
    sys.stdout.write(cgi.escape(s))
  print("</pre>")
