# This example uses an internal table, with 530 space group settings,
# that is based on Table 4.3.1 in the International Tables for
# Crystallography, Volume A (1983). Via the web interface the user
# specifies a space group symbol. This script determines the space
# group number corresponding to the given symbol, and then lists all
# tabulated settings for that space group number. If no space group
# symbol is given, all 530 entries in the internal table are listed.

import traceback
import sys, urllib, urlparse

from cctbx import sgtbx

class empty: pass

def interpret_form_data(form):
  inp = empty()
  for key in (("sgsymbol", ""),
              ("convention", "")):
    if (form.has_key(key[0])):
      inp.__dict__[key[0]] = form[key[0]].value.strip()
    else:
      inp.__dict__[key[0]] = key[1]
  return inp

def run(cctbx_url, inp):
  print "Content-type: text/html"
  print
  print "<td><a href=\"%s\">%s</a>" % (
    urlparse.urlunsplit(cctbx_url),
    urlparse.urlunsplit(cctbx_url))
  print "<p>"
  url_cctbx_web = urlparse.urlunsplit(
      cctbx_url[:2]
    + [cctbx_url[2] + "/cctbx_web.cgi"]
    + cctbx_url[3:])
  try:
    sg_number = 0
    if (len(inp.sgsymbol.strip()) != 0):
      sg_number = sgtbx.space_group_info(
        symbol=inp.sgsymbol,
        table_id=inp.convention).type().number()
    n_settings = 0
    print "<table border=2 cellpadding=2>"
    print "<tr>"
    print "<th>Space group<br>No."
    print "<th>Schoenflies<br>symbol"
    print "<th>Hermann-Mauguin<br>symbol"
    print "<th>Hall<br>symbol"
    for symbols in sgtbx.space_group_symbol_iterator():
      if (sg_number == 0 or symbols.number() == sg_number):
        print "<tr>"
        print "<td>(%d)<td>%s" % (
          symbols.number(), symbols.schoenflies())
        print ("<td><a href=\"%s?target_module=explore_symmetry"
              +"&sgsymbol=%s\">%s</a>") % (
          url_cctbx_web,
          urllib.quote_plus(symbols.extended_hermann_mauguin()),
          symbols.extended_hermann_mauguin())
        print "<td>%s" % (symbols.hall(),)
        n_settings += 1
    print "</table>"
    if (sg_number == 0):
      print "<p>"
      print "Number of settings listed:", n_settings
    print "<p>"
    print "<td><a href=\"%s\">%s</a>" % (
      urlparse.urlunsplit(cctbx_url),
      urlparse.urlunsplit(cctbx_url))

  except RuntimeError, e:
    print e
  except AssertionError:
    ei = sys.exc_info()
    print traceback.format_exception_only(ei[0], ei[1])[0]
