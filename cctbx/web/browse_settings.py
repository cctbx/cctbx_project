# This example uses an internal table, with 530 space group settings,
# that is based on Table 4.3.1 in the International Tables for
# Crystallography, Volume A (1983). Via the web interface the user
# specifies a space group symbol. This script determines the space
# group number corresponding to the given symbol, and then lists all
# tabulated settings for that space group number. If no space group
# symbol is given, all 530 entries in the internal table are listed.

from cctbx import sgtbx
import urllib
from cctbx.web import cgi_utils

def interpret_form_data(form):
  inp = cgi_utils.inp_from_form(form,
    (("sgsymbol", ""),
     ("convention", "")))
  return inp

def run(server_info, inp, status):
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
      query = "target_module=explore_symmetry&sgsymbol=" \
            + urllib.quote_plus(symbols.universal_hermann_mauguin())
      print ("<td><a href=\"%s\">%s</a>") % (
        server_info.script(query),
        symbols.universal_hermann_mauguin())
      print "<td>%s" % (symbols.hall(),)
      n_settings += 1
  print "</table>"
  if (sg_number == 0):
    print "<p>"
    print "Number of settings listed:", n_settings
  print "<p>"
