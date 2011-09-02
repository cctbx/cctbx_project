from cctbx.web.asu_gallery import web_links
from cctbx.web.asu_gallery import html_head_title
from cctbx import sgtbx
import sys

def table_format_html(table, f, n_columns, serial_fmt="%03d"):
  print >> f, "<table border=2 cellpadding=2>"
  n_symbols = len(table)
  n_rows = n_symbols // n_columns
  if (n_rows * n_columns < n_symbols): n_rows += 1
  for i_row in xrange(n_rows):
    print >> f, "<tr>"
    for i_column in xrange(n_columns):
      i = i_column * n_rows + i_row
      if (i < len(table)):
        symbols = table[i]
        print >> f, (
          '<td><a href="asu_'+serial_fmt+'.html">%s (%d)</a></td>') % (
            symbols.number(),
            symbols.hermann_mauguin().replace(" ", ""),
            symbols.number())
    print >> f, "</tr>"
  print >> f, "</table>"
  print >> f, "<p>"

class point_group_symbols(object):

  def __init__(O, number, hermann_mauguin):
    O._number = number
    O._hermann_mauguin = hermann_mauguin

  def number(O): return O._number

  def hermann_mauguin(O): return O._hermann_mauguin

class plane_group_table(object):

  def __init__(O):
    O.table = []
    from cctbx.sgtbx import plane_groups
    for i,(hm,_) in enumerate(plane_groups.hermann_mauguin_hall_table):
      O.table.append(point_group_symbols(i+1, hm.replace("_"," ")))

  def format_html(O, f, n_columns):
    print >> f, "<h3>Plane groups</h3>"
    table_format_html(O.table, f, n_columns, serial_fmt="%02d")

class symbol_table(object):

  def __init__(self, point_group_type):
    self.point_group_type = point_group_type
    self.table = []

  def add(self, group_symbols):
    self.table.append(group_symbols)

  def format_html(self, f, n_columns):
    table_format_html(self.table, f, n_columns)

class point_group_table(object):

  def __init__(self, crystal_system):
    self.crystal_system = crystal_system
    self.symbol_table = []

  def add_point_group_type(self, point_group_type):
    self.current_symbol_table = symbol_table(point_group_type)
    self.symbol_table.append(self.current_symbol_table)

  def add(self, space_group_symbols):
    self.current_symbol_table.add(space_group_symbols)

  def format_html(self, f, n_columns):
    print >> f, "<h3>%s</h3>" % self.crystal_system
    for symbols in self.symbol_table:
      symbols.format_html(f, n_columns)

class crystal_system_table(object):

  def __init__(self):
    self.point_group_tables = []
    previous_crystal_system = None
    previous_point_group_type = None
    for space_group_number in xrange(1,231):
      space_group_symbols = sgtbx.space_group_symbols(space_group_number)
      space_group = sgtbx.space_group(space_group_symbols)
      crystal_system = space_group.crystal_system()
      point_group_type = space_group.point_group_type()
      if (point_group_type != previous_point_group_type):
        if (crystal_system != previous_crystal_system):
          current_crystal_system = point_group_table(crystal_system)
          self.point_group_tables.append(current_crystal_system)
        current_crystal_system.add_point_group_type(point_group_type)
      current_crystal_system.add(space_group_symbols)
      previous_point_group_type = point_group_type
      previous_crystal_system = crystal_system

  def format_html(self, f=None, n_columns=6):
    if (f is None): f = sys.stdout
    title = "Gallery of direct-space asymmetric units"
    iucrcompcomm_jul2003 = web_links.iucrcompcomm_jul2003
    print >> f, html_head_title(title=title)
    print >> f, """\
<body>
<hr>
<h2>%(title)s</h2>
<hr>
References:
<ul>
<li><a href="http://scripts.iucr.org/cgi-bin/paper?pz5088" target="external"
    >Acta Cryst. (2011). A67, 269-275</a>
<p>
<li><a href="%(iucrcompcomm_jul2003)s" target="external"
    >IUCr Computing Commission Newsletter No. 2, July 2003</a>
</ul>
<hr>""" % vars()
    plane_group_table().format_html(f,n_columns)
    for point_group in self.point_group_tables:
      point_group.format_html(f, n_columns)
    print >> f, """\
<hr>
<a href="http://cctbx.sourceforge.net/">[cctbx home]</a>
</body>
</html>"""

def write_html(f=None, n_columns=6):
  cs_table = crystal_system_table()
  cs_table.format_html(f, n_columns)

if (__name__ == "__main__"):
  write_html()
