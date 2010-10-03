from cctbx.web.asu_gallery import web_links
from cctbx.web.asu_gallery import html_head_title
from cctbx import sgtbx
import sys

class symbol_table(object):

  def __init__(self, point_group_type):
    self.point_group_type = point_group_type
    self.table = []

  def add(self, space_group_symbols):
    self.table.append(space_group_symbols)

  def format_html(self, f, n_columns):
    print >> f, "<table border=2 cellpadding=2>"
    n_symbols = len(self.table)
    n_rows = n_symbols / n_columns
    if (n_rows * n_columns < n_symbols): n_rows += 1
    for i_row in xrange(n_rows):
      print >> f, "<tr>"
      for i_column in xrange(n_columns):
        i = i_column * n_rows + i_row
        if (i < len(self.table)):
          symbols = self.table[i]
          print >> f, '<td><a href="asu_%03d.html">%s (%d)</a></td>' % (
            symbols.number(),
            symbols.hermann_mauguin().replace(" ", ""),
            symbols.number())
      print >> f, "</tr>"
    print >> f, "</table>"
    print >> f, "<p>"

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
Reference:
<a href="%(iucrcompcomm_jul2003)s"
>IUCr Computing Commission Newsletter No. 2, July 2003</a>
<hr>""" % vars()
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
