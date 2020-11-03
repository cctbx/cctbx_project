
"""
Tools for handling plottable data, usually similar to CCP4's loggraph format
(which may be parsed and output by this module).
"""

from __future__ import absolute_import, division, print_function
from libtbx import adopt_init_args
from libtbx.utils import Sorry
import os.path
import math
import re
from six.moves import range
from six.moves import zip

class plot_data(object):
  def __init__(self,
               plot_title,
               x_label,
               y_label,
               x_data,
               y_data,
               y_legend,
               comments):
    self.plot_title = plot_title

    self.x_label = x_label
    self.y_label = y_label
    self.x_data = x_data
    self.comments = comments
    self.domain_flag = 'A'

    ## The x_data is an (flex) array of numbers
    ## The Y_data should be an list of (flex) arrays
    ## The legends should be an array of strings
    self.y_legend = []
    self.y_data = []
    if y_data is not None:
      assert ( len(y_data)==len(self.x_data) )
      self.y_data.append(y_data)
      self.y_legend.append(y_legend)

  def add_data(self,
               y_data,
               y_legend):
    assert ( len(y_data)==len(self.x_data) )
    self.y_data.append(y_data)
    self.y_legend.append(y_legend)


def plot_data_loggraph(plot_data,output):
  ## First we need to print the header information
  print(file=output)
  print(file=output)
  print('$TABLE: %s:'%(plot_data.plot_title), file=output)
  print('$GRAPHS', file=output)
  print(':%s' %(plot_data.comments), end=' ', file=output)
  index_string = ''
  for ii in range(len(plot_data.y_data)+1):
    index_string += '%d,'%(ii+1)
  print(':%s:%s:'%(plot_data.domain_flag,index_string[:-1]), file=output)
  print('$$', file=output)
  ## replace spaces for loggraph with underscores
  tmp_legend = plot_data.x_label
  spaces = 0
  spaces = tmp_legend.find(' ')
  if spaces>0:
    tmp_legend = tmp_legend.replace(' ','_')
  label_string = '%s'%(tmp_legend)
  for ii in range(len(plot_data.y_data)):
    tmp_legend = plot_data.y_legend[ii]
    ## loggraph does not like spaces in the legend names
    ## lets replace them with underscores
    spaces = 0
    spaces = tmp_legend.find(' ')
    if spaces>0:
      tmp_legend = tmp_legend.replace(' ','_')
    label_string += '   %s'%( tmp_legend )
  print('%s   $$ '%(label_string), file=output)
  print('$$', file=output)
  for ii in range(len(plot_data.x_data)):
    data_string = '%f'%(plot_data.x_data[ii])
    for jj in range(len(plot_data.y_data)):
      data_string +='   %f'%(plot_data.y_data[jj][ii])
    print('%s'%(data_string), file=output)
  print('$$', file=output)

#-----------------------------------------------------------------------
# Nat's utilities for plottable data
def flip_table(table):
  if table is None or len(table) == 0 :
    return []
  new_table = []
  for elem in table[0] :
    new_table.append([elem])
  if len(table) > 1 :
    for row in table[1:] :
      assert len(row) == len(new_table)
      for i, elem in enumerate(row):
        new_table[i].append(elem)
  return new_table

class table_data(object):
  """
  Container for tabular data.  Originally this was solely used as a container
  for CCP4 'loggraph' plots, but it can also be used to display tables that
  are not intended to be plotted.  A variety of output formats are supported,
  including export to JSON for eventual HTML display.  If graphs are defined,
  these objects can be passed to the loggraph frontend in wxtbx.plots.
  """
  def __init__(self,
      title,
      column_names=None,
      column_types=None,
      column_labels=None,
      column_formats=None,
      graph_names=None,
      graph_types=None,
      graph_labels=None,
      graph_columns=None,
      data=None,
      comments=None,
      x_is_inverse_d_min=False,
      first_two_columns_are_resolution=False,
      force_exact_x_labels=False,
      reference_marks=None):
    adopt_init_args(self, locals())
    if (reference_marks is not None):
      assert (len(reference_marks) == len(graph_columns) == 1)
    if (data is not None) : # check column size consistency
      lengths = { len(column) for column in data }
      assert len(lengths) == 1
    self._is_complete = False
    self._graphs = {}
    self._column_width = 10
    self.plot_type = "GRAPH"

  # Backwards compatibility
  def __setstate__(self, state):
    self.__dict__.update(state)
    if (not hasattr(self, "first_two_columns_are_resolution")):
      self.first_two_columns_are_resolution = None
    if (not hasattr(self, "force_exact_x_labels")):
      self.force_exact_x_labels = None
    if (not hasattr(self, "reference_marks")):
      self.reference_marks = None

  @property
  def n_rows(self):
    return len(self.data[0])

  @property
  def n_cols(self):
    return len(self.data)

  @property
  def n_graphs(self):
    return len(self.graph_names)

  def add_graph(self, name, type, columns):
    if self.graph_names is None :
      self.graph_names = [name.strip()]
    else :
      self.graph_names.append(name.strip())
    if self.graph_types is None :
      self.graph_types = [type]
    else :
      self.graph_types.append(type)
    if self.graph_columns is None :
      self.graph_columns = [columns]
    else :
      self.graph_columns.append(columns)

  def import_loggraph(self, lines):
    """
    Parse CCP4 loggraph format and populate internal data structures.
    Input may be provided as list of individual lines or as a string.
    """
    if isinstance(lines, list):
      lines = "\n".join(lines)

    blocks = [ b.strip() for b in lines.split('$$') ]

    # Loggraph format is defined by mandatory 4 blocks, separated by '$$', followed by nothing.
    # http://www.ccp4.ac.uk/html/loggraphformat.html
    assert len(blocks) == 5, 'input not in loggraph format (%d blocks found)' % len(blocks)
    header, columns, comment, data, remainder = blocks
    assert remainder == '', 'loggraph table has %d bytes following the table end' % len(remainder)

    if '$TABLE' in header:
      title = re.search('\\$TABLE\\s*:(.*?)(:|\n|$)', header, re.MULTILINE)
      if title:
        self.title = title.group(1).strip()

    graphs = re.search('\\$(GRAPHS|SCATTER)[\\s\n]*((:[^:\n]*:[^:\n]*:[^:\n]*(:|$)[\\s\n]*)*)($|\\$)', header, re.MULTILINE)
    if graphs:
      if graphs.group(1) == 'GRAPHS':
        self.plot_type = "GRAPH"
      elif graphs.group(1) == 'SCATTER':
        self.plot_type = "SCATTER"
      else:
        raise TypeError('Unknown graph type %s' % graphs.group(1))
      graphs = graphs.group(2)
      for graph in re.finditer(':([^:\n]*):([^:\n]*):([^:\n]*)(:|$)', graphs, re.MULTILINE):
        self.add_graph(name=graph.group(1),
                       type=graph.group(2),
                       columns=[ int(col)-1 for col in graph.group(3).split(',') ])

    # Newlines, spaces, tabs etc. are not allowed in column names.
    # Treat them as separators.
    columns = re.sub(r'\s+', ' ', columns).split(' ')
    data_width = len(columns)
    self.column_labels = [ lbl.replace('_', ' ') for lbl in columns ]
    if self.column_labels[0] in ("1/d^2", "1/d**2", "1/resol^2"):
      self.x_is_inverse_d_min = True

    self.data = [[] for x in range(data_width)]

    # Now load the data
    data = re.sub(r'\s+', ' ', data).split(' ')
    for entry, datum in enumerate(data):
      self.data[entry % data_width].append(_tolerant_float(datum))

    # Int-ify integer-like columns
    # It is unclear to me if this actually serves any purpose.
    for i, column in enumerate(self.data):
      column_is_ints = True
      for x in column:
        if (x is not None) and \
           (str(x).lower() not in ('nan', 'inf', '-inf')) and \
           (x != int(x)):
          column_is_ints = False
          break
      if column_is_ints:
        self.data[i] = [ int(x) if x and x == x else None
                         for x in column ]

  def add_row(self, row):
    '''Unclear if this is used from outside the class'''
    if self.data is None or len(self.data) == 0 :
      self.data = [ [x] for x in row ]
    else :
      assert len(self.data) == len(row), row
      for i, value in enumerate(row):
        self.data[i].append(value)

  def add_column(self, column, column_name=None, column_label=None):
    if self.data is None or len(self.data) == 0 :
      self.data = [ list(column) ]
    else :
      assert len(self.data[0]) == len(column)
      self.data.append(column)
    if column_name is not None :
      if self.column_names is None : self.column_names = [column_name]
      else :                         self.column_names.append(column_name)
    if column_label is not None :
      if self.column_labels is None : self.column_labels = [column_label]
      else :                          self.column_labels.append(column_label)

  def _max_column_width(self, precision):
    assert isinstance(precision, int)
    if self.column_labels is None : return precision
    label_widths = [ len(lab) for lab in self.column_labels ]
    cwidth = max(label_widths)
    if cwidth < precision :
      cwidth = precision
    return cwidth

  def get_value_as_resolution(self, x):
    if (self.x_is_inverse_d_min):
      return x ** -0.5
    return x

  def _format_num_row(self, row, precision=None, column_width=None):
    if row is None : return []
    if (self.column_formats is not None):
      if (self.first_two_columns_are_resolution):
        d_max = self.column_formats[0] % self.get_value_as_resolution(row[0])
        d_min = self.column_formats[1] % self.get_value_as_resolution(row[1])
        frow = [ "%s - %s" % (d_max, d_min) ]
        return frow + \
               [ (f % x) for f, x in zip(self.column_formats[2:], row[2:]) ]
      else :
        if (self.x_is_inverse_d_min):
          row = [ math.sqrt(1/row[0]) ] + row[1:]
        return [
          (f % x) if x is not None else "*" for f, x in zip(self.column_formats, row)
        ]
    else :
      f1 = "%-g"
      if (precision is not None):
        f1 = "%s-.%dg" % (r'%', precision)
      f2 = "%s"
      if (column_width is not None):
        f2 = "%s-%ds" % (r'%', column_width)
      return [ f2 % ftoa(x, f1) for x in row ]

  def _format_labels(self, labels, column_width):
    if labels is None : return []
    f1 = "%s-%ds" % (r'%', column_width)
    return [ f1 % lab for lab in labels ]

  def format_simple(self, precision=6, indent=0):
    data = self.data
    assert data is not None and len(data) > 0
    column_width = self._max_column_width(precision)
    column_headers = self._format_labels(self.column_labels, column_width)
    out = _formatting_buffer(indent)
    if self.title is not None :
      out += self.title
    trailing_spaces = re.compile(r"\ *$")
    labels = " ".join(self._format_labels(self.column_labels, column_width))
    out += trailing_spaces.sub("", labels)
    nrows = len(data[0])
    f1 = "%s-%ds" % (r'%', column_width)
    for j in range(nrows):
      row = [ col[j] for col in data ]
      frow = self._format_num_row(row, column_width, precision)
      frows = " ".join([ f1 % cv for cv in frow ])
      out += trailing_spaces.sub("", frows)
    return str(out)

  def format(self, precision=6, indent=0, equal_widths=True):
    """Formats the table for printout in a log file, with equal-width
    columns and cell boundaries, e.g.:
      -------------------------------
      | Title                       |
      |-----------------------------|
      | col1    | col1    | col3    |
      |-----------------------------|
      | 0.1     | 5       | *       |
      -------------------------------
    """
    data = self.data
    assert data is not None and len(data) > 0
    column_labels = self.column_labels
    # most of the code here deals with generic numerical values, but the
    # use of resolution ranges comes up often enough to merit special handling
    if (self.first_two_columns_are_resolution):
      column_labels = ["Res. range"] + column_labels[2:]
    formatted_rows = [ column_labels ]
    for j in range(self.n_rows):
      row = [ col[j] for col in data ]
      formatted_rows.append(self._format_num_row(row, precision=precision))
    column_widths = []
    for j in range(len(column_labels)):
      column_widths.append(max([ len(r[j]) for r in formatted_rows ]))
    if (equal_widths):
      # if the first column is the resolution range, we allow that to be
      # different in width than the remaining columns
      if (self.first_two_columns_are_resolution):
        max_w = max(column_widths[1:])
        column_widths = [column_widths[0]] + [ max_w ] * (len(column_widths)-1)
      else :
        max_w = max(column_widths)
        column_widths = [ max_w ] * len(column_widths)
    column_headers = []
    for cw, cl in zip(column_widths, column_labels):
      f1 = "%s-%ds" % (r'%', cw)
      column_headers.append(f1 % cl)
    column_headers = " | ".join(column_headers)
    f2 = "%s-%ds" % (r'%', len(column_headers))
    table_width = len(column_headers) + 4
    sep_line = "-" * table_width
    out = _formatting_buffer(indent)
    out += sep_line
    out += "| " + f2 % self.title + " |"
    out += "|" + sep_line[1:-1] + "|"
    out += "| " + column_headers + " |"
    out += "|" + sep_line[1:-1] + "|"
    for j in range(self.n_rows):
      frow = []
      for cw, cv in zip(column_widths, formatted_rows[j+1]):
        f1 = "%s-%ds" % (r'%', cw)
        frow.append(f1 % cv)
      frow = " | ".join(frow)
      out += "| " + frow + " |"
    out += sep_line
    return str(out)

  def format_loggraph(self, precision=6, column_width=None):
    """
    Create CCP4 loggraph format text for embedding in logfiles.
    """
    data = self.data
    graph_columns = self.graph_columns
    graph_names = self.graph_names
    graph_types = self.graph_types
    column_labels = self.column_labels
    assert data is not None and len(data) > 0
    assert column_width is None or isinstance(column_width, int)
    assert graph_types is None or len(graph_types) == len(graph_names)
    out = "$TABLE: %s\n" % self.title
    out += "$GRAPHS\n"
    if graph_types is None :
      graph_types = ["A" for graph in graph_names ]
    for i, graph_name in enumerate(graph_names):
      out += ":%s" % graph_name
      out += ":%s:%s:\n" % (graph_types[i],
                           ",".join([ "%d"%(x+1) for x in graph_columns[i] ]))
    out += "$$\n"
    re_spaces = re.compile(r"[\ ]{1,}")
    labels = [ re_spaces.sub("_", lab) for lab in column_labels ]
    if column_width is None :
      column_width = max(self._max_column_width(precision),
                         max([ len(lab) for lab in labels ]))
    f1 = "%s.%dg" % (r'%', precision)
    f2 = "%s-%ds" % (r'%', column_width)
    labels = [ re_spaces.sub("_", lab) for lab in column_labels ]
    out += "%s $$\n" % "  ".join([ f2 % lab for lab in labels ])
    out += "$$\n"
    trailing_spaces = re.compile(r"\ *$")
    for j in range(len(data[0])):
      row = [ col[j] for col in data ]
      frow = self._format_num_row(row, column_width, precision)
      frow = [ f2 % cv for cv in frow ]
      out += trailing_spaces.sub("", "  ".join(frow))
      out += "\n"
    out += "$$\n"
    return out

  def export_json(self, convert=True):
    import json
    graph_types = self.graph_types
    if graph_types is None :
      graph_types = ["A" for graph in self.graph_names ]
    export_dict = {
      "title" : self.title,
      "graph_columns" : self.graph_columns,
      "graph_names" : self.graph_names,
      "graph_types" : graph_types,
      "column_labels" : self.column_labels,
      "x_is_inverse_d_min" : self.x_is_inverse_d_min,
      "data" : self.data,
    }
    if (convert):
      return json.dumps(export_dict)
    return export_dict

  def export_rows(self, convert=True, precision=6):
    column_labels = self.column_labels
    # most of the code here deals with generic numerical values, but the
    # use of resolution ranges comes up often enough to merit special handling
    if (self.first_two_columns_are_resolution):
      column_labels = ["Res. range"] + column_labels[2:]
    formatted_rows = [ column_labels ]
    for j in range(self.n_rows):
      row = [ col[j] for col in self.data ]
      formatted_rows.append(self._format_num_row(row, precision=precision))
    return formatted_rows

  def export_json_table(self, convert=True):
    import json
    export_dict = {
      "title" : self.title,
      "rows" : self.export_rows(),
    }
    if (convert):
      return json.dumps(export_dict)
    return export_dict

  def as_rows(self):
    return flip_table(self.data)

  def only_plot(self):
    assert (len(self.graph_names) == 1)
    return self.graph_names[0]

  def get_reference_marks(self):
    if (self.reference_marks is not None):
      return self.reference_marks[0]
    return None

  def get_graph(self, graph_name=None, column_list=[]):
    graph_names = self.graph_names
    column_labels = self.column_labels
    graph_columns = self.graph_columns
    data = self.data
    _graphs = self._graphs
    if graph_name is not None :
      if not graph_name in _graphs :
        if not graph_name in graph_names :
          return None
        n = graph_names.index(graph_name)
        gdata = self._extract_data_column(graph_columns[n])
        if len(column_labels) == len(data):
          labels = [column_labels[i] for i in graph_columns[n]]
        else :
          labels = []
        x_axis = None
        y_axis = None
        if self.graph_labels is not None :
          (x_axis, y_axis) = self.graph_labels[n]
        _graphs[graph_name] = graph_data(graph_name, gdata, "plot", labels,
          x_axis, y_axis)
      return _graphs[graph_name]
    elif len(column_list) > 1 :
      if not column_list in self._graphs :
        data = None
        if isinstance(column_list[0], int):
          data = self._extract_data_column(column_list)
          labels = [ column_labels[i] for i in column_list ]
        elif isinstance(column_list[0], str):
          n_list = []
          for col in column_list :
            n_list.append(column_labels.index(col))
          gdata = self._extract_data_column(n_list)
          labels = [ column_labels[i] for i in n_list ]
        if gdata is None :
          return None
        _graphs[column_list] = graph_data(None, gdata, "plot", labels)
      return _graphs[column_list]

  def get_column_by_label(self, column_label):
    if not column_label in self.column_labels :
      raise RuntimeError(
        "Couldn't find column %s in this table.  (Valid columns: %s)" %
        (column_label, ",".join(self.column_labels)))
    i = self.column_labels.index(column_label)
    return self.data[i]

  def _extract_data_column(self, column_list):
    assert len(column_list) <= len(self.data)
    data = self.data
    new_data = [ [ x for x in data[i] ] for i in column_list ]
    return new_data

  def __str__(self):
    return self.format_simple()

  def get_x_values(self):
    return self.data[0]

  def get_x_as_resolution(self):
    assert self.x_is_inverse_d_min
    oldx = self.data[0]
    newx = []
    for x in oldx :
      newx.append(math.sqrt(1.0/x))
    return newx

class graph_data(object):
  def __init__(self, name, data, type="plot", data_labels=None, x_axis=None,
      y_axis=None):
    self.name = name
    self.data = data
    self.type = type
    if data_labels is None or len(data_labels) == 0 :
      self.x_label = "X"
      self.y_labels = [ "Y" for i in range(1, len(data)) ]
    else :
      self.x_label = data_labels[0]
      self.y_labels = [ data_labels[i] for i in range(1, len(data)) ]
    self.x_axis_label = x_axis
    self.y_axis_label = y_axis

  def get_plots(self, fill_in_missing_y=None):
    plots = []
    data = self.data
    for i in range(1, len(data)):
      plot_x = []
      plot_y = []
      for j in range(0, len(data[i])):
        if data[0][j] is not None :
          if data[i][j] is not None :
            plot_x.append(data[0][j])
            plot_y.append(data[i][j])
          elif fill_in_missing_y is not None :
            plot_x.append(data[0][j])
            plot_y.append(fill_in_missing_y)
      plots.append((plot_x, plot_y))
    return plots

class histogram_data(object):
  pass

def _tolerant_float(string):
  try:
    return float(string)
  except ValueError:
    return None

# backwards-atof
def ftoa(val, format_string='%.6g'):
  if val is None :
    return '*'
  else :
    return format_string % val

class _formatting_buffer(object):
  def __init__(self, indent=0):
    self._initial_space = " " * indent
    self._buffer = []

  def write(self, strdata):
    self._buffer.append(self._initial_space + strdata)

  def append(self, strdata):
    self.write(strdata)

  def __add__(self, strdata):
    self.write(strdata)
    return self

  def __str__(self):
    out = "\n".join(self._buffer) + "\n"
    return out

def import_ccp4i_logfile(file_name=None, log_lines=None):
  assert file_name is not None or log_lines is not None
  if not log_lines :
    with open(file_name) as f:
      log_lines = f.readlines()
  current_lines = None
  tables_raw = []
  sections_read = 0
  for line in log_lines :
    line = line.strip()
    if re.match(r"\$TABLE\s*:", line):
      current_lines = [ line ]
      sections_read = 0
    elif line[-2:] == "$$" and current_lines is not None :
      current_lines.append(line)
      sections_read += line.count("$$")
      if sections_read == 4 :
        tables_raw.append(current_lines)
        current_lines = None
    elif sections_read < 4 and current_lines is not None :
      current_lines.append(line)
  tables = []
  for loggraph in tables_raw :
    t = table_data(None)
    t.import_loggraph(loggraph)
    tables.append(t)
  return tables

class simple_matplotlib_plot(object):
  """
  Class for writing a Matplotlib plot to a static image file without a GUI.
  This should be subclassed and combined with whatever mixin is used to
  actually responsible for the plotting.
  """
  def __init__(self,
                figure_size=(8,8),
                font_size=12,
                title_font_size=12,
                facecolor='white',
                transparent=False):
    adopt_init_args(self, locals())
    try :
      import matplotlib
      import matplotlib.figure
      from matplotlib.backends.backend_agg import FigureCanvasAgg
    except ImportError as e :
      print(e)
      raise Sorry("Plotting requires that matplotlib be installed.")
    self.figure = matplotlib.figure.Figure(figure_size, 72, linewidth=0,
      facecolor=facecolor)
    if transparent :
      self.figure.figurePatch.set_alpha(0.0)
    self.canvas = FigureCanvasAgg(self.figure)
    #self.canvas.toolbar = oop.null()
    #self.figmgr = FigureManager(self.canvas, 1, self)

  def save_image(self, file_name, dpi):
    assert (file_name is not None) and (file_name != "")
    base, ext = os.path.splitext(file_name)
    if (ext == ".pdf"):
      self.figure.savefig(file_name, orientation="landscape", format="pdf")
    elif (ext == ".ps"):
      self.figure.savefig(file_name, orientation="landscape", format="ps")
    elif (ext == ".png"):
      self.figure.savefig(file_name, format="png", dpi=dpi)
    else :
      raise RuntimeError("Extension %s not supported" % s)

#---end
