from cctbx.array_family import flex
from libtbx import adopt_init_args
import string, re, math

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
  print >> output
  print >> output
  print >> output, '$TABLE: %s:'%(plot_data.plot_title)
  print >> output, '$GRAPHS'
  print >> output, ':%s' %(plot_data.comments),
  index_string = ''
  for ii in range(len(plot_data.y_data)+1):
    index_string += '%d,'%(ii+1)
  print >> output, ':%s:%s:'%(plot_data.domain_flag,index_string[:-1])
  print >> output, '$$'
  ## replace spaces for loggraph with underscores
  tmp_legend = plot_data.x_label
  spaces = 0
  spaces = string.find(tmp_legend,' ')
  if spaces>0:
    tmp_legend = tmp_legend.replace(' ','_')
  label_string = '%s'%(tmp_legend)
  for ii in range(len(plot_data.y_data)):
    tmp_legend = plot_data.y_legend[ii]
    ## loggraph does not like spaces in the legend names
    ## lets replace them with underscores
    spaces = 0
    spaces = string.find(tmp_legend,' ')
    if spaces>0:
      tmp_legend = tmp_legend.replace(' ','_')
    label_string += '   %s'%( tmp_legend )
  print >> output, '%s   $$ '%(label_string)
  print >> output, '$$'
  for ii in range(len(plot_data.x_data)):
    data_string = '%f'%(plot_data.x_data[ii])
    for jj in range(len(plot_data.y_data)):
      data_string +='   %f'%(plot_data.y_data[jj][ii])
    print >> output, '%s'%(data_string)
  print >> output, '$$'

#-----------------------------------------------------------------------
# Nat's utilities for plottable data
def flip_table (table) :
  if table is None or len(table) == 0 :
    return []
  new_table = []
  for elem in table[0] :
    new_table.append([elem])
  if len(table) > 1 :
    for row in table[1:] :
      assert len(row) == len(new_table)
      for i, elem in enumerate(row) :
        new_table[i].append(elem)
  return new_table

class table_data (object) :
  def __init__ (self,
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
      x_is_inverse_d_min=False) :
    adopt_init_args(self, locals())
    self._is_complete = False
    self._graphs = {}
    self._column_width = 10
    self.plot_type = "GRAPH"

  def add_graph (self, name, type, columns) :
    if self.graph_names is None :
      self.graph_names = [name]
    else :
      self.graph_names.append(name)
    if self.graph_types is None :
      self.graph_types = [type]
    else :
      self.graph_types.append(type)
    if self.graph_columns is None :
      self.graph_columns = [columns]
    else :
      self.graph_columns.append(columns)

  def import_loggraph (self, loggraph_lines) :
    sections_passed = 0
    initial_spaces = re.compile("^\s*")
    trailing_spaces = re.compile("\s*$")
    trailing_dollars = re.compile("\$\$\.*")
    trailing_dollars_phaser = re.compile("\$\$ loggraph \$\$.*")
    graph_lines = None
    if isinstance(loggraph_lines, str) :
      lines = loggraph_lines.split("\n")
    else :
      lines = loggraph_lines
    for raw_line in lines :
      line = raw_line.strip()
      #line = initial_spaces.sub("", raw_line)
      if line.startswith("$$") :
        sections_passed += 1 #line.count("$$")
        line = re.sub("^\$\$\ *", "", line)
      if line == "" :
        pass
      elif line.startswith("$TABLE") :
        self.title = initial_spaces.sub("", line.split(":")[1])
        #re.sub("\ *:\ *$", "", line[8:]))
      elif line.startswith("$GRAPHS") :
        graph_lines = line[8:]
        self.plot_type = "GRAPH"
      elif line.startswith("$SCATTER") :
        graph_lines = line[9:]
        self.plot_type = "SCATTER"
      elif graph_lines is not None and sections_passed == 0 :
        graph_lines += line
      elif sections_passed == 1 :
        clean_line = trailing_dollars.sub("",
          trailing_dollars_phaser.sub("", line))
        if self.column_labels is None :
          column_labels = clean_line.split()
          self.column_labels = [ re.sub("_"," ",lbl) for lbl in column_labels ]
      elif sections_passed == 3 and line[0:2] != "$$" :
        fields = [ _atof(x) for x in line.split() ]
        self.add_row(fields)
      if line.endswith("$$") : #trailing_spaces.sub("", line).endswith("$$") :
        sections_passed += line.count("$$")
        if sections_passed == 1 :
          if graph_lines is None : graph_lines = line[:-2]
          else                   : graph_lines += line[:-2]
          graph_string = re.sub(":\s*$", "", re.sub("^\s*:", "", graph_lines))
          fields = graph_string.split(":")
          i = 0
          while i < len(fields) :
            col_strs = fields[i+2].split(",")
            self.add_graph(name=fields[i],
                           type=fields[i+1],
                           columns=[ (_atoi(x_str)-1) for x_str in col_strs ])
            i += 4
        elif (sections_passed == 4) :
          break
    for i, column in enumerate(self.data) :
      column_is_ints = [ x is None or int(x)==x for x in column ]
      if not False in column_is_ints :
        newcol = []
        for x in column :
          if x is None : newcol.append(x)
          else :         newcol.append(int(x))
        self.data[i] = newcol
    if self.column_labels[0] in ["1/d^2","1/d**2","1/resol^2"] :
      self.x_is_inverse_d_min = True

  def add_row (self, row) :
    if self.data is None or len(self.data) == 0 :
      self.data = [ [x] for x in row ]
    else :
      assert len(self.data) == len(row)
      for i, value in enumerate(row) :
        self.data[i].append(value)

  def add_column (self, column, column_name=None, column_label=None) :
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

  def _max_column_width (self, precision) :
    assert isinstance(precision, int)
    if self.column_labels is None : return precision
    label_widths = [ len(lab) for lab in self.column_labels ]
    cwidth = max(label_widths)
    if cwidth < precision :
      cwidth = precision
    return cwidth

  def _format_num_row (self, row, column_width, precision) :
    if row is None : return []
    f1 = "%s-%dg" % (r'%', precision)
    f2 = "%s-%ds" % (r'%', column_width)
    return [ f2 % ftoa(x, f1) for x in row ]

  def _format_labels (self, labels, column_width) :
    if labels is None : return []
    f1 = "%s-%ds" % (r'%', column_width)
    return [ f1 % lab for lab in labels ]

  def format_simple (self, precision=6, indent=0) :
    data = self.data
    assert data is not None and len(data) > 0
    column_width = self._max_column_width(precision)
    column_headers = self._format_labels(self.column_labels, column_width)
    out = _formatting_buffer(indent)
    if self.title is not None :
      out += self.title
    trailing_spaces = re.compile("\ *$")
    labels = " ".join(self._format_labels(self.column_labels, column_width))
    out += trailing_spaces.sub("", labels)
    nrows = len(data[0])
    for j in xrange(nrows) :
      row = [ col[j] for col in data ]
      frows = " ".join(self._format_num_row(row, column_width, precision))
      out += trailing_spaces.sub("", frows)
    return str(out)

  def format (self, precision=6, indent=0) :
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
    column_width = self._max_column_width(precision)
    if column_width > precision :
      precision = column_width
    f1 = "%s-%ds" % (r'%', column_width)
    column_headers = " | ".join(self._format_labels(self.column_labels,
      column_width))
    f2 = "%s-%ds" % (r'%', len(column_headers))
    table_width = len(column_headers) + 4
    sep_line = "-" * table_width
    out = _formatting_buffer(indent)
    out += sep_line
    out += "| " + f2 % self.title + " |"
    out += "|" + sep_line[1:-1] + "|"
    out += "| " + column_headers + " |"
    out += "|" + sep_line[1:-1] + "|"
    for j in xrange(len(data[0])) :
      row = [ col[j] for col in data ]
      frow = " | ".join(self._format_num_row(row, column_width, precision))
      out += "| " + frow + " |"
    out += sep_line
    return str(out)

  def format_html (self) :
    pass

  def format_loggraph (self, precision=6, column_width=None) :
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
    for i, graph_name in enumerate(graph_names) :
      out += ":%s\n" % graph_name
      out += ":%s:%s:\n" % (graph_types[i],
                           ",".join([ "%d"%(x+1) for x in graph_columns[i] ]))
    out += "$$\n"
    if column_width is None :
      column_width = self._max_column_width(precision)
    f1 = "%s.%dg" % (r'%', precision)
    f2 = "%s-%ds" % (r'%', column_width)
    re_spaces = re.compile("[\ ]{1,}")
    labels = [ re_spaces.sub("_", lab) for lab in column_labels ]
    out += "%s $$\n" % "  ".join([ f2 % lab for lab in labels ])
    out += "$$\n"
    trailing_spaces = re.compile("\ *$")
    for j in xrange(len(data[0])) :
      row = [ col[j] for col in data ]
      frow = "  ".join(self._format_num_row(row, column_width, precision))
      out += trailing_spaces.sub("", frow)
      out += "\n"
    out += "$$\n"
    return out

  def as_rows (self) :
    return flip_table(self.data)

  def get_graph (self, graph_name=None, column_list=[]) :
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
        if len(column_labels) == len(data) :
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
        if isinstance(column_list[0], int) :
          data = self._extract_data_column(column_list)
          labels = [ column_labels[i] for i in column_list ]
        elif isinstance(column_list[0], str) :
          n_list = []
          for col in column_list :
            n_list.append(column_labels.index(col))
          gdata = self._extract_data_column(n_list)
          labels = [ column_labels[i] for i in n_list ]
        if gdata is None :
          return None
        _graphs[column_list] = graph_data(None, gdata, "plot", labels)
      return _graphs[column_list]

  def get_column_by_label (self, column_label) :
    if not column_label in self.column_labels :
      raise RuntimeError(
        "Couldn't find column %s in this table.  (Valid columns: %s)" %
        (column_label, ",".join(self.column_labels)))
    i = self.column_labels.index(column_label)
    return self.data[i]

  def _extract_data_column (self, column_list) :
    assert len(column_list) <= len(self.data)
    data = self.data
    new_data = [ [ x for x in data[i] ] for i in column_list ]
    return new_data

  def __str__ (self) :
    return self.format_simple()

  def get_x_as_resolution (self) :
    assert self.x_is_inverse_d_min
    oldx = self.data[0]
    newx = []
    for x in oldx :
      newx.append(math.sqrt(1.0/x))
    return newx

class graph_data (object) :
  def __init__ (self, name, data, type="plot", data_labels=None, x_axis=None,
      y_axis=None) :
    self.name = name
    self.data = data
    self.type = type
    if data_labels is None or len(data_labels) == 0 :
      self.x_label = "X"
      self.y_labels = [ "Y" for i in xrange(1, len(data)) ]
    else :
      self.x_label = data_labels[0]
      self.y_labels = [ data_labels[i] for i in xrange(1, len(data)) ]
    self.x_axis_label = x_axis
    self.y_axis_label = y_axis

  def get_plots (self, fill_in_missing_y=None) :
    plots = []
    data = self.data
    for i in xrange(1, len(data)) :
      plot_x = []
      plot_y = []
      for j in xrange(0, len(data[i])) :
        if data[0][j] is not None :
          if data[i][j] is not None :
            plot_x.append(data[0][j])
            plot_y.append(data[i][j])
          elif fill_in_missing_y is not None :
            plot_x.append(data[0][j])
            plot_y.append(fill_in_missing_y)
      plots.append((plot_x, plot_y))
    return plots

class histogram_data (object) :
  pass

def _atof (fstring) :
  try :
    val = string.atof(fstring)
  except :
    val = None
  return val

def _atoi (istring) :
  try :
    val = string.atoi(istring)
  except :
    val = None
  return val

# backwards-atof
def ftoa (val, format_string='%.6g') :
  if val is None :
    return '*'
  else :
    return format_string % val

class _formatting_buffer (object) :
  def __init__ (self, indent=0) :
    self._initial_space = " " * indent
    self._buffer = []

  def write (self, strdata) :
    self._buffer.append(self._initial_space + strdata)

  def append (self, strdata) :
    self.write(strdata)

  def __add__ (self, strdata) :
    self.write(strdata)
    return self

  def __str__ (self) :
    out = "\n".join(self._buffer) + "\n"
    return out

def import_ccp4i_logfile (file_name=None, log_lines=None) :
  assert file_name is not None or log_lines is not None
  if not log_lines :
    log_lines = open(file_name).readlines()
  current_lines = None
  tables_raw = []
  sections_read = 0
  for line in log_lines :
    line = line.strip()
    if re.match("\$TABLE\s*:", line) :
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

#---end
