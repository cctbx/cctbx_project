from cctbx.array_family import flex
from cStringIO import StringIO
import sys
import os
import string

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
  print >> output, ':%s ' %(plot_data.comments)
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
class table_data (object) :
  def __init__ (self,
      title,
      column_names=None,
      column_types=None,
      column_labels=None,
      graph_names=None,
      graph_columns=None,
      data=None) :
    self.title = title
    self._is_complete = False
    self.column_names = column_names
    self.column_types = column_types
    self.column_labels = column_labels
    self.graph_names = graph_names
    self.graph_types = []
    self.graph_columns = graph_columns
    self.data = data
    self._graphs = {}

  def add_graph (self, name, type, columns) :
    if self.graph_names is None :
      self.graph_names = [name]
    else :
      self.graph_names.append(name)
    self.graph_types.append(type)
    if self.graph_columns is None :
      self.graph_columns = [columns]
    else :
      self.graph_columns.append(columns)

  def import_loggraph (self, loggraph_lines) :
    import re
    from string import atof, atoi
    sections_passed = 0
    initial_spaces = re.compile("^\s*")
    trailing_spaces = re.compile("\s*$")
    trailing_dollars = re.compile("\$\$\ *$")
    graph_lines = None
    for raw_line in loggraph_lines :
      line = initial_spaces.sub("", raw_line)
      if line == "" :
        pass
      elif line[0:7] == "$TABLE:" :
        self.title = initial_spaces.sub("", re.sub("\ *:\ *$", "", line[8:-1]))
      elif line[0:7] == "$GRAPHS" :
        graph_lines = line[8:]
      elif graph_lines is not None and sections_passed == 0 :
        graph_lines += line
      elif sections_passed == 1 :
        clean_line = trailing_dollars.sub("", line)
        column_labels = clean_line.split()
        self.column_labels = [ re.sub("_", " ", lbl) for lbl in column_labels ]
      elif sections_passed == 3 and line[0:2] != "$$" :
        fields = [ atof(x) for x in line.split() ]
        self.add_row(fields)
      if trailing_spaces.sub("", line)[-2:] == "$$" :
        sections_passed += 1
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
                           columns=[ atoi(x_str) for x_str in col_strs ])
            i += 4
        elif sections_passed == 4 :
          break
    self.column_types = [ flex.double for column in self.data ]

  def add_row (self, row) :
    if self.data is None or len(self.data) == 0 :
      self.data = [ [x] for x in row ]
    else :
      assert len(self.data) == len(row)
      for i, value in enumerate(row) :
        self.data[i].append(value)

  def add_column (self, column) :
    self.add_column(column)

  def add_column (self, column) :
    if self.data is None or len(self.data) == 0 :
      self.data = [ column ]
    else :
      assert len(self.data[0]) == len(column)
      self.data.append(column)

  def formatted_simple (self) :
    pass

  def formatted (self) :
    pass

  def formatted_html (self) :
    pass

  def formatted_loggraph (self) :
    pass

  def get_graph (self, graph_name=None, column_list=[]) :
    if graph_name is not None :
      if not graph_name in self._graphs :
        if not graph_name in self.graph_names :
          return None
        n = self.graph_names.index(graph_name)
        data = self._extract_data_column(self.graph_column[n])
        if len(self.column_names) == len(self.data) :
          labels = [self.column_names[i] for i in self.graph_column[n]]
        else :
          labels = []
        self._graphs[graph_name] = graph_data(graph_name, data, "plot", labels)
      return self._graphs[graph_name]
    elif len(column_list) > 1 :
      if not column_list in self._graphs :
        data = None
        if isinstance(column_list[0], int) :
          data = self._extract_data_column(column_list)
          labels = [ self.column_labels[i] for i in column_list ]
        elif isinstance(column_list[0], str) :
          n_list = []
          for col in column_list :
            n_list.append(self.column_names.index(col))
          data = self._extract_data_column(n_list)
          labels = [ self.column_labels[i] for i in n_list ]
        if data is None :
          return None
        self._graphs[column_list] = graph_data(None, data, "plot", labels)
      return self._graphs[column_list]

  def _extract_data_column (self, column_list) :
    assert len(self.data) > 0
    assert len(column_list) <= len(self.data[0])
    data = self.data
    new_data = [ [ x for x in data[i] ] for i in column_list ]
    return new_data

  def __str__ (self) :
    return str(self.data)

class graph_data (object) :
  def __init__ (self, name, data, type="plot", data_labels=None) :
    self.name = name
    self.data = data
    self.type = type
    if data_labels is None or len(data_labels) == 0 :
      self.x_label = "X"
      self.y_label = [ "Y" for i in xrange(1, len(data)) ]
    else :
      self.x_label = data_labels[0]
      self.y_labels = [ data_labels[i] for i in xrange(1, len(data)) ]

  def get_plots (self, fill_in_missing_y=None) :
    plots = []
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

  def as_loggraph (self) :
    pass

class histogram_data (object) :
  pass


def import_ccp4i_logfile (file_name=None, log_lines=None) :
  import re
  assert file_name is not None or log_lines is not None
  if not log_lines :
    log_lines = open(file_name).readlines()
  current_lines = None
  tables_raw = []
  sections_read = 0
  for line in log_lines :
    if line[0:7] == "$TABLE:" :
      current_lines = [line.strip()]
      sections_read = 0
    elif re.sub("\s*$", "", line)[-2:] == "$$" and current_lines is not None :
      current_lines.append(line.strip())
      sections_read += 1
      if sections_read == 4 :
        tables_raw.append(current_lines)
        current_lines = None
    elif sections_read < 4 and current_lines is not None :
      current_lines.append(line.strip())
  tables = []
  for loggraph in tables_raw :
    t = table_data(None)
    t.import_loggraph(loggraph)
    tables.append(t)
  return tables

#---end
