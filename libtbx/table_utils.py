## Some functionality for formatted printing,
## indentation and simple tables
##
## Take from
##   http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/267662
## slightly modifed functionality.
##
import libtbx.forward_compatibility # for sum
import cStringIO,operator

def format(rows,
           comments=None,
           has_header=False,
           header_char='-',
           delim=' | ',
           justify='left',
           separate_rows=False,
           leading_and_terminal_separator=True ,
           prefix='',
           postfix='',
           wrapfunc=lambda x:x):
    """Indents a table by column.
       - rows: A sequence of sequences of items, one sequence per row.
       - hasHeader: True if the first row consists of the columns' names.
       - headerChar: Character to be used for the row separator line
         (if hasHeader==True or separateRows==True).
       - delim: The column delimiter.
       - justify: Determines how are data justified in their column.
         Valid values are 'left','right' and 'center'.
       - separateRows: True if rows are to be separated by a line
         of 'headerChar's.
       - prefix: A string prepended to each printed row.
       - postfix: A string appended to each printed row.
       - wrapfunc: A function f(text) for wrapping text; each element in
         the table is first wrapped by this function."""
    # closure for breaking logical rows to physical, using wrapfunc
    def row_wrapper(row):
        new_rows = [wrapfunc(item).split('\n') for item in row]
        return [[substr or '' for substr in item] for item in map(None,*new_rows)]
    # break each logical row into one or more physical ones
    logical_rows = [row_wrapper(row) for row in rows]
    # columns of physical rows
    columns = map(None,*reduce(operator.add,logical_rows))
    # get the maximum of each column by the string length of its items
    max_widths = [max([len(str(item))
      for item in column]) for column in columns]
    row_separator = header_char * (
                      len(prefix)
                    + len(postfix)
                    + sum(max_widths)
                    + len(delim)*(len(max_widths)-1))
    # select the appropriate justify method
    justify = {'center':str.center, 'right':str.rjust, 'left':str.ljust}[justify.lower()]
    output=cStringIO.StringIO()

    total_row_width=0
    line=None
    # Printing the table with values
    if leading_and_terminal_separator: print >> output, row_separator
    for physical_rows in logical_rows:
        for row in physical_rows:
            line =  prefix \
                + delim.join([justify(str(item),width) for (item,width) in zip(row,max_widths)]) \
                + postfix
            print >> output, line
        if separate_rows or has_header: print >> output, row_separator; has_header=False
    if not separate_rows:
        if leading_and_terminal_separator:
            print >> output, row_separator
    if comments is not None:
        print >> output, wrap_onspace(comments,len(line))
    return "\n".join(
      [line.rstrip() for line in output.getvalue().splitlines()])

# written by Mike Brown
# http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/148061
def wrap_onspace(text, width):
    """
    A word-wrap function that preserves existing line breaks
    and most spaces in the text. Expects that existing line
    breaks are posix newlines (\n).
    """
    return reduce(lambda line, word, width=width: '%s%s%s' %
                  (line,
                   ' \n'[(len(line[line.rfind('\n')+1:])
                         + len(word.split('\n',1)[0]
                              ) >= width)],
                   word),
                  text.split(' ')
                 )

import re
def wrap_onspace_strict(text, width):
    """Similar to wrap_onspace, but enforces the width constraint:
       words longer than width are split."""
    wordRegex = re.compile(r'\S{'+str(width)+r',}')
    return wrap_onspace(wordRegex.sub(lambda m: wrap_always(m.group(),width),text),width)

import math
def wrap_always(text, width):
    """A simple word-wrap function that wraps text on exactly width characters.
       It doesn't split the text in words."""
    return '\n'.join([ text[width*i:width*(i+1)] \
                       for i in xrange(int(math.ceil(1.*len(text)/width))) ])

# Nat's utilities for plottable data
class table_data (object) :
  def __init__ (self, title, series_names=[], series_types=[],
      series_labels=[], graph_names=[], graph_series=[], data=[]) :
    self.title = title
    self._is_complete = False
    self.series_names = series_names
    self.series_types = series_types
    self.series_labels = series_labels
    self.graph_names = graph_names
    self.graph_columns = graph_series
    self.data = data
    self._graphs = {}

  def formatted_simple (self) :
    pass

  def formatted (self) :
    pass

  def formatted_html (self) :
    pass

  def formatted_loggraph (self) :
    pass

  def get_graph (self, graph_name=None, series_list=[]) :
    if graph_name is not None :
      if not graph_name in self._graphs :
        if not graph_name in self.graph_names :
          return None
        n = self.graph_names.index(graph_name)
        x = self.data
        data = [ [ x[i][j] for i in xrange(len(x)) ] \
                    for j in self.graph_seriess[n] ]
        if len(self.series_names) == len(self.data) :
          labels = [self.series_names[i] for i in self.graph_columns[n]]
        else :
          labels = []
        self._graphs[graph_name] = graph_data(graph_name, data, "plot", labels)
      return self._graphs[graph_name]
    elif len(series_list) > 1 :
      if not series_list in self._graphs :
        data = None
        if isinstance(series_list[0], int) :
          data = self._extract_data_series(series_list)
          labels = [ self.series_labels[i] for i in column_list ]
        elif isinstance(series_list[0], str) :
          n_list = []
          for col in series_list :
            n_list.append(self.series_names.index(col))
          data = self._extract_data_series(n_list)
          labels = [ self.series_labels[i] for i in n_list ]
        if data is None :
          return None
        self._graphs[series_list] = graph_data(None, data, "plot", labels)
      return self._graphs[series_list]

  def _extract_data_series (self, series_list) :
    assert len(self.data) > 0
    assert len(series_list) <= len(self.data[0])
    data = self.data
    new_data = [ [ data[i][j] for j in xrange(len(x)) ] for i in series_list ]
    return new_data

  def __str__ (self) :
    return str(self.data)

class graph_data (object) :
  def __init__ (self, name, data, type="plot", data_labels=[]) :
    self.name = name
    self.data = data
    self.type = type
    if len(data_labels) == 0 :
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

#---
#--- tests
#---
def exercise():
  formatted_table = format(
    rows=[
      ["alpha","beta","gamma"],
      ["10","20","30"],
      ["40","50","60"]],
    comments="comments here",
    has_header=True)
  assert formatted_table == """\
--------------------
alpha | beta | gamma
--------------------
10    | 20   | 30
40    | 50   | 60
--------------------
comments here"""
  print "OK"

if (__name__ == "__main__"):
  exercise()
