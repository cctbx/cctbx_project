from __future__ import absolute_import, division, print_function
## Some functionality for formatted printing,
## indentation and simple tables
##
## Take from
##   http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/267662
## slightly modifed functionality.
##
from six.moves import range, zip, map
try: # Python 3
    from itertools import zip_longest
except ImportError: # Python 2
    from itertools import izip_longest as zip_longest
from functools import reduce
from six.moves import cStringIO as StringIO
import libtbx.forward_compatibility # for sum
import operator

def format(rows,
           comments=None,
           has_header=0,
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
    if has_header is True: has_header=1
    if has_header is False: has_header=0
    # closure for breaking logical rows to physical, using wrapfunc
    def row_wrapper(row):
        new_rows = [wrapfunc(item).split('\n') for item in row]
        return [[substr or '' for substr in item] for item in
                map(lambda *a: a, *list(zip(*zip_longest(*new_rows))))]
    # break each logical row into one or more physical ones
    logical_rows = [row_wrapper(row) for row in rows]
    # columns of physical rows
    columns = list(map(lambda *a: a, *list(zip(*zip_longest(*reduce(operator.add, logical_rows))))))
    # get the maximum of each column by the string length of its items
    max_widths = [max([len(str(item))
      for item in column]) for column in columns]
    prefix_minus_leading_spaces = prefix.lstrip()
    indent = len(prefix) - len(prefix_minus_leading_spaces)
    row_separator = " "*indent + header_char * (
                      len(prefix_minus_leading_spaces)
                    + len(postfix)
                    + sum(max_widths)
                    + len(delim)*(len(max_widths)-1))
    # select the appropriate justify method
    justify = {'center':str.center, 'right':str.rjust, 'left':str.ljust}[justify.lower()]
    output=StringIO()

    total_row_width=0
    line=None
    # Printing the table with values
    if leading_and_terminal_separator: print(row_separator, file=output)
    for physical_rows in logical_rows:
        for row in physical_rows:
            line =  prefix \
                + delim.join([justify(str(item),width) for (item,width) in zip(row,max_widths)]) \
                + postfix
            print(line, file=output)
        if separate_rows or has_header == 1: print(row_separator, file=output)
        if has_header : has_header-=1
    if not separate_rows:
        if leading_and_terminal_separator:
            print(row_separator, file=output)
    if comments is not None:
        print(wrap_onspace(comments,len(line)), file=output)
    return "\n".join(
      [line.rstrip() for line in output.getvalue().splitlines()])

def manage_columns(table_data, include_columns):
  for row in table_data:
    assert len(row) == len(include_columns)
  new_table = []
  for row in table_data:
    new_row = [row[idx] for idx in range(len(row)) if include_columns[idx]]
    new_table.append(new_row)
  return new_table

class simple_table(object):
  """
  Container for generic table contents, used in Xtriage and elsewhere.  The
  table cells are assumed to be pre-formatted as strings (but not fixed-width).
  This is designed to be easily fed into the format() function, but it can also
  be displayed in GUIs or HTML tables.
  """
  def __init__(self,
      table_rows,
      column_headers=(),
      comments=None):
    self.column_headers = list(column_headers)
    self._rows = table_rows

  def format(self, indent=0, equal_widths=None):
    # FIXME equal_widths is a placeholder, does not currently work
    prefix = " " * indent
    return format(
      rows=[ self.column_headers ]+self._rows,
      has_header=(len(self.column_headers) != 0),
      prefix=prefix+'| ',
      postfix=' |')

  def export(self):
    if (len(self.column_headers) != 0):
      return [ self.column_headers ] + self._rows
    else :
      return self._rows

  def export_rows(self) : return self.export()

  @property
  def n_rows(self):
    return len(self._rows)

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
                       for i in range(int(math.ceil(1.*len(text)/width))) ])

## Additional functionality; go beyond formatted printing; emulate spreadsheet
## Typical use case to include:
##   column headers, column formatting specifiers, default column values,
##   choice of which rows to print, optional summation or averaging of column values,
##   spreadsheet-style cell formulae
##
## This is original source code developed for the LABELIT project

class Spreadsheet:

  def __init__(self,rows=0):
    self.S_table_rows = rows
    self.summaries = []

  def addColumn(self,label,default_value = None):
    setattr(self,label,SpreadsheetColumn(self,self.S_table_rows,default_value))

  def addRow(self):
    self.S_table_rows += 1
    for item in self.__dict__:
      if isinstance(self.__dict__[item],SpreadsheetColumn):
        self.__dict__[item].append()

  def printTable(self,columns,printed_rows=None,out=None):
    if out==None: import sys; out=sys.stdout
    if printed_rows is None:  printed_rows=[True]*self.S_table_rows
    master = {}
    padding = {}
    for column in columns:
      padding[column] = max(
          len(column),
          len(self.__dict__[column].format%self.__dict__[column][0])
        )
      master[column] = "%%%ds "%(padding[column])
      print(master[column]%column, end=' ', file=out)
    print(file=out)
    for i in range(self.S_table_rows):
      if not printed_rows[i]: continue
      for column in columns:
        print(master[column]%(
          self.__dict__[column].format%self.__dict__[column][i] ,), end=' ', file=out)
      print(file=out)
    #print summations at the bottom
    for column in columns:
      if column in self.summaries:
        print(" "+"-"*(padding[column]), end=' ', file=out)
      else:
        print(" "*(1+padding[column]), end=' ', file=out)
    print(file=out)
    for column in columns:
      if column in self.summaries:
        print(master[column]%(
          self.__dict__[column].format%self.__dict__[column].sum(),), end=' ', file=out)
      elif hasattr(self,"wtmean") and column in self.wtmean:
        normalizer = self.weights.sum()
        summation = 0
        for x in range(self.S_table_rows):
            summation += self.weights[x] * self.__dict__[column][x]
        mean = summation/normalizer
        print((master[column]%(
          self.__dict__[column].format%mean,) ), end=' ', file=out)
      else:
        print(" "*(1+padding[column]), end=' ', file=out)
    print(file=out)

  def rows(self):
    return self.S_table_rows

  def eval(self,expression):
    return eval(expression)

class SpreadsheetColumn:

  def __init__(self,parent,rows,default_value):
    self.parent = parent
    self.default = default_value
    self.C_data = [self.default]*rows

  def __setitem__(self,index,value):
    #general case, normal integer index
    if isinstance(index, type(0)):
      self.C_data[index]=value

    #special case, consult parent to get pointer to index[0], then +=index[1]
    #this will be used to refer to the correct image frame for sublattice calc
    if isinstance(index, type((0,1))):
      pointer = self.parent.pointer(index[0])
      self.C_data[pointer + index[1]]=value

  def __getitem__(self,index):
    #general case, normal integer index with specializations for Formulae
   if isinstance(index, type(0)):
    V = self.C_data[index]
    if isinstance(V,Formula):
      if V.expression=='self.Population[%row] / self.Fract[%row]':
        if self.parent.Fract[index]==0:return 0.0 # save from Zero Division
        return self.parent.Population[index] / self.parent.Fract[index]
      else:
        E = V.expression.replace('%row','%d'%index)
        return self.parent.eval(E)
    else:
      return V
   if isinstance(index, type((0,1))):
      pointer = self.parent.pointer(index[0])
      return self.C_data[pointer + index[1]]

  def __len__(self):
    return len(self.C_data)

  def append(self):
    self.C_data.append(self.default)

  def setFormat(self,format):
    self.format = format

  def sum(self):
    def typesum(x,y): return(x+y)
    return reduce(typesum,[self[i] for i in range(len(self))])

class Formula:
  def __init__(self,expression):
    self.expression = expression

#---
#--- tests
#---
def excercise_spreadsheet():
  class typical_use_case_1(Spreadsheet): #Case 1 is required for spotfinder spot statistics
    def __init__(self):

      Spreadsheet.__init__(self,rows=6) # make 6 rows in the Table
      self.addColumn('Limit')
      self.addColumn('MeanI')
      self.addColumn('Population',0) # 0 is the default value

      for xrow in range(0,self.S_table_rows):
        self.Limit[xrow] = [12.543,9.959,8.700,7.90,7.34,6.90][xrow]
        self.MeanI[xrow] = [4348,461,313,378,376,0][xrow]

      population_data = [926,121,8,2,6]
      for c in range(len(population_data)): #reconcile inconsistent bin counts
          self.Population[c]=population_data[c]

      self.Limit.format = "%.2f"
      self.MeanI.format = "%.0f"
      self.Population.format = "%7d"

    def show(self,message,out,columns_to_print=['Limit','Population','MeanI']):
      self.summaries=['Population'] #display the sum of the Population column
      self.wtmean=['MeanI']         #display the weighted mean of MeanI column,
      self.weights=self.Population  #...weighted by the Population
      legend = """Detailed explanation of each column here"""
      print(message+":\n", file=out)
      to_print = [True,True,True,True,True,True]
      self.printTable(columns_to_print,printed_rows=to_print,out=out)
      print(legend, file=out)

  OC = StringIO()
  derived_1 = typical_use_case_1()
  derived_1.show(message="Analysis of signals after noise suppression",out=OC)
  assert "\n".join([i.rstrip() for i in OC.getvalue().split("\n")]) == \
"""Analysis of signals after noise suppression:

Limit  Population  MeanI
12.54         926   4348
 9.96         121    461
 8.70           8    313
 7.90           2    378
 7.34           6    376
 6.90           0      0
        ----------
             1063   3845
Detailed explanation of each column here
"""
  class typical_use_case_2(Spreadsheet): #Case 2 is required for spotfinder ice detection
    def __init__(self):
      #Total rows by formula
      Total_rows = 4
      Spreadsheet.__init__(self,rows=Total_rows)
      self.addColumn('Limit')
      self.addColumn('Fract')
      self.addColumn('Population',0)
      self.addColumn('adjustPop',
                     Formula('self.Population[%row] / self.Fract[%row]'))

      for xrow in range(Total_rows):
        self.Limit[xrow] = [39.34,31.22,27.28,24.78][xrow]
        self.Fract[xrow] = [1.0,2.0,3.0,4.0][xrow]
        self.Population[xrow]=[53,102,107,86][xrow]

      self.Limit.format = "%.2f"
      self.Fract.format = "%.2f"
      self.Population.format = "%7d"
      self.adjustPop.format = "%7.1f"

    def show(self,default=['Limit','Population','Fract','adjustPop'],out=None):
      self.printTable(default,out=out)
      #print >>out,"0"

  OC = StringIO()
  derived_2 = typical_use_case_2()
  derived_2.show(out=OC)

  assert "\n".join([i.rstrip() for i in OC.getvalue().split("\n")]) == \
"""Limit  Population  Fract  adjustPop
39.34          53   1.00       53.0
31.22         102   2.00       51.0
27.28         107   3.00       35.7
24.78          86   4.00       21.5


"""

def exercise_flat_table():
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
  print("OK")

def exercise_general():
  labels = ('First Name', 'Last Name', 'Age', 'Position')
  data = \
      '''John,Smith,24,Software Engineer
      Mary,Brohowski,23,Sales Manager
      Aristidis,Papageorgopoulos,28,Senior Reseacher'''
  rows = [row.strip().split(',')  for row in data.splitlines()]

  table = format([labels]+rows, has_header=True)
  assert table == '''\
-------------------------------------------------------
First Name | Last Name        | Age | Position
-------------------------------------------------------
John       | Smith            | 24  | Software Engineer
Mary       | Brohowski        | 23  | Sales Manager
Aristidis  | Papageorgopoulos | 28  | Senior Reseacher
-------------------------------------------------------'''

  # test indent with different wrapping functions
  width = 10
  wrapped_output = [
      '''\
----------------------------------------------
| First Name | Last Name  | Age | Position   |
----------------------------------------------
| John       | Smith      | 24  | Software E |
|            |            |     | ngineer    |
----------------------------------------------
| Mary       | Brohowski  | 23  | Sales Mana |
|            |            |     | ger        |
----------------------------------------------
| Aristidis  | Papageorgo | 28  | Senior Res |
|            | poulos     |     | eacher     |
----------------------------------------------''',
      '''\
---------------------------------------------------
| First Name | Last Name        | Age | Position  |
---------------------------------------------------
| John       | Smith            | 24  | Software  |
|            |                  |     | Engineer  |
---------------------------------------------------
| Mary       | Brohowski        | 23  | Sales     |
|            |                  |     | Manager   |
---------------------------------------------------
| Aristidis  | Papageorgopoulos | 28  | Senior    |
|            |                  |     | Reseacher |
---------------------------------------------------''',
      '''\
---------------------------------------------
| First Name | Last Name  | Age | Position  |
---------------------------------------------
| John       | Smith      | 24  | Software  |
|            |            |     | Engineer  |
---------------------------------------------
| Mary       | Brohowski  | 23  | Sales     |
|            |            |     | Manager   |
---------------------------------------------
| Aristidis  | Papageorgo | 28  | Senior    |
|            | poulos     |     | Reseacher |
---------------------------------------------'''
  ]
  for wrapper, output in zip((wrap_always,wrap_onspace,wrap_onspace_strict),
                             wrapped_output):
      table = format([labels]+rows, has_header=True, separate_rows=True,
                     prefix='| ', postfix=' |',
                     wrapfunc=lambda x: wrapper(x,width))
      assert table == output

if (__name__ == "__main__"):
  excercise_spreadsheet()
  exercise_flat_table()
  exercise_general()
