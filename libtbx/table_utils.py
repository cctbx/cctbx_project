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
