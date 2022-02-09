from __future__ import absolute_import, division, print_function

from functools import cmp_to_key
from six import string_types
from six.moves import cStringIO

import sys
from six.moves import zip

def format_none(format, null_value=0, replace_with="None"):
  assert isinstance(replace_with, str)
  value_width = len(replace_with)
  return " " * max(0, len(format % null_value) - value_width) + replace_with

def format_value(format, value, null_value=0, replace_none_with="None"):
  if (format is None): return str(value)
  if (value is None or str(value).upper()=="NONE"):
    return format_none(format=format, null_value=null_value,
      replace_with=replace_none_with)
  return format % value

def round_2_for_cif(value):
  """
  shortcut for cif output
  """
  if value is None or str(value).upper()=="NONE":
    return '?'
  return format_value(format="%.2f", value=value, null_value=0, replace_none_with="?")

def round_4_for_cif(value):
  """
  shortcut for cif output
  """
  if value is None or str(value).upper()=="NONE":
    return '?'
  return format_value(format="%.4f", value=value, null_value=0, replace_none_with="?")

def pad_string(line, width=72, border="|"):
  n_spaces = width - (len(border) * 2) - len(line)
  return border + line + (" " * n_spaces) + border

def py_string_representation(string, preferred_quote, alternative_quote):
  "based on stringobject.c, Python 2.7.2"
  quote = preferred_quote
  if (    alternative_quote != preferred_quote
      and string.find(preferred_quote) >= 0
      and string.find(alternative_quote) < 0):
    quote = alternative_quote
  result = [quote]
  rapp = result.append
  for c in string:
    if (c == quote or c == '\\'):
      rapp('\\')
      rapp(c)
    elif (c == '\t'):
      rapp('\\t')
    elif (c == '\n'):
      rapp('\\n')
    elif (c == '\r'):
      rapp('\\r')
    elif (c < ' ' or c >= chr(0x7f)):
      rapp("\\x%02x" % (ord(c) & 0xff))
    else:
      rapp(c)
  rapp(quote)
  return "".join(result)

_have_string_representation = False
def string_representation(string, preferred_quote, alternative_quote):
  global string_representation
  global _have_string_representation
  if (not _have_string_representation):
    _have_string_representation = True
    try:
      from boost_adaptbx.boost.python import ext as _
      string_representation = _.string_representation
    except Exception:
      string_representation = py_string_representation
  return string_representation(string, preferred_quote, alternative_quote)

def split_keeping_spaces(s):
  result = []
  field = []
  prev = " "
  for c in s:
    if ((prev == " ") != (c == " ")):
      if (len(field) != 0):
        result.append("".join(field))
        field = []
    field.append(c)
    prev = c
  if (len(field) != 0):
    result.append("".join(field))
  return result

def size_as_string_with_commas(sz):
  if (sz is None): return "unknown"
  if (sz < 0):
    sz = -sz
    sign = "-"
  else:
    sign = ""
  result = []
  while True:
    if (sz >= 1000):
      result.insert(0, "%03d" % (sz % 1000))
      sz //= 1000
      if (sz == 0): break
    else:
      result.insert(0, "%d" % sz)
      break
  return sign + ",".join(result)

def show_string(s):
  if (s is None): return None
  if (s.find('"') < 0): return '"'+s+'"'
  if (s.find("'") < 0): return "'"+s+"'"
  return '"'+s.replace('"','\\"')+'"'

def prefix_each_line_suffix(prefix, lines_as_one_string, suffix, rstrip=True):
  def do_rstip(s): return s.rstrip()
  def do_nothing(s): return s
  if (rstrip): do = do_rstip
  else:        do = do_nothing
  return "\n".join([do(prefix+line+suffix)
    for line in lines_as_one_string.splitlines()])

def prefix_each_line(prefix, lines_as_one_string, rstrip=True):
  return prefix_each_line_suffix(
    prefix=prefix,
    lines_as_one_string=lines_as_one_string,
    suffix="",
    rstrip=rstrip)

def show_text(line, prefix="", out=None):
  if(out is None): out = sys.stdout
  print(file=out)
  print("%s%s"%(prefix, line), file=out)
  print(file=out)

def make_header(line, out=None, header_len=80, repeat=1, new_line_start=False,
                new_line_end=False):
  if (out is None): out = sys.stdout
  #if(new_line_start):
  #  print("\n",file=out)
  def one():
    line_len = len(line)
    #assert line_len <= header_len
    fill_len = header_len - line_len
    fill_rl = fill_len//2
    fill_r = fill_rl
    fill_l = fill_rl
    if (fill_rl*2 != fill_len):
      fill_r +=1
    out_string = "\n"+"="*(fill_l-1)+" "+line+" "+"="*(fill_r-1)+"\n"
    if(len(out_string) > 80):
      out_string = "\n"+"="*(fill_l-1)+" "+line+" "+"="*(fill_r-2)+"\n"
    print(out_string, file=out)
  for it in range(0,repeat):
    one()
  #if(new_line_end):
  #  print("\n",file=out)
  out.flush()

def make_sub_header(text, out=None, header_len=80, sep='-'):
  """
  Prints a subheader to an output location.

  Parameters
  ----------
  text : str
  out : file, optional
  header_len : int, optional
  sep : str, optional

  Examples
  --------
  >>> from libtbx.str_utils import make_sub_header
  >>> make_sub_header("Predicted ions", header_len=20)

     ----------Predicted ions----------

  >>>
  """
  if (out is None): out = sys.stdout
  assert isinstance(sep, string_types)
  border = sep*10
  line = border+text+border
  line_len = len(line)
  #assert line_len <= header_len
  fill_len = header_len - line_len
  fill_rl = fill_len//2
  fill_r = fill_rl
  fill_l = fill_rl
  if (fill_rl*2 != fill_len):
    fill_r +=1
  out_string = "\n"+" "*(fill_l-1)+" "+line+" "+" "*(fill_r-1)+"\n"
  if(len(out_string) > 80):
    out_string = "\n"+" "*(fill_l-1)+" "+line+" "+" "*(fill_r-2)+"\n"
  print(out_string, file=out)
  out.flush()

def wordwrap(text, max_chars=80):
  words = text.split()
  lines = []
  current_line = ""
  for word in words :
    if ((len(current_line) + len(word)) >= max_chars):
      lines.append(current_line)
      current_line = word
    else :
      current_line += " " + word
  lines.append(current_line)
  return "\n".join(lines)

def reformat_terminal_text(text):
  text.strip()
  lines = text.splitlines()
  newtext = ""
  lastline = ""
  for i, line in enumerate(lines):
    if line == "" and i != 0 :
      newtext += "\n"
    else :
      if lastline != "" :
        newtext += " "
      newtext += line.strip()
    lastline = line
  return newtext

def rstrip_lines(text):
  new = []
  for line in text.splitlines():
    new.append(line.rstrip())
  return "\n".join(new)

def strip_lines(text):
  new = []
  for line in text.splitlines():
    new.append(line.strip())
  return "\n".join(new)

def show_sorted_by_counts(
      label_count_pairs,
      reverse=True,
      out=None,
      prefix="",
      annotations=None):
  assert annotations is None or len(annotations) == len(label_count_pairs)
  if (out is None): out = sys.stdout
  if (len(label_count_pairs) == 0): return False
  def sort_function(a, b):
    if (reverse):
      if (a[1] > b[1]): return -1
      if (a[1] < b[1]): return  1
    else:
      if (a[1] < b[1]): return -1
      if (a[1] > b[1]): return  1
    return (a[0] > b[0]) - (a[0] < b[0]) # cmp(a[0], b[0])
  if (annotations is None):
    annotations = [None]*len(label_count_pairs)
  lca = [(show_string(l), c, a)
    for (l,c),a in zip(label_count_pairs, annotations)]
  lca.sort(key=cmp_to_key(sort_function))
  fmt = "%%-%ds %%%dd" % (
    max([len(l) for l,c,a in lca]),
    max([len(str(lca[i][1])) for i in [0,-1]]))
  for l,c,a in lca:
    if (a is not None and len(a) > 0):
      print(prefix + fmt % (l,c), end=' ', file=out)
      print(a, end='', file=out)
    else:
      print(prefix + fmt % (l,c), end='', file=out)
    print(file=out)
  return True

def overwrite_at(s, offset, replacement):
  return s[:offset] + replacement + s[offset+len(replacement):]

def contains_one_of(label, patterns):
  for pattern in patterns:
    if (label.find(pattern) >= 0):
      return True
  return False

def line_breaker(string, width):
  if (width <= 0 or len(string) == 0):
    yield string
  else:
    i_block_start = 0
    i_last_space = None
    for i,c in enumerate(string):
      if (c == " "):
        i_last_space = i
      if (i-i_block_start >= width and i_last_space is not None):
        yield string[i_block_start:i_last_space]
        i_block_start = i_last_space + 1
        i_last_space = None
    if (i_block_start < len(string)):
      yield string[i_block_start:]

class line_feeder(object):

  def __init__(self, f):
    self.f = iter(f)
    self.eof = False

  def __iter__(self):
    return self

  def next(self):
    if (not self.eof):
      try:
        return self.f.next()[:-1]
      except StopIteration:
        self.eof = True
    return ""

  def next_non_empty(self):
    while True:
      result = next(self)
      if (self.eof or len(result.strip()) != 0):
        return result

# StringIO with pickling support and extra read-only semantics
class StringIO(object):
  def __init__(self, *args, **kwds):
    self._buffer = cStringIO(*args, **kwds)
    self.read_only = bool(args)

  def __getattr__(self, *args, **kwds):
    if self.read_only and args[0] in ('softspace', 'truncate', 'write', 'writelines'):
      raise AttributeError("'libtbx.str_utils.StringIO' has no attribute '%s'" % args[0])
    return getattr(self._buffer, *args, **kwds)

  def __getstate__(self):
    is_output_object = (type(self._buffer).__name__ == 'StringO')
    return (self._buffer.getvalue(), is_output_object)

  def __setstate__(self, state):
    (buffer_value, is_output_object) = state
    if is_output_object :
      self.__init__()
      self._buffer.write(buffer_value)
      self.read_only = True
    else :
      self.__init__(buffer_value)

def expandtabs_track_columns(s, tabsize=8):
  result_e = []
  result_j = []
  j = 0
  for i,c in enumerate(s):
    if (c == "\t"):
      if (tabsize > 0):
        n = tabsize - (j % tabsize)
        result_e.extend([" "]*n)
        result_j.append(j)
        j += n
    else:
      result_e.append(c)
      result_j.append(j)
      if (c == "\n" or c == "\r"):
        j = 0
      else:
        j += 1
  return "".join(result_e), result_j

class framed_output(object):
  """
  Pseudo-output file wrapper for drawing a box around arbitrary content, for
  example:

    |----------Refinement stats--------|
    | r_free = 0.1234  r_work = 0.1567 |
    |----------------------------------|

  The actual content will be buffered until the close() method is called or
  the object is deleted, at which point the formatted text will be printed to
  the actual filehandle.

  :param out: actual output filehandle
  """
  def __init__(self,
      out,
      title=None,
      width=None,
      frame='-',
      center=False,
      center_title=None,
      #center_frame=None,
      max_width=80,
      prefix=None):
    self.out = out
    if (title is not None):
      title = title.strip()
    self.title = title
    if (center_title is None):
      center_title = center
    self.center_title = center_title
    self.width = width
    self.frame = frame
    self.center = center
    #self.center_frame = center_frame
    self.max_width = max_width
    self.prefix = prefix
    self.content = []
    self._current_line = None
    self._closed = False

  def get_best_text_width(self):
    return self.width - 4

  def write(self, text):
    current_text = ""
    if (self._current_line is not None):
      current_text += self._current_line
      self._current_line = None
    for char in text :
      if (char == "\n"):
        self.content.append(current_text)
        current_text = ""
      else :
        current_text += char
    if (current_text != ""):
      self._current_line = current_text

  def flush(self):
    pass

  def add_separator(self):
    if (self._current_line is not None):
      self.content.append(self._current_line)
      self._current_line = None
    self.content.append(None)

  def close(self):
    assert (not self._closed)
    self._closed = True
    if (not self._current_line in [None, ""]):
      self.content.append(self._current_line)
    out = self.out
    lines = self.content
    # misc setup
    text_lines = [ s for s in lines if (s is not None) ]
    side_frame = self.frame
    if (side_frame in ['-', '_']):
      side_frame = "|"
    width = self.width
    if (width is None):
      width = min(self.max_width, 4 + max([ len(s) for s in text_lines ]))
    def get_padding(text, margin=2, center=self.center):
      from libtbx.math_utils import ifloor, iceil
      fill = max(0, width - len(text) - (margin * 2))
      if (center):
        rfill = ifloor(fill / 2)
        lfill = iceil(fill / 2)
      else :
        rfill = 0
        lfill = fill
      return (rfill, lfill)
    def write_corner():
      if (self.frame != '_'):
        out.write(side_frame)
      else :
        out.write("_")
    def write_sep():
      if (self.prefix is not None):
        out.write(self.prefix)
      out.write(side_frame)
      out.write(self.frame * (width - 2))
      out.write(side_frame + "\n")
    # header
    out.write("\n")
    if (self.prefix is not None):
      out.write(self.prefix)
    write_corner()
    if (self.title is not None):
      rf, lf = get_padding(self.title, margin=1, center=self.center_title)
      rf += 1
      lf -= 1
      out.write(self.frame * rf)
      out.write(self.title)
      if (lf > 0):
        out.write(self.frame * lf)
    else :
      out.write(self.frame * (width - 2))
    write_corner()
    out.write("\n")
    # main output
    for line in lines :
      if (line is None):
        write_sep()
      else :
        if (self.prefix is not None):
          out.write(self.prefix)
        out.write(side_frame + " ")
        rf, lf = get_padding(line)
        if (rf > 0):
          out.write(" " * rf)
        out.write(line)
        if (lf > 0):
          out.write(" " * lf)
        out.write(" " + side_frame)
        out.write("\n")
    # footer
    if (self.prefix is not None):
      out.write(self.prefix)
    write_corner()
    self.out.write(self.frame * (width - 2))
    write_corner()
    out.write("\n")

  def __del__(self):
    if (not self._closed):
      self.close()

def print_message_in_box(message, **kwds):
  box = framed_output(**kwds)
  for line in line_breaker(message, box.get_best_text_width()):
    print(line, file=box)
  del box

def make_big_header(line, out=None, header_len=80, border_char="#"):
  """
  Alternate API for print_message_in_box, for compatibility with make_header
  """
  if (out is None) : out = sys.stdout
  print_message_in_box(
    out=out,
    message=line,
    width=header_len,
    frame=border_char,
    center=True)

class find_matching_closing_symbol(object):
  """ This functor deals with a pair of symbol, an opening symbol and a closing one,
      the archetypical example being parentheses. Given the position of an opening
      symbol, it returns the position of the matching closing symbol, or -1 if it does not
      exist.
  """

  def __init__(self, opening, closing):
    """ Initialise the functor to work with the given pair. Each of `opening` and `closing`
        may be several character long.
    """
    import re
    self.opening, self.closing = opening, closing
    self._regex = re.compile("(%s)|(%s)" % (re.escape(opening), re.escape(closing)))

  def __call__(self, string, pos):
    """ An opening symbol shall start at position `pos` of the given `string`.
    """
    level = 1
    delta = (None, +1, -1)
    for m in self._regex.finditer(string, pos + len(self.opening)):
      if not m:
        return -1
      level += delta[m.lastindex]
      if level == 0:
        return m.start()
    return -1
