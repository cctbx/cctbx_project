import cStringIO
import sys

def format_none(format, null_value=0, replace_with="None"):
  assert isinstance(replace_with, str)
  value_width = len(replace_with)
  return " " * max(0, len(format % null_value) - value_width) + replace_with

def format_value(format, value, null_value=0, replace_none_with="None"):
  if (format is None): return str(value)
  if (value is None):
    return format_none(format=format, null_value=null_value,
      replace_with=replace_none_with)
  return format % value

def pad_string (line, width=72, border="|") :
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
      from boost.python import ext as _
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

def make_header (line, out=None, header_len=80):
  if (out is None): out = sys.stdout
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
  print >> out, out_string
  out.flush()

def make_sub_header(text, out=None, header_len=80):
  if (out is None): out = sys.stdout
  line = "----------"+text+"----------"
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
  print >> out, out_string
  out.flush()

def wordwrap (text, max_chars=80) :
  words = text.split()
  lines = []
  current_line = ""
  for word in words :
    if ((len(current_line) + len(word)) >= max_chars) :
      lines.append(current_line)
      current_line = word
    else :
      current_line += " " + word
  lines.append(current_line)
  return "\n".join(lines)

def reformat_terminal_text (text) :
  text.strip()
  lines = text.splitlines()
  newtext = ""
  lastline = ""
  for i, line in enumerate(lines) :
    if line == "" and i != 0 :
      newtext += "\n"
    else :
      if lastline != "" :
        newtext += " "
      newtext += line.strip()
    lastline = line
  return newtext

def rstrip_lines (text) :
  new = []
  for line in text.splitlines() :
    new.append(line.rstrip())
  return "\n".join(new)

def strip_lines (text) :
  new = []
  for line in text.splitlines() :
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
      if (a[1] > b[1]): return -1
      if (a[1] < b[1]): return  1
    return cmp(a[0], b[0])
  if (annotations is None):
    annotations = [None]*len(label_count_pairs)
  lca = [(show_string(l), c, a)
    for (l,c),a in zip(label_count_pairs, annotations)]
  lca.sort(sort_function)
  fmt = "%%-%ds %%%dd" % (
    max([len(l) for l,c,a in lca]),
    max([len(str(lca[i][1])) for i in [0,-1]]))
  for l,c,a in lca:
    print >> out, prefix+fmt % (l,c),
    if (a is not None and len(a) > 0): print >> out, a,
    print >> out
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
    while 1:
      result = self.next()
      if (self.eof or len(result.strip()) != 0):
        return result

# cStringIO with pickling support
class StringIO (object) :
  def __init__ (self, *args, **kwds) :
    self._buffer = cStringIO.StringIO(*args, **kwds)

  def __getattr__ (self, *args, **kwds) :
    return getattr(self._buffer, *args, **kwds)

  def __getstate__ (self) :
    is_output_object = (type(self._buffer).__name__ == 'StringO')
    return (self._buffer.getvalue(), is_output_object)

  def __setstate__ (self, state) :
    (buffer_value, is_output_object) = state
    if is_output_object :
      self.__init__()
      self._buffer.write(buffer_value)
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
