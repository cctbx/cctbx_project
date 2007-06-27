from libtbx.utils import Sorry

def remove_redundant_spaces(s):
  s = s.split()
  s = " ".join([i for i in s])
  return s

def format_string(object, ind, log, scope=False, allowed_line_lenght=79):
  name = str(object.name).strip()
  help = str(object.help).strip()
  help = remove_redundant_spaces(help)
  if(help == "None"): help = ""
  if(scope):
     if(len(" "*ind+name) >= allowed_line_lenght):
        raise Sorry(
          "Cannot create doc file: scope is too deep or its name is too long.")
     fmt = "%s<b>%s</b> <b><FONT color=blue>%s</FONT></b>"
     elements = (" "*ind, name, help)
     line = fmt % elements
     line_to_appear = "%s%s %s"%elements
     if(len(line_to_appear) <= allowed_line_lenght):
        print >> log, line
     else:
        help_start = len(" "*ind+name)+0
        help_offset = " "*help_start
        n_cut = partial_line(line, allowed_line_lenght, line_to_appear)
        print >> log, line[:n_cut]
        line = help_offset+line[n_cut:]
        while True:
          if(n_cut is None): break
          n_cut = partial_line(line, allowed_line_lenght, line_to_appear)
          if(len(line[:n_cut].strip())==0):
             print >> log, line[:allowed_line_lenght]
             line = help_offset+line[allowed_line_lenght:]
          else:
             print >> log, line[:n_cut]
             line = help_offset+line[n_cut:]
  else:
     values = (" ".join([str(i) for i in object.words])).strip()
     if(len(" "*ind+name) >= allowed_line_lenght):
        raise Sorry(
          "Cannot create doc file: scope is too deep or its name is too long.")
     elements = (" "*ind, name, values, help)
     fmt = "<!--ANCHOR--> %s%s= <FONT color=CC0000>%s</FONT> <FONT color=blue>%s</FONT>"
     line = fmt % elements
     line_length = len(" "*ind+str(object.name)+str(object.help)+values)
     line_to_appear = "%s%s= %s %s"%elements
     if(len(line_to_appear) <= allowed_line_lenght):
        print >> log, line
     else:
        help_start = len(" "*ind+name)+1
        help_offset = " "*help_start
        n_cut = partial_line(line, allowed_line_lenght, line_to_appear)
        print >> log, line[:n_cut]
        line = help_offset+line[n_cut:]
        while True:
          if(n_cut is None): break
          n_cut = partial_line(line, allowed_line_lenght, line_to_appear)
          if(len(line[:n_cut].strip())==0):
             print >> log, line[:allowed_line_lenght]
             line = help_offset+line[allowed_line_lenght:]
          else:
             print >> log, line[:n_cut]
             line = help_offset+line[n_cut:]

def partial_line(line, allowed_line_lenght, line_to_appear):
  assert len(line) > 0
  brakets = []
  found_l, found_r = False, False
  n_l, n_r = None, None
  for i_seq, item in enumerate(line):
    if(item == "<"):
       found_l = True
       n_l = i_seq
       if(not (not found_r and n_r is None)):
          raise Sorry("Bad line: %s"%line)
    if(item == ">"):
       found_r = True
       n_r = i_seq
       if(not (found_l and n_l is not None)):
          raise Sorry("Bad line: %s"%line)
    if(found_l and found_r):
       brakets.append([n_l, n_r])
       found_l, found_r = False, False
       n_l, n_r = None, None
  counter = 0
  for n in xrange(len(line)):
    in_braket = False
    for braket in brakets:
      if(n>=braket[0] and n<=braket[len(braket)-1]):
         in_braket = True
         break
    if(not in_braket):
       counter += 1
       in_braket = False
    if(counter >= allowed_line_lenght):
       i = counter
       j = n
       while i > 0:
         if(line[j] == " "):
            return j
         else:
            in_braket = False
            for braket in brakets:
              if(j>=braket[0] and j<=braket[len(braket)-1]):
                in_braket = True
                return braket[0]
            i -= 1
            j -= 1

def scope_walk(p, ind, log):
  values = []
  objects = []
  for obj in p.objects:
    try:
      x = obj.name
      y = obj.type
      values.append(obj)
    except:
      objects.append(obj)
  for object in values:
    format_string(object, ind, log)
  for object in objects:
    format_string(object, ind, log, scope = True)
    scope_walk(object, ind+3, log)

def header(log):
  legend = """\
Legend: black <b>bold</b> - scope names
        black - parameter names
        <FONT color=CC0000>red</FONT> - parameter values
        <FONT color=blue>blue</FONT> - parameter help
        <FONT color=blue><b>blue bold</b></FONT> - parameter help
        Parameter values:
          <FONT color=CC0000>*</FONT> means selected parameter (where multiple choices are available)
          <FONT color=CC0000>False</FONT> is No
          <FONT color=CC0000>True</FONT> is Yes
          <FONT color=CC0000>None</FONT> means not provided or not predefined or left up to program
          <FONT color=CC0000>"%3d"</FONT> is Python style formatting descriptor """
  print >> log, "<PRE><FONT face=courier>"
  print >> log, "<b>%s"%("-"*79),"</b>"
  print >> log, legend
  print >> log, "<b>%s"%("-"*79),"</b>"

def run(phil_object, log,skip_head=False):
  if not skip_head:
    print >> log, """<META http-equiv=Content-Type content="text/html; charset=utf-8">"""
    print >> log, """<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">"""

  header(log)
  scope_walk(p = phil_object, ind = 0, log = log)
  print >> log, "</FONT></PRE>"
