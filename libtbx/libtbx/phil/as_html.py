from libtbx.utils import Sorry

def remove_redundant_spaces(s):
  s = s.split()
  s = " ".join([i for i in s])
  return s

def format_string(object, ind, log, scope=False, allowed_line_length=79):
  name = str(object.name).strip()
  help = str(object.help).strip()
  help = remove_redundant_spaces(help)
  if(help == "None"): help = ""
  if(scope):
     if (ind + len(name) >= allowed_line_length):
        raise Sorry(
          "Cannot create doc file: scope is too deep or its name is too long.")
     line = '%s<b>%s</b>' % (" "*ind, name)
     line_to_appear = " "*ind + name
     if (len(help) != 0):
        line += ' <b><SPAN STYLE="color: blue">%s</SPAN></b>' % help
        line_to_appear += " " + help
     if (len(line_to_appear) <= allowed_line_length):
        print >> log, line
     else:
        help_offset = " "*(ind+len(name))
        n_cut = partial_line(line, allowed_line_length, line_to_appear)
        assert n_cut is not None
        print >> log, line[:n_cut]
        line = help_offset+line[n_cut:]
        while True:
          n_cut = partial_line(line, allowed_line_length, line_to_appear)
          if(n_cut is None): break
          if(len(line[:n_cut].strip())==0):
             print >> log, line[:allowed_line_length]
             line = help_offset+line[allowed_line_length:]
          else:
             print >> log, line[:n_cut]
             line = help_offset+line[n_cut:]
  else:
     values = (" ".join([str(i) for i in object.words])).strip()
     if (ind + len(name) >= allowed_line_length):
        raise Sorry(
          "Cannot create doc file: scope is too deep or its name is too long.")
     elements = [" "*ind, name, values]
     fmt = '<!--ANCHOR-->%s%s= <SPAN STYLE="color: #CC0000">%s</SPAN>'
     fmt_plain = '%s%s= %s'
     if (len(help) != 0):
        fmt += ' <SPAN STYLE="color: blue">%s</SPAN>'
        fmt_plain += ' %s'
        elements.append(help)
     elements = tuple(elements)
     line = fmt % elements
     line_to_appear = fmt_plain % elements
     if(len(line_to_appear) <= allowed_line_length):
        print >> log, line
     else:
        help_offset = " "*(ind+len(name)+1)
        n_cut = partial_line(line, allowed_line_length, line_to_appear)
        assert n_cut is not None
        print >> log, line[:n_cut]
        line = help_offset+line[n_cut:]
        while True:
          n_cut = partial_line(line, allowed_line_length, line_to_appear)
          if(n_cut is None): break
          if(len(line[:n_cut].strip())==0):
             print >> log, line[:allowed_line_length]
             line = help_offset+line[allowed_line_length:]
          else:
             print >> log, line[:n_cut]
             line = help_offset+line[n_cut:]

def partial_line(line, allowed_line_length, line_to_appear):
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
    if(counter >= allowed_line_length):
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

def write_legend(log):
  text = """\
Legend: black <B>bold</B> - scope names
        black - parameter names
        <SPAN STYLE="color: #CC0000">red</SPAN> - parameter values
        <SPAN STYLE="color: blue">blue</SPAN> - parameter help
        <SPAN STYLE="color: blue"><B>blue bold</B></SPAN> - scope help
        Parameter values:
          <SPAN STYLE="color: #CC0000">*</SPAN> means selected parameter (where multiple choices are available)
          <SPAN STYLE="color: #CC0000">False</SPAN> is No
          <SPAN STYLE="color: #CC0000">True</SPAN> is Yes
          <SPAN STYLE="color: #CC0000">None</SPAN> means not provided, not predefined, or left up to the program
          <SPAN STYLE="color: #CC0000">"%3d"</SPAN> is a Python style formatting descriptor"""
  print >> log, '<PRE style="face=courier">'
  print >> log, "<B>%s"%("-"*79),"</B>"
  print >> log, text
  print >> log, "<B>%s"%("-"*79),"</B>"

def run(phil_object, log, skip_head=False, head_title="libtbx.phil.as_html"):
  if (not skip_head):
    print >> log, """\
<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
   "http://www.w3.org/TR/html4/loose.dtd">
<HEAD>
<META http-equiv=Content-Type content="text/html; charset=utf-8">
<TITLE>%s</TITLE>
</HEAD>
<BODY>
""" % head_title
  write_legend(log=log)
  scope_walk(p=phil_object, ind=0, log=log)
  print >> log, "</PRE>"
  if (not skip_head):
    print >> log
    print >> log, "</BODY>"
