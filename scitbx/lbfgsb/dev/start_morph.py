#! /usr/local_cci/Python-2.2.1/bin/python

# Simple script used for converting L-BFGS-B 2.1 routines.f.
# The rest of the conversion was done manually.
# WARNING:
# This simple script inserts a stray semicolon in
# subroutine subsm which leads to wrong results:
#           d(i) = d(i) + wy(k,pointr)*wv(jy)/theta;
#                                                  ^
#                                                  here

from fileinput import input
import re
import sys

def process_do(line, do_line_dict):
  s = line.lower().replace(" ", "")
  m = re.match(r"do(\d+)([^\=]+)\=([^,]+),([^,]+),([^,]+)", s)
  if (m == None):
    m = re.match(r"do(\d+)([^\=]+)\=([^,]+),([^,]+)", s)
    assert m != None
    do_line_dict[int(m.group(1))] = line
    s = "for(int %s=%s;%s<=%s;%s++) {" % (
      m.group(2), m.group(3),
      m.group(2), m.group(4),
      m.group(2))
  else:
    do_line_dict[int(m.group(1))] = line
    s = "for(int %s=%s;%s<=%s;%s+=%s) {" % (
      m.group(2), m.group(3),
      m.group(2), m.group(4),
      m.group(2), m.group(5))
  return re.sub("do.*", s, line)

def process_continue(line, do_line_dict):
  m = re.match(r"\s*(\d+)\s*continue", line.lower())
  label = m.group(1)
  do_line = do_line_dict.get(int(label), None)
  if (do_line is not None):
    i = do_line.lower().find("do")
    return do_line[:i] + "}"
  i = line.lower().find("continue")
  return " "*i + "lbl_%s:" % label

do_line_dict = {}
for line in input():
  line = line.rstrip()
  if (line.lower().strip() == ""):
    line = ""
  elif (line.lower().strip() == "end"):
    line = ""
    do_line_dict = {}
  elif (line.lower().startswith("c")):
    assert len(line) == 1 or line[1] in [" ", "=", "-", "c"], line
    line = "//" + line[2:]
  elif (line.lower().strip().startswith("do ")):
    line = process_do(line, do_line_dict)
  elif (re.match(r"\s*\d+\s*continue", line.lower())):
    line = process_continue(line, do_line_dict)
  else:
    line = (line.lower()
      .replace(".eq.", "==")
      .replace(".ne.", "!=")
      .replace(".lt.", "<")
      .replace(".gt.", ">")
      .replace(".le.", "<=")
      .replace(".ge.", ">=")
      .replace(".not.", "!")
      .replace(".and.", "&&")
      .replace(".or.", "||")
      .replace(".false.", "false")
      .replace(".true.", "true")
      .replace("logical", "bool")
      .replace("double precision", "double")
      .replace("integer", "int")
      .replace("go to", "goto")
      .replace("goto ", "goto lbl_")
      .replace("call ", "")
      )
    if (len(line) > 6 and line[5] != " "):
      line = line[:5] + " " + line[6:]
    if (line.endswith("then")):
      line = line[:-4] + "{"
    elif (line.strip() == "endif"):
      line = line.replace("endif", "}")
    elif (line[-1] != ","):
      line += ";"
  print line
