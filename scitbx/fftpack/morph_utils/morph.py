import re

def process_if(line):
  s = line.replace("THEN", "{")
  for rel_f, rel_c in ((".EQ.", "=="), (".NE.", "!="),
                       (".LT.", "<"),  (".LE.", "<="),
                       (".GT.", ">"),  (".GE.", ">=")):
    s = s.replace(rel_f, rel_c)
  return s

def process_do(line):
  s = line.replace(" ", "")
  m = re.match(r"DO([^\=]+)\=([^,]+),([^,]+),([^,]+)", s)
  if (m == None):
    m = re.match(r"DO([^\=]+)\=([^,]+),([^,]+)", s)
    assert m != None
    s = "for (std::size_t %s = %s; %s <= %s; %s++) {" % (
      m.group(1), m.group(2),
      m.group(1), m.group(3),
      m.group(1))
  else:
    s = "for (std::size_t %s = %s; %s <= %s; %s += %s) {" % (
      m.group(1), m.group(2),
      m.group(1), m.group(3),
      m.group(1), m.group(4))
  return re.sub("DO.*", s, line)

def process_WA(line):
  return re.sub(r"([ =+-]WA[1-4]?)\(([^\)]+)\)", r"\1[\2-1]", line)

def morph_line(fortran_line):
  stripped_line = fortran_line.strip()
  if   (stripped_line.startswith("IF")):
    s = fortran_line.replace("IF", "if");
    s = process_if(s)
  elif (stripped_line.startswith("ELSE IF")):
    s1 = re.sub("ELSE IF.*", "}", fortran_line)
    s2 = fortran_line.replace("ELSE IF", "else if");
    s2 = process_if(s2)
    s = s1 + "\n" + s2
  elif (stripped_line == "ELSE"):
    s1 = fortran_line.replace("ELSE", "}")
    s2 = fortran_line.replace("ELSE", "else {");
    s = s1 + "\n" + s2
  elif (stripped_line.startswith("DO ")):
    s = process_do(fortran_line)
  elif (stripped_line.startswith("CALL")):
    s = re.sub(r"CALL\s*", "", fortran_line)
  elif (stripped_line in ("ENDIF", "ENDDO")):
    s = re.sub("END..", "}", fortran_line)
  elif (stripped_line.startswith("DOUBLE PRECISION")):
    s = re.sub(r"DOUBLE PRECISION\s*", "FloatType ", fortran_line)
  elif (stripped_line.startswith("INTEGER")):
    s = re.sub(r"INTEGER\s*", "std::size_t ", fortran_line)
  else:
    s = fortran_line
  if (fortran_line[0] != " "): return s
  if (stripped_line in ("IMPLICIT NONE", "END")): return s
  if (stripped_line.startswith("SUBROUTINE")): return s
  if (s[-1] in "{}"): return s
  s = process_WA(s)
  return s.replace("RETURN", "return") + ";"

def morph_file(fortran_lines):
  out = []
  for l in fortran_lines:
    ul = l.upper()
    out.append(morph_line(ul[:-1]))
  return out

if (__name__ == "__main__"):
  import sys
  f = open(sys.argv[1], "r")
  fortran_lines = f.readlines()
  f.close()
  morphed_lines = morph_file(fortran_lines)
  for lines in morphed_lines:
    for l in lines.split("\n"):
      if (l[:2] == "  "):
        print l[2:]
      else:
        print l
