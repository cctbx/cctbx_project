import sys, re

def process_group(groups):
  result = ""
  for elem in groups:
    if (elem == "I"):
      elem = "I1"
    elif (elem == "I-1"):
      elem = "I0"
    else:
      try: i = int(elem)
      except: pass
      else:
        elem = str(i-1)
    result += "," + elem
  return "(" + result[1:] + ")"

def offs0(file):
  ix3d = re.compile(r"\(([0-9A-Za-z-]+),([0-9A-Za-z-]+),([0-9A-Za-z-]+)\)")
  ix2d = re.compile(r"\(([0-9A-Za-z-]+),([0-9A-Za-z-]+)\)")
  for line in file:
    line = line[:-1]
    next_start = 0
    mod_line = ""
    while 1:
      m = ix3d.search(line, next_start)
      if (m == None): m = ix2d.search(line, next_start)
      if (m == None): break
      repl = process_group(m.groups())
      #repl = m.group()
      mod_line += line[next_start:m.start()] + repl
      next_start = m.end()
    mod_line += line[next_start:]
    print mod_line

if (__name__ == "__main__"):
  file = sys.stdin.readlines()
  offs0(file)
