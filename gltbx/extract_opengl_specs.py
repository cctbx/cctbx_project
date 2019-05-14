"""
To update opengl_specs.txt:

  1. In an empty temporary directory run:
       wget_opengl_specs.csh

  2. In $GLTBX_DIST run:
       python extract_opengl_specs.py /tmpdir/html/*/*.html > opengl_specs.txt
"""
from __future__ import absolute_import, division, print_function

import sys, os

sequential_defines = ["GL_AUX", "GL_CLIP_PLANE", "GL_LIGHT"]

def extract_defines(html_string, all_defines):
  for word in html_string.split():
    if (   word.startswith("<STRONG>GL_")
        or word.startswith("<STRONG>GLU_")):
      i = word.find("</STRONG>")
      if (i < 0): continue
      define = word[8:i]
      if (define.endswith("_")): continue
      if (define.upper() != define): continue
      if (not define.replace("_","").isalnum()): continue
      keep = True
      for sequential_define in sequential_defines:
        if (define == sequential_define):
          keep = False
          break
        elif (define.startswith(sequential_define)):
          num = define[len(sequential_define):]
          try: num = int(num)
          except ValueError: pass
          else:
            keep = False
            break
      if (not keep): continue
      all_defines.add(define)

def extract_signatures(html_string, all_signatures):
  signature_block = []
  active_block = False
  current_line = None
  c_specification = "<STRONG>C</STRONG> <STRONG>SPECIFICATION</STRONG>"
  for line in html_string.splitlines():
    if (line.strip() == c_specification):
      active_block = True
    elif (line.strip() in [
            "<STRONG>PARAMETERS</STRONG>",
            "<STRONG>DESCRIPTION</STRONG>"]):
      active_block = False
      if (current_line is not None):
        current_line = current_line.strip()
        if (len(current_line) > 0):
          all_signatures.append(current_line)
        current_line = None
    elif (active_block):
      line = line.expandtabs()
      line = line.replace("<STRONG>", "").replace("</STRONG>", "")
      line = line.replace("<EM>", "").replace("</EM>", "")
      line = line.replace("GLvoid (*CallBackFunc)(", "glu_function_pointer fn")
      line = line.replace("(", " ( ")
      line = line.replace(")", " ) ")
      line = " ".join(line.split())
      if (current_line is None):
        current_line = line
      elif (current_line.endswith(",")):
        current_line += " " + line
      else:
        current_line = current_line.strip()
        if (len(current_line) > 0):
          all_signatures.append(current_line)
        current_line = line

def run(args):
  all_defines = set()
  for arg in args:
    extract_defines(html_string=open(arg).read(), all_defines=all_defines)
  for define in sorted(all_defines):
    print(define)
  #
  all_signatures = []
  for arg in args:
    if (os.path.basename(arg).lower() == "index.html"): continue
    prev_len = len(all_signatures)
    extract_signatures(
      html_string=open(arg).read(),
      all_signatures=all_signatures)
    assert len(all_signatures) > prev_len
  for signature in all_signatures:
    print(signature)

if (__name__ == "__main__"):
  run(sys.argv[1:])
