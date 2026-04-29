"""
Edit index.html files to simplify the path for display
"""
from __future__ import absolute_import, division, print_function
import os

def run(args):
  """
  Expect path to file and name of indexing directory relative to top level
  """
  assert len(args) == 3
  fn = args[0]
  indexing_dir = args[1]
  prefix = args[2]

  text = open(fn).read()
  text_to_find = get_text_to_replace(text)
  if text_to_find:
    ct = len(text)
    lines = text.splitlines()
    new_lines = []
    for line in lines:
      if line.find("<code") > -1:  # edit lines with code marking (not pre)
        line = line.replace(text_to_find,"")
      new_lines.append(line.rstrip())
    text = "\n".join(new_lines)
    new_ct = len(text)
    print("Text to find: '%s' (%s chars removed)" %(text_to_find, ct - new_ct))

  # Figure out path to top level from file name
  spl = fn.split(os.path.sep)
  paths_to_top_level = []
  for i in range(len(spl) - 1):
    paths_to_top_level.append("..")
  path_to_index = os.path.join(os.path.sep.join(paths_to_top_level),indexing_dir,"index.html")
  text = add_super_module_text(text, path_to_index, prefix)

  f = open(fn,'w')
  print(text, file = f)
  f.close()
  print("Wrote edited text to %s" %(fn))

def add_super_module_text(text, path_to_index, prefix):
  """
   Find the text:
   "<li><h3><a href="#header-submodules">Sub-modules</a></h3>"

   and replace it with:
   '''
<h3><a href="../../index_files/index.html">CCTBX API Index</a></h3>

<li><h3><a href="#header-supermodules">Super-module</a></h3>
<li><code><a title="cctbx_project" href="../index.html">cctbx_project</a></code></li>
</ul>

<ul id="index">
<li><h3><a href="#header-submodules">Sub-modules</a></h3>
   '''

  If Super-module is already present just add index above it
  """

  text1 = """<li><h3><a href="#header-submodules">Sub-modules</a></h3>"""
  text1a = """<li><h3>Super-module</h3>"""

  text2 = """
<h3><a href="%s ">%s API Index</a></h3>

<li><h3><a href="#header-supermodules">Super-module</a></h3>
<li><code><a title="cctbx_project" href="../index.html">cctbx_project</a></code></li>
</ul>

<ul id="index">
<li><h3><a href="#header-submodules">Sub-modules</a></h3>
""" %(path_to_index, prefix)

  text2a = """
<h3><a href="%s ">%s API Index</a></h3>

<li><h3>Super-module</h3>

""" %(path_to_index, prefix)


  if text.find(text1a) > -1:
     text = text.replace(text1a, text2a)
  elif text.find(text1) > -1:
     text = text.replace(text1, text2)
  return text

def get_text_to_replace(text):
  """Text to replace looks like iotbx.bioinformatics. , where this comes from:
    ">Module <code>iotbx.bioinformatics<"
    "iotbx.bioinformatics." To be replaced with "" everywhere
    Same for Package
    Hardwired for pdoc3 output html
  """
  spl = text.split(">Module <code>")
  if len(spl) < 2:
     spl = text.split(">Package <code>")
  if len(spl) < 2:
     return None
  if not "<" in spl[1]:
     return None
  spl2 = spl[1].split("<")[0]
  text_to_find = "%s." %(spl2)
  return text_to_find

if __name__ == "__main__":
  import sys
  run(sys.argv[1:])
