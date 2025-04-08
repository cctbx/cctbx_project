"""
Edit index.html files to simplify the path for display
"""
from __future__ import absolute_import, division, print_function

def run(args):
  assert len(args) == 1
  fn = args[0]
  text = open(fn).read()
  text_to_find = get_text_to_replace(text)
  if text_to_find:
    print("Text to find: '%s'" %(text_to_find))
    text = text.replace(text_to_find,"")
  text = add_super_module_text(text)

  f = open(fn,'w')
  print(text, file = f)
  f.close()
  print("Wrote edited text to %s" %(fn))

def add_super_module_text(text):
  """
   Find the text: 
   "<li><h3><a href="#header-submodules">Sub-modules</a></h3>"

   and replace it with:
   '''
<li><h3><a href="#header-supermodules">Super-module</a></h3>
<li><code><a title="cctbx_project" href="../index.html">cctbx_project</a></code></li>
</ul>

<ul id="index">
<li><h3><a href="#header-submodules">Sub-modules</a></h3>
   '''
  """

  text1 = """<li><h3><a href="#header-submodules">Sub-modules</a></h3>"""
  text2 = """<li><h3><a href="#header-supermodules">Super-module</a></h3>
<li><code><a title="cctbx_project" href="../index.html">cctbx_project</a></code></li>
</ul>

<ul id="index">
<li><h3><a href="#header-submodules">Sub-modules</a></h3>
"""
  if text.find(text1) > -1:
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
