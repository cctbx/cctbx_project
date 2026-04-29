"""
Edit python file that contains boost instructions to make it useable
by pdoc3.

Replace code block like:

bp.inject(ext.root, __hash_eq_mixin)
@bp.inject_into(ext.root)
class _():

with

class root:

"""
from __future__ import absolute_import, division, print_function

def get_boost_names(text):
  """find code blocks like:

@bp.inject_into(ext.root)
class _():

and save "root"
  """

  key_text = "@bp.inject_into(ext."
  boost_names = []
  for line in text.splitlines():
    if line.rstrip().startswith(key_text):
      boost_names.append(line.replace(key_text,"").replace(")",""))
  return boost_names


def edit_boost_code_blocks(text, boost_names):
  """Edit the boost code blocks"""

  for boost_name in boost_names:
    search_text = """bp.inject(ext.%s, __hash_eq_mixin)
@bp.inject_into(ext.%s)
class _():""" %(boost_name, boost_name)

    search_text2 = """@bp.inject_into(ext.%s)
class _():""" %(boost_name)

    replacement_text ="class %s:" %(boost_name)
    start_text = str(text)
    text = text.replace(search_text, replacement_text)

    if text == start_text: # try alternate version
      text = text.replace(search_text2, replacement_text)

  return text

def modify_boost_text(text):
  boost_names = get_boost_names(text)
  print("Editing the following boost_names: %s" %(" ".join(boost_names)))
  text = edit_boost_code_blocks(text, boost_names)
  return text

def run(args):
  """
  Expect path to file
  """
  assert len(args) == 1
  fn = args[0]
  text = open(fn).read()
  text = modify_boost_text(text)
  f = open(fn, 'w')
  print(text, file = f)
  f.close()
  print("Wrote edited text to %s" %(fn))

if __name__ == "__main__":
  import sys
  run(sys.argv[1:])
