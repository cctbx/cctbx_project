from __future__ import absolute_import, division, print_function
from io import open
from html.parser import HTMLParser
import os
from libtbx.utils import to_unicode

# ------------------------------------------------------------------------------

class MyHTMLParser(HTMLParser):
  '''
  This is a customized class to extract code from an html file

  Anything that is between <pre></pre> tags will be considered code.
  So only code in the html file that is to be executed in the script
  should be between these tags
  '''

  def __init__(self):
    HTMLParser.__init__(self)
    self.is_code = False
    #self.code_str = list()
    self.code_str = ''

  def handle_starttag(self, tag, attrs):
    #print("Encountered a start tag:", tag)
    if tag == 'pre' and attrs == [(u'class', u'codeDL')]:
      self.is_code = True

  def handle_endtag(self, tag):
    #print("Encountered an end tag :", tag)
    if tag == 'pre':
      self.is_code = False

  def handle_data(self, data):
    #print("Encountered some data  :", data)
    if self.is_code == True:
      self.code_str = self.code_str + data
      #self.code_str.extend(data.strip().splitlines())

  def return_result(self):
    return self.code_str

# ------------------------------------------------------------------------------

def run(parent_dir):
  '''
  Walk through directory cctbx_project/cctbx_website/

  Extract code snippets from all .html files in the directory and store them
  in a script with the same basename

  For example:
  template.html --> template.py

  The file template.py will be tested as part of cctbx tests.
  '''
  # this is one directory up: /cctbx_project/cctbx_website/
  #parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

  msg = '''#********************************************************************
# This script is automatically generated when running libtbx.refresh (or bootstrap.py)
# It is not part of the GitHub repository
# So if this script is manually changed, the changes will be lost when updating
#********************************************************************\n'''

  directory = os.path.join(parent_dir, 'html_files')

  for filename in os.listdir(directory):
    # look at all .html files
    if filename.endswith(".html"):
      fn = os.path.join(directory, filename)
      with open(fn, 'r', encoding='utf-8') as html_file:
        # get code in html file
        data = html_file.read()
      parser = MyHTMLParser()
      parser.feed(data)
      code_str = parser.return_result()
      # get filename
      base = os.path.splitext(filename)[0]
      # save code in script and put it in folder cctbx_website/examples/
      dest_dir = os.path.join(parent_dir, 'examples')
      if (not os.path.isdir(dest_dir)):
        os.makedirs(dest_dir)
      script_filename = os.path.join(dest_dir, base+'.py')
      with open(script_filename, 'w', encoding='utf-8') as file:
        file.write(to_unicode('from __future__ import absolute_import, division, print_function\n'))
        file.write(to_unicode(msg))
        file.write(to_unicode(code_str.strip()))
        file.write(to_unicode('\n'))
    else:
      continue


if __name__ == '__main__':
  import libtbx.load_env
  for parent_dir in [libtbx.env.dist_path("cctbx_website", default=None),
    libtbx.env.dist_path("phenix_dev_doc", default=None)]:
    if parent_dir is not None:
      run(parent_dir)
