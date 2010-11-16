""" Tools to statically analyse source code in various languages """

import os
import re
from glob import glob


class comments(object):
  """ Information about comment lines """

  def __init__(self, files, debug=False):
    """ Construct the information for the given files """
    self.debug = debug
    self.commented_lines = 0
    self.lines = 0
    for f in files:
      for p in glob(f):
        p = os.path.expanduser(p)
        root, ext = os.path.splitext(p)
        ext = ext.lower()
        if ext in ('.f', '.fpp', '.inc', '.f77', '.f90'):
          self.process_fortran(p)
        else:
          raise NotImplemented(ext)

  fortran_start_of_line = re.compile(r'^ (?: C | \* | !) \W* \w', re.X)
  fortran_bang_anywhere = re.compile(r'! \W* \w', re.X)

  def process_fortran(self, file_path):
    n = 0
    for line in open(file_path):
      self.lines += 1
      if self.fortran_start_of_line.search(line):
        n += 1
        if self.debug: print self.lines, line
        continue
      bang = line.rfind('!')
      if bang < 0: continue
      single_quote, double_quote = line.rfind("'"), line.rfind('"')
      if single_quote >= 0 or double_quote >= 0:
        if single_quote >= 0 and double_quote >= 0:
          quote = max(single_quote, double_quote)
        elif single_quote >= 0:
          quote = single_quote
        else:
          quote = double_quote
        if quote < bang: continue
      if self.fortran_bang_anywhere.search(line, pos=bang):
        n += 1
        if self.debug: print self.lines, line
    self.commented_lines += n
