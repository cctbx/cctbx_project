from __future__ import absolute_import, division, print_function
import boost.python
ext = boost.python.import_ext("prime_ext")
from prime_ext import *
