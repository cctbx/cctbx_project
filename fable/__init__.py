import boost.python
ext = boost.python.import_ext("fable_ext")
from fable_ext import *

class SemanticError(Exception): pass
