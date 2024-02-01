from __future__ import absolute_import, division, print_function

from builtins import object

import ast

class imported_name(object):

  __slots__ = ('name', 'lineno')

  def __init__(self, name, lineno):
    self.name = tuple(name.split('.'))
    self.lineno = lineno

  def __eq__(self, other):
    return (self.name == other.name
            and self.lineno == other.lineno)

  def __ne__(self, other):
    return not self.__eq__(other)

  def __hash__(self):
    return hash((self.name, self.lineno))

  def __repr__(self):
    return '%s imported at line %i' % (self.name_as_str(), self.lineno)

  def name_as_str(self):
    return '.'.join(self.name)


class unused_imports(ast.NodeVisitor, object):
  """ Unused import's in a module.
  This finds more of them than the algorithm in
  libtbx.find_unused_imports_crude
  Work in progress: careful checking of the outcome mandatory
  """

  @classmethod
  def is_subpath_of(cls, a, b):
    """ Whether a is a subpath of b """
    if len(a) > len(b): return False
    return a == b[:len(a)]

  def __init__(self, python_source_code=None, python_source_filename=None,
               ignored_imports=(), ignored_imports_from=(),
               ignore_imports_flagged_by_comments=()):
    assert (python_source_code, python_source_filename).count(None) == 1
    if python_source_code is None:
      python_source_code = open(python_source_filename).read()
    super(unused_imports, self).__init__()
    self.comment_flags = ignore_imports_flagged_by_comments
    self.python_source_line = python_source_code.splitlines()
    self.ignored_imports = ignored_imports
    self.ignored_imports_from = ignored_imports_from
    self.current_context = () # start at module level
    self.imported_in_context = {}
    self.imported_from_full_name_in_context = {}
    self.used_in_context = {}
    tree = ast.parse(python_source_code)
    self._used = set()
    self.visit(tree)
    self._unused = set()
    for imported in self.imported_in_context.values():
      self._unused.update(imported)
    self._unused -= self._used

  def __repr__(self):
    return '\n'.join( str(imp) for imp in self )

  def __iter__(self):
    return iter(self._unused)

  def __bool__(self):
    return bool(self._unused)

  @property
  def names(self):
    return set( imp.name_as_str() for imp in self )

  def visit_Import(self, imp):
    for comment in self.comment_flags:
      if self.python_source_line[imp.lineno - 1].endswith(comment): return
    imported = self.imported_in_context.setdefault(self.current_context, set())
    imported.update(
      imported_name(name.name if name.asname is None else name.asname,
                    imp.lineno)
      for name in imp.names if name.name not in self.ignored_imports)
    self.consolidate_imports_info()

  def visit_ImportFrom(self, imp):
    for comment in self.comment_flags:
      if self.python_source_line[imp.lineno - 1].endswith(comment): return
    imported = self.imported_in_context.setdefault(self.current_context, set())
    imported.update(
      imported_name(name.name if name.asname is None else name.asname,
                    imp.lineno)
      for name in imp.names
      if name.name != '*' and imp.module not in self.ignored_imports_from)
    imported = self.imported_from_full_name_in_context.setdefault(
      self.current_context, set())
    imported.update( tuple(imp.module.split('.')) + (name.name,)
                     for name in imp.names if name.name != '*')
    self.consolidate_imports_info()

  def consolidate_imports_info(self):
    for ctx, imported in self.imported_in_context.items():
      discarded = set()
      for imp1 in self.imported_from_full_name_in_context.get(ctx, set()):
        for imp in imported:
          if self.is_subpath_of(imp1, imp.name):
            discarded.add(imp)
      imported -= discarded

  def _process_namespace(self, namespace, lineno):
    for import_ctx, imported in self.imported_in_context.items():
      if self.is_subpath_of(import_ctx, self.current_context):
        for imp in imported:
          if lineno < imp.lineno: continue
          if not self.is_subpath_of(imp.name, namespace): continue
          self._used.add(imp)

  def visit_Name(self, name):
    self._process_namespace((name.id,), name.lineno)

  def visit_Attribute(self, attr):
    namespace = []
    x = attr.value
    while not isinstance(x, ast.Name):
      if isinstance(x, ast.Attribute):
        namespace.append(x.attr)
        x = x.value
      elif isinstance(x, ast.Call):
        for arg in x.args:
          self.visit(arg)
        for keyword in x.keywords:
          self.visit(keyword)
        namespace = []
        x = x.func
      elif isinstance(x, ast.Subscript):
        self.visit(x.slice)
        namespace = []
        x = x.value
      else:
        return
    namespace.append(x.id) # got an instance of Name now
    namespace = tuple(reversed(namespace))
    self._process_namespace(namespace, attr.lineno)
    self._process_namespace(namespace + (attr.attr,), attr.lineno)

  def visit_FunctionDef(self, func):
    self.current_context += (func.name,)
    for default in func.args.defaults:
      self.visit(default)
    for stmt in func.body: self.visit(stmt)
    self.current_context = self.current_context[:-1]



class old_style_class(object):

  __slots__ = ('name', 'lineno')

  def __init__(self, name, lineno):
    self.name = tuple(name.split('.'))
    self.lineno = lineno

  def __eq__(self, other):
    return (self.name == other.name
            and self.lineno == other.lineno)

  def __ne__(self, other):
    return not self.__eq__(other)

  def __hash__(self):
    return hash((self.name, self.lineno))

  def __repr__(self):
    return 'class %s defined at line %i' % (self.name_as_str(), self.lineno)

  def name_as_str(self):
    return '.'.join(self.name)


class find_old_style_classes(ast.NodeVisitor, object):
  """ Finds old-style classes (i.e. ones that don't inherit from object)
  """

  @classmethod
  def is_subpath_of(cls, a, b):
    """ Whether a is a subpath of b """
    if len(a) > len(b): return False
    return a == b[:len(a)]

  def __init__(self, python_source_code=None, python_source_filename=None):
    assert (python_source_code, python_source_filename).count(None) == 1
    if python_source_code is None:
      python_source_code = file(python_source_filename).read()
    super(find_old_style_classes, self).__init__()
    self.python_source_line = python_source_code.splitlines()
    self.current_context = () # start at module level
    tree = ast.parse(python_source_code)
    self._old_style_classes = set()
    self.visit(tree)

  @property
  def names(self):
    return set( imp.name_as_str() for imp in self )

  def __repr__(self):
    return '\n'.join( str(imp) for imp in self )

  def __iter__(self):
    return iter(self._old_style_classes)

  def __bool__(self):
    return bool(self._old_style_classes)

  def visit_ClassDef(self, node):
    #print node.name, [n.id for n in node.bases]
    if len(node.bases) == 0:
      self._old_style_classes.add(old_style_class(node.name, node.lineno))
      #print node.name
