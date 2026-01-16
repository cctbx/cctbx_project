# -*- coding: utf-8 -*-

from __future__ import division, print_function
from pdoc import Doc, Context, _pep224_docstrings, _is_public
from pdoc import _is_whitelisted, _is_blacklisted, _is_function, Class
from pdoc import Function, link_inheritance, External, _render_template
from pdoc import Variable, _filter_type

from pathlib import Path

################# MODIFIED CODE FROM pdoc3/__init__.py #################
"""
Python package `pdoc` provides types, functions, and a command-line
interface for accessing public documentation of Python modules, and
for presenting it in a user-friendly, industry-standard open format.
It is best suited for small- to medium-sized projects with tidy,
hierarchical APIs.

.. include:: ./documentation.md
"""
import importlib
import inspect
import os
import os.path as path
import re
import sys
import typing
from contextlib import contextmanager
from copy import copy
from functools import lru_cache
from types import ModuleType
from typing import (  # noqa: F401
    cast, Any, Callable, Dict, Generator, Iterable, List, Literal, Mapping, NewType,
    Optional, Set, Tuple, Type, TypeVar, Union,
)
from warnings import warn

from mako.lookup import TemplateLookup

try:
    from pdoc._version import version as __version__  # noqa: F401
except ImportError:
    __version__ = '???'  # Package not installed


class dummy_module:
  def __init__():
    pass

_get_type_hints = lru_cache()(typing.get_type_hints)

_URL_MODULE_SUFFIX = '.html'
_URL_INDEX_MODULE_SUFFIX = '.m.html'  # For modules named literal 'index'
_URL_PACKAGE_SUFFIX = '/index.html'

# type.__module__ can be None by the Python spec. In those cases, use this value
_UNKNOWN_MODULE = '?'

T = TypeVar('T', 'Module', 'Class', 'Function', 'Variable')

__pdoc__: Dict[str, Union[bool, str]] = {}

tpl_lookup = TemplateLookup(
    cache_args=dict(cached=True,
                    cache_type='memory'),
    input_encoding='utf-8',
    directories=[path.join(path.dirname(__file__), "templates")],
)
"""
A `mako.lookup.TemplateLookup` object that knows how to load templates
from the file system. You may add additional paths by modifying the
object's `directories` attribute.
"""
if os.getenv("XDG_CONFIG_HOME"):
    tpl_lookup.directories.insert(0, path.join(os.getenv("XDG_CONFIG_HOME", ''), "pdoc"))


def import_module(module: Union[str, ModuleType],
                  *, reload: bool = False,
                   convert_comments_to_docstring: bool = True,
                   files_to_edit_for_boost: list[str] = None) -> ModuleType:
    """
    Return module object matching `module` specification (either a python
    module path or a filesystem path to file/directory).
    """
    @contextmanager
    def _module_path(module):
        from os.path import abspath, dirname, isfile, isdir, split
        path = '_pdoc_dummy_nonexistent'
        module_name = inspect.getmodulename(module)
        if isdir(module):
            path, module = split(abspath(module))
        elif isfile(module) and module_name:
            path, module = dirname(abspath(module)), module_name
        try:
            sys.path.insert(0, path)
            yield module
        finally:
            sys.path.remove(path)

    if isinstance(module, Module):
        module = module.obj
    if isinstance(module, str):
        with _module_path(module) as module_path:
            original_module = module
            got_it = False
            if convert_comments_to_docstring:
              # Convert comments at top of classes and methods to docstrings
              #   if no docstrings present
              try:
                # module = importlib.import_module(module_path)
                fn = Path(module)
                module = import_and_edit(module_path, fn)
                got_it = True
              except Exception as e:
                got_it = False

            if not got_it: # try without editing
              try:
                module = importlib.import_module(module_path)
              except Exception as e:
                print("FAILED TO IMPORT ",original_module,module_path,"::",str(e),"::")
                from copy import deepcopy
                module = deepcopy(dummy_module) # skipping it and marking
                module.__name__= "%s_dummy_module" %(original_module)
                print("FAILED_TO_IMPORT_MODULE:",original_module)
                failed_file = 'pdoc.failed'
                if os.path.isfile(failed_file):
                  f = open(failed_file,'a')
                else:
                  f = open(failed_file,'w')
                print(original_module, file = f)
                f.close()

    assert inspect.ismodule(module)
    # If this is pdoc itself, return without reloading. Otherwise later
    # `isinstance(..., pdoc.Doc)` calls won't work correctly.
    if reload and not module.__name__.startswith(__name__):
        module = importlib.reload(module)
        # We recursively reload all submodules, in case __all_ is used - cf. issue #264
        for mod_key, mod in list(sys.modules.items()):
            if mod_key.startswith(module.__name__):
                importlib.reload(mod)
    return module


class Module(Doc):
    """
    Representation of a module's documentation.
    """
    __pdoc__["Module.name"] = """
        The name of this module with respect to the context/path in which
        it was imported from. It is always an absolute import path.
        """

    __slots__ = ('supermodule', 'doc', '_context', '_is_inheritance_linked',
                 '_skipped_submodules')

    def __init__(self, module: Union[ModuleType, str], *,
                 docfilter: Optional[Callable[[Doc], bool]] = None,
                 supermodule: Optional['Module'] = None,
                 context: Optional[Context] = None,
                 skip_errors: bool = False,
                 convert_comments_to_docstring: bool = True,
                 files_to_edit_for_boost: list[str] = None):
        """
        Creates a `Module` documentation object given the actual
        module Python object.

        `docfilter` is an optional predicate that controls which
        sub-objects are documentated (see also: `pdoc.html()`).

        `supermodule` is the parent `pdoc.Module` this module is
        a submodule of.

        `context` is an instance of `pdoc.Context`. If `None` a
        global context object will be used.

        If `skip_errors` is `True` and an unimportable, erroneous
        submodule is encountered, a warning will be issued instead
        of raising an exception.
        """
        if isinstance(module, str):
          original_module = module
          try:
            module = import_module(module,
              convert_comments_to_docstring = convert_comments_to_docstring,
              files_to_edit_for_boost = files_to_edit_for_boost)
          except Exception as e:
            from copy import deepcopy
            module = deepcopy(dummy_module) # skipping it and marking name
            module.__name__= "%s_dummy_module" %(original_module)
            print("FAILED_IMPORT_MODULE:",original_module)
            failed_file = 'pdoc.failed'
            if os.path.isfile(failed_file):
              f = open(failed_file,'a')
            else:
              f = open(failed_file,'w')
            print(original_module, file = f)
            f.close()

        super().__init__(module.__name__, self, module)
        if self.name.endswith('.__init__') and not self.is_package:
            self.name = self.name[:-len('.__init__')]

        self._context = _global_context if context is None else context
        """
        A lookup table for ALL doc objects of all modules that share this context,
        mainly used in `Module.find_ident()`.
        """
        #assert isinstance(self._context, Context), \ 'pdoc.Module(context=) should be a pdoc.Context instance'

        self.supermodule = supermodule
        """
        The parent `pdoc.Module` this module is a submodule of, or `None`.
        """

        self.doc: Dict[str, Union[Module, Class, Function, Variable]] = {}
        """A mapping from identifier name to a documentation object."""

        self._is_inheritance_linked = False
        """Re-entry guard for `pdoc.Module._link_inheritance()`."""

        self._skipped_submodules = set()

        var_docstrings, _ = _pep224_docstrings(self)

        # Populate self.doc with this module's public members
        public_objs = []
        if hasattr(self.obj, '__all__'):
            for name in self.obj.__all__:
                try:
                    obj = getattr(self.obj, name)
                except AttributeError:
                    warn(f"Module {self.module!r} doesn't contain identifier `{name}` "
                         "exported in `__all__`")
                else:
                    if not _is_blacklisted(name, self):
                        obj = inspect.unwrap(obj)
                    public_objs.append((name, obj))
        else:
            def is_from_this_module(obj):
                mod = inspect.getmodule(inspect.unwrap(obj))
                return mod is None or mod.__name__ == self.obj.__name__

            for name, obj in inspect.getmembers(self.obj):
                if ((_is_public(name) or
                     _is_whitelisted(name, self)) and
                        (_is_blacklisted(name, self) or  # skips unwrapping that follows
                         is_from_this_module(obj) or
                         name in var_docstrings)):

                    if _is_blacklisted(name, self):
                        self._context.blacklisted.add(f'{self.refname}.{name}')
                        continue

                    obj = inspect.unwrap(obj)
                    public_objs.append((name, obj))

            index = list(self.obj.__dict__).index
            public_objs.sort(key=lambda i: index(i[0]))

        for name, obj in public_objs:
            if _is_function(obj):
                self.doc[name] = Function(name, self, obj)
            elif inspect.isclass(obj):
                self.doc[name] = Class(name, self, obj)
            elif name in var_docstrings:
                self.doc[name] = Variable(name, self, var_docstrings[name], obj=obj)

        # If the module is a package, scan the directory for submodules
        if self.is_package:

            def iter_modules(paths):
                """
                Custom implementation of `pkgutil.iter_modules()`
                because that one doesn't play well with namespace packages.
                See: https://github.com/pypa/setuptools/issues/83
                """
                from os.path import isdir, join
                for pth in paths:
                    if pth.startswith("__editable__."):
                        # See https://github.com/pypa/pip/issues/11380
                        continue
                    for file in os.listdir(pth):
                        if file.startswith(('.', '__pycache__', '__init__.py')):
                            continue
                        module_name = inspect.getmodulename(file)
                        if module_name:
                            yield module_name
                        if isdir(join(pth, file)) and '.' not in file:
                            yield file

            for root in iter_modules(self.obj.__path__):
                # Ignore if this module was already doc'd.
                if root in self.doc:
                    continue

                # Ignore if it isn't exported
                if not _is_public(root) and not _is_whitelisted(root, self):
                    continue
                if _is_blacklisted(root, self):
                    self._skipped_submodules.add(root)
                    continue

                assert self.refname == self.name
                fullname = f"{self.name}.{root}"
                try:
                    m = Module(import_module(fullname,
                                 convert_comments_to_docstring =
                                   convert_comments_to_docstring),
                               docfilter=docfilter, supermodule=self,
                               context=self._context, skip_errors=skip_errors,
                               convert_comments_to_docstring =
                              convert_comments_to_docstring,
                               files_to_edit_for_boost =
                                 files_to_edit_for_boost)
                except Exception as ex:
                    if skip_errors:
                        warn(str(ex), Module.ImportWarning)
                        continue
                    raise

                self.doc[root] = m
                # Skip empty namespace packages because they may
                # as well be other auxiliary directories
                if m.is_namespace and not m.doc:
                    del self.doc[root]
                    self._context.pop(m.refname)

        # Apply docfilter
        if docfilter:
            for name, dobj in self.doc.copy().items():
                if not docfilter(dobj):
                    self.doc.pop(name)
                    self._context.pop(dobj.refname, None)

        # Build the reference name dictionary of the module
        self._context[self.refname] = self
        for docobj in self.doc.values():
            self._context[docobj.refname] = docobj
            if isinstance(docobj, Class):
                self._context.update((obj.refname, obj)
                                     for obj in docobj.doc.values())

    class ImportWarning(UserWarning):
        """
        Our custom import warning because the builtin is ignored by default.
        https://docs.python.org/3/library/warnings.html#default-warning-filter
        """

    __pdoc__['Module.ImportWarning'] = False

    @property
    def __pdoc__(self) -> dict:
        """This module's __pdoc__ dict, or an empty dict if none."""
        return getattr(self.obj, '__pdoc__', {})

    def _link_inheritance(self):
        # Inherited members are already in place since
        # `Class._fill_inheritance()` has been called from
        # `pdoc.fill_inheritance()`.
        # Now look for docstrings in the module's __pdoc__ override.

        if self._is_inheritance_linked:
            # Prevent re-linking inheritance for modules which have already
            # had done so. Otherwise, this would raise "does not exist"
            # errors if `pdoc.link_inheritance()` is called multiple times.
            return

        # Apply __pdoc__ overrides
        for name, docstring in self.__pdoc__.items():
            # In case of whitelisting with "True", there's nothing to do
            if docstring is True:
                continue

            refname = f"{self.refname}.{name}"
            if docstring in (False, None):
                if docstring is None:
                    warn('Setting `__pdoc__[key] = None` is deprecated; '
                         'use `__pdoc__[key] = False` '
                         f'(key: {name!r}, module: {self.name!r}).')

                if name in self._skipped_submodules:
                    continue

                if (not name.endswith('.__init__') and
                        name not in self.doc and
                        refname not in self._context and
                        refname not in self._context.blacklisted):
                    warn(f'__pdoc__-overriden key {name!r} does not exist '
                         f'in module {self.name!r}')

                obj = self.find_ident(name)
                cls = getattr(obj, 'cls', None)
                if cls:
                    del cls.doc[obj.name]
                self.doc.pop(name, None)
                self._context.pop(refname, None)

                # Pop also all that startwith refname
                for key in list(self._context.keys()):
                    if key.startswith(refname + '.'):
                        del self._context[key]

                continue

            dobj = self.find_ident(refname)
            if isinstance(dobj, External):
                continue
            if not isinstance(docstring, str):
                raise ValueError('__pdoc__ dict values must be strings; '
                                 f'__pdoc__[{name!r}] is of type {type(docstring)}')
            dobj.docstring = inspect.cleandoc(docstring)

        # Now after docstrings are set correctly, continue the
        # inheritance routine, marking members inherited or not
        for c in _filter_type(Class, self.doc):
            c._link_inheritance()

        self._is_inheritance_linked = True

    def text(self, **kwargs) -> str:
        """
        Returns the documentation for this module as plain text.
        """
        txt = _render_template('/text.mako', module=self, **kwargs)
        return re.sub("\n\n\n+", "\n\n", txt)

    def html(self, minify=True, **kwargs) -> str:
        """
        Returns the documentation for this module as
        self-contained HTML.

        If `minify` is `True`, the resulting HTML is minified.

        For explanation of other arguments, see `pdoc.html()`.

        `kwargs` is passed to the `mako` render function.
        """
        html = _render_template('/html.mako', module=self, **kwargs)
        if minify:
            from pdoc.html_helpers import minify_html
            html = minify_html(html)
        if not html.endswith('\n'):
            html = html + '\n'
        return html

    @property
    def is_package(self) -> bool:
        """
        `True` if this module is a package.

        Works by checking whether the module has a `__path__` attribute.
        """
        return hasattr(self.obj, "__path__")

    @property
    def is_namespace(self) -> bool:
        """
        `True` if this module is a namespace package.
        """
        try:
            return self.obj.__spec__.origin in (None, 'namespace')  # None in Py3.7+
        except AttributeError:
            return False

    def find_class(self, cls: type) -> Doc:
        """
        Given a Python `cls` object, try to find it in this module
        or in any of the exported identifiers of the submodules.
        """
        # XXX: Is this corrent? Does it always match
        # `Class.module.name + Class.qualname`?. Especially now?
        # If not, see what was here before.
        return self.find_ident(f'{cls.__module__ or _UNKNOWN_MODULE}.{cls.__qualname__}')

    def find_ident(self, name: str) -> Doc:
        """
        Searches this module and **all** other public modules
        for an identifier with name `name` in its list of
        exported identifiers.

        The documentation object corresponding to the identifier is
        returned. If one cannot be found, then an instance of
        `External` is returned populated with the given identifier.
        """
        _name = name.rstrip('()')  # Function specified with parentheses

        if _name.endswith('.__init__'):  # Ref to class' init is ref to class itself
            _name = _name[:-len('.__init__')]

        return (self.doc.get(_name) or
                self._context.get(_name) or
                self._context.get(f'{self.name}.{_name}') or
                External(name))

    def _filter_doc_objs(self, type: Type[T], sort=True) -> List[T]:
        result = _filter_type(type, self.doc)
        return sorted(result) if sort else result

    def variables(self, sort=True) -> List['Variable']:
        """
        Returns all documented module-level variables in the module,
        optionally sorted alphabetically, as a list of `pdoc.Variable`.
        """
        return self._filter_doc_objs(Variable, sort)

    def classes(self, sort=True) -> List['Class']:
        """
        Returns all documented module-level classes in the module,
        optionally sorted alphabetically, as a list of `pdoc.Class`.
        """
        return self._filter_doc_objs(Class, sort)

    def functions(self, sort=True) -> List['Function']:
        """
        Returns all documented module-level functions in the module,
        optionally sorted alphabetically, as a list of `pdoc.Function`.
        """
        return self._filter_doc_objs(Function, sort)

    def submodules(self) -> List['Module']:
        """
        Returns all documented sub-modules of the module sorted
        alphabetically as a list of `pdoc.Module`.
        """
        return self._filter_doc_objs(Module)

    def _url(self):
        url = self.module.name.replace('.', '/')
        if self.is_package:
            return url + _URL_PACKAGE_SUFFIX
        elif url.endswith('/index'):
            return url + _URL_INDEX_MODULE_SUFFIX
        return url + _URL_MODULE_SUFFIX
################# END MODIFIED CODE FROM pdoc3/__init__.py #################

def run(args, top_level = None):

  # Catch failed imports in 'pdoc.failed'
  failed_file = 'pdoc.failed'
  if os.path.isfile (failed_file):
     try:
       os.remove(failed_file)
     except Exception as e:
       pass # another job removed it...

  # Choose whether to preprocess all files and convert
  #  comments to docstrings where docstrings are missing

  if 'convert_comments_to_docstring' in args:
    convert_comments_to_docstring = True
    args.remove('convert_comments_to_docstring')
  else:
    convert_comments_to_docstring = False

  # Identify any files that will need to be edited for boost:
  files_to_edit_for_boost = []
  key = "files_to_edit_for_boost="
  for arg in args:
    if arg.startswith(key):
      files_to_edit_for_boost = arg.replace(key,"").split()
      args.remove(arg)
      break

  modules = args
  context = Context()

  new_modules = []
  for mod in modules:
    x = Module(mod, context=context, skip_errors = True,
      convert_comments_to_docstring = convert_comments_to_docstring,
      files_to_edit_for_boost = files_to_edit_for_boost)
    new_modules.append(x)
  modules = new_modules
  if not modules:
    print("Nothing to do")
    return

  link_inheritance(context)

  def recursive_htmls(mod):
      yield mod, mod.name, mod.html()
      for submod in mod.submodules():
          yield from recursive_htmls(submod)
  def top_level_recursive_htmls(mod):
      yield mod, mod.name, mod.html()
      for submod in mod.submodules():
          yield submod, submod.name, submod.html()

  if top_level:
    get_htmls = top_level_recursive_htmls
    print("\nUsing top-level only")
  else:
    get_htmls = recursive_htmls
    print("\nUsing all levels")

  if convert_comments_to_docstring:
    print("Converting comments at top of classes and functions to doc strings")
  if files_to_edit_for_boost:
    print("Files to edit for boost: %s" %(" ".join(files_to_edit_for_boost)))


  ok_modules = []
  failed_modules = []

  for mod in modules:
      for m,  module_name, html in get_htmls(mod):
          if module_name.find('_dummy_module')> -1:
             failed_modules.append(module_name.replace('_dummy_module',''))
             continue
          else:
             ok_modules.append(module_name)
          # Determine if this is a directory or a file
          path = m.url()
          # Edit html to remove eg "qttbx.command_line." everywhere
          #remove_text = "%s." %(module_name)
          #html = html.replace(remove_text,"")
          working_path = "."
          for p in path.split(os.path.sep)[:-1]:
            working_path = os.path.join(working_path, p)
            if not os.path.isdir(working_path):
              os.mkdir(working_path)
          f = open('%s' %(path), 'w')
          print(html, file = f)
          f.close()
          print("MODULE: %s lines: %s" %(
              path, len(html.splitlines())))

  if os.path.isfile(failed_file):
    for x in open(failed_file).read().split():
      if not x in failed_modules:
        failed_modules.append(x)

  print("List of failed modules:")
  for m in failed_modules:
    print(m)
  print("\nTotal of %s ok modules and %s failed modules" %(
      len(ok_modules), len(failed_modules)))

def edit_file_for_boost(modified_source):
  """
   Edits source to convert usage of boost into something that pdoc3 can
   recognize
  """
  from cctbx_web_site.command_line.edit_for_boost import modify_boost_text
  new_source = modify_boost_text(modified_source)
  return new_source

def import_and_edit(module_name: str, file_path: Path):
    """
    Loads a module from a file path, modifying its source code before execution.

    Args:
        module_name: The name for the new module (e.g., 'file').
        file_path: A Path object pointing to the Python source file.

    Returns:
        The newly loaded module object.
    """
    import importlib.util

    # Step 1: Get the module's "specification" from the file path.
    # This spec contains metadata but doesn't execute the module.
    spec = importlib.util.spec_from_file_location(module_name, file_path)

    if not spec or not spec.loader:
      spec = importlib.util.find_spec(module_name)

    if not spec or not spec.loader:
        raise ImportError(f"Could not load spec for module {module_name} at {file_path}")

    # Step 2: Get the source code from the spec's loader.
    source_code = spec.loader.get_source(module_name)
    if source_code is None:
        raise ImportError(f"Could not get source for module {module_name}")

    # Step 3: Modify the source code in memory.
    from cctbx_website.command_line.comment_to_docstring import convert_comments_to_docstrings
    modified_source = convert_comments_to_docstrings(source_code)

    if files_to_edit_for_boost and file_path in files_to_edit_for_boost:
      print("Editing %s for boost" %(file_path))
      modified_source = edit_file_for_boost(modified_source)

    # Step 4: Create a new module object based on the original spec.
    # This ensures __name__, __file__, __spec__, etc., are set correctly.
    module = importlib.util.module_from_spec(spec)

    # Step 5: Compile and execute the *modified* source code.
    # The second argument to compile() is the filename, used for tracebacks.
    # We use spec.origin (the original file path) to make tracebacks accurate.
    code_obj = compile(modified_source, spec.origin, 'exec')
    exec(code_obj, module.__dict__)

    # Step 6: Add the new module to sys.modules to make it "officially" imported.
    # This allows this module to be imported by other parts of the program.
    sys.modules[module_name] = module

    return module

if __name__=="__main__":
  import sys
  args = sys.argv[1:]
  if 'top_level' in args:
     top_level = True
     args.remove('top_level')
  else:
     top_level = False

  run(args, top_level = top_level)
