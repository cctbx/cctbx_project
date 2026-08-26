import ast
import importlib
import inspect
import sys
from pathlib import Path

def ensure_modules_on_syspath(file_path: Path):
  parts = file_path.resolve().parts
  if "modules" not in parts:
    raise RuntimeError(f"Path does not contain 'modules': {file_path}")
  idx = parts.index("modules")
  modules_dir = Path(*parts[: idx + 1])  # .../modules
  modules_dir_str = str(modules_dir)
  if modules_dir_str not in sys.path:
    sys.path.insert(0, modules_dir_str)

def package_from_realpath(file_path: Path) -> str:
  parts = file_path.resolve().parts
  idx = parts.index("modules")
  pkg_parts = parts[idx + 1 : -1]  # after modules, exclude filename
  return ".".join(pkg_parts)

def import_test_module(file_path: Path):
  stem = file_path.stem  # correct for .py
  pkg = package_from_realpath(file_path)
  module_name = f"{pkg}.{stem}" if pkg else stem
  return importlib.import_module(module_name)

def is_if_name_equals_main(test_expr: ast.AST) -> bool:
  if not isinstance(test_expr, ast.Compare):
    return False
  if len(test_expr.ops) != 1 or len(test_expr.comparators) != 1:
    return False
  if not isinstance(test_expr.ops[0], ast.Eq):
    return False
  left = test_expr.left
  right = test_expr.comparators[0]
  def is_name(node):
    return isinstance(node, ast.Name) and node.id == "__name__"
  def is_main_str(node):
    return isinstance(node, ast.Constant) and node.value == "__main__"
  return (is_name(left) and is_main_str(right)) or \
         (is_main_str(left) and is_name(right))

def count_run_calls_in_main_block(file_path: Path) -> int:
  tree = ast.parse(
    file_path.read_text(encoding="utf-8"), filename=str(file_path))
  count = 0
  for node in ast.walk(tree):
    if isinstance(node, ast.If) and is_if_name_equals_main(node.test):
      for inner in ast.walk(ast.Module(body=node.body, type_ignores=[])):
        if (isinstance(inner, ast.Call) and isinstance(inner.func, ast.Name) and
            inner.func.id == "run"):
          count += 1
  return count

def run(root="."):
  root = Path(root).expanduser().resolve()
  if not root.is_dir():
    raise RuntimeError(f"Not a directory: {root}")
  for file_path in sorted(root.glob("tst_*.py")): #only this folder (non-recursive)
    ensure_modules_on_syspath(file_path)
    print(file_path.name)
    mod = import_test_module(file_path)
    if not hasattr(mod, "run") or not callable(mod.run):
      raise RuntimeError(f"{file_path.name}: missing callable run()")
    doc = mod.run.__doc__
    if not doc or not doc.strip():
      raise RuntimeError(f"{file_path.name}: run() has no docstring")
    sig = inspect.signature(mod.run)
    if "prefix" not in sig.parameters:
      raise RuntimeError(f"{file_path.name}: run() must accept prefix")
    if sig.parameters["prefix"].default is inspect._empty:
      raise RuntimeError(
        f"{file_path.name}: run(prefix=...) must have a default")
    expected_default = Path(mod.__file__).name.replace(".py", "")
    if sig.parameters["prefix"].default != expected_default:
      raise RuntimeError(
        f"{file_path.name}: expected run(prefix=...) default '{expected_default}', "
        f"got '{sig.parameters['prefix'].default}'"
      )
    n = count_run_calls_in_main_block(file_path)
    if n != 1:
      raise RuntimeError(
        f"{file_path.name}: expected 1 run(...) call in __main__, found {n}")
    print(doc)

if __name__ == "__main__":
  import sys
  run(sys.argv[1] if len(sys.argv) > 1 else ".")
