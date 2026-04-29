from __future__ import absolute_import, division, print_function
import os
import shutil
import sys
import tempfile

from libtbx.command_line.run_tests_parallel import (
  _derive_unique_dir_name,
  _build_unique_dir_mapping,
)


def exercise_derive_quoted_python_command():
  used = set()
  name = _derive_unique_dir_name(
    'libtbx.python "/abs/path/tst_foo.py"', used)
  assert name == "tst_foo_py", name
  assert "tst_foo_py" in used


def exercise_derive_quoted_with_args():
  used = set()
  name = _derive_unique_dir_name(
    'libtbx.python "/abs/path/tst_foo.py" 10000', used)
  assert name == "tst_foo_py_10000", name


def exercise_derive_quoted_with_keyword_text():
  used = set()
  name = _derive_unique_dir_name(
    'libtbx.python --some_kw "/abs/path/tst_foo.py"', used)
  assert name == "tst_foo_py", name


def exercise_derive_bare_shell_script():
  used = set()
  name = _derive_unique_dir_name('/abs/path/tst_bar.sh', used)
  assert name == "tst_bar_sh", name


def exercise_derive_bare_csh_script():
  used = set()
  name = _derive_unique_dir_name('/abs/path/tst_bar.csh', used)
  assert name == "tst_bar_csh", name


def exercise_derive_sanitizes_argument_characters():
  used = set()
  name = _derive_unique_dir_name(
    'libtbx.python "/abs/path/tst_foo.py" a.b/c d', used)
  # non-[A-Za-z0-9_-] chars collapse to single underscores; stripped
  assert name == "tst_foo_py_a_b_c_d", name


def exercise_derive_collision_suffix():
  used = set()
  a = _derive_unique_dir_name('libtbx.python "/a/tst_foo.py"', used)
  b = _derive_unique_dir_name('libtbx.python "/b/tst_foo.py"', used)
  c = _derive_unique_dir_name('libtbx.python "/c/tst_foo.py"', used)
  assert a == "tst_foo_py", a
  assert b == "tst_foo_py_1", b
  assert c == "tst_foo_py_2", c


def exercise_derive_empty_basename_fallback():
  used = set()
  name = _derive_unique_dir_name('', used)
  # Falls back to a stable non-empty name so mkdir still works
  assert name == "test", name


def exercise_build_mapping_creates_directories():
  base = tempfile.mkdtemp()
  try:
    cmds = [
      'libtbx.python "/a/tst_foo.py"',
      'libtbx.python "/a/tst_bar.py" 10',
    ]
    mapping = _build_unique_dir_mapping(cmds, base)
    assert set(mapping.keys()) == set(cmds)
    assert os.path.isdir(os.path.join(base, "tst_foo_py"))
    assert os.path.isdir(os.path.join(base, "tst_bar_py_10"))
  finally:
    shutil.rmtree(base)


def exercise_build_mapping_returns_absolute_dir_paths():
  base = tempfile.mkdtemp()
  try:
    cmd = 'libtbx.python "/a/tst_foo.py"'
    mapping = _build_unique_dir_mapping([cmd], base)
    dir_path = mapping[cmd]
    assert os.path.isabs(dir_path), dir_path
    assert os.path.basename(dir_path) == "tst_foo_py", dir_path
    assert os.path.isdir(dir_path)
  finally:
    shutil.rmtree(base)


def exercise_build_mapping_deduplicates():
  base = tempfile.mkdtemp()
  try:
    cmd = 'libtbx.python "/a/tst_foo.py"'
    # Same command repeated should only create one directory
    mapping = _build_unique_dir_mapping([cmd, cmd], base)
    assert len(mapping) == 1
    assert os.path.isdir(os.path.join(base, "tst_foo_py"))
    # Second occurrence should not produce tst_foo_py_1
    assert not os.path.isdir(os.path.join(base, "tst_foo_py_1"))
  finally:
    shutil.rmtree(base)


def exercise_build_mapping_collides_different_commands():
  base = tempfile.mkdtemp()
  try:
    cmds = [
      'libtbx.python "/a/tst_foo.py"',
      'libtbx.python "/b/tst_foo.py"',
    ]
    mapping = _build_unique_dir_mapping(cmds, base)
    # Both commands kept, different directories
    assert len(mapping) == 2
    dirs = sorted(os.path.basename(p) for p in mapping.values())
    assert dirs == ["tst_foo_py", "tst_foo_py_1"], dirs
  finally:
    shutil.rmtree(base)


def run(args):
  assert len(args) == 0
  exercise_derive_quoted_python_command()
  exercise_derive_quoted_with_args()
  exercise_derive_quoted_with_keyword_text()
  exercise_derive_bare_shell_script()
  exercise_derive_bare_csh_script()
  exercise_derive_sanitizes_argument_characters()
  exercise_derive_collision_suffix()
  exercise_derive_empty_basename_fallback()
  exercise_build_mapping_creates_directories()
  exercise_build_mapping_returns_absolute_dir_paths()
  exercise_build_mapping_deduplicates()
  exercise_build_mapping_collides_different_commands()
  print("OK")


if __name__ == "__main__":
  run(args=sys.argv[1:])
