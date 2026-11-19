"""qttbx.viewers.coot.coot_launcher_path: where the coot conda package puts
the command that starts Coot. On Linux and macOS that is <prefix>/bin/coot.
On Windows the package installs no coot.exe, only Library\\bin\\coot.bat,
which runs Library\\bin\\WinCoot.exe."""

import os
import sys
from pathlib import Path

from libtbx.utils import format_cpu_times
from qttbx.viewers.coot import coot_environment, coot_launcher_path


def exercise_layouts():
  """bin/coot on Linux and macOS, Library\\bin\\coot.bat on Windows."""
  prefix = os.path.join("opt", "env")
  for platform in ("linux", "darwin"):
    assert coot_launcher_path(prefix, platform=platform) == \
      os.path.join(prefix, "bin", "coot")
  assert coot_launcher_path(prefix, platform="win32") == \
    os.path.join(prefix, "Library", "bin", "coot.bat")


def exercise_path_prefix_and_default_platform():
  """A pathlib prefix gives the same str as a str prefix, and the default
  platform is sys.platform read at call time, so callers that patch
  sys.platform (as the Phenix tests do) get the matching layout."""
  prefix = Path("opt") / "env"
  result = coot_launcher_path(prefix)
  assert isinstance(result, str), type(result)
  assert result == coot_launcher_path(str(prefix), platform=sys.platform)
  orig = sys.platform
  sys.platform = "win32"
  try:
    assert coot_launcher_path(prefix) == \
      os.path.join(str(prefix), "Library", "bin", "coot.bat")
  finally:
    sys.platform = orig


def exercise_environment():
  """coot_environment: a copy of the environment for starting the Coot in a
  given conda env. The caller's own conda activation exported COOT_PREFIX
  and COOT_DATA_DIR for its env, which is wrong for a Coot from another one.
  Linux/macOS: point both at the target env, as its coot-activate.sh does
  ($CONDA_PREFIX and $CONDA_PREFIX/share/coot); bin/coot sets SYMINFO
  itself. Windows: drop them and SYMINFO, which coot.bat fills in from its
  own location when unset. PYTHONHOME is dropped everywhere, PYTHONPATH
  for a Coot from another env (see exercise_environment_own_env). The input
  is not modified, and unrelated variables pass through."""
  inherited = {"COOT_PREFIX": "/phenix", "COOT_DATA_DIR": "/phenix/share/coot",
               "SYMINFO": "/phenix/syminfo.lib", "PYTHONHOME": "/phenix",
               "PYTHONPATH": "/phenix/modules", "PATH": "/usr/bin"}
  before = dict(inherited)
  prefix = os.path.join("opt", "other")
  for platform in ("linux", "darwin"):
    env = coot_environment(prefix, environ=inherited, platform=platform)
    assert env["COOT_PREFIX"] == prefix, env
    assert env["COOT_DATA_DIR"] == os.path.join(prefix, "share", "coot"), env
    assert env["SYMINFO"] == "/phenix/syminfo.lib", env
    assert "PYTHONHOME" not in env and "PYTHONPATH" not in env, env
    assert env["PATH"] == "/usr/bin", env
  env = coot_environment(prefix, environ=inherited, platform="win32")
  for name in ("COOT_PREFIX", "COOT_DATA_DIR", "SYMINFO", "PYTHONHOME",
               "PYTHONPATH"):
    assert name not in env, (name, env)
  assert env["PATH"] == "/usr/bin", env
  assert inherited == before, inherited
  # A pathlib prefix gives the same result as a str prefix.
  assert coot_environment(Path(prefix), environ=inherited,
                          platform="linux") == \
    coot_environment(prefix, environ=inherited, platform="linux")
  # No prefix, or an empty one: COOT_* stay as they are on Linux/macOS
  # rather than becoming a relative share/coot.
  for empty in (None, ""):
    env = coot_environment(empty, environ=inherited, platform="linux")
    assert env["COOT_PREFIX"] == "/phenix", (empty, env)
    assert env["COOT_DATA_DIR"] == "/phenix/share/coot", (empty, env)
    assert "PYTHONHOME" not in env, (empty, env)


def exercise_environment_launcher():
  """On Windows COOT_PREFIX, COOT_DATA_DIR and SYMINFO are removed only when
  the launcher is a batch file (the package's coot.bat, which fills them in
  again) or unspecified. WinCoot.exe run directly has nothing to refill
  them, so they are kept; PYTHONHOME and PYTHONPATH go either way. The
  launcher does not matter on Linux/macOS."""
  inherited = {"COOT_PREFIX": "/phenix", "COOT_DATA_DIR": "/phenix/share/coot",
               "SYMINFO": "/phenix/syminfo.lib", "PYTHONHOME": "/phenix",
               "PYTHONPATH": "/phenix/modules"}
  for launcher in (None, r"C:\env\Library\bin\coot.bat",
                   r"C:\WinCoot\WINCOOT.CMD"):
    env = coot_environment(None, environ=inherited, platform="win32",
                           launcher=launcher)
    for name in ("COOT_PREFIX", "COOT_DATA_DIR", "SYMINFO"):
      assert name not in env, (launcher, name, env)
  env = coot_environment(None, environ=inherited, platform="win32",
                         launcher=r"C:\env\Library\bin\WinCoot.exe")
  for name in ("COOT_PREFIX", "COOT_DATA_DIR", "SYMINFO"):
    assert env[name] == inherited[name], (name, env)
  assert "PYTHONHOME" not in env, env
  prefix = os.path.join("opt", "other")
  env = coot_environment(prefix, environ=inherited, platform="linux",
                         launcher="/x/WinCoot.exe")
  assert env["COOT_PREFIX"] == prefix, env


def exercise_environment_own_env():
  """PYTHONPATH is dropped only for a Coot from another env. For the Coot in
  the caller's own env (prefix is sys.prefix) or an unknown one (prefix
  None) it stays: in a development build it is how Coot-side scripts such
  as cootbx's view_refinement.py import libtbx."""
  inherited = {"PYTHONPATH": "/phenix/modules", "PYTHONHOME": "/phenix"}
  for platform in ("linux", "win32"):
    for prefix in (sys.prefix, Path(sys.prefix), None):
      env = coot_environment(prefix, environ=inherited, platform=platform)
      assert env.get("PYTHONPATH") == "/phenix/modules", \
        (platform, prefix, env)
      assert "PYTHONHOME" not in env, (platform, prefix, env)
    env = coot_environment(os.path.join("opt", "other"), environ=inherited,
                           platform=platform)
    assert "PYTHONPATH" not in env, (platform, env)


def exercise_environment_defaults():
  """With no environ, coot_environment copies os.environ and leaves it
  untouched; with no platform it reads sys.platform at call time."""
  saved = dict(os.environ)
  os.environ["COOT_PREFIX"] = "/phenix"
  os.environ["PYTHONHOME"] = "/phenix"
  orig = sys.platform
  try:
    snapshot = dict(os.environ)
    sys.platform = "win32"
    env = coot_environment(os.path.join("opt", "other"))
    assert "COOT_PREFIX" not in env and "PYTHONHOME" not in env, env
    sys.platform = "linux"
    env = coot_environment(os.path.join("opt", "other"))
    assert env["COOT_PREFIX"] == os.path.join("opt", "other"), env
    assert dict(os.environ) == snapshot, "os.environ was modified"
  finally:
    sys.platform = orig
    os.environ.clear()
    os.environ.update(saved)


def exercise():
  exercise_layouts()
  exercise_path_prefix_and_default_platform()
  exercise_environment()
  exercise_environment_launcher()
  exercise_environment_own_env()
  exercise_environment_defaults()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
