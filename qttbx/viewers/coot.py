"""
Where the coot conda package puts the command that starts Coot, and the
environment that command needs.

The package comes from the coot channel on anaconda.org. Phenix bundles it
in its own conda environment, and a Coot can also be installed in a separate
conda environment.
"""

import os
import sys


def coot_launcher_path(prefix, platform=None):
  """
  Path of the command that starts Coot in a conda environment.

  On Linux and macOS this is the ``bin/coot`` script. On Windows the package
  installs no ``coot.exe``: ``Library\\bin\\coot.bat`` puts its environment
  on PATH, sets LC_ALL and LANG to C, fills in ``COOT_DATA_DIR``,
  ``COOT_PREFIX``, ``SYMINFO`` and ``PYTHONHOME`` where they are unset, and
  runs ``Library\\bin\\WinCoot.exe``, so the batch file is the command to
  run. The path is not checked and may not exist.

  Parameters
  ----------
  prefix : str or pathlib.Path
      The conda environment prefix.
  platform : str, optional
      A ``sys.platform`` value. When omitted, ``sys.platform`` is read at
      call time.

  Returns
  -------
  str
      The path of the launcher under ``prefix``.
  """
  if platform is None:
    platform = sys.platform
  if platform == "win32":
    return os.path.join(str(prefix), "Library", "bin", "coot.bat")
  return os.path.join(str(prefix), "bin", "coot")


def coot_environment(prefix, environ=None, platform=None, launcher=None):
  """
  Environment for starting the Coot installed in a conda environment.

  A program running in another conda environment (Phenix with its own
  bundled Coot, for example) has COOT_PREFIX and COOT_DATA_DIR exported for
  that environment by its activation, and they are wrong for a Coot from a
  different one. On Linux and macOS both are set for ``prefix`` the way the
  package's coot-activate.sh sets them ($CONDA_PREFIX and
  $CONDA_PREFIX/share/coot); ``bin/coot`` sets SYMINFO itself. On Windows
  they are removed together with SYMINFO when the launcher is a batch file
  such as the package's coot.bat, which fills all three in from its own
  location when they are unset; WinCoot.exe run directly has nothing to
  fill them in, so they are kept for it. PYTHONHOME is always removed.
  PYTHONPATH is removed for a Coot from another environment, whose Python
  must not import the caller's packages; for the Coot in the caller's own
  environment (``prefix`` is ``sys.prefix``) or an unknown one it is kept,
  because in a development build it is how Coot-side scripts import
  libtbx.

  Only these variables are handled. Others that the caller's activation
  exported, such as GSETTINGS_SCHEMA_DIR and XML_CATALOG_FILES, still
  describe the caller's environment.

  Parameters
  ----------
  prefix : str or pathlib.Path or None
      The conda environment holding the Coot. When None or empty,
      COOT_PREFIX and COOT_DATA_DIR are left as they are on Linux and macOS.
  environ : mapping, optional
      The environment to start from; ``os.environ`` when omitted. It is not
      modified.
  platform : str, optional
      A ``sys.platform`` value. When omitted, ``sys.platform`` is read at
      call time.
  launcher : str or pathlib.Path, optional
      The command that will start Coot. Only used on Windows; when omitted
      it is taken to be the package's coot.bat.

  Returns
  -------
  dict
      A new environment for the Coot process.
  """
  if platform is None:
    platform = sys.platform
  env = dict(os.environ if environ is None else environ)
  env.pop("PYTHONHOME", None)
  if prefix is not None and len(str(prefix)) != 0 and \
     os.path.realpath(str(prefix)) != os.path.realpath(sys.prefix):
    env.pop("PYTHONPATH", None)
  if platform == "win32":
    if launcher is None or \
       str(launcher).lower().endswith((".bat", ".cmd")):
      for name in ("COOT_PREFIX", "COOT_DATA_DIR", "SYMINFO"):
        env.pop(name, None)
  elif prefix is not None and len(str(prefix)) != 0:
    env["COOT_PREFIX"] = str(prefix)
    env["COOT_DATA_DIR"] = os.path.join(str(prefix), "share", "coot")
  return env
