"""Environment hygiene shared by every dispatcher writer.

A dispatcher runs in a copy of the caller's environment. Several variables in
that environment make Python or the dynamic loader pick up code from outside
the installation: PYTHONPATH, the per-user site-packages
directory, the loader search paths and DYLD_INSERT_LIBRARIES, the module
paths of GTK, GIO, GDK-pixbuf, fontconfig and Qt that a desktop session or
another conda environment exports. The guard emitted here clears them and
disables user site-packages before the interpreter starts. Setting either
name in ``dispatcher_environment_opt_outs`` skips the guard, for callers who
deliberately combine environments.

PYTHONHOME is different: a stray value stops the interpreter from starting,
and no caller needs to keep one, so ``unconditional_lines`` clears it outside
the guard.

Left alone on purpose: LD_PRELOAD, because VirtualGL's vglrun works by
setting it; MPLBACKEND and MATPLOTLIBRC, because headless users set
MPLBACKEND=Agg and a clean sys.path already stops a foreign matplotlib
backend from loading; and GTK_IM_MODULE, which only chooses an input method
among the modules found through paths the guard clears.

Known reasons to opt out (set one of the names to a non-empty value): a
development build whose CUDA runtime comes from ``module load`` through
LD_LIBRARY_PATH (XFEL and DIALS builds link cudart without an rpath), and
dxtbx format classes installed with ``dxtbx.install_format --user``, which
live in the user site-packages the guard disables.

This module must stay importable with only the standard library: bootstrap
runs the development-build include generator with a Python that has no
cctbx_project on sys.path.
"""
from __future__ import absolute_import, division, print_function

dispatcher_environment_opt_outs = (
  "CCTBX_CONDA_USE_ENVIRONMENT_VARIABLES",
  "PHENIX_TRUST_OTHER_ENV",
)

# Cleared by the POSIX guard.
dispatcher_cleared_variables = (
  # Python (PYTHONHOME is cleared unconditionally, see unconditional_lines)
  "PYTHONPATH",
  # dynamic loader
  # (LD_PRELOAD is kept: vglrun works by setting it)
  "LD_LIBRARY_PATH",
  "DYLD_LIBRARY_PATH", "DYLD_FALLBACK_LIBRARY_PATH",
  "DYLD_FRAMEWORK_PATH", "DYLD_FALLBACK_FRAMEWORK_PATH",
  "DYLD_INSERT_LIBRARIES",
  # GTK, GIO, GDK-pixbuf, GObject introspection, fontconfig (wxPython on
  # Linux, bundled Coot)
  "GTK_MODULES", "GTK_PATH", "GTK_IM_MODULE_FILE",
  "GTK2_RC_FILES", "GTK_EXE_PREFIX", "GTK_DATA_PREFIX",
  "GDK_PIXBUF_MODULE_FILE", "GDK_PIXBUF_MODULEDIR",
  "GIO_MODULE_DIR", "GIO_EXTRA_MODULES", "GI_TYPELIB_PATH",
  "FONTCONFIG_FILE", "FONTCONFIG_PATH",
  # Qt
  "QT_PLUGIN_PATH", "QT_QPA_PLATFORM_PLUGIN_PATH", "QT_QPA_PLATFORMTHEME",
  "QT_STYLE_OVERRIDE",
)

# Cleared by the Windows guard. DLLs resolve through PATH there, which the
# dispatcher already puts the installation first on.
dispatcher_cleared_variables_win32 = (
  "PYTHONPATH",
  "QT_PLUGIN_PATH", "QT_QPA_PLATFORM_PLUGIN_PATH", "QT_QPA_PLATFORMTHEME",
  "QT_STYLE_OVERRIDE",
)

# Set by the guard.
dispatcher_set_variables = (
  ("PYTHONNOUSERSITE", "1"),
)
dispatcher_set_variables_win32 = (
  ("PYTHONNOUSERSITE", "1"),
)

def _guard_comment():
  """The guard's comment, built from the current opt-out names."""
  return (
    "Isolate this run from the caller's environment: PYTHONPATH, user",
    "site-packages, library search paths and GUI-toolkit module paths from",
    "another Python or desktop session would otherwise shadow the packages",
    "shipped here. To keep the caller's environment instead, set",
    " or ".join(dispatcher_environment_opt_outs) + ".",
  )


def unconditional_lines(shell):
  """Lines every dispatcher runs whether or not the caller opts out of the
  guard: the interpreter must use its own standard library.

  Parameters
  ----------
  shell : str
      "sh" or "bat".
  """
  if shell == "sh":
    return ["unset PYTHONHOME"]
  if shell == "bat":
    return ['@set "PYTHONHOME="']
  raise ValueError("shell must be 'sh' or 'bat', not %r" % (shell,))


def environment_guard_lines(shell):
  """Lines that clear the caller's environment before Python starts.

  Parameters
  ----------
  shell : str
      "sh" or "bat".
  """
  if shell == "sh":
    lines = ["# " + c for c in _guard_comment()]
    lines.append("if %s; then" % " && ".join(
      '[ -z "${%s}" ]' % v for v in dispatcher_environment_opt_outs))
    for v in dispatcher_cleared_variables:
      lines.append("  unset %s" % v)
    for name, value in dispatcher_set_variables:
      lines.append("  %s=%s" % (name, value))
      lines.append("  export %s" % name)
    lines.append("fi")
    return lines
  if shell == "bat":
    lines = ["@rem " + c for c in _guard_comment()]
    lines.append("@" + " ".join(
      "if not defined %s" % v for v in dispatcher_environment_opt_outs) + " (")
    for v in dispatcher_cleared_variables_win32:
      lines.append('  @set "%s="' % v)
    for name, value in dispatcher_set_variables_win32:
      lines.append('  @set "%s=%s"' % (name, value))
    lines.append(")")
    return lines
  raise ValueError("shell must be 'sh' or 'bat', not %r" % (shell,))


def external_path_lines(shell):
  """Opt-in search paths for the DivCon (QuantumBio) integration, appended
  after the guard so they are kept whether or not the caller's environment is
  trusted. Only the POSIX dispatchers carry them."""
  if shell == "bat":
    return []
  if shell != "sh":
    raise ValueError("shell must be 'sh' or 'bat', not %r" % (shell,))
  lines = [
    "# Opt-in paths for the DivCon (QuantumBio) integration; appended after",
    "# the guard so they are kept either way.",
  ]
  for var in ("PYTHONPATH", "LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH"):
    qb = "QB_" + var
    lines += [
      'if [ -n "${%s}" ]; then' % qb,
      '  %s="${%s:+${%s}:}${%s}"' % (var, var, var, qb),
      '  export %s' % var,
      'fi',
    ]
  return lines
