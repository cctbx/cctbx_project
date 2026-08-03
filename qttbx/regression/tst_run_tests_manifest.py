"""Every regression file in this directory must be listed in run_tests.py.

qttbx's manifest is hand-maintained and has no glob, so a file that is never
added simply never runs -- silently, and indefinitely. The manifest's own
comment claimed alphabetization would make that visible; it did not. Two files
sat unlisted, one of them a 339-line suite that was the only coverage of the
subagent-view and subject-digest machinery, passing in 0.07s whenever anyone
ran it by hand and never once in CI.

Alphabetization is a reading aid. This is the check.
"""

import re
from pathlib import Path

from libtbx.utils import format_cpu_times

_HERE = Path(__file__).resolve().parent
_MANIFEST = _HERE.parent / "run_tests.py"

#: Listed by a different mechanism, or deliberately not run as a unit test.
#: Keep this empty unless there is a real reason, and state the reason.
_EXEMPT = frozenset()

#: A path inside an ACTIVE tst_list entry. Anchored to the call, so a path in
#: prose or a docstring is not mistaken for a registration.
_ENTRY = re.compile(r"tst_list[A-Za-z0-9_]*\.(?:append|extend)\b[^\n]*?"
                    r"regression/(tst_[A-Za-z0-9_]+\.py)")


def _listed_in_manifest():
  """Test files the manifest ACTIVELY registers.

  Comment-stripped, because the check exists for the file that never runs and
  a commented-out ``tst_list.append(...)`` is precisely that: the path is
  still in the text, so a raw scan reads it as listed and the guard passes
  over the one case it was written to catch. Commenting a test out to get
  around a failure is banned by both CLAUDE.md files, which is what makes
  this the likely form of the mistake rather than an exotic one.

  Deliberately NOT implemented by importing ``run_tests`` and reading
  ``tst_list``: two ``try: import qttbx.qt`` blocks guard the Qt-dependent
  entries, so on a Qt-less host the list is 17 entries short and every one of
  them would be reported as unlisted.
  """
  lines = [ln for ln in _MANIFEST.read_text(encoding="utf-8").splitlines()
           if not ln.lstrip().startswith("#")]
  return set(_ENTRY.findall("\n".join(lines)))


def exercise_every_regression_file_is_listed():
  """A file on disk that the manifest does not name never runs."""
  listed = _listed_in_manifest()
  on_disk = {p.name for p in _HERE.glob("tst_*.py")}
  missing = sorted(on_disk - listed - _EXEMPT)
  assert not missing, (
    "these regression files exist but run_tests.py does not list them, so "
    "they run nowhere:\n  " + "\n  ".join(missing))


def exercise_the_manifest_names_no_file_that_is_gone():
  """The mirror case: a renamed or deleted test leaves an entry that fails the
  whole module's run for a file nobody meant to keep."""
  listed = _listed_in_manifest()
  on_disk = {p.name for p in _HERE.glob("tst_*.py")}
  dangling = sorted(listed - on_disk)
  assert not dangling, (
    "run_tests.py lists regression files that do not exist:\n  "
    + "\n  ".join(dangling))


def exercise():
  exercise_every_regression_file_is_listed()
  exercise_the_manifest_names_no_file_that_is_gone()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
