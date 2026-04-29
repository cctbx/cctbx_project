from __future__ import absolute_import, division, print_function

import errno
import os
import shutil
import stat
import sys
import tempfile
import time

from six.moves import cStringIO as StringIO

from libtbx.utils import format_cpu_times
from libtbx.test_utils import parallel as parallel_mod


def exercise_chmod_then_retry_handles_readonly_file():
  """_chmod_then_retry should chmod a read-only file then remove it."""
  tmpdir = tempfile.mkdtemp(prefix="tst_retry_chmod_")
  try:
    victim = os.path.join(tmpdir, "locked.txt")
    with open(victim, "w") as f:
      f.write("payload")
    os.chmod(victim, stat.S_IREAD)
    # Fabricate an exc_info as if os.remove had raised EACCES.
    exc_info = (OSError, OSError(errno.EACCES, "simulated"), None)
    parallel_mod._chmod_then_retry(os.remove, victim, exc_info)
    assert not os.path.exists(victim), \
      "_chmod_then_retry did not remove %s" % victim
  finally:
    # Best-effort cleanup in case the assertion fired.
    shutil.rmtree(tmpdir, ignore_errors=True)


def exercise_chmod_then_retry_ignores_unrelated_errors():
  """_chmod_then_retry should be a no-op for non-EACCES/EPERM errors."""
  tmpdir = tempfile.mkdtemp(prefix="tst_retry_chmod_")
  try:
    victim = os.path.join(tmpdir, "present.txt")
    with open(victim, "w") as f:
      f.write("payload")
    # Fabricate an ENOENT — helper should do nothing.
    exc_info = (OSError, OSError(errno.ENOENT, "simulated"), None)
    parallel_mod._chmod_then_retry(os.remove, victim, exc_info)
    assert os.path.exists(victim), \
      "_chmod_then_retry removed file on unrelated error"
  finally:
    shutil.rmtree(tmpdir, ignore_errors=True)


def exercise_reset_cwd_removes_top_level_files_and_dirs():
  """_reset_cwd_for_retry should empty the directory in place."""
  tmpdir = tempfile.mkdtemp(prefix="tst_retry_reset_")
  try:
    # Top-level file.
    with open(os.path.join(tmpdir, "a.txt"), "w") as f:
      f.write("x")
    # Top-level directory with nested content.
    nested = os.path.join(tmpdir, "sub")
    os.makedirs(os.path.join(nested, "deeper"))
    with open(os.path.join(nested, "deeper", "b.txt"), "w") as f:
      f.write("y")
    # Symlink (skip on platforms where os.symlink is unprivileged).
    link_src = os.path.join(tmpdir, "link")
    try:
      os.symlink(nested, link_src)
    except (OSError, NotImplementedError, AttributeError):
      link_src = None  # Non-POSIX-or-unprivileged platform.

    parallel_mod._reset_cwd_for_retry(tmpdir)

    assert os.path.isdir(tmpdir), "directory itself should survive"
    assert os.listdir(tmpdir) == [], \
      "directory not emptied, still contains: %r" % os.listdir(tmpdir)
    if link_src is not None:
      # The symlink target must NOT have been touched.
      # (We did not create the target outside tmpdir, but we assert the
      # symlink itself is gone and did not cause recursion outside tmpdir.)
      assert not os.path.lexists(link_src)
  finally:
    shutil.rmtree(tmpdir, ignore_errors=True)


def exercise_reset_cwd_handles_readonly_file():
  """_reset_cwd_for_retry should remove read-only files at top level."""
  tmpdir = tempfile.mkdtemp(prefix="tst_retry_reset_ro_")
  try:
    locked = os.path.join(tmpdir, "locked.txt")
    with open(locked, "w") as f:
      f.write("x")
    os.chmod(locked, stat.S_IREAD)

    parallel_mod._reset_cwd_for_retry(tmpdir)

    assert os.path.isdir(tmpdir)
    assert os.listdir(tmpdir) == [], \
      "read-only file not removed: %r" % os.listdir(tmpdir)
  finally:
    shutil.rmtree(tmpdir, ignore_errors=True)


def exercise_reset_cwd_is_noop_on_missing_or_none():
  """_reset_cwd_for_retry should accept None / missing paths silently."""
  parallel_mod._reset_cwd_for_retry(None)
  parallel_mod._reset_cwd_for_retry("")
  parallel_mod._reset_cwd_for_retry("/nonexistent_path_for_tst_retry")


def _write_flaky_script(dir_path, pass_on_attempt):
  """Write a Python script that passes only on attempt N.

  The script reads/writes ``counter.txt`` inside ``dir_path`` using
  an absolute path baked into the generated source. This means the
  counter survives ``cwd`` cleanup between retries, which lets callers
  freely pass any ``cwd`` to ``run_command`` while still accumulating
  the attempt count. Exits 0 once the counter reaches
  ``pass_on_attempt``; otherwise exits 1.
  """
  script = os.path.join(dir_path, "flaky.py")
  counter_path = os.path.join(dir_path, "counter.txt").replace("\\", "/")
  source = (
    "from __future__ import print_function\n"
    "import os, sys\n"
    "counter_path = r'''%s'''\n"
    "try:\n"
    "  with open(counter_path) as f:\n"
    "    n = int(f.read().strip() or '0')\n"
    "except (IOError, OSError, ValueError):\n"
    "  n = 0\n"
    "n += 1\n"
    "with open(counter_path, 'w') as f:\n"
    "  f.write(str(n))\n"
    "print('attempt', n)\n"
    "if n >= %d:\n"
    "  print('OK')\n"
    "  sys.exit(0)\n"
    "sys.exit(1)\n"
  ) % (counter_path, pass_on_attempt)
  with open(script, "w") as f:
    f.write(source)
  return script


def exercise_run_command_retries_until_success():
  """run_command(max_retries=2) passes a script that succeeds on attempt 3."""
  tmpdir = tempfile.mkdtemp(prefix="tst_retry_run_")
  try:
    script = _write_flaky_script(tmpdir, pass_on_attempt=3)
    # NB: cwd is NOT tmpdir (that's where the script + counter live).
    # We use a separate subdir as the cwd so cleanup between retries
    # does not touch the counter.
    cwd = os.path.join(tmpdir, "run_cwd")
    os.makedirs(cwd)
    cmd = 'libtbx.python "%s"' % script
    result = parallel_mod.run_command(
      cmd, cwd=cwd, max_retries=2)
    assert result is not None
    assert result.return_code == 0, \
      "expected eventual success, got %r" % result.return_code
    assert result.attempt == 3, \
      "expected attempt=3, got %r" % result.attempt
    assert result.attempts_total == 3
  finally:
    shutil.rmtree(tmpdir, ignore_errors=True)


def exercise_run_command_no_retries_fails_once():
  """run_command(max_retries=0) never retries."""
  tmpdir = tempfile.mkdtemp(prefix="tst_retry_run_")
  try:
    script = _write_flaky_script(tmpdir, pass_on_attempt=2)
    cwd = os.path.join(tmpdir, "run_cwd")
    os.makedirs(cwd)
    cmd = 'libtbx.python "%s"' % script
    result = parallel_mod.run_command(cmd, cwd=cwd, max_retries=0)
    assert result is not None
    assert result.return_code != 0
    assert result.attempt == 1
    assert result.attempts_total == 1
  finally:
    shutil.rmtree(tmpdir, ignore_errors=True)


def exercise_run_command_skip_retry_overrides_max_retries():
  """skip_retry=True forces one attempt even with max_retries=5."""
  tmpdir = tempfile.mkdtemp(prefix="tst_retry_run_")
  try:
    script = _write_flaky_script(tmpdir, pass_on_attempt=99)
    cwd = os.path.join(tmpdir, "run_cwd")
    os.makedirs(cwd)
    cmd = 'libtbx.python "%s"' % script
    result = parallel_mod.run_command(
      cmd, cwd=cwd, max_retries=5, skip_retry=True)
    assert result is not None
    assert result.return_code != 0
    assert result.attempt == 1
    assert result.attempts_total == 1
  finally:
    shutil.rmtree(tmpdir, ignore_errors=True)


def exercise_run_command_cleans_cwd_between_attempts():
  """Cleanup runs between retries: files left by attempt 1 are gone by attempt 2."""
  tmpdir = tempfile.mkdtemp(prefix="tst_retry_run_")
  try:
    script_path = os.path.join(tmpdir, "cleanup.py")
    counter_path = os.path.join(tmpdir, "external_counter.txt").replace("\\", "/")
    source = (
      "from __future__ import print_function\n"
      "import os, sys\n"
      "counter = r'''%s'''\n"
      "try:\n"
      "  with open(counter) as f:\n"
      "    n = int(f.read().strip() or '0')\n"
      "except (IOError, OSError, ValueError):\n"
      "  n = 0\n"
      "n += 1\n"
      "with open(counter, 'w') as f:\n"
      "  f.write(str(n))\n"
      "# On non-first attempts, marker.txt from prior attempts must be gone.\n"
      "if n > 1:\n"
      "  assert not os.path.exists('marker.txt'), \\\n"
      "    'cleanup failed: marker.txt still present on attempt %%d' %% n\n"
      "with open('marker.txt', 'w') as f:\n"
      "  f.write('from attempt %%d' %% n)\n"
      "if n < 2:\n"
      "  sys.exit(1)\n"
      "print('OK')\n"
      "sys.exit(0)\n"
    ) % counter_path
    with open(script_path, "w") as f:
      f.write(source)
    cwd = os.path.join(tmpdir, "run_cwd")
    os.makedirs(cwd)
    cmd = 'libtbx.python "%s"' % script_path
    result = parallel_mod.run_command(cmd, cwd=cwd, max_retries=2)
    assert result.return_code == 0, \
      "script failed; stderr=%r" % result.stderr_lines
    assert result.attempt == 2
  finally:
    shutil.rmtree(tmpdir, ignore_errors=True)


def exercise_run_command_backoff_timing():
  """A fully-failing 3-attempt run sleeps ~3s total (1 + 2)."""
  tmpdir = tempfile.mkdtemp(prefix="tst_retry_run_")
  try:
    script = _write_flaky_script(tmpdir, pass_on_attempt=99)
    cwd = os.path.join(tmpdir, "run_cwd")
    os.makedirs(cwd)
    cmd = 'libtbx.python "%s"' % script
    t0 = time.time()
    result = parallel_mod.run_command(cmd, cwd=cwd, max_retries=2)
    elapsed = time.time() - t0
    assert result.return_code != 0
    assert result.attempt == 3
    # Sleeps of 1s + 2s = 3s, plus 3 process spawns. Upper bound is
    # generous to avoid CI flakiness.
    assert elapsed >= 3.0, \
      "expected >= 3s elapsed (backoff), got %.2fs" % elapsed
    assert elapsed <= 15.0, \
      "expected <= 15s elapsed, got %.2fs" % elapsed
  finally:
    shutil.rmtree(tmpdir, ignore_errors=True)


def exercise_run_command_list_respects_max_retries():
  """run_command_list forwards max_retries to each run_command."""
  from libtbx.utils import null_out
  tmpdir = tempfile.mkdtemp(prefix="tst_retry_rcl_")
  try:
    script = _write_flaky_script(tmpdir, pass_on_attempt=2)
    cwd = os.path.join(tmpdir, "cwd_a")
    os.makedirs(cwd)
    cmd = 'libtbx.python "%s"' % script
    log = null_out()
    out = StringIO()
    rcl = parallel_mod.run_command_list(
      cmd_list=[cmd],
      nprocs=1,
      out=out,
      log=log,
      verbosity=0,
      cwd_map={cmd: cwd},
      max_retries=2)
    # Final state: the single test passed on attempt 2.
    assert len(rcl.results) == 1
    assert rcl.results[0].return_code == 0
    assert rcl.results[0].attempt == 2
    assert rcl.failure == 0
  finally:
    shutil.rmtree(tmpdir, ignore_errors=True)


def exercise_run_command_list_skips_retry_for_expected_failures():
  """Commands in expected_failure_list run once even with max_retries>0."""
  from libtbx.utils import null_out
  tmpdir = tempfile.mkdtemp(prefix="tst_retry_rcl_")
  try:
    script = _write_flaky_script(tmpdir, pass_on_attempt=999)
    cwd = os.path.join(tmpdir, "cwd_a")
    os.makedirs(cwd)
    cmd = 'libtbx.python "%s"' % script
    out = StringIO()
    rcl = parallel_mod.run_command_list(
      cmd_list=[cmd],
      expected_failure_list=[cmd],
      nprocs=1,
      out=out,
      log=null_out(),
      verbosity=0,
      cwd_map={cmd: cwd},
      max_retries=5)
    assert len(rcl.results) == 1
    assert rcl.results[0].attempt == 1, \
      "expected_failure should not retry; attempt=%r" \
      % rcl.results[0].attempt
  finally:
    shutil.rmtree(tmpdir, ignore_errors=True)


def exercise_cli_auto_enables_unique_dirs_when_max_retries_positive():
  """CLI run() flips run_in_unique_dirs=True when max_retries>0."""
  from libtbx.command_line import run_tests_parallel as rtp
  captured = StringIO()
  old_stdout = sys.stdout
  sys.stdout = captured
  try:
    # return_list_of_tests=True short-circuits execution so we can test
    # the PHIL parsing + auto-enable logic without actually running
    # tests.
    rtp.run(
      ["module=libtbx", "max_retries=1", "run_in_unique_dirs=False"],
      return_list_of_tests=True)
  finally:
    sys.stdout = old_stdout
  out = captured.getvalue()
  assert "enabled run_in_unique_dirs=True" in out, \
    "expected auto-enable notice, got:\n%s" % out


def exercise_cli_leaves_flag_alone_when_max_retries_is_zero():
  """CLI run() with max_retries=0 does not auto-enable."""
  from libtbx.command_line import run_tests_parallel as rtp
  captured = StringIO()
  old_stdout = sys.stdout
  sys.stdout = captured
  try:
    rtp.run(
      ["module=libtbx", "max_retries=0", "run_in_unique_dirs=False"],
      return_list_of_tests=True)
  finally:
    sys.stdout = old_stdout
  out = captured.getvalue()
  assert "enabled run_in_unique_dirs=True" not in out, \
    "auto-enable fired with max_retries=0:\n%s" % out


def exercise_display_result_shows_retry_count():
  """run_command_list output includes 'passed on attempt N of M' when retried."""
  from libtbx.utils import null_out
  tmpdir = tempfile.mkdtemp(prefix="tst_retry_disp_")
  try:
    script = _write_flaky_script(tmpdir, pass_on_attempt=2)
    cwd = os.path.join(tmpdir, "cwd")
    os.makedirs(cwd)
    cmd = 'libtbx.python "%s"' % script
    out = StringIO()
    parallel_mod.run_command_list(
      cmd_list=[cmd],
      nprocs=1,
      out=out,
      log=null_out(),
      verbosity=1,
      cwd_map={cmd: cwd},
      max_retries=2)
    captured = out.getvalue()
    assert "passed on attempt 2 of 3" in captured, \
      "missing per-test retry annotation in output:\n%s" % captured
    assert "Retries used" in captured, \
      "missing summary retry footer in output:\n%s" % captured
  finally:
    shutil.rmtree(tmpdir, ignore_errors=True)


def exercise_display_result_omits_retry_note_on_clean_pass():
  """No retry annotation when a test passes on the first attempt."""
  from libtbx.utils import null_out
  tmpdir = tempfile.mkdtemp(prefix="tst_retry_disp_")
  try:
    script = _write_flaky_script(tmpdir, pass_on_attempt=1)
    cwd = os.path.join(tmpdir, "cwd")
    os.makedirs(cwd)
    cmd = 'libtbx.python "%s"' % script
    out = StringIO()
    parallel_mod.run_command_list(
      cmd_list=[cmd],
      nprocs=1,
      out=out,
      log=null_out(),
      verbosity=1,
      cwd_map={cmd: cwd},
      max_retries=2)
    captured = out.getvalue()
    assert "passed on attempt" not in captured, \
      "unexpected retry annotation on clean pass:\n%s" % captured
    assert "Retries used" not in captured, \
      "unexpected retry footer on clean run:\n%s" % captured
  finally:
    shutil.rmtree(tmpdir, ignore_errors=True)


def exercise_run_command_rejects_negative_max_retries():
  """run_command raises ValueError on negative max_retries."""
  cmd = 'libtbx.python -c "pass"'
  try:
    parallel_mod.run_command(cmd, max_retries=-1)
  except ValueError as e:
    assert "max_retries must be >= 0" in str(e), \
      "unexpected ValueError message: %r" % str(e)
  else:
    raise AssertionError(
      "run_command did not raise ValueError on max_retries=-1")


def exercise_cli_phil_rejects_negative_max_retries():
  """The PHIL parser rejects max_retries<0 at .extract() time."""
  from libtbx.command_line import run_tests_parallel as rtp
  try:
    rtp.run(
      ["module=libtbx", "max_retries=-1"],
      return_list_of_tests=True)
  except RuntimeError as e:
    # PHIL raises RuntimeError from _check_value when value_min is
    # violated; the message mentions the parameter name and the bound.
    msg = str(e)
    assert "max_retries" in msg and "minimum" in msg, \
      "unexpected PHIL RuntimeError message: %r" % msg
  else:
    raise AssertionError(
      "CLI run() did not reject max_retries=-1")


def run():
  exercise_chmod_then_retry_handles_readonly_file()
  exercise_chmod_then_retry_ignores_unrelated_errors()
  exercise_reset_cwd_removes_top_level_files_and_dirs()
  exercise_reset_cwd_handles_readonly_file()
  exercise_reset_cwd_is_noop_on_missing_or_none()
  exercise_run_command_retries_until_success()
  exercise_run_command_no_retries_fails_once()
  exercise_run_command_skip_retry_overrides_max_retries()
  exercise_run_command_cleans_cwd_between_attempts()
  exercise_run_command_backoff_timing()
  exercise_run_command_list_respects_max_retries()
  exercise_run_command_list_skips_retry_for_expected_failures()
  exercise_cli_auto_enables_unique_dirs_when_max_retries_positive()
  exercise_cli_leaves_flag_alone_when_max_retries_is_zero()
  exercise_display_result_shows_retry_count()
  exercise_display_result_omits_retry_note_on_clean_pass()
  exercise_run_command_rejects_negative_max_retries()
  exercise_cli_phil_rejects_negative_max_retries()


if __name__ == "__main__":
  run()
  print(format_cpu_times())
  print("OK")
