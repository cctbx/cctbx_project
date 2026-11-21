import os
import shutil
import sys
import tempfile
from io import StringIO

from libtbx.utils import Sorry, format_cpu_times
from libtbx.test_utils import parallel, show_diff
from libtbx.command_line import run_tests_parallel


def _sorry_from(function, *args, **keyword_args):
  """Return the message of the Sorry that the call has to raise."""
  try:
    function(*args, **keyword_args)
  except Sorry as e:
    return str(e)
  raise AssertionError("no Sorry from %s%s" % (function.__name__, args))


def exercise_test_names():
  """A test is named by its path below the directories of the install, or of
  another install by the modules and build in its path, plus arguments."""
  build = "/opt/modules/phenix-2.0/build_py311"
  repositories = ["/opt/modules/phenix-2.0/modules",
                  "/opt/modules/phenix-2.0/modules/cctbx_project"]
  windows = ["C:\\ph\\src", "C:\\ph\\src\\cctbx_project"]
  for command, build_path, repository_paths, expected in [
      # not in this install, as in a log from another machine
      ('libtbx.python "/a/b/modules/cctbx_project/iotbx/tst_x.py"', None, None,
        "cctbx_project/iotbx/tst_x.py"),
      # arguments are part of the name, interpreter options are not
      ('libtbx.python "/r/modules/cctbx_project/mmtbx/tst_link.py" 1-5',
        None, None, "cctbx_project/mmtbx/tst_link.py 1-5"),
      ('libtbx.python -m coverage run "/r/modules/cctbx_project/tst_z.py"',
        None, None, "cctbx_project/tst_z.py"),
      # compiled tests, a source directory called build, a shell script
      ('"/r/build/scitbx/array_family/tst_af_1"', None, None,
        "build/scitbx/array_family/tst_af_1"),
      ('libtbx.python "/r/modules/cctbx_project/mmtbx/ss/build/tst_1.py"',
        None, None, "cctbx_project/mmtbx/ss/build/tst_1.py"),
      ("/r/modules/phenix_regression/misc/tst_y.csh", None, None,
        "phenix_regression/misc/tst_y.csh"),
      # only the first modules directory counts; .. and // are resolved
      ('libtbx.python "/r/modules/repo/modules/tst_x.py"', None, None,
        "repo/modules/tst_x.py"),
      ('libtbx.python "/r/modules//repo/sub/../tst_x.py"', None, None,
        "repo/tst_x.py"),
      # in the install: whatever its directories are called or lie under
      ('libtbx.python '
        '"/opt/modules/phenix-2.0/modules/cctbx_project/libtbx/tst_x.py" 3',
        build, repositories, "cctbx_project/libtbx/tst_x.py 3"),
      ('libtbx.python '
        '"/opt/modules/phenix-2.0/modules/cctbx_project/libtbx/tst_x.py" 3',
        build, list(reversed(repositories)), "cctbx_project/libtbx/tst_x.py 3"),
      ('"/opt/modules/phenix-2.0/build_py311/scitbx/tst_af_1"',
        build, repositories, "build/scitbx/tst_af_1"),
      ('libtbx.python "C:\\ph\\src\\cctbx_project\\libtbx\\tst_x.py"',
        "C:\\ph\\bld", windows, "cctbx_project/libtbx/tst_x.py"),
      ('"C:\\ph\\bld\\scitbx\\tst_af_1"', "C:\\ph\\bld", windows,
        "build/scitbx/tst_af_1")]:
    got = parallel.test_name_from_command(
      command, build_path, repository_paths)
    assert got == expected, (command, got)


def exercise_test_name_of_a_file_does_not_depend_on_how_it_is_written():
  """A test that exists gets the same name whether its path is absolute,
  relative to the current directory, or has .. and doubled separators."""
  tmpdir = os.path.realpath(tempfile.mkdtemp(prefix="tst_timings_"))
  cwd = os.getcwd()
  try:
    src = os.path.join(tmpdir, "src")
    os.makedirs(os.path.join(src, "repo", "sub"))
    with open(os.path.join(src, "repo", "sub", "tst_x.py"), "w") as f:
      f.write("")
    os.chdir(os.path.join(src, "repo"))
    for path in [
        os.path.join(src, "repo", "sub", "tst_x.py"),
        os.path.join("sub", "tst_x.py"),
        os.path.join(src, "repo", "sub", "..", "sub", "tst_x.py"),
        os.path.join(src, "repo") + os.sep + os.sep
          + os.path.join("sub", "tst_x.py")]:
      got = parallel.test_name_from_command(
        'libtbx.python "%s" 7' % path, os.path.join(tmpdir, "bld"), [src])
      assert got == "repo/sub/tst_x.py 7", (path, got)
  finally:
    os.chdir(cwd)
    shutil.rmtree(tmpdir, ignore_errors=True)


def exercise_timings_file():
  """Rows are name,t1,t2,... in UTF-8, sorted by name. A missing file reads
  as empty, a failed write leaves the file as it was, and a file that is
  not a timings file is reported with a Sorry."""
  tmpdir = tempfile.mkdtemp(prefix="tst_timings_")
  try:
    file_name = os.path.join(tmpdir, "timings.csv")
    assert parallel.read_test_timings(file_name) == {}
    timings = {"b/tst_2.py": [1.0, 2.5], "a/tst_\u00e9.py x,y": [0.25]}
    parallel.write_test_timings(file_name, timings)
    with open(file_name, "rb") as f:
      text = f.read()
    assert text == b'"a/tst_\xc3\xa9.py x,y",0.25\nb/tst_2.py,1.00,2.50\n', text
    assert parallel.read_test_timings(file_name) == timings
    # a name that cannot be written stops the write part way
    try:
      parallel.write_test_timings(file_name,
        {"a/tst_1.py": [2.0], "b/tst_\udce9.py": [3.0]})
    except UnicodeError:
      pass
    else:
      raise AssertionError("wrote a name that cannot be encoded")
    assert parallel.read_test_timings(file_name) == timings
    assert os.listdir(tmpdir) == ["timings.csv"], os.listdir(tmpdir)
    # empty rows are skipped
    with open(file_name, "w") as f:
      f.write("a/tst_1.py,1.00\n\nb/tst_2.py,2.00\n")
    got = parallel.read_test_timings(file_name)
    assert got == {"a/tst_1.py": [1.0], "b/tst_2.py": [2.0]}, got
    for bad_value in ["fast", "nan", "inf", "-1.00"]:
      with open(file_name, "w") as f:
        f.write("a/tst_1.py,1.00\nb/tst_2.py,%s\n" % bad_value)
      message = _sorry_from(parallel.read_test_timings, file_name)
      assert "timings.csv" in message and "line 2" in message, message
    # bytes that are not UTF-8; a quote that is never closed, which would
    # turn the rest of the file into one name
    for bad_content in [b"caf\xe9/tst_1.py,1.00\n",
                        b'"a/tst_1.py,1.00\nb/tst_2.py,2.00\n']:
      with open(file_name, "wb") as f:
        f.write(bad_content)
      message = _sorry_from(parallel.read_test_timings, file_name)
      assert "timings.csv" in message, message
  finally:
    shutil.rmtree(tmpdir, ignore_errors=True)


def exercise_add_test_timings():
  """New times go last and the five newest of a test are kept. A negative
  time, of a clock set back, is not added: the file would not read back."""
  timings = {"r/tst_a.py": [1.0, 2.0, 3.0, 4.0, 5.0]}
  parallel.add_test_timings(timings, [
    ('libtbx.python "/x/modules/r/tst_a.py"', 6.0),
    ('libtbx.python "/x/modules/r/tst_b.py"', 0.5),
    ('libtbx.python "/x/modules/r/tst_b.py"', -0.3)])
  assert timings == {
    "r/tst_a.py": [2.0, 3.0, 4.0, 5.0, 6.0],
    "r/tst_b.py": [0.5]}, timings


def exercise_check_timings_file():
  """The path of a timings file is made absolute, a link is followed, and a
  directory, or a file to be written where it cannot be, is a Sorry."""
  tmpdir = os.path.realpath(tempfile.mkdtemp(prefix="tst_timings_"))
  cwd = os.getcwd()
  try:
    os.chdir(tmpdir)
    got = parallel.check_timings_file("timings.csv")
    assert got == os.path.join(tmpdir, "timings.csv"), got
    assert tmpdir in _sorry_from(parallel.check_timings_file, tmpdir)
    missing = os.path.join(tmpdir, "no_such_dir", "timings.csv")
    assert "no_such_dir" in _sorry_from(parallel.check_timings_file, missing)
    # a file that is only read may be one that cannot exist yet
    assert parallel.check_timings_file(missing, will_write=False) == missing
    os.mkdir("shared")
    try:
      os.symlink(os.path.join("shared", "timings.csv"), "link.csv")
    except (OSError, NotImplementedError, AttributeError):
      print("Cannot make links here, skipping the link check")
    else:
      got = parallel.check_timings_file("link.csv")
      assert got == os.path.join(tmpdir, "shared", "timings.csv"), got
  finally:
    os.chdir(cwd)
    shutil.rmtree(tmpdir, ignore_errors=True)


def _commands(root, names):
  result = []
  for name in names:
    result.append('libtbx.python "%s/modules/repo/%s.py"' % (root, name))
  return result


def _one_time_each(names_and_seconds):
  timings = {}
  for name in names_and_seconds:
    timings["repo/%s.py" % name] = [names_and_seconds[name]]
  return timings


def exercise_chunks_by_time():
  """Each test goes to the chunk that is lightest so far, longest test
  first. Dealing 10, 9, 8, 7, 6, 5 s into two chunks gives 10+7+6 and 9+8+5.
  The time of a test is the median of its times, which ignores an outlier;
  a test without times gets the median of the others."""
  cmds = _commands("/r", ["a", "b", "c", "d", "e", "f"])
  timings = _one_time_each(
    {"a": 10.0, "b": 9.0, "c": 8.0, "d": 7.0, "e": 6.0, "f": 5.0})
  a, b, c, d, e, f = cmds
  assert parallel.select_chunk(cmds, 2, 1, timings) == [a, d, e]
  assert parallel.select_chunk(cmds, 2, 2, timings) == [b, c, f]
  # medians 3, 10 and 30; d has no times and gets 10, after b by its name
  timings = {"repo/a.py": [1.0, 100.0, 3.0], "repo/b.py": [10.0],
    "repo/c.py": [20.0, 40.0], "repo/not_in_this_run.py": [999.0]}
  got = parallel.select_chunk([d, c, b, a], 1, 1, timings)
  assert got == [c, b, d, a], got


def exercise_chunks_without_times():
  """Without recorded times, or with times of zero, the chunks get the same
  number of tests, dealt in the order of their names, whatever the commands
  look like."""
  cmds = _commands("/r", ["n1", "n2", "n3", "n4", "n5"])
  n1, n2, n3, n4, n5 = cmds
  assert parallel.select_chunk(cmds, 2, 1, {}) == [n1, n3, n5]
  assert parallel.select_chunk(cmds, 2, 2, {}) == [n2, n4]
  zeros = _one_time_each({"n1": 0.0, "n2": 0.0, "n3": 0.0, "n4": 0.0})
  assert parallel.select_chunk(cmds[:4], 2, 1, zeros) == [n1, n3]
  assert parallel.select_chunk(cmds[:4], 2, 2, zeros) == [n2, n4]
  # repo/a.py comes before repo/z.csh
  script = "/r/modules/repo/z.csh"
  python = 'libtbx.python "/r/modules/repo/a.py"'
  assert parallel.select_chunk([script, python], 2, 1, {}) == [python]
  assert parallel.select_chunk([script, python], 2, 2, {}) == [script]


def exercise_chunks_cover_every_test_once_on_every_machine():
  """Whatever the number of chunks, also more chunks than tests, the chunks
  together hold each test once; and two machines that list the tests in a
  different order, below a different root, select the same chunks."""
  names = []
  seconds = {}
  for i in range(23):
    names.append("t%02d" % i)
    if i % 4 != 0:  # every fourth test has no recorded time
      seconds["t%02d" % i] = float((i * 37) % 11)
  timings = _one_time_each(seconds)
  here = _commands("/Users/me/phenix", names)
  there = _commands("/dev/shm/other", list(reversed(names)))
  for n_chunks in [1, 2, 3, 5, 7, 30]:
    union = []
    for chunk in range(1, n_chunks + 1):
      chunk_here = parallel.select_chunk(here, n_chunks, chunk, timings)
      chunk_there = parallel.select_chunk(there, n_chunks, chunk, timings)
      names_here = []
      for command in chunk_here:
        names_here.append(parallel.test_name_from_command(command))
      names_there = []
      for command in chunk_there:
        names_there.append(parallel.test_name_from_command(command))
      assert names_here == names_there, (n_chunks, chunk)
      union.extend(chunk_here)
    assert sorted(union) == sorted(here), (n_chunks, len(union))


def _make_tests(tmpdir, scripts):
  """Write {file name: source} into <tmpdir>/repo and return that directory."""
  repo = os.path.join(tmpdir, "repo")
  os.makedirs(repo)
  for file_name in scripts:
    with open(os.path.join(repo, file_name), "w") as f:
      f.write(scripts[file_name])
  return repo


def _name(repo, file_name):
  """The name of a test made by _make_tests. It depends on where the
  temporary directory is, so the tests that run the runner ask for it; what
  a name looks like is tested with literal paths above."""
  return parallel.test_name_from_command(
    'libtbx.python "%s"' % os.path.join(repo, file_name))


def _basenames(commands):
  result = []
  for command in commands:
    result.append(os.path.basename(command.split('"')[1]))
  return result


def _run_quietly(args, **keyword_args):
  """Call the runner and keep what it prints out of the output of this test:
  a result line of a run in here, such as the [FAIL] of a test made to fail,
  would read as a result of the run this test is part of. Returns what the
  runner returns and the text."""
  run_command_list = parallel.run_command_list
  def quiet_run_command_list(**arguments):
    arguments["out"] = StringIO()
    return run_command_list(**arguments)
  printed = StringIO()
  stdout = sys.stdout
  parallel.run_command_list = quiet_run_command_list
  sys.stdout = printed
  try:
    result = run_tests_parallel.run(args, **keyword_args)
  finally:
    sys.stdout = stdout
    parallel.run_command_list = run_command_list
  return result, printed.getvalue()


def _list_tests(args):
  return _run_quietly(args, return_list_of_tests=True)


def _run_in_empty_directory(tmpdir, dir_name, args, **keyword_args):
  """Run the tests for real in a new directory under tmpdir. Returns the
  number of failures and the printed text."""
  run_dir = os.path.join(tmpdir, dir_name)
  os.makedirs(run_dir)
  cwd = os.getcwd()
  os.chdir(run_dir)
  try:
    return _run_quietly(args, **keyword_args)
  finally:
    os.chdir(cwd)


def exercise_run_without_the_new_parameters_keeps_its_order():
  """Without timings= and n_chunks= the tests stay in the order the
  directories were given, here not the order of their names."""
  tmpdir = tempfile.mkdtemp(prefix="tst_timings_")
  try:
    repo = _make_tests(tmpdir, {})
    for dir_name, file_name in [("z_dir", "tst_z.py"), ("a_dir", "tst_a.py")]:
      os.makedirs(os.path.join(repo, dir_name))
      with open(os.path.join(repo, dir_name, file_name), "w") as f:
        f.write("")
    tests = _list_tests(["directory=%s" % os.path.join(repo, "z_dir"),
                         "directory=%s" % os.path.join(repo, "a_dir")])[0]
    assert _basenames(tests) == ["tst_z.py", "tst_a.py"], tests
  finally:
    shutil.rmtree(tmpdir, ignore_errors=True)


def exercise_run_lists_the_tests_longest_first():
  """With timings= the list comes back longest first, and listing the tests
  does not touch the file. Two files give two orders, so the order of the
  directory listing cannot satisfy both."""
  tmpdir = tempfile.mkdtemp(prefix="tst_timings_")
  try:
    repo = _make_tests(tmpdir, {"tst_a.py": "", "tst_b.py": "", "tst_c.py": ""})
    file_name = os.path.join(tmpdir, "timings.csv")
    a = _name(repo, "tst_a.py")
    c = _name(repo, "tst_c.py")
    for timings, expected in [
        ({a: [5.0], c: [30.0]}, ["tst_c.py", "tst_b.py", "tst_a.py"]),
        ({a: [30.0], c: [5.0]}, ["tst_a.py", "tst_b.py", "tst_c.py"])]:
      parallel.write_test_timings(file_name, timings)
      tests = _list_tests(
        ["directory=%s" % repo, "timings=%s" % file_name])[0]
      assert _basenames(tests) == expected, tests
      assert parallel.read_test_timings(file_name) == timings
  finally:
    shutil.rmtree(tmpdir, ignore_errors=True)


def exercise_run_records_the_times_of_passed_tests():
  """A run with timings= adds the run time of each passed test on every run,
  records nothing for a test that fails, and keeps a row that appeared in
  the file while the tests were running."""
  tmpdir = tempfile.mkdtemp(prefix="tst_timings_")
  try:
    file_name = os.path.join(tmpdir, "timings.csv")
    # the passing test takes half a second and, on the first run, plays
    # another run that writes the timings file in the meantime
    repo = _make_tests(tmpdir, {
      "tst_ok.py": "import os, time\n"
        "time.sleep(0.5)\n"
        "if not os.path.exists(%r):\n"
        "  with open(%r, 'w') as f:\n"
        "    f.write('other/tst_x.py,7.00\\n')\n" % (file_name, file_name),
      "tst_bad.py": "import sys\nsys.exit(1)\n"})
    args = ["directory=%s" % repo, "timings=%s" % file_name, "nproc=1"]
    ok = _name(repo, "tst_ok.py")
    for n_times in [1, 2]:
      n_failed = _run_in_empty_directory(tmpdir, "run_%d" % n_times, args)[0]
      assert n_failed == 1, n_failed
      timings = parallel.read_test_timings(file_name)
      assert sorted(timings.keys()) == sorted(["other/tst_x.py", ok]), timings
      assert timings["other/tst_x.py"] == [7.0], timings
      assert len(timings[ok]) == n_times, timings
      for seconds in timings[ok]:
        assert 0.5 <= seconds < 300, timings
  finally:
    shutil.rmtree(tmpdir, ignore_errors=True)


def exercise_run_records_a_test_that_passes_on_retry():
  """With max_retries a test that fails and then passes is recorded, with
  the time of the attempt that passed."""
  tmpdir = tempfile.mkdtemp(prefix="tst_timings_")
  try:
    marker = os.path.join(tmpdir, "second_attempt")
    # fails at once the first time, takes 0.4 s and passes the second time
    repo = _make_tests(tmpdir, {
      "tst_flaky.py": "import os, sys, time\n"
        "if not os.path.exists(%r):\n"
        "  open(%r, 'w').close()\n"
        "  sys.exit(1)\n"
        "time.sleep(0.4)\n" % (marker, marker)})
    file_name = os.path.join(tmpdir, "timings.csv")
    n_failed = _run_in_empty_directory(tmpdir, "run", ["directory=%s" % repo,
      "timings=%s" % file_name, "nproc=1", "max_retries=1"])[0]
    assert n_failed == 0, n_failed
    timings = parallel.read_test_timings(file_name)
    flaky = _name(repo, "tst_flaky.py")
    assert list(timings.keys()) == [flaky], timings
    assert len(timings[flaky]) == 1 and timings[flaky][0] >= 0.4, timings
  finally:
    shutil.rmtree(tmpdir, ignore_errors=True)


def exercise_run_reports_its_result_when_the_timings_update_fails():
  """If the timings file cannot be updated at the end, the run says so and
  still returns its number of failures. Here a test damages the file while
  it runs, or puts a directory in its place."""
  damage = {
    "content": "with open(%r, 'w') as f:\n"
      "  f.write('other/tst_x.py,not a time\\n')\n",
    "directory": "import os\nos.mkdir(%r)\n"}
  for kind in sorted(damage):
    tmpdir = tempfile.mkdtemp(prefix="tst_timings_")
    try:
      file_name = os.path.join(tmpdir, "timings.csv")
      repo = _make_tests(tmpdir, {"tst_ok.py": damage[kind] % file_name})
      n_failed, printed = _run_in_empty_directory(tmpdir, "run", [
        "directory=%s" % repo, "timings=%s" % file_name, "nproc=1"])
      assert n_failed == 0, n_failed
      assert "was not updated" in printed, printed
      assert os.path.isdir(file_name) == (kind == "directory")
    finally:
      shutil.rmtree(tmpdir, ignore_errors=True)


def exercise_run_checks_the_timings_path_where_it_starts():
  """A relative timings= names a file where the run was started, also when
  run_in_tmp_dir=True moves the run elsewhere. A place where the file
  cannot be written is reported before any test runs, but not for a chunk,
  which only reads the file."""
  import libtbx.test_utils
  tmpdir = tempfile.mkdtemp(prefix="tst_timings_")
  cwd = os.getcwd()
  open_tmp_directory = libtbx.test_utils.open_tmp_directory
  try:
    repo = _make_tests(tmpdir, {"tst_ok.py": "print('OK')\n"})
    start = os.path.join(tmpdir, "start")
    elsewhere = os.path.join(tmpdir, "elsewhere")
    os.makedirs(start)
    os.makedirs(elsewhere)
    libtbx.test_utils.open_tmp_directory = lambda: elsewhere
    os.chdir(start)
    n_failed = _run_quietly(["directory=%s" % repo, "timings=timings.csv",
      "run_in_tmp_dir=True", "nproc=1"])[0]
    assert n_failed == 0, n_failed
    assert os.path.isfile(os.path.join(start, "timings.csv")), \
      (os.listdir(start), os.listdir(elsewhere))
    args = ["directory=%s" % repo,
      "timings=%s" % os.path.join(tmpdir, "no_such_dir", "timings.csv")]
    assert "no_such_dir" in _sorry_from(_list_tests, args)
    tests = _list_tests(args + ["n_chunks=2", "chunk=1"])[0]
    assert _basenames(tests) == ["tst_ok.py"], tests
  finally:
    libtbx.test_utils.open_tmp_directory = open_tmp_directory
    os.chdir(cwd)
    shutil.rmtree(tmpdir, ignore_errors=True)


def exercise_run_lists_one_chunk():
  """n_chunks= and chunk= select complementary parts of the list, by the
  timings file or, without one, by name; the chunk has to exist, and a
  parameter that does not exist is an error."""
  tmpdir = tempfile.mkdtemp(prefix="tst_timings_")
  try:
    repo = _make_tests(tmpdir, {
      "tst_a.py": "", "tst_b.py": "", "tst_c.py": "", "tst_d.py": ""})
    file_name = os.path.join(tmpdir, "timings.csv")
    parallel.write_test_timings(file_name, {
      _name(repo, "tst_a.py"): [9.0], _name(repo, "tst_b.py"): [1.0],
      _name(repo, "tst_c.py"): [2.0]})
    # 9 s alone balances 2 + 2 + 1 s, tst_d.py getting the median, 2 s; the
    # run says that three of the four tests have times
    args = ["directory=%s" % repo, "n_chunks=2"]
    timed = args + ["timings=%s" % file_name]
    for extra, expected in [
        (timed + ["chunk=1"], ["tst_a.py"]),
        (timed + ["chunk=2"], ["tst_c.py", "tst_d.py", "tst_b.py"]),
        (args + ["chunk=1"], ["tst_a.py", "tst_c.py"]),
        (args + ["chunk=2"], ["tst_b.py", "tst_d.py"])]:
      tests, printed = _list_tests(extra)
      assert _basenames(tests) == expected, (extra, tests)
      assert ("known for 3." in printed) == ("timings=" in extra[-2]), printed
    for extra in [[], ["chunk=3"], ["chunk=0"]]:
      assert "chunk" in _sorry_from(_list_tests, args + extra)
    message = _sorry_from(_list_tests, ["directory=%s" % repo, "chunk=2"])
    assert "chunk" in message, message
    # a misspelled parameter would otherwise make every chunk run all tests
    message = _sorry_from(_list_tests, args[:1] + ["nchunks=2", "chunk=1"])
    assert "nchunks" in message, message
  finally:
    shutil.rmtree(tmpdir, ignore_errors=True)


def _returns(value):
  def function(*args, **keyword_args):
    return list(value)
  return function


def exercise_run_chunks_the_one_at_a_time_tests_of_a_module():
  """The tst_list_parallel tests of a module are divided among the chunks
  with its other tests. Dealt by name, a c p2 make chunk 1 and b p1 chunk 2;
  run() returns the tests of the pool and prints the size of the chunk."""
  a, b, c = _commands("/r", ["a", "b", "c"])
  module_lists = {
    "get_module_tests": [a, b, c],
    "get_module_expected_test_failures": [],
    "get_module_expected_unstable_tests": [],
    "get_module_parallel_tests": _commands("/r", ["p1", "p2"])}
  saved = {}
  tmpdir = tempfile.mkdtemp(prefix="tst_timings_")
  try:
    for function_name in module_lists:
      saved[function_name] = getattr(parallel, function_name)
      setattr(parallel, function_name, _returns(module_lists[function_name]))
    # one known time, of a one-at-a-time test; the others get the same
    file_name = os.path.join(tmpdir, "timings.csv")
    parallel.write_test_timings(file_name, {"repo/p1.py": [5.0]})
    args = ["module=any", "timings=%s" % file_name, "n_chunks=2"]
    chunk_1, printed_1 = _list_tests(args + ["chunk=1"])
    chunk_2, printed_2 = _list_tests(args + ["chunk=2"])
    assert chunk_1 == [a, c], chunk_1
    assert chunk_2 == [b], chunk_2
    assert "Chunk 1 of 2: 3 of 5 tests; run times known for 1." in printed_1
    assert "Chunk 2 of 2: 2 of 5 tests; run times known for 1." in printed_2
  finally:
    for function_name in saved:
      setattr(parallel, function_name, saved[function_name])
    shutil.rmtree(tmpdir, ignore_errors=True)


def exercise_chunked_run():
  """A chunk runs its tests, in directories of their own only for these,
  says in its log which chunk it is, and does not create the timings file.
  A chunk without tests is not an error and leaves a log as well."""
  tmpdir = tempfile.mkdtemp(prefix="tst_timings_")
  try:
    repo = _make_tests(tmpdir, {"tst_a.py": "print('OK')\n", "tst_b.py": ""})
    file_name = os.path.join(tmpdir, "timings.csv")
    args = ["directory=%s" % repo, "timings=%s" % file_name, "nproc=1",
      "n_chunks=3", "run_in_unique_dirs=True"]
    logs = []
    # tst_a.py is chunk 1 and tst_b.py, an expected failure, is chunk 2
    for chunk, expected in [(1, ["tst_a_py"]), (3, [])]:
      run_dir = "chunk_%d" % chunk
      n_failed = _run_in_empty_directory(tmpdir, run_dir,
        args + ["chunk=%d" % chunk],
        expected_failures_from_phenix_regression=["tst_b.py"])[0]
      assert n_failed == 0, n_failed
      made = sorted(os.listdir(os.path.join(tmpdir, run_dir)))
      if "output.xml" in made:
        made.remove("output.xml")
      assert made == ["run_tests_parallel_zlog"] + expected, made
      logs.append(os.path.join(tmpdir, run_dir, "run_tests_parallel_zlog"))
      with open(logs[-1]) as f:
        text = f.read()
      assert parallel.chunk_of_log(text) == (chunk, 3, len(expected), 2), text
      assert ("tst_a.py" in text) == (chunk == 1), text
    assert not os.path.exists(file_name)
    problems = parallel.chunk_log_problems(logs)
    assert len(problems) == 1 and "Chunk 2 of 3" in problems[0], problems
  finally:
    shutil.rmtree(tmpdir, ignore_errors=True)


def _write_log(tmpdir, file_name, names, chunk=None, n_chunks=None,
               n_tests=None, quiet=False):
  """Write the log of a run of the tests repo/<name>.py, all passed, as the
  runner prints it; with chunk given, of that chunk of n_chunks chunks of
  n_tests tests; with quiet, as quiet=True logs it, without the list of
  tests at the top. Returns the path."""
  lines = []
  if chunk is not None:
    lines.append("Chunk %d of %d: %d of %d tests; run times known for 0."
      % (chunk, n_chunks, len(names), n_tests))
  if not quiet:
    lines.append("Running %d tests on 2 processors:" % len(names))
    for name in names:
      lines.append('  libtbx.python "/r/modules/repo/%s.py"' % name)
    lines.append("")
  for name in names:
    lines.append('libtbx.python "/r/modules/repo/%s.py" [OK] 2.0s' % name)
  lines.append("=" * 80)
  path = os.path.join(tmpdir, file_name)
  with open(path, "w") as f:
    f.write("\n".join(lines) + "\n")
  return path


def exercise_chunk_logs_are_checked():
  """chunk_log_problems is silent for the logs of all chunks of one divided
  run and names a missing chunk, a chunk given twice, logs of different
  divisions, a log that is not of a chunk, chunks that do not add up, and a
  test in two logs; timings_from_log= prints what it finds."""
  tmpdir = tempfile.mkdtemp(prefix="tst_timings_")
  try:
    one = _write_log(tmpdir, "one", ["a", "c"], 1, 2, 3)
    two = _write_log(tmpdir, "two", ["b"], 2, 2, 3)
    assert parallel.chunk_log_problems([one, two]) == []
    assert parallel.chunk_log_problems([two, one]) == []
    # a chunk is missing
    problems = parallel.chunk_log_problems([one])
    assert len(problems) == 1, problems
    assert "Chunk 2 of 2" in problems[0] and "missing" in problems[0], problems
    # a chunk is given twice, so its tests are in two logs
    problems = parallel.chunk_log_problems([one, two, one])
    assert len(problems) == 2, problems
    assert "Chunk 1 of 2" in problems[0] and "2 times" in problems[0], problems
    assert "2 tests" in problems[1] and "repo/a.py" in problems[1], problems
    # chunks of a division into three do not belong with these
    other = _write_log(tmpdir, "other", ["d"], 3, 3, 3)
    problems = parallel.chunk_log_problems([one, two, other])
    assert len(problems) == 1 and "different" in problems[0], problems
    # the same chunk numbers, but the tests were divided differently, as
    # with two timings files: a test is in both chunks, or in neither
    two_c = _write_log(tmpdir, "two_c", ["b", "c"], 2, 2, 3)
    problems = parallel.chunk_log_problems([one, two_c])
    assert len(problems) == 2, problems
    assert "4 tests" in problems[0] and "of 3" in problems[0], problems
    assert "1 test" in problems[1] and "repo/c.py" in problems[1], problems
    one_of_4 = _write_log(tmpdir, "one_of_4", ["a", "c"], 1, 2, 4)
    two_of_4 = _write_log(tmpdir, "two_of_4", ["b"], 2, 2, 4)
    problems = parallel.chunk_log_problems([one_of_4, two_of_4])
    assert len(problems) == 1, problems
    assert "3 tests" in problems[0] and "of 4" in problems[0], problems
    # logs of runs that were not divided: only repeated tests matter
    whole = _write_log(tmpdir, "whole", ["a", "b"])
    more = _write_log(tmpdir, "more", ["e"])
    assert parallel.chunk_log_problems([whole, more]) == []
    problems = parallel.chunk_log_problems([one, two, more])
    assert len(problems) == 1 and "more" in problems[0], problems
    # a quiet log has no list of tests at the top; its results count
    hushed = _write_log(tmpdir, "hushed", ["a", "b"], quiet=True)
    assert parallel.chunk_log_problems([hushed]) == []
    problems = parallel.chunk_log_problems([hushed, hushed])
    assert len(problems) == 1 and "2 tests" in problems[0], problems
    # adding the times of logs says if they are not one whole run
    args = ["timings=%s" % os.path.join(tmpdir, "timings.csv")]
    printed = _run_quietly(args + ["timings_from_log=%s" % one,
                                   "timings_from_log=%s" % two])[1]
    assert "Warning" not in printed, printed
    printed = _run_quietly(args + ["timings_from_log=%s" % one])[1]
    assert "Warning: Chunk 2 of 2 is missing." in printed, printed
  finally:
    shutil.rmtree(tmpdir, ignore_errors=True)


run_log = """\
Running 4 tests on 2 processors:
  libtbx.python "/r/modules/cctbx_project/libtbx/tst_a.py"
  libtbx.python "/r/modules/cctbx_project/libtbx/tst_b.py"
  libtbx.python "/r/modules/cctbx_project/libtbx/tst_c.py" 1-5
  libtbx.python "/r/modules/cctbx_project/libtbx/tst_d.py"

libtbx.python "/r/modules/cctbx_project/libtbx/tst_a.py" [OK] 1.5s
libtbx.python "/r/modules/cctbx_project/libtbx/tst_d.py" [WARNING] 4.0s
  Time:  4.00
  Return code: 0
  OKs: 1
libtbx.python "/r/modules/cctbx_project/libtbx/tst_b.py" [FAIL] 2.0s
  Time:  2.00
  Return code: 1
  OKs: 0
libtbx.python "/r/modules/cctbx_project/libtbx/tst_c.py" 1-5 [OK] 30.2s  (passed on attempt 2 of 3)
================================================================================

Tests finished. Elapsed time: 31.70s
"""


def exercise_timings_from_log():
  """timings_from_log= adds one time for each test that exited 0 in a run
  log and runs nothing, so it works in a directory that has files. Over
  several logs a test gets one time, from the first log it passed in."""
  tmpdir = tempfile.mkdtemp(prefix="tst_timings_")
  cwd = os.getcwd()
  try:
    os.chdir(tmpdir)
    with open("run.log", "wb") as f:
      # captured test output need not be UTF-8
      f.write(run_log.replace("OKs: 0", "OKs: 0\n    \xff\xfe").encode(
        "latin-1"))
    args = ["timings=timings.csv", "timings_from_log=run.log"]
    assert _run_quietly(args)[0] == 0
    with open("timings.csv") as f:
      assert not show_diff(f.read(),
        "cctbx_project/libtbx/tst_a.py,1.50\n"
        "cctbx_project/libtbx/tst_c.py 1-5,30.20\n"
        "cctbx_project/libtbx/tst_d.py,4.00\n")
    assert sorted(os.listdir(".")) == ["run.log", "timings.csv"]
    # in a second log tst_b.py passes too, and tst_a.py takes longer
    with open("run_2.log", "w") as f:
      f.write(run_log.replace("[FAIL] 2.0s", "[OK] 7.0s").replace(
        "[OK] 1.5s", "[OK] 9.9s"))
    n_added = parallel.add_test_timings_from_logs(
      "both.csv", ["run.log", "run_2.log"])
    assert n_added == [3, 1], n_added
    with open("both.csv") as f:
      assert not show_diff(f.read(),
        "cctbx_project/libtbx/tst_a.py,1.50\n"
        "cctbx_project/libtbx/tst_b.py,7.00\n"
        "cctbx_project/libtbx/tst_c.py 1-5,30.20\n"
        "cctbx_project/libtbx/tst_d.py,4.00\n")
    # what cannot work is a Sorry: no file to update, a directory as the
    # file, a log that does not exist
    os.mkdir("a_directory")
    _sorry_from(parallel.add_test_timings_from_logs, "a_directory", ["run.log"])
    for bad_args in [["timings_from_log=run.log"],
                     ["timings=timings.csv", "timings_from_log=no_such.log"]]:
      _sorry_from(_run_quietly, bad_args)
  finally:
    os.chdir(cwd)
    shutil.rmtree(tmpdir, ignore_errors=True)


def run(args):
  assert len(args) == 0
  exercise_test_names()
  exercise_test_name_of_a_file_does_not_depend_on_how_it_is_written()
  exercise_timings_file()
  exercise_add_test_timings()
  exercise_check_timings_file()
  exercise_chunks_by_time()
  exercise_chunks_without_times()
  exercise_chunks_cover_every_test_once_on_every_machine()
  exercise_run_without_the_new_parameters_keeps_its_order()
  exercise_run_lists_the_tests_longest_first()
  exercise_run_records_the_times_of_passed_tests()
  exercise_run_records_a_test_that_passes_on_retry()
  exercise_run_reports_its_result_when_the_timings_update_fails()
  exercise_run_checks_the_timings_path_where_it_starts()
  exercise_run_lists_one_chunk()
  exercise_run_chunks_the_one_at_a_time_tests_of_a_module()
  exercise_chunked_run()
  exercise_chunk_logs_are_checked()
  exercise_timings_from_log()


if __name__ == "__main__":
  run(args=sys.argv[1:])
  print(format_cpu_times())
  print("OK")
