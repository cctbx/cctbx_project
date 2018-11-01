from __future__ import absolute_import, division, print_function

import os
from libtbx import easy_run

def find_new_python3_incompatible_code(module_under_test):
  '''
    Check source code to see if any files violate Python 3 syntax that
    previously did not. Example call:

    def test_find_python3_violations():
      import xia2
      import pytest
      import libtbx.test_utils.python3_regression as py3test
      result = py3test.find_new_python3_incompatible_code(xia2)
      if result is None:
        pytest.skip('No python3 interpreter available')
      elif result:
        pytest.fail(result)

    Known violations are kept in file .known-python3-violations in the
    module directory.
  '''

  # Check for python 3
  try:
    result = easy_run.fully_buffered('python3').raise_if_errors()
  except RuntimeError as e:
    return None

  # File containing list of excluded files
  allowed_broken_files_list = '.known-python3-violations'

  module_path = module_under_test.__path__[0]
  try:
    result = easy_run.fully_buffered(['PYTHONPATH=""', 'python3', '-m', 'compileall', '-x', '\.git', '-q', module_path])
    if len(result.stderr_lines) > 0:
      result.show_stderr()
      result.raise_if_errors()
  except RuntimeError as e:
    return 'Python3 compilation exited with unexpected STDERR output'

  if not result.return_code: # No compilation errors
    return False

  errors = map(lambda x: x.replace(module_path + os.path.sep, '').strip(), '\n'.join(result.stdout_lines).split('***'))
  errors = filter(lambda x: "'" in x, errors)
  broken_files = { error.split("'")[1]: error for error in errors }

  exclusion_file = os.path.join(module_path, allowed_broken_files_list)
  with open(exclusion_file + '.log', 'w') as fh:
    fh.write("\n".join(sorted(broken_files)))
  if os.path.exists(exclusion_file):
    with open(exclusion_file, 'r') as fh:
      excluded_files = fh.read().splitlines()
    broken_files = { filename: broken_files[filename] for filename in broken_files if filename not in excluded_files }

  if not broken_files: # No syntax violations in new files
    return False

  for filename in sorted(broken_files):
    print(broken_files[filename], end="\n\n")

  return "{} file[s] contain newly introduced Python3 syntax errors".format(len(broken_files))
