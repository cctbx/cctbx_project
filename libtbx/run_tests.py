from __future__ import absolute_import, division, print_function
from libtbx import test_utils
import sys
import libtbx.load_env

tst_list_base = [
  "$D/metric_prefixes.py",
  "$D/tst_utils.py",
  "$D/tst_word_index_generator.py",
  "$D/test_utils/__init__.py",
  "$D/queuing_system_utils/pbs_utils.py",
  "$D/queuing_system_utils/sge_utils.py",
  "$D/introspection.py",
  "$D/tst_thread_utils.py",
  "$D/tst_easy_mp.py",
  "$D/tst_easy_mp_state.py",
  "$D/tst_easy_pickle.py",
  "$D/tst_fully_buffered_timeout.py",
  "$D/tst_add_docstrings_with_ai.py",
  "$D/tst_scheduling.py",
  "$D/easy_run.py",
  "$D/tst_containers.py",
  "$D/tst_path.py",
  "$D/tst_math_utils.py",
  "$D/assert_utils.py",
  "$D/tst_str_utils.py",
  "$D/table_utils.py",
  "$D/tst_dlite.py",
  "$D/phil/tst_tokenizer.py",
  "$D/phil/tst.py",
  "$D/phil/tst_experimental.py",
  "$D/phil/tst_interface.py",
  "$D/tst_object_oriented_patterns.py",
  "$D/find_reference_cycles.py",
  "$D/tst_binary_search.py",
  "$D/tst_topological_sort.py",
  "$D/clusterTests.py",
  "$D/tst_citations.py",
  "$D/tst_python_code_parsing.py",
  "$D/tst_representation.py",
  "$D/tst_find_unused_imports.py",
  "$D/tst_program_template.py",
  "$D/tst_version.py",
  '$D/tst_easy_mp_multicore.py',
  ]

# generally failing tests
tst_list_fail = [
  "$D/tst_xmlrpc_utils.py",
  ]
# failing tests on Windows, Python 2.7
tst_list_windows_fail = [
  "$D/tst_runtime_utils.py",
  ]
if sys.platform == 'win32':
  tst_list_fail += tst_list_windows_fail
else:
  tst_list_base += tst_list_windows_fail

# final lists
tst_list = tst_list_base
tst_list_expected_failures = tst_list_fail

def run():
  build_dir = None
  dist_dir = libtbx.env.dist_path("libtbx")

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
