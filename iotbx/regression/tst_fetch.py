from __future__ import absolute_import, division, print_function

from libtbx.utils import format_cpu_times
import requests
from iotbx.pdb.fetch import fetch, get_link
from libtbx.test_utils import assert_lines_in_text

def exercise_1():
  string_1yjp = fetch(id = "1yjp").read().decode()
  # print("%s" % string_1yjp)
  assert_lines_in_text(string_1yjp, """\
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
""")

def exercise_get_link():
  r = []
  for ft in ['model_pdb', 'model_cif', 'sequence', 'sf', 'map']:
    r.append(get_link('rcsb', ft, '1ab2', emdb_number="1111"))
  # print (r)
  assert r == ['https://files.rcsb.org/pub/pdb/data/structures/divided/pdb/ab/pdb1ab2.ent.gz',
               'https://files.rcsb.org/pub/pdb/data/structures/divided/mmCIF/ab/1ab2.cif.gz',
               'https://www.rcsb.org/fasta/entry/1ab2',
               'https://files.rcsb.org/download/1ab2-sf.cif.gz',
               'https://files.rcsb.org/pub/emdb/structures/EMD-1111/map/emd_1111.map.gz'], r

if (__name__ == "__main__"):
  exercise_get_link()
  exception_occured = False
  try:
    r = requests.get('https://search.rcsb.org/')
  except Exception:
    print("OK but exception.")
    exception_occured = True
  if not exception_occured and r.ok and len(r.text) > 100:
    exercise_1()
    print("OK")
  else:
    print("OK but skipped.")
  print(format_cpu_times())
