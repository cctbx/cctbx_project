from __future__ import absolute_import, division, print_function

from libtbx.test_utils import Exception_expected
import sys


def exercise_origin_id_registered():
  """TDD step 1: 'reference hydrogen bonds' must be registered as an origin_id.

  The new origin_id is added to cctbx/geometry_restraints/auto_linking_types.py
  so that bond/angle proxies derived from reference-model H-bonds can be tagged
  distinctly from secondary-structure H-bonds.
  """
  from cctbx.geometry_restraints.linking_class import linking_class
  lc = linking_class()
  oid = lc.get_origin_id('reference hydrogen bonds')
  assert isinstance(oid, int), \
    "expected int origin_id, got %r" % (oid,)
  assert oid > 0, "origin_id must be > 0 (0 is reserved for 'covalent geometry')"
  # distinct from the existing SS 'hydrogen bonds' origin_id
  ss_oid = lc.get_origin_id('hydrogen bonds')
  assert oid != ss_oid, \
    "new origin_id must differ from SS 'hydrogen bonds' (%d)" % ss_oid


def run(args):
  assert not args, args
  exercise_origin_id_registered()
  print("OK")


if __name__ == "__main__":
  run(sys.argv[1:])
