"""Parity test: iotbx.cif.reader with engine='ucif' vs engine='xcif'.

Parses a corpus of CIF fixtures with both engines and asserts the resulting
iotbx.cif.model.cif objects are structurally equivalent (pair items, loops,
save frames). Order-insensitive at the block/loop level; order-preserving
inside loop columns.
"""
from __future__ import absolute_import, division, print_function

import os
import sys

import libtbx.load_env
import iotbx.cif


# ---------------------------------------------------------------------------
# Inline fixtures — cover corner cases the real files may not exercise.
# ---------------------------------------------------------------------------

CIF_SIMPLE_PAIRS = """\
data_pairs_only
_one 1
_two 'two words'
_three 3.14
_four .
_five ?
"""

CIF_LOOP_ONLY = """\
data_loop_only
loop_
_atom.label
_atom.x
_atom.y
_atom.z
 A  1.0  2.0  3.0
 B  4.0  5.0  6.0
 C  7.0  8.0  9.0
"""

CIF_INTERLEAVED = """\
data_interleaved
_leading_key 'before loop'
loop_
_site.id
_site.occ
 1  1.0
 2  0.5
_middle_key 'between loops'
loop_
_other.id
_other.val
 1  a
 2  b
_trailing_key 'after loops'
"""

CIF_MULTIBLOCK = """\
data_alpha
_id alpha
_flag y

data_beta
_id beta
loop_
_beta.i
_beta.v
 1 one
 2 two
"""

CIF_SAVE_FRAMES = """\
data_with_saves
_top_level 'outer'
save_frame1
_inner.key1 value1
_inner.key2 value2
loop_
_inner_loop.a
_inner_loop.b
 1 x
 2 y
save_

save_frame2
_other.key other_value
save_
"""

CIF_UNDERSCORE_IN_NAME = """\
data_foo_bar_baz
_key val
"""

CIF_QUOTED_AND_SEMI = """\
data_quotes
_single_quoted 'with spaces'
_double_quoted "also spaces"
_semicolon_text
;
line one
line two
;
_normal normal
"""


INLINE_FIXTURES = [
  ("simple_pairs", CIF_SIMPLE_PAIRS),
  ("loop_only", CIF_LOOP_ONLY),
  ("interleaved", CIF_INTERLEAVED),
  ("multiblock", CIF_MULTIBLOCK),
  ("save_frames", CIF_SAVE_FRAMES),
  ("underscore_in_name", CIF_UNDERSCORE_IN_NAME),
  ("quoted_and_semi", CIF_QUOTED_AND_SEMI),
]


# ---------------------------------------------------------------------------
# File fixtures from xcif/regression/
# ---------------------------------------------------------------------------

def _xcif_regression_dir():
  path = libtbx.env.dist_path("xcif")
  return os.path.join(path, "regression")

def _file_fixtures():
  d = _xcif_regression_dir()
  out = []
  for name in ("example.cif", "1yjp.cif", "1yjp-sf.cif"):
    p = os.path.join(d, name)
    if os.path.isfile(p):
      out.append((name, p))
  return out


# ---------------------------------------------------------------------------
# Structural comparator
# ---------------------------------------------------------------------------

def _fail(context, msg):
  return "parity mismatch at %s: %s" % (context, msg)

def _compare_loop(ctx, lp_u, lp_x):
  tags_u = list(lp_u.keys())
  tags_x = list(lp_x.keys())
  # Tag set must match (order within a loop is observed but allow either
  # engine to enumerate in a different cif_model_builder order — compare as
  # set first, then per-tag column content).
  if set(tags_u) != set(tags_x):
    return _fail(ctx, "loop tags differ: ucif=%s xcif=%s" % (
      sorted(tags_u), sorted(tags_x)))
  for tag in tags_u:
    col_u = list(lp_u[tag])
    col_x = list(lp_x[tag])
    if col_u != col_x:
      return _fail(ctx + "/" + tag,
                   "column differs: ucif[:5]=%s xcif[:5]=%s lengths=%d,%d" %
                   (col_u[:5], col_x[:5], len(col_u), len(col_x)))
  return None

def _compare_block_like(ctx, blk_u, blk_x):
  # Pair items + loop names: check ORDER first via block._set.
  # _set is the OrderedSet that preserves source-order interleave
  # across pair items and loops within a block; it drives
  # str(model) / block.show() output. Comparing _items and .loops
  # as separate dicts cannot catch pair/loop interleave mismatches.
  set_u = list(blk_u._set)
  set_x = list(blk_x._set)
  if set_u != set_x:
    return _fail(ctx,
      "block key order differs (source-order mismatch):\n"
      "  ucif=%s\n  xcif=%s" % (set_u, set_x))
  # Pair items (content; order already verified via _set above)
  items_u = dict(blk_u._items)
  items_x = dict(blk_x._items)
  if set(items_u.keys()) != set(items_x.keys()):
    only_u = set(items_u) - set(items_x)
    only_x = set(items_x) - set(items_u)
    return _fail(ctx, "pair keys differ: only-in-ucif=%s only-in-xcif=%s" % (
      sorted(only_u), sorted(only_x)))
  for k, v_u in items_u.items():
    v_x = items_x[k]
    if v_u != v_x:
      return _fail(ctx + "/" + k,
                   "pair value differs: ucif=%r xcif=%r" % (v_u, v_x))
  # Loops (content)
  loops_u = blk_u.loops
  loops_x = blk_x.loops
  if set(loops_u.keys()) != set(loops_x.keys()):
    only_u = set(loops_u) - set(loops_x)
    only_x = set(loops_x) - set(loops_u)
    return _fail(ctx, "loop names differ: only-in-ucif=%s only-in-xcif=%s" % (
      sorted(only_u), sorted(only_x)))
  for name in loops_u:
    err = _compare_loop(ctx + "/loop:" + name, loops_u[name], loops_x[name])
    if err: return err
  return None

def _compare_cif(ctx, cif_u, cif_x):
  names_u = list(cif_u.keys())
  names_x = list(cif_x.keys())
  if set(names_u) != set(names_x):
    only_u = set(names_u) - set(names_x)
    only_x = set(names_x) - set(names_u)
    return _fail(ctx, "block names differ: only-in-ucif=%s only-in-xcif=%s" %
                 (sorted(only_u), sorted(only_x)))
  for name in names_u:
    blk_u = cif_u[name]
    blk_x = cif_x[name]
    err = _compare_block_like("%s/block:%s" % (ctx, name), blk_u, blk_x)
    if err: return err
    # Save frames (block has .saves, save does not)
    saves_u = getattr(blk_u, "saves", {})
    saves_x = getattr(blk_x, "saves", {})
    if set(saves_u.keys()) != set(saves_x.keys()):
      only_u = set(saves_u) - set(saves_x)
      only_x = set(saves_x) - set(saves_u)
      return _fail(
        "%s/block:%s" % (ctx, name),
        "save frames differ: only-in-ucif=%s only-in-xcif=%s" % (
          sorted(only_u), sorted(only_x)))
    for sf_name in saves_u:
      err = _compare_block_like(
        "%s/block:%s/save:%s" % (ctx, name, sf_name),
        saves_u[sf_name], saves_x[sf_name])
      if err: return err
  return None


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

def _parse_both(input_string=None, file_path=None):
  if file_path is not None:
    m_u = iotbx.cif.reader(file_path=file_path, engine="ucif").model()
    m_x = iotbx.cif.reader(file_path=file_path, engine="xcif").model()
  else:
    m_u = iotbx.cif.reader(input_string=input_string, engine="ucif").model()
    m_x = iotbx.cif.reader(input_string=input_string, engine="xcif").model()
  return m_u, m_x

def _run_case(label, *, input_string=None, file_path=None):
  try:
    m_u, m_x = _parse_both(input_string=input_string, file_path=file_path)
  except Exception as e:
    return "%s: parse raised: %s" % (label, e)
  err = _compare_cif(label, m_u, m_x)
  return err

def run():
  failures = []
  for label, text in INLINE_FIXTURES:
    err = _run_case(label, input_string=text)
    if err:
      failures.append(err)
    else:
      print("OK:", label)

  for name, path in _file_fixtures():
    err = _run_case("file:" + name, file_path=path)
    if err:
      failures.append(err)
    else:
      print("OK: file:" + name)

  if failures:
    print()
    print("FAILURES:")
    for f in failures:
      print(" -", f)
    sys.exit(1)
  print()
  print("OK all parity cases")

if __name__ == "__main__":
  run()
