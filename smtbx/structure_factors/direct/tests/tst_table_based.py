""" A tabulated scattering table must survive its structure being refined.

Each entry of a .tsc carries a scatterer id built from the element, the part
and the fractional coordinate. The coordinate is encoded to about 1e-8 of a
cell edge, so a single refinement cycle invalidates every id in the file and
the reader has to fall back on matching the position the id carries. This
checks that the fallback recognises a structure which has moved by a refinement
without ever confusing one scatterer with another.
"""
from __future__ import absolute_import, division, print_function

from cctbx import crystal, sgtbx, xray
from cctbx.array_family import flex
import smtbx.development
import smtbx.structure_factors.direct as structure_factors
import os
import random
import struct

indices = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (2, 1, 1)]


def write_table(file_name, xs, columns=None):
  """ A table for the structure as it stands, entry i marked with i+1.

  Labelling the entries by construction means a table which reads back in the
  wrong order is caught, and not merely one which fails to read at all.

  columns: which scatterers to give a column, defaulting to all of them. A
  table is free to cover only part of a structure, and the entries keep the
  label of the scatterer they describe either way.
  """
  scatterers = xs.scatterers()
  if columns is None:
    columns = range(scatterers.size())
  with open(file_name, 'w') as out:
    out.write('Title: labelled table for the scatterer id test')
    out.write('\nScatterer_ids:')
    for i in columns:
      out.write(' %s' % scatterers[i].get_id_big(
        scatterers[i].get_part()).to_hex_string())
    out.write('\nData:')
    for h in indices:
      out.write('\n%d %d %d' % h)
      for i in columns:
        out.write(' %.6f,%.6f' % (i + 1, 0))
    out.write('\n')


def read_table(file_name, xs):
  """ What the table says each scatterer contributes, in structure order. """
  contribution = structure_factors.ext.table_based_scatterer_contribution.build(
    xs.unit_cell(), xs.scatterers(), file_name, xs.space_group(),
    not xs.space_group().is_origin_centric())
  return [[contribution.get(i, h).real for i in range(xs.scatterers().size())]
          for h in indices]


def resolves(file_name, xs):
  """ Whether the table reads back onto xs with every entry where it belongs.
  """
  expected = [float(i + 1) for i in range(xs.scatterers().size())]
  try:
    rows = read_table(file_name, xs)
  except RuntimeError as e:
    # an entry that cannot be placed leaves its scatterer uncovered, which
    # without a registry to fall back on is refused outright
    assert 'does not cover the whole structure' in str(e), e
    return False
  for row in rows:
    assert row == expected, (row, expected)
  return True


def moved_by(xs, distance, which=None):
  """ A copy of xs with the named scatterers displaced by distance Angstrom.
  """
  out = xs.deep_copy_scatterers()
  cell = out.unit_cell()
  for i, sc in enumerate(out.scatterers()):
    if which is not None and i not in which:
      continue
    direction = flex.random_double(3) - 0.5
    step = cell.fractionalize(tuple(direction / direction.norm() * distance))
    sc.site = tuple(a + b for a, b in zip(sc.site, step))
  return out


def exercise_table_survives_a_refinement(file_name):
  """ Shifts of the size a refinement applies shall not lose a scatterer. """
  xs = smtbx.development.random_xray_structure(
    sgtbx.space_group_info('P1'), elements=['C'] * 6 + ['H'] * 6,
    u_iso=0.05, random_u_iso=False)
  write_table(file_name, xs)
  assert resolves(file_name, xs), 'the structure it was written for'
  # every id is invalidated by any shift at all, however small
  assert not any(a == b for a, b in zip(
    [sc.get_id_big(sc.get_part()).to_hex_string() for sc in xs.scatterers()],
    [sc.get_id_big(sc.get_part()).to_hex_string()
     for sc in moved_by(xs, 1e-4).scatterers()]))
  for shift in (1e-4, 1e-3, 1e-2, 0.05, 0.1):
    assert resolves(file_name, moved_by(xs, shift)), \
      'a structure which moved %g A' % shift


def exercise_another_structure_is_refused(file_name):
  """ A table which does not belong to the structure shall not be used. """
  xs = smtbx.development.random_xray_structure(
    sgtbx.space_group_info('P1'), elements=['C'] * 6 + ['N'] * 6,
    u_iso=0.05, random_u_iso=False)
  write_table(file_name, xs)
  assert not resolves(file_name, moved_by(xs, 2.0)), \
    'a table for another structure was accepted'


def exercise_an_ambiguous_match_is_refused(file_name):
  """ Two candidates comparably close shall be an error, not a guess.

  Picking one of them would put a whole column of contributions on the wrong
  atom, and nothing downstream would notice: the table is the right shape and
  every number in it is plausible. Refusing is the only safe answer.

  A cell edge of 10 A makes the fractional offsets below read directly as
  hundredths of an Angstrom.
  """
  cs = crystal.symmetry(unit_cell=(10, 10, 10, 90, 90, 90),
                        space_group_symbol='P1')
  sites = [(0.1, 0.1, 0.1), (0.125, 0.1, 0.1), (0.6, 0.6, 0.6)]
  xs = xray.structure(crystal_symmetry=cs)
  for i, site in enumerate(sites):
    xs.add_scatterer(
      xray.scatterer(label='C%d' % i, site=site, u=0.05))
  write_table(file_name, xs)
  assert resolves(file_name, xs), 'the structure it was written for'

  # leave the first entry with two candidates 0.10 and 0.12 A away: near
  # enough to be a refined position, too near each other to tell apart
  moved = xs.deep_copy_scatterers()
  moved.scatterers()[0].site = (0.110, 0.1, 0.1)
  moved.scatterers()[1].site = (0.088, 0.1, 0.1)
  assert not resolves(file_name, moved), \
    'a scatterer was picked out of two comparably good candidates'


def exercise_the_generated_table_is_readable(file_name):
  """ The table the isotropic helper writes shall read back in.

  Writer and reader agree on the id format by convention alone, and nothing
  else in the tree calls the helper, so this is the only thing keeping the two
  from drifting apart. An expanded table is read through get_full, which is not
  wrapped, so building it is as far as this can check from Python.
  """
  for symbol in ('P1', 'P21/c'):
    xs = smtbx.development.random_xray_structure(
      sgtbx.space_group_info(symbol), elements=['C'] * 3 + ['O'] * 2,
      u_iso=0.05, random_u_iso=False)
    structure_factors.generate_isc_table_file(file_name, xs, indices)
    structure_factors.ext.table_based_scatterer_contribution.build(
      xs.unit_cell(), xs.scatterers(), file_name, xs.space_group(),
      not xs.space_group().is_origin_centric())


def exercise_a_partial_table_falls_back_on_spherical(file_name):
  """ Atoms the table misses shall get spherical form factors, not zeroes.

  A zero would be silent: the refinement would run and simply pretend those
  atoms are not there. So this checks the value each uncovered scatterer gets
  against what it would have had with no table at all, and checks that the
  covered ones still come from the table.
  """
  xs = smtbx.development.random_xray_structure(
    sgtbx.space_group_info('P1'), elements=['C'] * 4 + ['O'] * 2,
    u_iso=0.05, random_u_iso=False)
  n = xs.scatterers().size()
  missing = [1, 4]
  covered = [i for i in range(n) if i not in missing]
  write_table(file_name, xs, columns=covered)

  # without a registry there is nothing to fall back on, so it is an error
  try:
    read_table(file_name, xs)
    raise AssertionError('a partial table was accepted with no fallback')
  except RuntimeError as e:
    assert 'does not cover the whole structure' in str(e), e

  contribution = structure_factors.ext.table_based_scatterer_contribution.\
    build_with_fallback(
      xs.unit_cell(), xs.scatterers(), file_name, xs.space_group(),
      not xs.space_group().is_origin_centric(), xs.scattering_type_registry())
  assert list(contribution.scatterers_not_in_table()) == missing, \
    list(contribution.scatterers_not_in_table())

  spherical = structure_factors.ext.isotropic_scatterer_contribution(
    xs.scatterers(), xs.scattering_type_registry())
  for h in indices:
    d_star_sq = xs.unit_cell().d_star_sq(h)
    contribution.at_d_star_sq(d_star_sq)
    spherical.at_d_star_sq(d_star_sq)
    for i in range(n):
      got = contribution.get(i, h)
      if i in missing:
        expected = spherical.get(i, h)
        assert abs(got - expected) < 1e-12, (i, h, got, expected)
        assert abs(got) > 1e-6, (i, h, got)
      else:
        assert abs(got - complex(i + 1, 0)) < 1e-6, (i, h, got)


ENCODE_SCALE = 2147483647 / 16.0


def write_tscb(file_name, xs, columns):
  """ The same table in the binary format, which has its own reader.

  The two readers keep separate bookkeeping of which scatterers a table
  covered, and the binary one is what a NoSpherA2 run actually produces, so
  exercising only the text one would leave the used path untested.
  """
  scatterers = xs.scatterers()
  header = b'SCATTERER_IDS\nAD: FALSE\n'
  with open(file_name, 'wb') as out:
    out.write(struct.pack('<i', len(header)))
    out.write(header)
    out.write(struct.pack('<i', len(columns)))
    for i in columns:
      sc = scatterers[i]
      # the 16 bytes scatterer_id_big reads straight off the stream
      out.write(struct.pack(
        '<iiihBB',
        int(round(sc.site[0] * ENCODE_SCALE)),
        int(round(sc.site[1] * ENCODE_SCALE)),
        int(round(sc.site[2] * ENCODE_SCALE)),
        sc.get_part(), sc.electron_count(), 0))
    out.write(struct.pack('<i', len(indices)))
    for h in indices:
      out.write(struct.pack('<iii', *h))
      for i in columns:
        out.write(struct.pack('<dd', i + 1, 0.0))


def exercise_a_partial_tscb_falls_back_on_spherical(file_name):
  """ As above, for the binary format. Values are unreachable from Python
  here -- a binary table is always expanded, so it is read through get_full,
  which is not wrapped -- so this checks the coverage bookkeeping only.
  """
  xs = smtbx.development.random_xray_structure(
    sgtbx.space_group_info('P1'), elements=['C'] * 4 + ['O'] * 2,
    u_iso=0.05, random_u_iso=False)
  n = xs.scatterers().size()
  ext = structure_factors.ext.table_based_scatterer_contribution

  write_tscb(file_name, xs, list(range(n)))
  full = ext.build(xs.unit_cell(), xs.scatterers(), file_name,
                   xs.space_group(), not xs.space_group().is_origin_centric())
  assert list(full.scatterers_not_in_table()) == []

  missing = [1, 4]
  write_tscb(file_name, xs, [i for i in range(n) if i not in missing])
  try:
    ext.build(xs.unit_cell(), xs.scatterers(), file_name, xs.space_group(),
              not xs.space_group().is_origin_centric())
    raise AssertionError('a partial .tscb was accepted with no fallback')
  except RuntimeError as e:
    assert 'does not cover the whole structure' in str(e), e

  partial = ext.build_with_fallback(
    xs.unit_cell(), xs.scatterers(), file_name, xs.space_group(),
    not xs.space_group().is_origin_centric(), xs.scattering_type_registry())
  assert list(partial.scatterers_not_in_table()) == missing, \
    list(partial.scatterers_not_in_table())


def run():
  flex.set_random_seed(0)
  random.seed(0)
  file_name = 'tst_table_based_tmp.tsc'
  binary_file_name = 'tst_table_based_tmp.tscb'
  try:
    exercise_table_survives_a_refinement(file_name)
    exercise_another_structure_is_refused(file_name)
    exercise_an_ambiguous_match_is_refused(file_name)
    exercise_the_generated_table_is_readable(file_name)
    exercise_a_partial_table_falls_back_on_spherical(file_name)
    exercise_a_partial_tscb_falls_back_on_spherical(binary_file_name)
  finally:
    for name in (file_name, binary_file_name):
      if os.path.isfile(name):
        os.remove(name)
  print('OK')


if __name__ == '__main__':
  run()
