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

indices = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (2, 1, 1)]


def write_table(file_name, xs):
  """ A table for the structure as it stands, entry i marked with i+1.

  Labelling the entries by construction means a table which reads back in the
  wrong order is caught, and not merely one which fails to read at all.
  """
  scatterers = xs.scatterers()
  with open(file_name, 'w') as out:
    out.write('Title: labelled table for the scatterer id test')
    out.write('\nScatterer_ids:')
    for sc in scatterers:
      out.write(' %s' % sc.get_id_big(sc.get_part()).to_hex_string())
    out.write('\nData:')
    for h in indices:
      out.write('\n%d %d %d' % h)
      for i in range(scatterers.size()):
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
    assert 'Could not locate' in str(e), e
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


def run():
  flex.set_random_seed(0)
  random.seed(0)
  file_name = 'tst_table_based_tmp.tsc'
  try:
    exercise_table_survives_a_refinement(file_name)
    exercise_another_structure_is_refused(file_name)
    exercise_an_ambiguous_match_is_refused(file_name)
    exercise_the_generated_table_is_readable(file_name)
  finally:
    if os.path.isfile(file_name):
      os.remove(file_name)
  print('OK')


if __name__ == '__main__':
  run()
