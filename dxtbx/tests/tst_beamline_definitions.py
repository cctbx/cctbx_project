#!/usr/bin/env python
#
#   Copyright (C) 2016 Diamond Light Source, Markus Gerstel
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Tests for the beamline definition database

from __future__ import absolute_import, division
from __future__ import print_function
from builtins import str
import dxtbx.data.beamline_defs as dxbd

def test_lookup_unknown_detector():
  n = dxbd.get_beamline_definition('This is a detector serial number that does not exist')
  assert 'Dummy CIF' in str(n), n
  assert str(n.CIF_block()) == ''
  assert n.CIF_block().__module__ == 'iotbx.cif.model'
  assert str(n.mmCIF_block()) == ''
  assert n.mmCIF_block().__module__ == 'iotbx.cif.model'
  print('OK')

def test_lookup_known_detector():
  n = dxbd.get_beamline_definition('PILATUS 2M, S/N 24-0107 Diamond')
  assert 'Dummy' not in str(n), n
  assert str(n) != ''
  cif = n.CIF_block()
  assert cif.__module__ == 'iotbx.cif.model'
  assert '_diffrn_radiation_type' in str(cif), cif
  cif = n.mmCIF_block()
  assert cif.__module__ == 'iotbx.cif.model'
  assert '_diffrn_radiation.type' in str(cif), cif
  print('OK')

if __name__ == '__main__':
  test_lookup_unknown_detector()
  test_lookup_known_detector()

