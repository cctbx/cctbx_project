"""
Test Phenix/Schrodinger interface

For now it only checks if Schrodinger is installed and the interface requested.
"""
from __future__ import absolute_import, division, print_function
import tempfile
import unittest

from mmtbx.geometry_restraints import external


class TestSchrodingerInterface(unittest.TestCase):

  def test_is_schrodinger_installed(self):
    env = dict(
        SCHRODINGER=tempfile.gettempdir(),
    )
    assert not external.is_schrodinger_installed(env)
    env['PHENIX_SCHRODINGER'] = True
    assert external.is_schrodinger_installed(env)
    del env['SCHRODINGER']
    assert not external.is_schrodinger_installed(env)
    env['SCHRODINGER'] = None
    assert not external.is_schrodinger_installed(env)


if __name__ == '__main__':
  unittest.main(verbosity=0)
