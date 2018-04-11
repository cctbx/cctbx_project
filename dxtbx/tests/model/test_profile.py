from __future__ import absolute_import, division, print_function

from dxtbx.model import ProfileModelFactory
import pytest

def test_profile_modelling():
  grs = pytest.importorskip('dials.algorithms.profile_model.gaussian_rs')

  profile1 = grs.Model(None, 3, 0.1, 0.2, deg=True)
  dictionary = profile1.to_dict()
  profile2 = ProfileModelFactory.from_dict(dictionary)
  assert profile1.sigma_b() == pytest.approx(profile2.sigma_b(), abs=1e-7)
  assert profile1.sigma_m() == pytest.approx(profile2.sigma_m(), abs=1e-7)
