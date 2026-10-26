from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from six.moves import cStringIO as StringIO
import mmtbx.refinement.rigid_body as rigid_body

""" Rigid-body refinement with the LLGI target: one zone at the full
resolution limit (unless llgi_resolution_zones is True), llgi at every
resolution; and fmodel.resolution_filter keeping the tNCS arrays of
llgi_data aligned with the reflections. """

def zones(target_name, d_min=2.0, auto_switch=6.0, **kwds):
  params = rigid_body.master_params.extract()
  for k, v in kwds.items():
    setattr(params, k, v)
  # spread d-spacings uniformly in 1/d^2 between 30 A and d_min
  s_lo, s_hi = 1 / 30.**2, 1 / d_min**2
  d = flex.double([(s_lo + (s_hi - s_lo) * i / 1999.)**-0.5
    for i in range(2000)])
  params.target_auto_switch_resolution = auto_switch
  log = StringIO()
  d_mins, targets = rigid_body.zone_setup(d_spacings=d, n_bodies=1,
    params=params, original_target_name=target_name, log=log)
  return d_mins, targets

def exercise_llgi_single_zone():
  params = rigid_body.master_params.extract()
  assert params.number_of_zones == 5 and params.high_resolution == 3.0
  assert params.llgi_resolution_zones is False
  d_mins, targets = zones("llgi")
  assert len(d_mins) == 1 and abs(d_mins[0] - 2.0) < 1.e-6, d_mins
  assert targets == ["llgi"], targets
  # the old behaviour
  d_mins, targets = zones("llgi", llgi_resolution_zones=True)
  assert len(d_mins) == 5, d_mins
  assert targets[0] == "ls_wunit_k1" and targets[-1] == "llgi", targets
  assert d_mins[-1] >= 3.0 - 1.e-6, d_mins
  # other targets are unchanged
  d_mins, targets = zones("ml")
  assert len(d_mins) == 5, d_mins
  assert targets[0] == "ls_wunit_k1" and targets[-1] == "ml", targets
  assert d_mins == zones("ml", llgi_resolution_zones=True)[0]
  # an explicit target is still honoured with the single zone
  d_mins, targets = zones("llgi", target="ls_wunit_k1")
  assert len(d_mins) == 1 and targets == ["ls_wunit_k1"], targets

def exercise_resolution_filter_keeps_tncs():
  import tst_llgi_tncs as t
  h, xrs = t.build_model()
  fmodel = t.cctbx_fmodel(h, xrs)
  ld = fmodel.llgi_data()
  assert ld.tncs_state is not None
  n = fmodel.f_obs().size()
  for f in [fmodel.resolution_filter(d_min=4.0), fmodel.deep_copy()]:
    new = f.llgi_data()
    assert new.tncs_state is not None
    assert new.tncs_u.indices().all_eq(f.f_obs().indices())
    assert new.tncs_t.indices().all_eq(f.f_obs().indices())
  f = fmodel.resolution_filter(d_min=4.0)
  keep = fmodel.f_obs().d_spacings().data() >= 4.0 - 1.e-9
  assert 0 < f.f_obs().size() < n
  assert f.f_obs().size() == keep.count(True)
  new = f.llgi_data()
  assert new.tncs_u.data().all_eq(ld.tncs_u.data().select(keep))
  assert new.tncs_t.data().all_eq(ld.tncs_t.data().select(keep))

if __name__ == "__main__":
  exercise_llgi_single_zone()
  exercise_resolution_filter_keeps_tncs()
  print("OK")
