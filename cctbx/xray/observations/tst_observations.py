from cctbx import sgtbx, uctbx, crystal, xray
from cctbx.xray import observations
from iotbx.shelx import hklf
from cStringIO import StringIO

def excersise():
  s = """\
  -1  -6   9624.8650 54.9190   1
  -1   6   9749.5660 52.1870  -2
   1   6  -9749.5660 52.1870   1
   1  -6  -8695.6020 53.8100  -2
  -1  -6   8695.6020 53.8100   1
  -1   6   8746.7970 51.3980  -2
   1   6  -8746.7970 51.3980   1
  -1  -6 -12281.3590 49.5060   1
   1   6  12185.9370 47.3950   1
  -1  -6 -1316.0700044.01400   1
   1   6  13 4.7500043.54900   1
  -1  -6 -1479.6380048.66400   2
   1   6  1432.7830051.51700   2
  """
  sg = sgtbx.space_group("P 1")
  uc = uctbx.unit_cell((11,11,13,90,90,120))
  cs = crystal.symmetry(unit_cell=uc, space_group=sg)
  ma = hklf.reader(file_object=StringIO(s))\
         .as_miller_arrays(crystal_symmetry=cs, merge_equivalents=False)
  fo_sq = ma[0]
  batch_numbers = ma[1]
  obs = fo_sq.as_xray_observations(
          scale_indices=batch_numbers.data(),
          twin_fractions=(xray.twin_fraction(0.4,True),))
  measured_cnt = 0
  measured_scale_indices = obs.measured_scale_indices
  for bn in batch_numbers.data():
    if bn > 0:
      assert(measured_scale_indices[measured_cnt]==bn)
      measured_cnt = measured_cnt+1
  assert(measured_cnt == obs.indices.size())

  itr = obs.iterator(0)
  assert(not itr.has_next())
  itr = obs.iterator(1)
  assert(itr.next().h==(-1,6,9))
  itr = obs.iterator(2)
  assert(itr.next().h==(1,-6,-8))

  obs = observations.customized_copy(obs,
          twin_fractions=(xray.twin_fraction(0.7,True),),
          twin_components=(xray.twin_component(
              sgtbx.rot_mx((-1,0,0,0,-1,0,0,0,-1)), 0.25, True),))

  itr = obs.iterator(0)
  assert(itr.has_next())
  assert(itr.next().h==(1,6,-9))
  assert(not itr.has_next())
  itr = obs.iterator(1)
  assert(itr.next().h==(-1,-6,9))
  assert(itr.next().h==(-1,6,9))
  assert(itr.next().h==(1,-6,-9))
  assert(not itr.has_next())

  ts = 1-obs.ref_twin_components[0].value
  ps = 1-obs.ref_twin_fractions[0].value
  itr = obs.iterator(0)
  assert obs.scale(0) == ts*ps
  nv = itr.next()
  assert nv.scale == obs.ref_twin_components[0].value*ps


def run():
  excersise()
  print "OK"

if (__name__ == "__main__"):
  run()
