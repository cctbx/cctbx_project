# Note: some related tests are in
#       cctbx/regression/tst_sgtbx_special_op_simplifier.py

def exercise():
  from cctbx import xray
  from cctbx import crystal
  from cctbx.array_family import flex
  sc = xray.scatterer
  xs = xray.structure(
    # COD entry 2103456, selected atoms
    crystal_symmetry=crystal.symmetry(
      unit_cell=(10.1053, 10.1053, 10.1053, 90, 90, 90),
      space_group_symbol="P 21 3"),
    scatterers=flex.xray_scatterer((
      sc(label='ER1', site=(-0.16698, 0.66698, 0.33302)),
      sc(label='ER2', site=(0.39957, 0.60043, 0.10043)))))
  #
  from iotbx.shelx import fvar_encoding
  xsm = fvar_encoding.move_sites_if_necessary_for_shelx_fvar_encoding(xs)
  from cStringIO import StringIO
  sio = StringIO()
  xsm.show_scatterers(f=sio)
  from libtbx.test_utils import show_diff
  assert not show_diff(sio.getvalue(), """\
Label, Scattering, Multiplicity, Coordinates, Occupancy, Uiso, Ustar as Uiso
ER1  Er     4 (-0.1670 -0.3330  0.3330) 1.00 0.0000 [ - ]
ER2  Er     4 ( 0.3996 -0.3996  0.1004) 1.00 0.0000 [ - ]
""")
  #
  from iotbx.shelx.errors import error
  from libtbx.test_utils import Exception_expected, approx_equal
  try: fvar_encoding.dev_build_shelx76_fvars(xs)
  except error: pass
  else: raise Exception_expected
  #
  fvars, encoded_sites = fvar_encoding.dev_build_shelx76_fvars(xsm)
  assert approx_equal(fvars, [1, 0.33396, 0.79914])
  assert approx_equal(encoded_sites,
    [[19.5, -19.5, -20.5], [30.5, 29.5, -30.5]])

def run(args):
  assert len(args) == 0
  exercise()
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
