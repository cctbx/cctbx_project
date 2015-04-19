
from __future__ import division
from mmtbx.ringer import em_rscc
from cctbx import crystal

def exercise () :
  import mmtbx.regression
  from iotbx import file_reader
  from cStringIO import StringIO
  pdb_file = "tmp_em_rscc.pdb"
  map_file = "tmp_em_rscc.map"
  f = open(pdb_file, "w")
  for line in mmtbx.regression.model_1yjp.splitlines() :
    if line.startswith("ATOM") :
      f.write(line + "\n")
  f.close()
  pdb_in = file_reader.any_file(pdb_file).file_object
  symm = crystal.symmetry(
    space_group_symbol="P1",
    unit_cell=(30, 30, 30, 90, 90, 90))
  xrs = pdb_in.input.xray_structure_simple(crystal_symmetry=symm)
  xrs.scattering_type_registry(
    d_min=3.0,
    table="electron")
  fc = xrs.structure_factors(d_min=3.0).f_calc()
  fft_map = fc.fft_map(resolution_factor=1/3).apply_sigma_scaling()
  assert (fft_map.n_real() == (32,32,32))
  fft_map.as_ccp4_map(
    file_name=map_file,
    gridding_first=(-16,-16,-16),
    gridding_last=(15,15,15))
  out = StringIO()
  em_rscc.run(args=[pdb_file, map_file], out=out)
  assert ("""\
PER-RESIDUE CORRELATION:
 A   1  1.0
 A   2  1.0
 A   3  1.0
 A   4  1.0
 A   5  1.0
 A   6  1.0
 A   7  1.0
""" in out.getvalue()), out.getvalue()

if (__name__ == "__main__") :
  exercise()
  print "OK"
