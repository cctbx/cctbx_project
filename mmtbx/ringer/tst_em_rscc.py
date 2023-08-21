
from __future__ import absolute_import, division, print_function
from mmtbx.ringer import em_rscc
from cctbx import crystal

def exercise():
  import mmtbx.regression
  import iotbx.pdb
  from six.moves import cStringIO as StringIO
  pdb_file = "tmp_em_rscc.pdb"
  map_file = "tmp_em_rscc.map"
  f = open(pdb_file, "w")
  for line in mmtbx.regression.model_1yjp.splitlines():
    if line.startswith("ATOM"):
      f.write(line + "\n")
  f.close()
  pdb_in = iotbx.pdb.input(pdb_file)
  symm = crystal.symmetry(
    space_group_symbol="P1",
    unit_cell=(30, 30, 30, 90, 90, 90))
  xrs = pdb_in.xray_structure_simple(crystal_symmetry=symm)
  xrs.scattering_type_registry(
    d_min=3.0,
    table="electron")
  fc = xrs.structure_factors(d_min=3.0).f_calc()
  fft_map = fc.fft_map(resolution_factor=1/3).apply_sigma_scaling()
  i,j,k = fft_map.n_real()
  s = i//2
  f = i//2-1
  print(i,j,k,s,f)
  fft_map.as_ccp4_map(
    file_name=map_file,
    gridding_first=(-s,-s,-s),
    gridding_last=(f,f,f))
  out = StringIO()
  em_rscc.run(args=[pdb_file, map_file], out=out)
  for line in out.getvalue().splitlines():
    if line.find(" A  ")==-1: continue
    assert abs(float(line.split()[2])-1)<0.1

if (__name__ == "__main__"):
  exercise()
  print("OK")
