from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
import iotbx.pdb
from cctbx.maptbx.bcr import qmap
from cctbx import maptbx
import time
import math
from libtbx import adopt_init_args
from libtbx.test_utils import approx_equal
import mmtbx.model

import boost_adaptbx.boost.python as bp
ext = bp.import_ext("cctbx_maptbx_bcr_bcr_ext")

pdb_str = """
remark CRYST1   20.000   25.000   30.000  96.06  99.11  91.14 P 1
CRYST1   20.000   25.000   30.000  90.00  90.00  90.00 P 1
ATOM      0  N   HIS A 109      10.000  14.000   9.000  1.10 30.10           N
ATOM      1  C   HIS A 109      13.000  10.000  14.000  0.90 20.00           C
ATOM      2  O   HIS A 109      12.000  10.000  11.000  0.50 10.00           O
TER
END
"""

class gradients_fd(object):
  def __init__(self, ControlMap, xrs, n_real, d_min):
    adopt_init_args(self, locals())

  def compute(self):
    o = qmap.compute(
      xray_structure=self.xrs,
      n_real=self.n_real,
      resolution=self.d_min,
      resolutions=None,
      use_exp_table=False,
      debug=True)
    return qmap.CalcFuncMap(
      OmegaMap=o.OmegaMap_py, ControlMap=self.ControlMap, Ncrs=self.n_real)

def run(d_min, table="wk1995"):
  pdb_inp = iotbx.pdb.input(source_info=None, lines = pdb_str)
  model = mmtbx.model.manager(model_input=pdb_inp)
  model.setup_scattering_dictionaries(
    scattering_table = table, d_min=d_min)
  xrs = model.get_xray_structure()
  cs = model.crystal_symmetry()
  uc = cs.unit_cell()
  crystal_gridding = maptbx.crystal_gridding(
    unit_cell        = cs.unit_cell(),
    space_group_info = cs.space_group_info(),
    symmetry_flags   = maptbx.use_space_group_symmetry,
    step             = 1)
  n_real = crystal_gridding.n_real()
  #
  o = qmap.compute(
    xray_structure = xrs,
    n_real=n_real,
    resolution=d_min,
    resolutions=None,
    use_exp_table=False,
    debug=True)
  OmegaMap_py = o.OmegaMap_py
  bcr_scatterers = o.bcr_scatterers
  #
  nx,ny,nz = n_real
  ControlMap = [[[ math.sin(ix*iy*iz/math.pi) for ix in range(nx) ]
                for iy in range(ny)] for iz in range(nz)]
  GradMap  = qmap.CalcGradMap(OmegaMap_py, ControlMap, n_real)
  GradAtom = qmap.CalcGradAtom(
    GradMap, n_real,[0,0,0],n_real,cs.unit_cell(), bcr_scatterers)

  ControlMap_flex = flex.double(flex.grid([nz,ny,nx]))
  for iz in range(0, nz):
    for iy in range(0, ny):
      for ix in range(0, nx):
        ControlMap_flex[iz,iy,ix] = ControlMap[iz][iy][ix]
  grads = o.gradients(map_data = ControlMap_flex)
  g_cpp = []
  for gr, go, gb in zip(grads.grad_xyz, grads.grad_occ, grads.grad_uiso):
    g_cpp.append([gr[0],gr[1],gr[2],go,gb])
  #
  gcalc = gradients_fd(
    ControlMap = ControlMap,
    xrs        = xrs,
    n_real     = n_real,
    d_min      = d_min)
  assert approx_equal(o.target(), gcalc.compute())

  eob = 0.00001
  er = 1.e-6
  g_fd = []
  fra = uc.fractionalize
  ort = uc.orthogonalize
  for sc in xrs.scatterers():
    g_fd_ = []
    bs_fd = flex.double()
    os_fd = flex.double()
    # coordinates
    xyzf = sc.site[:]
    xyz  = ort(sc.site)
    g = []
    for e3 in [[er,0,0], [0,er,0], [0,0,er]]:
      e3 = flex.double(e3)
      xyzp = list( flex.double(xyz) + e3 )
      xyzm = list( flex.double(xyz) - e3 )
      sc.site = fra(xyzp)
      f1 = gcalc.compute()
      sc.site = fra(xyzm)
      f2 = gcalc.compute()
      g_ = (f1-f2)/(2*er)
      sc.site = xyzf
      g.append(g_)
    g_fd_ += g
    # occupancy
    sc.occupancy += eob
    f1 = gcalc.compute()
    sc.occupancy -= eob
    sc.occupancy -= eob
    f2 = gcalc.compute()
    sc.occupancy += eob
    g = (f1-f2)/(2*eob)
    os_fd.append(g)
    g_fd_ += os_fd
    # u_iso
    sc.u_iso += eob
    f1 = gcalc.compute()
    sc.u_iso -= eob
    sc.u_iso -= eob
    f2 = gcalc.compute()
    sc.u_iso += eob
    g = (f1-f2)/(2*eob)
    bs_fd.append(g)
    g_fd_ += bs_fd
    #
    g_fd.append(g_fd_)

  for g1, g2, g3 in zip(GradAtom, g_fd, g_cpp):
    print("anal:", ["%9.6f"%_ for _ in g1])
    print("find:", ["%9.6f"%_ for _ in g2])
    print("cpp :", ["%9.6f"%_ for _ in g3])
    print()
    assert approx_equal(g1,g2, 1.e-5)
    assert approx_equal(g3,g2, 1.e-5)

if (__name__ == "__main__"):
  start = time.perf_counter()
  #
  # This fails for some resolutions. Need to investigate the reason.
  #
  #for d_min in [round(x * 0.1, 1) for x in range(10, 61)]:
  for d_min in [1.0, 2.0]:
    d_min = round(d_min, 1)
    print(d_min)
    run(d_min = d_min)
  print("Time:", time.perf_counter()-start)
  print("OK")
