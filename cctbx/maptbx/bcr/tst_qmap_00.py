from __future__ import absolute_import, division, print_function
import math
import time
from scitbx.array_family import flex
import iotbx.map_manager
from libtbx.test_utils import approx_equal
import iotbx.pdb
from cctbx import adptbx
from cctbx.maptbx.bcr import qmap

import boost_adaptbx.boost.python as bp
ext = bp.import_ext("cctbx_maptbx_bcr_bcr_ext")

# N, Res=2.5
R = [   # mu, mu_max = atomic_resolution * (RadFact+RadMu)
  0.000000000000000,
  2.137901575982925,
  3.573558167112186,
  4.910453205604299,
  6.197723641551286,
  7.463266240825266,
  8.722743932451982,
  9.985192083454228,
 11.231510113425726,
 12.526176570560663]
B = [ #sigma
  72.054181932630286,
  35.262725859069029,
  20.965775916519476,
  14.404958530433785,
  11.340349195234730,
   9.713633940197088,
   8.449135265753416,
   8.350351017987917,
   8.441267229500561,
   8.246840051545909]
C = [ # kappa
   20.108947809148834,
  -19.705151510498212,
   11.106663690650764,
   -8.298529895167423,
    7.221275369283975,
   -6.725763928844316,
    6.354013560352876,
   -6.330524081979416,
    6.317387204802556,
   -6.358094541306267]

nterms = 5 # because mu_max = atomic_resolution * (RadFact+RadMu)
ScaleB = 1.0 / (8.0 * math.pi**2)
_X_mu     = R[:nterms]
_X_kappa  = C[:nterms]
_X_nu     = list(flex.double(B)  * ScaleB)[:nterms]
_X_musq   = list(flex.double(R)*flex.double(R))[:nterms]
_X_kappi  = list(flex.double(C)/(math.pi**1.5))[:nterms]

# three atoms
X_mu    = [_X_mu    , _X_mu    , _X_mu    ]
X_kappa = [_X_kappa , _X_kappa , _X_kappa ]
X_nu    = [_X_nu    , _X_nu    , _X_nu    ]
X_musq  = [_X_musq  , _X_musq  , _X_musq  ]
X_kappi = [_X_kappi , _X_kappi , _X_kappi ]




####################################################
#
#  MAIN PROGRAM

#==================== input information =====================

time_t01 = time.time()

RadFact = 2.0
RadAdd  = 0.5

FileMap   = 'Map_Orth_d25_B40.mrc'
FileOut   = 'NewMap1.mrc'

mm = iotbx.map_manager.map_manager(file_name = FileMap)
Ncrs = mm.map_data().all()
Scrs = [0,0,0]
Nxyz = mm.map_data().all()
unit_cell = mm.crystal_symmetry().unit_cell()

Mx,My,Mz = Nxyz
ControlMap = [[[ mm.map_data()[ix,iy,iz] for ix in range(Mx) ] for iy in range(My)] for iz in range(Mz)]

atoms = iotbx.pdb.input(file_name = "Atom2_orth.pdb").construct_hierarchy().atoms()
sites_cart = atoms.extract_xyz()
adp_as_u = atoms.extract_b()*adptbx.b_as_u(1.)
occupancy= atoms.extract_occ()

ModelValues = []
ModelTypes = []
bcr_scatterers = []
for s, u, o, e in zip(sites_cart, adp_as_u, occupancy, atoms.extract_element()):
  bcr_scatterer = ext.bcr_scatterer(
     site_cart = s,
     u_iso     = u,
     occ       = o,
     radius    = 2.5*RadFact, # atomic radius = atomic_resolution * RadFact
     resolution=2.5,
     mu        = _X_mu,
     kappa     = _X_kappa,
     nu        = _X_nu,
     musq      = _X_musq,
     kappi     = _X_kappi)
  bcr_scatterers.append(bcr_scatterer)

#-------------- calculation --------

OmegaMap = qmap.CalcOmegaMap(Ncrs,Scrs,Nxyz, unit_cell, bcr_scatterers)
print("OmegaMap", OmegaMap)
FuncMap  = qmap.CalcFuncMap(OmegaMap, ControlMap, Ncrs)
GradMap  = qmap.CalcGradMap(OmegaMap, ControlMap, Ncrs)
GradAtom = qmap.CalcGradAtom(GradMap, Ncrs,Scrs,Nxyz,unit_cell, bcr_scatterers)

#=============== output results ======================================

assert approx_equal(GradAtom[0], [0.6210971239199201, -1.2420777490312254, -2.484579744482787, 6.2528792787377645, -7.192527753707328])
assert approx_equal(GradAtom[1], (0,0,0,0,0), 1.e-4 )
assert approx_equal(GradAtom[2], (0,0,0,0,0), 1.e-4 )
assert approx_equal(FuncMap, 8.219926081295183)
#
nx,ny,nz = Nxyz
m0 = flex.double(flex.grid(Nxyz))
for iz in range(0, nz):
  for iy in range(0, ny):
    for ix in range(0, nx):
      v = OmegaMap[iz][iy][ix]
      m0[ix,iy,iz] = v

mm = iotbx.map_manager.map_manager(file_name = FileOut)
m1 = mm.map_data()
assert m1.all() == Nxyz

diff = flex.abs(m1-m0)
print("max:", flex.max(diff))
assert flex.max(diff) < 0.0045
for d in diff:
  assert d < 0.0045

#
mm2 = iotbx.map_manager.map_manager(
  map_data                   = m1,
  unit_cell_grid             = m1.all(),
  unit_cell_crystal_symmetry = mm.crystal_symmetry(),
  wrapping                   = True)

time_t08 = time.time()

print('')
print('time total                  ...... ',f'{(time_t08-time_t01)/60.:12.2f} min')
