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

# N, Res=2.
R = [   # mu, mu_max = atomic_resolution * (RadFact+RadMu)
  0.000000000000000,
  1.712360593307695,
  2.860580166643273,
  3.930302726197112,
  4.962251374908396,
  5.970647512379903,
  6.952758683715003,
  7.963220822471991,
  8.984906194981603,
  9.993120074344716,
 10.992451032022949,
 11.983697479405516,]
B = [ #sigma
  48.244553782807202,
  24.628243632189346,
  15.810470552464865,
  12.067116501902749,
  10.574205805298980,
  10.283495265305211,
   8.826162387624654,
   6.294058820812079,
   5.068141952638741,
   5.301406174543284,
   5.191481939931816,
   5.211623691024553]
C = [ # kappa
     19.271321843662683,
    -18.541138399974106,
     10.683508403513310,
     -8.236978149345655,
      7.479991534297030,
     -7.187545771994462,
      6.619893613856681,
     -5.708064789716665,
      5.210018761718937,
     -5.248008819482612,
      5.288269347413617,
     -5.302056075769035]

nterms = 5 # because mu_max = atomic_resolution * (RadFact+RadMu)
ScaleB = 1.0 / (8.0 * math.pi**2)
N_X_mu     = R[:nterms]
N_X_kappa  = C[:nterms]
N_X_nu     = list(flex.double(B)  * ScaleB)[:nterms]
N_X_musq   = list(flex.double(R)*flex.double(R))[:nterms]
N_X_kappi  = list(flex.double(C)/(math.pi**1.5))[:nterms]

# O, Res=3.
R = [   # mu, mu_max = atomic_resolution * (RadFact+RadMu)
   0.000000000000000,
   2.560903848918213,
   4.285732715914404,
   5.883984166277236,
   7.429375952487164,
   8.946299027954700,
  10.461160399797066,
  11.981469849074973]
B = [ #sigma
  100.849360060054153,
   49.035068535643454,
   29.701338976588403,
   20.750319103360798,
   16.347331326197232,
   14.096707421660756,
   11.234748486659131,
   12.486313907705881]
C = [ # kappa
    24.517074406860615,
   -25.156073018660404,
    14.748022516225689,
   -11.254350447134822,
     9.770456163765779,
    -9.079011754959023,
     8.269105203895405,
    -8.690867480475717]

nterms = 5 # because mu_max = atomic_resolution * (RadFact+RadMu)
ScaleB = 1.0 / (8.0 * math.pi**2)
O_X_mu     = R[:nterms]
O_X_kappa  = C[:nterms]
O_X_nu     = list(flex.double(B)  * ScaleB)[:nterms]
O_X_musq   = list(flex.double(R)*flex.double(R))[:nterms]
O_X_kappi  = list(flex.double(C)/(math.pi**1.5))[:nterms]


# C, Res=2.5
R = [   # mu, mu_max = atomic_resolution * (RadFact+RadMu)
   0.000000000000000,
   2.149314576062981,
   3.581569289873087,
   4.909977154934223,
   6.213154151570343,
   7.460190603722058,
   8.694650323766439,
   9.975498448979685,
  11.231418942888796,
  12.470010423176785]
B = [ #sigma
  72.307716144351417,
  35.251770397789826,
  22.221470328689239,
  15.692910656671938,
  13.502235006234359,
  13.350937276893909,
  10.486834228262911,
   6.693423840357337,
   8.391914576270764,
   8.280936413782010]
C = [ # kappa
     16.220333626106022,
    -15.450048791374753,
      8.946243578041553,
     -6.834266037834003,
      6.209643355061319,
     -5.984150779865546,
      5.432810395348175,
     -4.602147663635787,
      4.943594694422395,
     -4.967984493030292]

nterms = 5 # because mu_max = atomic_resolution * (RadFact+RadMu)
ScaleB = 1.0 / (8.0 * math.pi**2)
C_X_mu     = R[:nterms]
C_X_kappa  = C[:nterms]
C_X_nu     = list(flex.double(B)  * ScaleB)[:nterms]
C_X_musq   = list(flex.double(R)*flex.double(R))[:nterms]
C_X_kappi  = list(flex.double(C)/(math.pi**1.5))[:nterms]

# three atoms
X_mu    = [N_X_mu    , C_X_mu    , O_X_mu    ]
X_kappa = [N_X_kappa , C_X_kappa , O_X_kappa ]
X_nu    = [N_X_nu    , C_X_nu    , O_X_nu    ]
X_musq  = [N_X_musq  , C_X_musq  , O_X_musq  ]
X_kappi = [N_X_kappi , C_X_kappi , O_X_kappi ]




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

atoms = iotbx.pdb.input(file_name = "Atom3_orth.pdb").construct_hierarchy().atoms()
sites_cart = atoms.extract_xyz()
adp_as_u = atoms.extract_b()*adptbx.b_as_u(1.)
occupancy= atoms.extract_occ()
resolutions = [2, 2.5, 3.0]

bcr_scatterers = []
cntr=0
for s, u, o, e, r in zip(sites_cart, adp_as_u, occupancy, atoms.extract_element(),
                         resolutions):
  print(e, u)
  bcr_scatterer = ext.bcr_scatterer(
     site_cart = s,
     u_iso     = u,
     occ       = o,
     radius    = r*RadFact, # atomic radius = atomic_resolution * RadFact
     resolution=r,
     mu        = X_mu[cntr],
     kappa     = X_kappa[cntr],
     nu        = X_nu[cntr],
     musq      = X_musq[cntr],
     kappi     = X_kappi[cntr])
  bcr_scatterers.append(bcr_scatterer)
  cntr+=1

#-------------- calculation --------

OmegaMap = qmap.CalcOmegaMap(Ncrs,Scrs,Nxyz, unit_cell, bcr_scatterers)
FuncMap  = qmap.CalcFuncMap(OmegaMap, ControlMap, Ncrs)
GradMap  = qmap.CalcGradMap(OmegaMap, ControlMap, Ncrs)
GradAtom = qmap.CalcGradAtom(GradMap, Ncrs,Scrs,Nxyz,unit_cell, bcr_scatterers)

#=============== output results ======================================

assert approx_equal(GradAtom[0], [0.6122411344818136, -1.224368700760162, -2.4486010941878025, 6.355750038594599, -7.611067428662894])
assert approx_equal(GradAtom[1], (0.12988093751281296, -0.09073256875785853, -0.09230546853093913, -0.01424746771152637, -0.17724454937327658), 1.e-3 )
assert approx_equal(GradAtom[2], (-0.2157388505239968, -0.2281456322261232, -0.23658988439512613, -1.621959968943133, 0.7862790931118283))
assert approx_equal(FuncMap, 8.706029526196591, 1.e-4)
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
#assert flex.max(diff) < 0.0045
#for d in diff:
#  assert d < 0.0045
#
##
#mm2 = iotbx.map_manager.map_manager(
#  map_data                   = m1,
#  unit_cell_grid             = m1.all(),
#  unit_cell_crystal_symmetry = mm.crystal_symmetry(),
#  wrapping                   = True)

time_t08 = time.time()

print('')
print('time total                  ...... ',f'{(time_t08-time_t01)/60.:12.2f} min')
