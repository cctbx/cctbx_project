from scitbx.array_family import flex
from scitbx import math

def makexyz():
  xyz=flex.vec3_double()
  xyz.append( ( 25, 0, 0) )
  xyz.append( (-25, 0, 0) )
  xyz.append( ( 0, 10, 0) )
  xyz.append( ( 0,-10, 0) )
  xyz.append( ( 0,0, 15) )
  xyz.append( ( 0,0,-15) )

  return xyz

def tst_voxel(nmax, np):
  splat_range = 1
  fraction = 0.9
  default_dx = 0.7
  external_rmax = -1.0
  uniform = True
  adjust_dx = False
  centering = False

  xyz = makexyz()
  density=flex.double(xyz.size(),1.0)
  voxel_obj = math.sphere_voxel(np,splat_range,uniform,adjust_dx,external_rmax,default_dx, fraction,xyz,density)

  assert voxel_obj.np() == np
  assert abs(voxel_obj.rmax() - 25)<1e-4
  assert voxel_obj.occupied_sites() == 42

  info = voxel_obj.status()
  expect_info = """number of grid point is:    68921
rmax is                : 25.00000000
max fraction one 1-d is: 0.90000000
non-empty grid point is:       42
non-empty grid fract is: 0.00060939
"""
  assert info == expect_info

  new_xyz = voxel_obj.xyz()
  for xxx, yyy in zip(xyz, new_xyz):
    for x1, x2 in zip(xxx, yyy):
      assert x1 == x2

  return voxel_obj

def tst_grid(nmax,np):
  voxel_obj = tst_voxel( nmax, np )

  grid_obj = math.sphere_grid(np, nmax)
  pdb_out = False
  grid_obj.clean_space(voxel_obj, pdb_out)
  grid_obj.construct_space_sum()

  assert grid_obj.occupied_sites() == voxel_obj.occupied_sites()

  return grid_obj

def tst_moments(nmax,np):
  grid_obj = tst_grid(nmax,np)
  mom_obj = math.zernike_moments(grid_obj, nmax)

  moments = mom_obj.moments()
  Fnl= mom_obj.fnl()
  Fnn= mom_obj.fnn()
  Fnnl= mom_obj.fnnl()
  eps = 1e-8

#  print list(Fnl.coefs())
  check_nl =[1.5708741323298389e-06, 1.4726944990592269e-08, 7.729857198833008e-07, 7.7652533267041947e-07, 1.2381980826509915e-08, 1.8344106269055982e-08, 4.6642346948942935e-07, 1.5390436682985499e-06, 2.8370187789679361e-06]
  for i_pre, i_cal in zip( check_nl, Fnl.coefs()):
    assert abs(i_pre-i_cal) < eps

  check_nn = [7.8543706616491945e-07, 7.3634724952961344e-09, 1.1019361469817733e-06, 7.7475552627686014e-07, 1.4997075816312812e-09, 1.5363043547782947e-08, 8.5597462750509169e-07, -2.9659150565642845e-07, 2.4212429583779575e-06]
  for i_pre, i_cal in zip( check_nn, Fnn.coefs()):
    assert abs(i_pre-i_cal) < eps
  check_nnl=[(7.8543706616491945e-07+0j), 0j, 0j, (7.3634724952961344e-09+0j), (-1.1019361469817733e-06+0j), 0j, 0j, (3.864928599416504e-07+0j), 0j, (3.8826266633520974e-07+0j), 0j, 0j, (-1.4997075816312812e-09+0j), 0j, 0j, 0j, 0j, (6.1909904132549573e-09+0j), 0j, (9.1720531345279909e-09+0j), (8.5597462750509169e-07+0j), 0j, 0j, (-6.0044873331014139e-07+0j), 0j, (8.9704023896656984e-07+0j), 0j, 0j, 0j, 0j, (2.3321173474471467e-07+0j), 0j, (7.6952183414927495e-07+0j), 0j, (1.4185093894839679e-06+0j)]
  for i_pre, i_cal in zip( check_nnl, Fnnl.coefs()):
    assert abs(i_pre-i_cal) < eps

  mom_0_0_0 = (0.00125334517685-0j)
  mom_1_1_0 = (7.00641253676e-05-0j)
  mom_1_1_1 = (-4.95428181654e-05+4.95428181654e-05j)
  mom_2_0_0 = (-0.000879196064529-0j)
  mom_2_2_1 = (-3.27694940287e-06+3.27694940287e-06j)
  mom_4_2_0 = (-0.000953571758906-0j)
  mom_4_4_0 = (0.000823613460613-0j)

  assert abs( moments.get_coef(0,0,0) - mom_0_0_0 ) < eps
  assert abs( moments.get_coef(1,1,0) - mom_1_1_0 ) < eps
  assert abs( moments.get_coef(1,1,1) - mom_1_1_1 ) < eps
  assert abs( moments.get_coef(2,0,0) - mom_2_0_0 ) < eps
  assert abs( moments.get_coef(2,2,1) - mom_2_2_1 ) < eps
  assert abs( moments.get_coef(4,2,0) - mom_4_2_0 ) < eps
  assert abs( moments.get_coef(4,4,0) - mom_4_4_0 ) < eps



  # do some alignment please
  from scitbx.math import zernike_align_fft as zafft
  fixed = mom_obj.moments()
  moving = mom_obj.moments()
  al_obj = zafft.align( fixed, moving)
  assert abs(al_obj.get_cc()-1) < 1e-3




  return True



if __name__ == "__main__":
  nmax = 4
  np = 20
  tst_moments(nmax, np)
  print "OK"
