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
  voxel_obj = math.sphere_voxel(np,splat_range,uniform,adjust_dx,external_rmax,default_dx, fraction,xyz)

  assert voxel_obj.np() == np
  assert abs(voxel_obj.rmax() - 25)<1e-4
  assert voxel_obj.occupied_sites() == 162

  info = voxel_obj.status()
  expect_info = """number of grid point is:    68921
rmax is                : 25.00000000
max fraction one 1-d is: 0.90000000
non-empty grid point is:      162
non-empty grid fract is: 0.00235052
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

  #print list(Fnl.coefs())
  check_nl = [2.3370760050376607e-05, 2.1910087547227873e-07, 1.1145167882067956e-05, 1.1552795255443502e-05, 1.8332528001413832e-07, 2.7291537694166531e-07, 7.2474025554586763e-06, 2.3917327527551719e-05, 4.2203807328102694e-05]
  for i_pre, i_cal in zip( check_nl, Fnl.coefs()):
    assert abs(i_pre-i_cal) < eps

  check_nn = [1.1685380025188304e-05, 1.0955043773613935e-07, 1.6139115350383189e-05, 1.134898156875573e-05, 1.7412731546555787e-08, 2.2812032847790184e-07, 1.3014503682895901e-05, -4.7928189112959946e-06, 3.668426870555654e-05]
  for i_pre, i_cal in zip( check_nn, Fnn.coefs()):
    assert abs(i_pre-i_cal) < eps

  check_nnl=[(1.1685380025188304e-05+0j), 0j, 0j, (1.0955043773613935e-07+0j), (-1.6139115350383189e-05+0j), 0j, 0j, (5.5725839410339778e-06+0j), 0j, (5.7763976277217503e-06+0j), 0j, 0j, (-1.7412731546555787e-08+0j), 0j, 0j, 0j, 0j, (9.1662640007069148e-08+0j), 0j, (1.3645768847083266e-07+0j), (1.3014503682895901e-05+0j), 0j, 0j, (-8.9874088696083742e-06+0j), 0j, (1.3780227780904368e-05+0j), 0j, 0j, 0j, 0j, (3.6237012777293381e-06+0j), 0j, (1.1958663763775861e-05+0j), 0j, (2.1101903664051344e-05+0j)]
  for i_pre, i_cal in zip( check_nnl, Fnnl.coefs()):
    assert abs(i_pre-i_cal) < eps

  mom_0_0_0 = (0.00483433139642-0j)
  mom_1_1_0 = (0.000270247340704-0j)
  mom_1_1_1 = (-0.000191093727209+0.000191093727209j)
  mom_2_0_0 = (-0.00333843794042-0j)
  mom_2_2_1 = (-1.26396619825e-05+1.26396619825e-05j)
  mom_4_2_0 = (-0.00371195771764-0j)
  mom_4_4_0 = (0.00317650549416-0j)


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
