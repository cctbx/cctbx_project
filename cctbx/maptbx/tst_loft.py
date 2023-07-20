from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
import math, time
from cctbx import miller
from cctbx import crystal
import iotbx.pdb
from cctbx import maptbx
import scitbx.math

sin = math.sin
cos = math.cos
pi  = math.pi

sin_cos_table = scitbx.math.sin_cos_table(n=10000)

pdb_str = """
ATOM      1  O   HOH A   2    %8.3f%8.3f%8.3f  1.00  1.00           O
END
"""

def write_ccp4_map(cs, file_name, map_data):
  from iotbx import mrcfile
  mrcfile.write_ccp4_map(
      file_name=file_name,
      unit_cell=cs.unit_cell(),
      space_group=cs.space_group(),
      #gridding_first=(0,0,0),# This causes a bug (map gets shifted)
      #gridding_last=n_real,  # This causes a bug (map gets shifted)
      map_data=map_data,
      labels=flex.std_string([""]))

def rfactor(x,y):
  x = abs(x).data()
  y = abs(y).data()
  k = flex.sum(x*y)/flex.sum(y*y)
  diff = x-y*k
  summ = x+y*k
  return flex.sum(flex.abs(diff))/flex.sum(summ)*100.*2

def G(s, R):
  """
  G-function
  """
  arg = 2*pi*s*R
  vol = 4*pi*R**3/3
  return vol*3*(sin(arg) - arg*cos(arg))/(arg)**3

def make_map(N, x,y,z, a,b,c, R):
  """
  Make sperical mask of radius R on a grid N in (a,b,c) box with the center
  at (x,y,z).
  """
  sx,sy,sz = a/N[0], b/N[1], c/N[2]
  mask = flex.int(flex.grid(N),0)
  i=0
  while i<N[0]:
    j=0
    while j<N[1]:
      k=0
      while k<N[2]:
        d = math.sqrt( (x - i*sx)**2 + (y - j*sy)**2 + (z - k*sz)**2 )
        if d<=R:
          mask[i,j,k]=1
        k+=1
      j+=1
    i+=1
  return mask

def sf_analytical(miller_array, x,y,z, R):
  """
  Compute SF from spherical mask by exact formula
  """
  s = 1./miller_array.d_spacings().data()
  xf,yf,zf = miller_array.unit_cell().fractionalize([x,y,z])
  fm = flex.complex_double()
  for hkl, si in zip(miller_array.indices(), s):
    h,k,l = hkl
    arg = 2*pi*(h*xf+k*yf+l*zf)
    fm.append( G(si, R) * complex(cos(arg), sin(arg)) )
  return miller_array.customized_copy(data = fm)

def sf_direct(miller_array, mask):
  """
  Compute structure factors using mask_FT_05.docx
  """
  # Extract region (sphere)
  nx,ny,nz = mask.all()
  jxyz = flex.vec3_int()
  jx=0
  while jx<nx:
    jy=0
    while jy<ny:
      jz=0
      while jz<nz:
        mv = mask[jx,jy,jz]
        if mv == 1: jxyz.append([jx,jy,jz])
        jz+=1
      jy+=1
    jx+=1
  # Prepare/precompute arrays
  indices = miller_array.indices()
  size = indices.size()
  f  = flex.complex_double(size)
  hi = flex.double(size)
  ki = flex.double(size)
  li = flex.double(size)
  tponx,tpony,tponz  = 2*pi/nx, 2*pi/ny, 2*pi/nz
  for i_hkl, hkl in enumerate(indices):
    hi[i_hkl] = tponx*hkl[0]
    ki[i_hkl] = tpony*hkl[1]
    li[i_hkl] = tponz*hkl[2]
  # Loop over hkl indices and over region grid points
  for i_hkl in range(size):
    h,k,l = hi[i_hkl],ki[i_hkl],li[i_hkl]
    arg = 0
    for it in jxyz:
      jx,jy,jz = it
      arg_ = h*jx + k*jy + l*jz
      arg += complex(cos(arg_), sin(arg_))
    f[i_hkl] = arg
  return miller_array.customized_copy(data = f)

def run(d_min=3):
  #
  # Make cubic box and compute spherical map (mask) of radius R in its center
  #
  sc=3
  cs = crystal.symmetry((20.0*sc , 24.0*sc, 32.0*sc, 90,90,90), "P1")
  uc = cs.unit_cell()
  a,b,c = uc.parameters()[:3]
  x,y,z = 11.0, 13.0, 15.0
  #with open("one.pdb","w") as fo:
  #  fo.write(iotbx.pdb.format_cryst1_record(crystal_symmetry = cs))
  #  fo.write(pdb_str%(x,y,z))
  R = 3
  N = [int(a/0.6),int(b/0.6),int(c/0.6)]
  mask = make_map(N=N, x=x,y=y,z=z, a=a,b=b,c=c, R=R)
  #write_ccp4_map(cs=cs, file_name="mask.map", map_data=mask)
  #
  # Make complete set of indices hkl up to specified resolution
  #
  ms = miller.build_set(crystal_symmetry = cs, anomalous_flag = False, d_min=d_min)
  #
  # Compute structure factors using mask_FT_05.docx
  #
  #f3 = sf_direct(miller_array=f1, mask=mask)
  t0 = time.time()
  o3 = maptbx.loft(
    miller_indices = ms.indices(),
    map_data       = mask,
    abc            = flex.double([a,b,c]),
    d_min          = d_min)
  t3 = time.time()-t0
  print("time f3:", t3)
  ms = miller.set(crystal_symmetry=cs, indices=o3.indices())
  f3 = ms.array(data = o3.structure_factors())
  del ms
  #
  # Compute structure factors using FT from CCTBX
  #
  t0 = time.time()
  f1 = f3.structure_factors_from_map(map=mask, use_scale=True)
  t1 = time.time()-t0
  print("time f1:", t1)
  #
  # Compute structure factors using exact formula
  #
  t0 = time.time()
  f2 = sf_analytical(miller_array=f1, x=x,y=y,z=z, R=R)
  t2 = time.time()-t0
  print("time f2:", t2)
  #
  # Write mask structure factors as MTZ
  #
  mtz_dataset = f1.as_mtz_dataset(column_root_label = "f1")
  mtz_dataset.add_miller_array(
    miller_array      = f2,
    column_root_label = "f2")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "maps.mtz")
  #
  # Check all three Fmask structure factors are the same
  #
  cc12 = f1.map_correlation(other=f2)
  cc13 = f1.map_correlation(other=f3)
  cc23 = f2.map_correlation(other=f3)
  print()
  print("t1/t3: ", t1/t3)
  print()
  print("Map CC:", cc12)
  print("Map CC:", cc13)
  print("Map CC:", cc23)
  #
  #print("t2/t3: ", t2/t3)
  print("R-factor:", rfactor(f1,f2))
  print("R-factor:", rfactor(f1,f3))
  print("R-factor:", rfactor(f2,f3))
  #
  assert cc12 > 0.99
  assert cc13 > 0.99
  assert cc23 > 0.99

if (__name__ == "__main__"):
  run()
