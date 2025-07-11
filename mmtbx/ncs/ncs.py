"""Tools to hold and manipulate local (NCS) symmetry
 There can be any number of ncs groups in a ncs symmetry object.
   Each group has a set of NCS operators and centers and may apply
      to a part of the structure.
  for group in ncs.ncs_groups(): returns list of groups
  id=group.chain_and_residue_id() returns id of where it applies
  for center in group.centers(): returns list of centers of ncs regions in group
  for rota_matr in group.rota_matrices(): returns rota matrices
  for trans_orth in group.translations_orth(): returns translation matrices

 NOTE: symmetry operators map NCS position i on to NCS position 0 (they are
  inverses of the operators mapping position 0 on to i).

"""

from __future__ import absolute_import, division, print_function
# from mmtbx.ncs.ncs_utils import convert_phil_format
import sys, os, string
from operator import itemgetter
from libtbx.utils import Sorry
from libtbx.utils import null_out
import scitbx.rigid_body
from six.moves import zip
from six.moves import range
from copy import deepcopy

# Defaults for tolerances:
# Set 2017-12-23 to match values in find_ncs.py; these are very relaxed...
# previous values were tol_z=0.01, tol_r=.01,abs_tol_t=.10,rel_tol_t=0.001

default_tol_z=0.01
default_tol_r=0.02
default_abs_tol_t=2.0
default_rel_tol_t=0.05

def abs(aa):
  """Return absolute value of aa"""
  if aa>=0:return aa
  return -aa

def is_in_range(z,z_min,z_max):
    """Return True if z is >= z_min and <= z_max"""
    if z<z_min or z>z_max: return False
    return True

def remove_quotes_from_chain_id(chain_residue_id):
  """Remove the quotes from the chain names in group"""
  if chain_residue_id is None:
     return

  [group,list_of_resseq_list]=chain_residue_id
  new_group=[]
  for chain_id in group:
    new_group.append(remove_single_quotes(chain_id))
  return [new_group,list_of_resseq_list]

def remove_single_quotes(text):
  """Remove single quotes from ends of a string if they occur on both ends."""
  if not text: return text
  if not type(text)==type("abc"): return text
  if text.startswith("'") and text.endswith("'"):
    text=text[1:-1]
  return text


def is_identity(r,t,tol=1.e-2):
  """Return True if r, t is the identity"""
  identity_r=[1,0,0,0,1,0,0,0,1]
  identity_t=[0,0,0]
  for i in range(9):
    if abs(r[i]-identity_r[i])>tol: return False
  for i in range(3):
    if abs(t[i]-identity_t[i])>tol: return False
  return True

def is_same_transform(r1,t1,r2,t2,
   tol_r=default_tol_r,
   abs_tol_t=default_abs_tol_t,
   rel_tol_t=default_rel_tol_t):
    """Return True if r1, t1 is the same as r2, t2 within tolerance"""
    # require everything to be very similar
    for i in range(9):
      if abs(r1[i]-r2[i])>tol_r: return False
    for i in range(3):
      dd=abs(t1[i]-t2[i])
      dd2=0.5*(abs(t1[i])+abs(t2[i]))
      if dd>abs_tol_t and dd>rel_tol_t*dd2: # definitely does not match
        return False
    return True
def crystal_symmetry_to_ncs(crystal_symmetry=None):

  """Convert r,t fractional to r_orth,t_orth  orthogonal"""
  #  r x + t = x'
  #   x_orth=Ax   A = orthogonalization_matrix
  #  r_orth  (Ax) + t_orth = Ax'
  #  A-1 r_orth A  x + A-1 t_orth = x'
  #  r= (A-1 r_orth A)     t=A-1 t_orth
  #  r_orth = A r A-1     t_orth = A t

  from scitbx import matrix
  a=matrix.sqr(crystal_symmetry.unit_cell().orthogonalization_matrix())
  a_inv=a.inverse()

  trans_orth=[]
  ncs_rota_matr=[]
  ncs_center_orth=[]
  center=matrix.col((0.41,0.43,0.39)) # just a point near center of cell
  center_orth=crystal_symmetry.unit_cell().orthogonalize(center)

  for rt_mx in crystal_symmetry.space_group().all_ops():
    r=matrix.sqr(rt_mx.r().as_rational().as_float())
    t=matrix.col(rt_mx.t().as_rational().as_float())

    # Try to put each center for this ncs operator inside the cell
    c=matrix.col(rt_mx * center) # Note: forward (not inverse)
    coordinate_offset=matrix.col(offset_inside_cell(c,
       unit_cell=crystal_symmetry.unit_cell(),orthogonalize=False))
    c_offset=c+coordinate_offset

    r_orth=matrix.sqr(a * r * a_inv)
    t_orth= matrix.col(a * (t + coordinate_offset))
    c_value=crystal_symmetry.unit_cell().orthogonalize(c_offset)

    r_orth_inv=r_orth.inverse()
    t_orth_inv=-r_orth_inv*t_orth

    ncs_rota_matr.append(r_orth_inv)
    trans_orth.append(t_orth_inv)
    ncs_center_orth.append(c_value)

    chain_residue_id= [
        len(r_orth)*[None],
        len(r_orth)*[[]]
        ]

  ncs_obj=ncs()
  ncs_obj.import_ncs_group(
       ncs_rota_matr=ncs_rota_matr,
       center_orth=ncs_center_orth,
       trans_orth=trans_orth,
       chain_residue_id=chain_residue_id)

  return ncs_obj

def offset_inside_zero_one(x):
    """Place x inside [0,1] by adding integers"""
    if x >=0.0:
      return -1.0*int(x)  # 2.1 gives -2 to place inside (0,1)
    else:
      return 1.0-int(x)   # -2.1 gives + 3 to place inside (0,1)

def offset_inside_cell(center,unit_cell,orthogonalize=True):
    """Put the center inside (0,1)"""
    from scitbx.math import  matrix
    c=matrix.col(center)
    if orthogonalize:
      c_frac=unit_cell.fractionalize(c)
    else:
      c_frac=c
    offset_frac=[]
    for x in c_frac:
     offset_frac.append(offset_inside_zero_one(x))
    if orthogonalize:
      return unit_cell.orthogonalize(matrix.col(offset_frac))
    else:
      return matrix.col(offset_frac)

def get_ncs_from_text(text=None,text_is_ncs_spec=None,rotate_about_z=None,
    rotate_about_y=None,rotate_about_new_y=None,ncs_name=None,out=sys.stdout):
  """Read a text file containing ncs information and get an ncs object.
   Allow rotation about y, z, and new_y"""
  from mmtbx.ncs.ncs import ncs
  import iotbx.pdb
  ncs_object=ncs()
  if text_is_ncs_spec:
    ncs_object.read_ncs(lines=text.splitlines())
  else: # read BIOMTR
    from cctbx.array_family import flex
    pdb_inp=iotbx.pdb.input(lines=flex.split_lines(text),source_info='string')
    ncs_object.ncs_from_pdb_input_BIOMT(pdb_inp=pdb_inp,log=out)

  if rotate_about_new_y:
    ncs_object.rotate_about_y(rot_deg=rotate_about_new_y,invert_matrices=True)
  if rotate_about_z:
    ncs_object.rotate_about_z(rot_deg=rotate_about_z,invert_matrices=True)
  if rotate_about_y:
    ncs_object.rotate_about_y(rot_deg=rotate_about_y,invert_matrices=True)
  if ncs_name:
    ncs_object.set_ncs_name(ncs_name)
  return ncs_object


def get_helical_symmetry(helical_rot_deg=None,
     helical_trans_z_angstrom=None,max_ops=None):
  """Get helical symmetry object from rotation and translation along z"""

  from scitbx import matrix
  rot=get_rot_z(rot_deg=helical_rot_deg)
  rot_inv=rot.inverse()
  trans_along_z=matrix.col((0,0,helical_trans_z_angstrom))

  n=max(2,int(360/helical_rot_deg))
  if max_ops and n> max(1,max_ops//2): n= max(1,max_ops//2)
  rots=[]
  trans=[]
  rots.append(matrix.sqr((1,0,0,0,1,0,0,0,1),))
  trans.append(matrix.col((0,0,0,)))

  for i in range(n):
    rots=[rot*rots[0]]+rots[:]
    trans=[trans[0]-trans_along_z]+trans[:]

  for i in range(n):
    rots.append(rot_inv*rots[-1])
    trans.append(trans[-1]+trans_along_z)

  from mmtbx.ncs.ncs import ncs
  ncs_object=ncs()
  ncs_name="Helical %5.2f deg  %6.2f Z-trans " %(
        helical_rot_deg,helical_trans_z_angstrom)
  ncs_object.ncs_from_import(rot_list=rots,trans_list=trans)
  ncs_object.set_ncs_name(ncs_name)


  return ncs_object

def get_d_symmetry(n=None,two_fold_along_x=True,ncs_name=None):
  """Get ncs object for D symmetry"""
  return get_c_symmetry(n=n,is_d=True,two_fold_along_x=two_fold_along_x,
     ncs_name=ncs_name)

def get_rot_z(rot_deg=None):
  """Get rotation about z for rot_deg"""
  import math
  theta=rot_deg*3.14159/180.
  cc=math.cos(theta)
  ss=math.sin(theta)
  from scitbx import matrix
  return matrix.sqr((cc,ss,0,-ss,cc,0,0,0,1,))

def get_rot_y(rot_deg=None):
  """Get rotation about y for rot_deg"""
  import math
  theta=rot_deg*3.14159/180.
  cc=math.cos(theta)
  ss=math.sin(theta)
  from scitbx import matrix
  return matrix.sqr((cc,0,ss,
                     0,1,0,
                    -ss,0,cc,
                     ))

def get_c_symmetry(n=None,is_d=False,two_fold_along_x=None,ncs_name=None):
  """Generate n-fold C symmetry"""
  oper=get_rot_z(rot_deg=360./n)
  oper_inv=oper.inverse()
  rots=[]
  trans=[]
  from scitbx import matrix
  rots.append(matrix.sqr((1,0,0,0,1,0,0,0,1),))
  trans.append(matrix.col((0,0,0,)))

  for i in range(n-1):
    rots.append(oper_inv*rots[-1])
    trans.append(matrix.col((0,0,0,)))

  if is_d:
    if two_fold_along_x:
      d_oper=matrix.sqr((1,0,0,0,-1,0,0,0,-1,))
    else: # along y
      d_oper=matrix.sqr((-1,0,0,0,1,0,0,0,-1,))
    new_rots=[]
    new_trans=[]
    for r,t in zip(rots,trans):
      new_rots.append(d_oper*r)
      new_trans.append(t)
    rots+=new_rots
    trans+=new_trans

  from mmtbx.ncs.ncs import ncs

  ncs_object=ncs()
  ncs_object.ncs_from_import(rot_list=rots,trans_list=trans)
  if ncs_name:
    ncs_object.set_ncs_name(ncs_name)

  return ncs_object

def remove_extra(text):
  """Keep text that is not alphabetical lower case or parentheses"""
  new_text=""
  for t in text:
    if not t.lower() in "()abcdefghijklmnopqrstuvwxyz":
      new_text+=t
  return new_text

def value(str):
  """Get the integer value of the numbers in str"""
  try:
    return int(remove_extra(str[1:]))
  except Exception as e:
    return 0


def generate_ncs_ops(symmetry=None,
   must_be_consistent_with_space_group_number = None,
   helical_rot_deg=None,
   helical_trans_z_angstrom=None,
   op_max=None,
   two_fold_along_x=None,
   include_helical_symmetry=None,
   max_helical_ops_to_check=None,
   require_helical_or_point_group_symmetry=None,
   out=sys.stdout):
  """Generate ncs objects corresponding to common point-group symmetries in
  conventional orientations such as C2 D7 etc."""

  ncs_list=[]
  all=False
  sym_type=None
  sym_n=None
  if symmetry.lower() in ['all','any']:
    all=True
  elif symmetry.lower() in ["i"]:
    sym_type='I'
  elif symmetry.lower() in ["o"]:
    sym_type='O'
  elif symmetry.lower() in ["t"]:
    sym_type='T'
  elif symmetry.lower().startswith("d"):
    sym_type='D'
    if len(symmetry)>1:
      sym_n=value(symmetry)
  elif symmetry.lower().startswith("c"):
    sym_type='C'
    if len(symmetry)>1:
      sym_n=value(symmetry)
  elif symmetry.lower() in ['helical','helix']:
    sym_type='helical'

  print("Sym type: %s  Sym N: %s" %(
     sym_type,sym_n), file=out)
  if sym_type is None and not all:
    from libtbx.utils import Sorry
    raise Sorry("Unknown symmetry type: '%s' " %(symmetry))

  if sym_n:
    i_start=sym_n
    i_end=sym_n
  else:
    i_start=2
    if op_max is None:
      i_end=14
    else:
      i_end=op_max

  # NOTE: the (a) (b) etc designations are arbitrary

  from mmtbx.ncs.ncs import get_ncs_from_text, \
      get_c_symmetry,get_d_symmetry,get_helical_symmetry
  if sym_type=='I' or all:
    if two_fold_along_x is None or two_fold_along_x==False:
      ncs_list.append(get_ncs_from_text(text=icosahedral_b,
        text_is_ncs_spec=True,ncs_name='I (b)'))
      ncs_list.append(get_ncs_from_text(text=icosahedral_d,
        text_is_ncs_spec=True,ncs_name='I (d)'))
      ncs_list.append(get_ncs_from_text(text=icosahedral_f,
         text_is_ncs_spec=True,
         ncs_name='I (f)'))
    if two_fold_along_x is None or two_fold_along_x==True:
      ncs_list.append(get_ncs_from_text(text=icosahedral_b,
          rotate_about_z=90,
          text_is_ncs_spec=True,ncs_name='I (a)'))
      ncs_list.append(get_ncs_from_text(text=icosahedral_d,
          text_is_ncs_spec=True,rotate_about_z=90,
           ncs_name='I (c)'))
      ncs_list.append(get_ncs_from_text(text=icosahedral_f,
          rotate_about_z=90,text_is_ncs_spec=True,
          ncs_name='I (e)'))
  if sym_type=='O' or all:
    if two_fold_along_x is None or two_fold_along_x==True:
      ncs_list.append(get_ncs_from_text(text=octahedral_a,
        text_is_ncs_spec=True,ncs_name='O (a)'))
    if two_fold_along_x is None or two_fold_along_x==False:
      ncs_list.append(get_ncs_from_text(text=octahedral_a,
        rotate_about_z=45,
        text_is_ncs_spec=True,ncs_name='O (b)'))
  if sym_type=='T' or all:
    if two_fold_along_x is None or two_fold_along_x==False:
      ncs_list.append(get_ncs_from_text(text=tetrahedral_a,
        text_is_ncs_spec=True,ncs_name='T (a)'))
    if two_fold_along_x is None or two_fold_along_x==True:
      ncs_list.append(get_ncs_from_text(text=tetrahedral_a,
        # rotate_about_y=54.73563863,
        # rotate_about_z=-45,
        rotate_about_y=45,
        text_is_ncs_spec=True,ncs_name='T (b)'))
    ncs_list.append(get_ncs_from_text(text=tetrahedral_b,
           text_is_ncs_spec=True,ncs_name='T (c)'))

  if sym_type=='C' or all:
    for i in range(i_start,i_end+1):
      ncs_list.append(get_c_symmetry(n=i,ncs_name='C%d ' %(i)))
  if sym_type=='D' or all:
    for i in range(i_start,i_end+1):
      if two_fold_along_x is None or two_fold_along_x==True:
        ncs_list.append(get_d_symmetry(n=i,two_fold_along_x=True,
        ncs_name='D%d (a)' %(i)))
      if two_fold_along_x is None or two_fold_along_x==False:
        ncs_list.append(get_d_symmetry(n=i,two_fold_along_x=False,
        ncs_name='D%d (b)' %(i)))

  #  Pare down ncs_list if must_be_consistent_with_space_group_number is set
  if must_be_consistent_with_space_group_number:
    ncs_list=remove_ncs_not_consistent_with_space_group_number(
       must_be_consistent_with_space_group_number,ncs_list)

  if sym_type=='helical' or (
      all and include_helical_symmetry):
     if helical_rot_deg is not None and helical_trans_z_angstrom is not None:
      ncs_list.append(get_helical_symmetry(
       helical_rot_deg=helical_rot_deg,
       helical_trans_z_angstrom=helical_trans_z_angstrom,
       max_ops=max_helical_ops_to_check))

  return ncs_list

def remove_ncs_not_consistent_with_space_group_number(sg_number,ncs_list):
  ''' remove ncs ojbect that are not consistent with space_group number
     sg_number
   NOTE: allows ncs that has higher symmetry than sg_number
  '''

  from cctbx import crystal
  from cctbx import sgtbx
  from cctbx import uctbx
  sg = sgtbx.space_group_info(sg_number)
  uc = uctbx.unit_cell((10.,10.,10.,90.,90.,90.))
  cs = crystal.symmetry(unit_cell=uc ,space_group_info = sg)
  ncs_from_cs = crystal_symmetry_to_ncs(cs)
  # Make sure all ops of ncs_from_cs are in each ncs used
  new_ncs_list =[]
  for ncs_obj in ncs_list:
    if ncs_from_cs.is_similar_ncs_object(ncs_obj,abs_tol_t=10000):
      new_ncs_list.append(ncs_obj)
  ncs_list = new_ncs_list
  return ncs_list

# Symmetry for icosahedron in 3 settings


icosahedral_f=\
"""

Summary of NCS information
Mon Jul  2 11:08:24 2018
/Users/terwill/Downloads/view/ncs_text

source_info icosahedral_f.dat




new_ncs_group
new_operator

rota_matrix    1.0000    0.0000    0.0000
rota_matrix    0.0000    1.0000    0.0000
rota_matrix    0.0000    0.0000    1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth  -14.5127   20.8478  102.7355
new_operator

rota_matrix    0.5000    0.8090   -0.3090
rota_matrix   -0.8090    0.3090   -0.5000
rota_matrix   -0.3090    0.5000    0.8090
tran_orth    -0.0000    0.0000   -0.0000

center_orth  -55.8696   46.0690   77.1756
new_operator

rota_matrix   -0.3090    0.5000   -0.8090
rota_matrix   -0.5000   -0.8090   -0.3090
rota_matrix   -0.8090    0.3090    0.5000
tran_orth     0.0000   -0.0000    0.0000

center_orth  -89.0540    7.6245   56.6665
new_operator

rota_matrix   -0.3090   -0.5000   -0.8090
rota_matrix    0.5000   -0.8090    0.3090
rota_matrix   -0.8090   -0.3090    0.5000
tran_orth    -0.0000    0.0000   -0.0000

center_orth  -68.2062  -41.3569   69.5511
new_operator

rota_matrix    0.5000   -0.8090   -0.3090
rota_matrix    0.8090    0.3090    0.5000
rota_matrix   -0.3090   -0.5000    0.8090
tran_orth    -0.0000   -0.0000   -0.0000

center_orth  -22.1372  -33.1844   98.0233
new_operator

rota_matrix   -1.0000    0.0000    0.0000
rota_matrix    0.0000   -1.0000    0.0000
rota_matrix    0.0000    0.0000    1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth   14.5127  -20.8478  102.7355
new_operator

rota_matrix   -0.5000   -0.8090   -0.3090
rota_matrix    0.8090   -0.3090   -0.5000
rota_matrix    0.3090   -0.5000    0.8090
tran_orth    -0.0000    0.0000   -0.0000

center_orth   55.8696  -46.0690   77.1756
new_operator

rota_matrix    0.3090   -0.5000   -0.8090
rota_matrix    0.5000    0.8090   -0.3090
rota_matrix    0.8090   -0.3090    0.5000
tran_orth     0.0000   -0.0000    0.0000

center_orth   89.0540   -7.6245   56.6665
new_operator

rota_matrix    0.3090    0.5000   -0.8090
rota_matrix   -0.5000    0.8090    0.3090
rota_matrix    0.8090    0.3090    0.5000
tran_orth    -0.0000    0.0000   -0.0000

center_orth   68.2062   41.3569   69.5511
new_operator

rota_matrix   -0.5000    0.8090   -0.3090
rota_matrix   -0.8090   -0.3090    0.5000
rota_matrix    0.3090    0.5000    0.8090
tran_orth    -0.0000   -0.0000   -0.0000

center_orth   22.1372   33.1844   98.0233
new_operator

rota_matrix    1.0000    0.0000    0.0000
rota_matrix    0.0000   -1.0000    0.0000
rota_matrix    0.0000    0.0000   -1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth  -14.5127  -20.8478  -102.7355
new_operator

rota_matrix    0.5000   -0.8090    0.3090
rota_matrix   -0.8090   -0.3090    0.5000
rota_matrix   -0.3090   -0.5000   -0.8090
tran_orth    -0.0000    0.0000   -0.0000

center_orth  -55.8696  -46.0690  -77.1756
new_operator

rota_matrix   -0.3090   -0.5000    0.8090
rota_matrix   -0.5000    0.8090    0.3090
rota_matrix   -0.8090   -0.3090   -0.5000
tran_orth     0.0000   -0.0000    0.0000

center_orth  -89.0540   -7.6245  -56.6665
new_operator

rota_matrix   -0.3090    0.5000    0.8090
rota_matrix    0.5000    0.8090   -0.3090
rota_matrix   -0.8090    0.3090   -0.5000
tran_orth    -0.0000    0.0000   -0.0000

center_orth  -68.2062   41.3569  -69.5511
new_operator

rota_matrix    0.5000    0.8090    0.3090
rota_matrix    0.8090   -0.3090   -0.5000
rota_matrix   -0.3090    0.5000   -0.8090
tran_orth    -0.0000   -0.0000   -0.0000

center_orth  -22.1372   33.1844  -98.0233
new_operator

rota_matrix   -1.0000    0.0000    0.0000
rota_matrix    0.0000    1.0000    0.0000
rota_matrix    0.0000    0.0000   -1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth   14.5127   20.8478  -102.7355
new_operator

rota_matrix   -0.5000    0.8090    0.3090
rota_matrix    0.8090    0.3090    0.5000
rota_matrix    0.3090    0.5000   -0.8090
tran_orth    -0.0000    0.0000   -0.0000

center_orth   55.8696   46.0690  -77.1756
new_operator

rota_matrix    0.3090    0.5000    0.8090
rota_matrix    0.5000   -0.8090    0.3090
rota_matrix    0.8090    0.3090   -0.5000
tran_orth     0.0000   -0.0000    0.0000

center_orth   89.0540    7.6245  -56.6665
new_operator

rota_matrix    0.3090   -0.5000    0.8090
rota_matrix   -0.5000   -0.8090   -0.3090
rota_matrix    0.8090   -0.3090   -0.5000
tran_orth    -0.0000    0.0000   -0.0000

center_orth   68.2062  -41.3569  -69.5511
new_operator

rota_matrix   -0.5000   -0.8090    0.3090
rota_matrix   -0.8090    0.3090   -0.5000
rota_matrix    0.3090   -0.5000   -0.8090
tran_orth    -0.0000   -0.0000   -0.0000

center_orth   22.1372  -33.1844  -98.0233
new_operator

rota_matrix   -0.0000   -0.0000   -1.0000
rota_matrix   -1.0000    0.0000   -0.0000
rota_matrix    0.0000    1.0000   -0.0000
tran_orth     0.0000   -0.0000    0.0000

center_orth  -20.8478  102.7355   14.5127
new_operator

rota_matrix   -0.8090   -0.3090   -0.5000
rota_matrix   -0.3090   -0.5000    0.8090
rota_matrix   -0.5000    0.8090    0.3090
tran_orth    -0.0000    0.0000   -0.0000

center_orth  -46.0690   77.1756   55.8696
new_operator

rota_matrix   -0.5000   -0.8090    0.3090
rota_matrix    0.8090   -0.3090    0.5000
rota_matrix   -0.3090    0.5000    0.8090
tran_orth     0.0000   -0.0000    0.0000

center_orth   -7.6245   56.6665   89.0540
new_operator

rota_matrix    0.5000   -0.8090    0.3090
rota_matrix    0.8090    0.3090   -0.5000
rota_matrix    0.3090    0.5000    0.8090
tran_orth    -0.0000    0.0000   -0.0000

center_orth   41.3569   69.5511   68.2062
new_operator

rota_matrix    0.8090   -0.3090   -0.5000
rota_matrix   -0.3090    0.5000   -0.8090
rota_matrix    0.5000    0.8090    0.3090
tran_orth    -0.0000   -0.0000   -0.0000

center_orth   33.1844   98.0233   22.1372
new_operator

rota_matrix    0.0000    0.0000    1.0000
rota_matrix    1.0000    0.0000    0.0000
rota_matrix   -0.0000    1.0000    0.0000
tran_orth     0.0000   -0.0000    0.0000

center_orth   20.8478  102.7355  -14.5127
new_operator

rota_matrix    0.8090   -0.3090    0.5000
rota_matrix    0.3090   -0.5000   -0.8090
rota_matrix    0.5000    0.8090   -0.3090
tran_orth    -0.0000    0.0000   -0.0000

center_orth   46.0690   77.1756  -55.8696
new_operator

rota_matrix    0.5000   -0.8090   -0.3090
rota_matrix   -0.8090   -0.3090   -0.5000
rota_matrix    0.3090    0.5000   -0.8090
tran_orth     0.0000   -0.0000    0.0000

center_orth    7.6245   56.6665  -89.0540
new_operator

rota_matrix   -0.5000   -0.8090   -0.3090
rota_matrix   -0.8090    0.3090    0.5000
rota_matrix   -0.3090    0.5000   -0.8090
tran_orth    -0.0000    0.0000   -0.0000

center_orth  -41.3569   69.5511  -68.2062
new_operator

rota_matrix   -0.8090   -0.3090    0.5000
rota_matrix    0.3090    0.5000    0.8090
rota_matrix   -0.5000    0.8090   -0.3090
tran_orth    -0.0000   -0.0000   -0.0000

center_orth  -33.1844   98.0233  -22.1372
new_operator

rota_matrix   -0.0000   -0.0000   -1.0000
rota_matrix    1.0000   -0.0000   -0.0000
rota_matrix   -0.0000   -1.0000    0.0000
tran_orth    -0.0000   -0.0000   -0.0000

center_orth   20.8478  -102.7355   14.5127
new_operator

rota_matrix    0.8090    0.3090   -0.5000
rota_matrix    0.3090    0.5000    0.8090
rota_matrix    0.5000   -0.8090    0.3090
tran_orth    -0.0000    0.0000   -0.0000

center_orth   46.0690  -77.1756   55.8696
new_operator

rota_matrix    0.5000    0.8090    0.3090
rota_matrix   -0.8090    0.3090    0.5000
rota_matrix    0.3090   -0.5000    0.8090
tran_orth     0.0000   -0.0000    0.0000

center_orth    7.6245  -56.6665   89.0540
new_operator

rota_matrix   -0.5000    0.8090    0.3090
rota_matrix   -0.8090   -0.3090   -0.5000
rota_matrix   -0.3090   -0.5000    0.8090
tran_orth    -0.0000    0.0000   -0.0000

center_orth  -41.3569  -69.5511   68.2062
new_operator

rota_matrix   -0.8090    0.3090   -0.5000
rota_matrix    0.3090   -0.5000   -0.8090
rota_matrix   -0.5000   -0.8090    0.3090
tran_orth    -0.0000   -0.0000   -0.0000

center_orth  -33.1844  -98.0233   22.1372
new_operator

rota_matrix   -0.0000    0.0000    1.0000
rota_matrix   -1.0000   -0.0000   -0.0000
rota_matrix    0.0000   -1.0000   -0.0000
tran_orth    -0.0000   -0.0000    0.0000

center_orth  -20.8478  -102.7355  -14.5127
new_operator

rota_matrix   -0.8090    0.3090    0.5000
rota_matrix   -0.3090    0.5000   -0.8090
rota_matrix   -0.5000   -0.8090   -0.3090
tran_orth    -0.0000    0.0000   -0.0000

center_orth  -46.0690  -77.1756  -55.8696
new_operator

rota_matrix   -0.5000    0.8090   -0.3090
rota_matrix    0.8090    0.3090   -0.5000
rota_matrix   -0.3090   -0.5000   -0.8090
tran_orth     0.0000   -0.0000    0.0000

center_orth   -7.6245  -56.6665  -89.0540
new_operator

rota_matrix    0.5000    0.8090   -0.3090
rota_matrix    0.8090   -0.3090    0.5000
rota_matrix    0.3090   -0.5000   -0.8090
tran_orth    -0.0000    0.0000   -0.0000

center_orth   41.3569  -69.5511  -68.2062
new_operator

rota_matrix    0.8090    0.3090    0.5000
rota_matrix   -0.3090   -0.5000    0.8090
rota_matrix    0.5000   -0.8090   -0.3090
tran_orth    -0.0000   -0.0000   -0.0000

center_orth   33.1844  -98.0233  -22.1372
new_operator

rota_matrix   -0.0000   -1.0000    0.0000
rota_matrix    0.0000    0.0000    1.0000
rota_matrix   -1.0000    0.0000    0.0000
tran_orth    -0.0000    0.0000   -0.0000

center_orth  -102.7355   14.5127   20.8478
new_operator

rota_matrix    0.3090   -0.5000    0.8090
rota_matrix    0.5000    0.8090    0.3090
rota_matrix   -0.8090    0.3090    0.5000
tran_orth    -0.0000    0.0000   -0.0000

center_orth  -77.1756   55.8696   46.0690
new_operator

rota_matrix    0.8090    0.3090    0.5000
rota_matrix    0.3090    0.5000   -0.8090
rota_matrix   -0.5000    0.8090    0.3090
tran_orth     0.0000   -0.0000    0.0000

center_orth  -56.6665   89.0540    7.6245
new_operator

rota_matrix    0.8090    0.3090   -0.5000
rota_matrix   -0.3090   -0.5000   -0.8090
rota_matrix   -0.5000    0.8090   -0.3090
tran_orth    -0.0000    0.0000   -0.0000

center_orth  -69.5511   68.2062  -41.3569
new_operator

rota_matrix    0.3090   -0.5000   -0.8090
rota_matrix   -0.5000   -0.8090    0.3090
rota_matrix   -0.8090    0.3090   -0.5000
tran_orth    -0.0000   -0.0000   -0.0000

center_orth  -98.0233   22.1372  -33.1844
new_operator

rota_matrix   -0.0000    1.0000    0.0000
rota_matrix    0.0000    0.0000   -1.0000
rota_matrix   -1.0000   -0.0000   -0.0000
tran_orth    -0.0000    0.0000   -0.0000

center_orth  -102.7355  -14.5127  -20.8478
new_operator

rota_matrix    0.3090    0.5000   -0.8090
rota_matrix    0.5000   -0.8090   -0.3090
rota_matrix   -0.8090   -0.3090   -0.5000
tran_orth    -0.0000    0.0000   -0.0000

center_orth  -77.1756  -55.8696  -46.0690
new_operator

rota_matrix    0.8090   -0.3090   -0.5000
rota_matrix    0.3090   -0.5000    0.8090
rota_matrix   -0.5000   -0.8090   -0.3090
tran_orth     0.0000   -0.0000    0.0000

center_orth  -56.6665  -89.0540   -7.6245
new_operator

rota_matrix    0.8090   -0.3090    0.5000
rota_matrix   -0.3090    0.5000    0.8090
rota_matrix   -0.5000   -0.8090    0.3090
tran_orth    -0.0000    0.0000   -0.0000

center_orth  -69.5511  -68.2062   41.3569
new_operator

rota_matrix    0.3090    0.5000    0.8090
rota_matrix   -0.5000    0.8090   -0.3090
rota_matrix   -0.8090   -0.3090    0.5000
tran_orth    -0.0000   -0.0000   -0.0000

center_orth  -98.0233  -22.1372   33.1844
new_operator

rota_matrix    0.0000   -1.0000    0.0000
rota_matrix   -0.0000   -0.0000   -1.0000
rota_matrix    1.0000    0.0000   -0.0000
tran_orth    -0.0000    0.0000   -0.0000

center_orth  102.7355   14.5127  -20.8478
new_operator

rota_matrix   -0.3090   -0.5000   -0.8090
rota_matrix   -0.5000    0.8090   -0.3090
rota_matrix    0.8090    0.3090   -0.5000
tran_orth    -0.0000    0.0000   -0.0000

center_orth   77.1756   55.8696  -46.0690
new_operator

rota_matrix   -0.8090    0.3090   -0.5000
rota_matrix   -0.3090    0.5000    0.8090
rota_matrix    0.5000    0.8090   -0.3090
tran_orth     0.0000   -0.0000    0.0000

center_orth   56.6665   89.0540   -7.6245
new_operator

rota_matrix   -0.8090    0.3090    0.5000
rota_matrix    0.3090   -0.5000    0.8090
rota_matrix    0.5000    0.8090    0.3090
tran_orth    -0.0000    0.0000   -0.0000

center_orth   69.5511   68.2062   41.3569
new_operator

rota_matrix   -0.3090   -0.5000    0.8090
rota_matrix    0.5000   -0.8090   -0.3090
rota_matrix    0.8090    0.3090    0.5000
tran_orth    -0.0000   -0.0000   -0.0000

center_orth   98.0233   22.1372   33.1844
new_operator

rota_matrix   -0.0000    1.0000    0.0000
rota_matrix   -0.0000   -0.0000    1.0000
rota_matrix    1.0000    0.0000    0.0000
tran_orth     0.0000    0.0000    0.0000

center_orth  102.7355  -14.5127   20.8478
new_operator

rota_matrix   -0.3090    0.5000    0.8090
rota_matrix   -0.5000   -0.8090    0.3090
rota_matrix    0.8090   -0.3090    0.5000
tran_orth    -0.0000    0.0000   -0.0000

center_orth   77.1756  -55.8696   46.0690
new_operator

rota_matrix   -0.8090   -0.3090    0.5000
rota_matrix   -0.3090   -0.5000   -0.8090
rota_matrix    0.5000   -0.8090    0.3090
tran_orth     0.0000   -0.0000    0.0000

center_orth   56.6665  -89.0540    7.6245
new_operator

rota_matrix   -0.8090   -0.3090   -0.5000
rota_matrix    0.3090    0.5000   -0.8090
rota_matrix    0.5000   -0.8090   -0.3090
tran_orth    -0.0000    0.0000   -0.0000

center_orth   69.5511  -68.2062  -41.3569
new_operator

rota_matrix   -0.3090    0.5000   -0.8090
rota_matrix    0.5000    0.8090    0.3090
rota_matrix    0.8090   -0.3090   -0.5000
tran_orth    -0.0000   -0.0000   -0.0000

center_orth   98.0233  -22.1372  -33.1844


"""
icosahedral_d=\
"""

Summary of NCS information
Mon Jul  2 11:08:15 2018
/Users/terwill/Downloads/view/ncs_text

source_info icosahedral_d.dat




new_ncs_group
new_operator

rota_matrix    1.0000    0.0000    0.0000
rota_matrix    0.0000    1.0000    0.0000
rota_matrix    0.0000    0.0000    1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.6708    0.1625    0.7236
rota_matrix   -0.6882    0.5000    0.5257
rota_matrix   -0.2764   -0.8507    0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.1382   -0.4253    0.8944
rota_matrix   -0.9511   -0.3090    0.0000
rota_matrix    0.2764   -0.8507   -0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.1382   -0.9511    0.2764
rota_matrix   -0.4253   -0.3090   -0.8507
rota_matrix    0.8944    0.0000   -0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.6708   -0.6882   -0.2764
rota_matrix    0.1625    0.5000   -0.8507
rota_matrix    0.7236    0.5257    0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.4472    0.0000    0.8944
rota_matrix    0.0000   -1.0000    0.0000
rota_matrix    0.8944    0.0000   -0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.9472   -0.1625    0.2764
rota_matrix    0.1625   -0.5000   -0.8507
rota_matrix    0.2764    0.8507   -0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.8618    0.4253   -0.2764
rota_matrix   -0.4253    0.3090   -0.8507
rota_matrix   -0.2764    0.8507    0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.3090    0.9511    0.0000
rota_matrix   -0.9511    0.3090    0.0000
rota_matrix    0.0000   -0.0000    1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.0528    0.6882    0.7236
rota_matrix   -0.6882   -0.5000    0.5257
rota_matrix    0.7236   -0.5257    0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -1.0000   -0.0000   -0.0000
rota_matrix   -0.0000    1.0000    0.0000
rota_matrix   -0.0000    0.0000   -1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.6708    0.1625   -0.7236
rota_matrix    0.6882    0.5000   -0.5257
rota_matrix    0.2764   -0.8507   -0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.1382   -0.4253   -0.8944
rota_matrix    0.9511   -0.3090   -0.0000
rota_matrix   -0.2764   -0.8507    0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.1382   -0.9511   -0.2764
rota_matrix    0.4253   -0.3090    0.8507
rota_matrix   -0.8944    0.0000    0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.6708   -0.6882    0.2764
rota_matrix   -0.1625    0.5000    0.8507
rota_matrix   -0.7236    0.5257   -0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.4472   -0.0000   -0.8944
rota_matrix   -0.0000   -1.0000    0.0000
rota_matrix   -0.8944    0.0000    0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.9472   -0.1625   -0.2764
rota_matrix   -0.1625   -0.5000    0.8507
rota_matrix   -0.2764    0.8507    0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.8618    0.4253    0.2764
rota_matrix    0.4253    0.3090    0.8507
rota_matrix    0.2764    0.8507   -0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.3090    0.9511   -0.0000
rota_matrix    0.9511    0.3090   -0.0000
rota_matrix   -0.0000   -0.0000   -1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.0528    0.6882   -0.7236
rota_matrix    0.6882   -0.5000   -0.5257
rota_matrix   -0.7236   -0.5257   -0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.4472    0.5257    0.7236
rota_matrix   -0.8507   -0.0000   -0.5257
rota_matrix   -0.2764   -0.8507    0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.6382   -0.2629    0.7236
rota_matrix   -0.2629   -0.8090   -0.5257
rota_matrix    0.7236   -0.5257    0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.0528   -0.6882    0.7236
rota_matrix    0.6882   -0.5000   -0.5257
rota_matrix    0.7236    0.5257    0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.6708   -0.1625    0.7236
rota_matrix    0.6882    0.5000   -0.5257
rota_matrix   -0.2764    0.8507    0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.3618    0.5878    0.7236
rota_matrix   -0.2629    0.8090   -0.5257
rota_matrix   -0.8944   -0.0000    0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.4472   -0.5257    0.7236
rota_matrix    0.8507    0.0000    0.5257
rota_matrix   -0.2764    0.8507    0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.3618    0.2629    0.8944
rota_matrix    0.5878    0.8090   -0.0000
rota_matrix   -0.7236    0.5257   -0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.6708    0.6882    0.2764
rota_matrix    0.1625    0.5000   -0.8507
rota_matrix   -0.7236   -0.5257   -0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.9472    0.1625   -0.2764
rota_matrix    0.1625   -0.5000   -0.8507
rota_matrix   -0.2764   -0.8507    0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.8090   -0.5878    0.0000
rota_matrix    0.5878   -0.8090   -0.0000
rota_matrix    0.0000    0.0000    1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.4472   -0.5257   -0.7236
rota_matrix   -0.8507   -0.0000   -0.5257
rota_matrix    0.2764    0.8507   -0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.3618    0.2629   -0.8944
rota_matrix   -0.5878    0.8090    0.0000
rota_matrix    0.7236    0.5257    0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.6708    0.6882   -0.2764
rota_matrix   -0.1625    0.5000    0.8507
rota_matrix    0.7236   -0.5257    0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.9472    0.1625    0.2764
rota_matrix   -0.1625   -0.5000    0.8507
rota_matrix    0.2764   -0.8507   -0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.8090   -0.5878   -0.0000
rota_matrix   -0.5878   -0.8090    0.0000
rota_matrix   -0.0000    0.0000   -1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.4472    0.5257   -0.7236
rota_matrix    0.8507    0.0000    0.5257
rota_matrix    0.2764   -0.8507   -0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.6382   -0.2629   -0.7236
rota_matrix    0.2629   -0.8090    0.5257
rota_matrix   -0.7236   -0.5257   -0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.0528   -0.6882   -0.7236
rota_matrix   -0.6882   -0.5000    0.5257
rota_matrix   -0.7236    0.5257   -0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.6708   -0.1625   -0.7236
rota_matrix   -0.6882    0.5000    0.5257
rota_matrix    0.2764    0.8507   -0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.3618    0.5878   -0.7236
rota_matrix    0.2629    0.8090    0.5257
rota_matrix    0.8944   -0.0000   -0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.4472   -0.8507   -0.2764
rota_matrix    0.5257   -0.0000   -0.8507
rota_matrix    0.7236   -0.5257    0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.3090   -0.9511    0.0000
rota_matrix    0.9511    0.3090   -0.0000
rota_matrix    0.0000    0.0000    1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.3618   -0.5878    0.7236
rota_matrix    0.2629    0.8090    0.5257
rota_matrix   -0.8944    0.0000    0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.3618   -0.2629    0.8944
rota_matrix   -0.5878    0.8090    0.0000
rota_matrix   -0.7236   -0.5257   -0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.8618   -0.4253    0.2764
rota_matrix   -0.4253    0.3090   -0.8507
rota_matrix    0.2764   -0.8507   -0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.4472   -0.8507    0.2764
rota_matrix   -0.5257   -0.0000    0.8507
rota_matrix   -0.7236   -0.5257   -0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.3090   -0.9511   -0.0000
rota_matrix   -0.9511    0.3090    0.0000
rota_matrix   -0.0000    0.0000   -1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.3618   -0.5878   -0.7236
rota_matrix   -0.2629    0.8090   -0.5257
rota_matrix    0.8944    0.0000   -0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.3618   -0.2629   -0.8944
rota_matrix    0.5878    0.8090   -0.0000
rota_matrix    0.7236   -0.5257    0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.8618   -0.4253   -0.2764
rota_matrix    0.4253    0.3090    0.8507
rota_matrix   -0.2764   -0.8507    0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.4472    0.8507    0.2764
rota_matrix    0.5257    0.0000   -0.8507
rota_matrix   -0.7236    0.5257   -0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.1382    0.9511   -0.2764
rota_matrix   -0.4253   -0.3090   -0.8507
rota_matrix   -0.8944   -0.0000    0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.8090    0.5878    0.0000
rota_matrix   -0.5878   -0.8090    0.0000
rota_matrix    0.0000   -0.0000    1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.6382    0.2629    0.7236
rota_matrix    0.2629   -0.8090    0.5257
rota_matrix    0.7236    0.5257    0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.1382    0.4253    0.8944
rota_matrix    0.9511   -0.3090   -0.0000
rota_matrix    0.2764    0.8507   -0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.4472    0.8507   -0.2764
rota_matrix   -0.5257    0.0000    0.8507
rota_matrix    0.7236    0.5257    0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.1382    0.9511    0.2764
rota_matrix    0.4253   -0.3090    0.8507
rota_matrix    0.8944   -0.0000   -0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.8090    0.5878   -0.0000
rota_matrix    0.5878   -0.8090   -0.0000
rota_matrix   -0.0000   -0.0000   -1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.6382    0.2629   -0.7236
rota_matrix   -0.2629   -0.8090   -0.5257
rota_matrix   -0.7236    0.5257   -0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.1382    0.4253   -0.8944
rota_matrix   -0.9511   -0.3090    0.0000
rota_matrix   -0.2764    0.8507    0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000


"""
tetrahedral_a=\
"""
new_operator
rota_matrix  1.000  0.000  0.000
rota_matrix  0.000  1.000  0.000
rota_matrix  0.000  0.000  1.000
tran_orth  0.000  0.000  0.000
center_orth  0.000  0.000  0.000
new_operator
rota_matrix  0.833  0.289  -0.471
rota_matrix  -0.288  -0.500  -0.816
rota_matrix  -0.472  0.816  -0.334
tran_orth  0.000  0.000  0.000
center_orth  0.000  0.000  0.000
new_operator
rota_matrix  -0.668  -0.577  -0.471
rota_matrix  -0.578  0.003  0.816
rota_matrix  -0.469  0.817  -0.335
tran_orth  0.000  0.000  0.000
center_orth  0.000  0.000  0.000
new_operator
rota_matrix  -0.166  0.290  0.943
rota_matrix  0.867  0.499  -0.001
rota_matrix  -0.471  0.817  -0.334
tran_orth  0.000  0.000  0.000
center_orth  0.000  0.000  0.000
new_operator
rota_matrix  -0.500  0.866  -0.000
rota_matrix  -0.866  -0.500  0.000
rota_matrix  0.000  0.001  1.000
tran_orth  0.000  0.000  0.000
center_orth  0.000  0.000  0.000
new_operator
rota_matrix  -0.667  0.577  -0.471
rota_matrix  0.578  0.000  -0.816
rota_matrix  -0.471  -0.817  -0.333
tran_orth  0.000  0.000  0.000
center_orth  0.000  0.000  0.000
new_operator
rota_matrix  0.833  -0.290  -0.471
rota_matrix  0.287  -0.501  0.816
rota_matrix  -0.473  -0.815  -0.335
tran_orth  0.000  0.000  0.000
center_orth  0.000  0.000  0.000
new_operator
rota_matrix  -0.168  -0.288  0.943
rota_matrix  -0.866  0.501  -0.001
rota_matrix  -0.472  -0.816  -0.334
tran_orth  0.000  0.000  0.000
center_orth  0.000  0.000  0.000
new_operator
rota_matrix  -0.500  -0.866  0.000
rota_matrix  0.866  -0.500  0.001
rota_matrix  -0.000  0.000  1.000
tran_orth  0.000  0.000  0.000
center_orth  0.000  0.000  0.000
new_operator
rota_matrix  -0.166  -0.866  -0.471
rota_matrix  -0.289  0.500  -0.817
rota_matrix  0.943  0.000  -0.333
tran_orth  0.000  0.000  0.000
center_orth  0.000  0.000  0.000
new_operator
rota_matrix  -0.165  0.866  -0.471
rota_matrix  0.291  0.499  0.816
rota_matrix  0.942  -0.002  -0.335
tran_orth  0.000  0.000  0.000
center_orth  0.000  0.000  0.000
new_operator
rota_matrix  0.334  -0.001  0.943
rota_matrix  -0.001  -1.000  -0.001
rota_matrix  0.943  -0.001  -0.334
tran_orth  0.000  0.000  0.000
center_orth  0.000  0.000  0.000
"""


tetrahedral_b=\
"""
new_operator
rota_matrix  1.000  0.000  0.000
rota_matrix  0.000  1.000  0.000
rota_matrix  0.000  0.000  1.000
tran_orth  0.000  0.000  0.000
center_orth  0.000  0.000  0.000
new_operator
rota_matrix  -1.000  0.001  0.001
rota_matrix  0.001  0.333  0.943
rota_matrix  0.000  0.943  -0.333
tran_orth  0.000  0.000  0.000
center_orth  0.000  0.000  0.000
new_operator
rota_matrix  0.499  0.288  0.818
rota_matrix  0.867  -0.166  -0.470
rota_matrix  0.001  0.943  -0.333
tran_orth  0.000  0.000  0.000
center_orth  0.000  0.000  0.000
new_operator
rota_matrix  0.498  -0.287  -0.818
rota_matrix  -0.867  -0.167  -0.469
rota_matrix  -0.002  0.943  -0.332
tran_orth  0.000  0.000  0.000
center_orth  0.000  0.000  0.000
new_operator
rota_matrix  -0.500  0.866  -0.000
rota_matrix  -0.866  -0.500  -0.000
rota_matrix  -0.000  0.000  1.000
tran_orth  0.000  0.000  0.000
center_orth  0.000  0.000  0.000
new_operator
rota_matrix  0.499  -0.866  0.001
rota_matrix  -0.289  -0.166  0.943
rota_matrix  -0.817  -0.471  -0.334
tran_orth  0.000  0.000  0.000
center_orth  0.000  0.000  0.000
new_operator
rota_matrix  -0.498  0.287  0.818
rota_matrix  -0.289  0.834  -0.469
rota_matrix  -0.817  -0.471  -0.332
tran_orth  0.000  0.000  0.000
center_orth  0.000  0.000  0.000
new_operator
rota_matrix  0.000  0.574  -0.819
rota_matrix  0.578  -0.668  -0.469
rota_matrix  -0.816  -0.473  -0.332
tran_orth  0.000  0.000  0.000
center_orth  0.000  0.000  0.000
new_operator
rota_matrix  -0.500  -0.866  0.000
rota_matrix  0.866  -0.500  -0.000
rota_matrix  0.000  -0.000  1.000
tran_orth  0.000  0.000  0.000
center_orth  0.000  0.000  0.000
new_operator
rota_matrix  0.501  0.866  0.001
rota_matrix  0.289  -0.168  0.943
rota_matrix  0.816  -0.472  -0.334
tran_orth  0.000  0.000  0.000
center_orth  0.000  0.000  0.000
new_operator
rota_matrix  -0.000  -0.577  0.817
rota_matrix  -0.578  -0.667  -0.470
rota_matrix  0.816  -0.472  -0.333
tran_orth  0.000  0.000  0.000
center_orth  0.000  0.000  0.000
new_operator
rota_matrix  -0.498  -0.288  -0.818
rota_matrix  0.289  0.834  -0.470
rota_matrix  0.818  -0.470  -0.332
tran_orth  0.000  0.000  0.000
center_orth  0.000  0.000  0.000
"""
octahedral_a=\
"""
new_operator
rota_matrix  1.0  0.0  0.0
rota_matrix  0.0  1.0  0.0
rota_matrix  0.0  0.0  1.0
tran_orth  0.0  0.0  0.0
center_orth  0.0  0.0  0.0
new_operator
rota_matrix  0.0  -1.0  -0.0
rota_matrix  1.0  0.0  0.0
rota_matrix  -0.0  -0.0  1.0
tran_orth  0.0  0.0  0.0
center_orth  0.0  0.0  0.0
new_operator
rota_matrix  -1.0  0.0  -0.0
rota_matrix  -0.0  -1.0  0.0
rota_matrix  -0.0  0.0  1.0
tran_orth  0.0  0.0  0.0
center_orth  0.0  0.0  0.0
new_operator
rota_matrix  0.0  1.0  -0.0
rota_matrix  -0.0  -0.0  -1.0
rota_matrix  -1.0  0.0  0.0
tran_orth  0.0  0.0  0.0
center_orth  0.0  0.0  0.0
new_operator
rota_matrix  0.0  -0.0  1.0
rota_matrix  -0.0  1.0  0.0
rota_matrix  -1.0  -0.0  0.0
tran_orth  0.0  0.0  0.0
center_orth  0.0  0.0  0.0
new_operator
rota_matrix  -1.0  0.0  -0.0
rota_matrix  0.0  -0.0  -1.0
rota_matrix  -0.0  -1.0  0.0
tran_orth  0.0  0.0  0.0
center_orth  0.0  0.0  0.0
new_operator
rota_matrix  -0.0  0.0  -1.0
rota_matrix  1.0  -0.0  -0.0
rota_matrix  -0.0  -1.0  -0.0
tran_orth  0.0  0.0  0.0
center_orth  0.0  0.0  0.0
new_operator
rota_matrix  -0.0  1.0  -0.0
rota_matrix  -1.0  -0.0  0.0
rota_matrix  0.0  0.0  1.0
tran_orth  0.0  0.0  0.0
center_orth  0.0  0.0  0.0
new_operator
rota_matrix  -0.0  -1.0  -0.0
rota_matrix  0.0  0.0  -1.0
rota_matrix  1.0  -0.0  0.0
tran_orth  0.0  0.0  0.0
center_orth  0.0  0.0  0.0
new_operator
rota_matrix  1.0  -0.0  -0.0
rota_matrix  -0.0  0.0  -1.0
rota_matrix  0.0  1.0  0.0
tran_orth  0.0  0.0  0.0
center_orth  0.0  0.0  0.0
new_operator
rota_matrix  -0.0  -0.0  1.0
rota_matrix  0.0  -1.0  -0.0
rota_matrix  1.0  0.0  0.0
tran_orth  0.0  0.0  0.0
center_orth  0.0  0.0  0.0
new_operator
rota_matrix  -0.0  -0.0  1.0
rota_matrix  1.0  0.0  0.0
rota_matrix  -0.0  1.0  0.0
tran_orth  0.0  0.0  0.0
center_orth  0.0  0.0  0.0
new_operator
rota_matrix  0.0  0.0  -1.0
rota_matrix  -0.0  -1.0  -0.0
rota_matrix  -1.0  0.0  -0.0
tran_orth  0.0  0.0  0.0
center_orth  0.0  0.0  0.0
new_operator
rota_matrix  0.0  -1.0  0.0
rota_matrix  -0.0  0.0  1.0
rota_matrix  -1.0  -0.0  -0.0
tran_orth  0.0  0.0  0.0
center_orth  0.0  0.0  0.0
new_operator
rota_matrix  -0.0  1.0  0.0
rota_matrix  0.0  -0.0  1.0
rota_matrix  1.0  0.0  -0.0
tran_orth  0.0  0.0  0.0
center_orth  0.0  0.0  0.0
new_operator
rota_matrix  0.0  -0.0  -1.0
rota_matrix  -1.0  0.0  -0.0
rota_matrix  0.0  1.0  -0.0
tran_orth  0.0  0.0  0.0
center_orth  0.0  0.0  0.0
new_operator
rota_matrix  -1.0  -0.0  0.0
rota_matrix  0.0  0.0  1.0
rota_matrix  -0.0  1.0  -0.0
tran_orth  0.0  0.0  0.0
center_orth  0.0  0.0  0.0
new_operator
rota_matrix  1.0  0.0  0.0
rota_matrix  0.0  -1.0  -0.0
rota_matrix  0.0  0.0  -1.0
tran_orth  0.0  0.0  0.0
center_orth  0.0  0.0  0.0
new_operator
rota_matrix  0.0  1.0  0.0
rota_matrix  1.0  -0.0  -0.0
rota_matrix  -0.0  0.0  -1.0
tran_orth  0.0  0.0  0.0
center_orth  0.0  0.0  0.0
new_operator
rota_matrix  -1.0  0.0  0.0
rota_matrix  0.0  1.0  -0.0
rota_matrix  -0.0  -0.0  -1.0
tran_orth  0.0  0.0  0.0
center_orth  0.0  0.0  0.0
new_operator
rota_matrix  0.0  -1.0  0.0
rota_matrix  -1.0  -0.0  -0.0
rota_matrix  0.0  -0.0  -1.0
tran_orth  0.0  0.0  0.0
center_orth  0.0  0.0  0.0
new_operator
rota_matrix  1.0  0.0  0.0
rota_matrix  -0.0  -0.0  1.0
rota_matrix  0.0  -1.0  -0.0
tran_orth  0.0  0.0  0.0
center_orth  0.0  0.0  0.0
new_operator
rota_matrix  -0.0  0.0  -1.0
rota_matrix  0.0  1.0  0.0
rota_matrix  1.0  -0.0  -0.0
tran_orth  0.0  0.0  0.0
center_orth  0.0  0.0  0.0
new_operator
rota_matrix  -0.0  0.0  1.0
rota_matrix  -1.0  -0.0  -0.0
rota_matrix  0.0  -1.0  0.0
tran_orth  0.0  0.0  0.0
center_orth  0.0  0.0  0.0
"""

icosahedral_b=\
"""

Summary of NCS information
Mon Jul  2 11:07:59 2018
/Users/terwill/Downloads/view/ncs_text

source_info icosahedral_b.dat




new_ncs_group
new_operator

rota_matrix    1.0000   -0.0000    0.0000
rota_matrix    0.0000    1.0000    0.0000
rota_matrix    0.0000   -0.0000    1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.3090    0.9511   -0.0001
rota_matrix   -0.9511    0.3090   -0.0002
rota_matrix   -0.0001    0.0002    1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.8090    0.5878   -0.0003
rota_matrix   -0.5878   -0.8090   -0.0001
rota_matrix   -0.0003    0.0001    1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.8090   -0.5878   -0.0003
rota_matrix    0.5878   -0.8090    0.0001
rota_matrix   -0.0003   -0.0001    1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.3090   -0.9511   -0.0001
rota_matrix    0.9511    0.3090    0.0002
rota_matrix   -0.0001   -0.0002    1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.9473   -0.1623    0.2761
rota_matrix   -0.1623   -0.5000   -0.8507
rota_matrix    0.2761   -0.8507    0.4473
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.4471   -0.5256   -0.7238
rota_matrix    0.8508    0.0000   -0.5256
rota_matrix    0.2762   -0.8508    0.4471
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.6709   -0.1623   -0.7236
rota_matrix    0.6881    0.5000    0.5259
rota_matrix    0.2764   -0.8507    0.4471
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.8617    0.4255    0.2765
rota_matrix   -0.4255    0.3090    0.8506
rota_matrix    0.2765   -0.8506    0.4473
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.1384    0.4255    0.8943
rota_matrix   -0.9511   -0.3090   -0.0002
rota_matrix    0.2763   -0.8506    0.4474
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.8617   -0.4255   -0.2765
rota_matrix   -0.4255    0.3090    0.8506
rota_matrix   -0.2765    0.8506   -0.4473
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.6709    0.1623    0.7236
rota_matrix    0.6881    0.5000    0.5259
rota_matrix   -0.2764    0.8507   -0.4471
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.4471    0.5256    0.7238
rota_matrix    0.8508    0.0000   -0.5256
rota_matrix   -0.2762    0.8508   -0.4471
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.9473    0.1623   -0.2761
rota_matrix   -0.1623   -0.5000   -0.8507
rota_matrix   -0.2761    0.8507   -0.4473
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.1384   -0.4255   -0.8943
rota_matrix   -0.9511   -0.3090   -0.0002
rota_matrix   -0.2763    0.8506   -0.4474
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.8090    0.5878    0.0003
rota_matrix    0.5878   -0.8090    0.0001
rota_matrix    0.0003    0.0001   -1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.8090   -0.5878    0.0003
rota_matrix   -0.5878   -0.8090   -0.0001
rota_matrix    0.0003   -0.0001   -1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.3090   -0.9511    0.0001
rota_matrix   -0.9511    0.3090   -0.0002
rota_matrix    0.0001   -0.0002   -1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -1.0000    0.0000    0.0000
rota_matrix    0.0000    1.0000    0.0000
rota_matrix    0.0000    0.0000   -1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.3090    0.9511    0.0001
rota_matrix    0.9511    0.3090    0.0002
rota_matrix    0.0001    0.0002   -1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.1384    0.9511    0.2763
rota_matrix   -0.4255   -0.3090    0.8506
rota_matrix    0.8943    0.0002    0.4474
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.4476    0.0000    0.8943
rota_matrix    0.0000   -1.0000   -0.0000
rota_matrix    0.8943   -0.0000    0.4476
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.1384   -0.9511    0.2763
rota_matrix    0.4255   -0.3090   -0.8506
rota_matrix    0.8943   -0.0002    0.4474
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.3618   -0.5878   -0.7236
rota_matrix    0.2630    0.8090   -0.5257
rota_matrix    0.8944   -0.0001    0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.3618    0.5878   -0.7236
rota_matrix   -0.2630    0.8090    0.5257
rota_matrix    0.8944    0.0001    0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.4471   -0.8508   -0.2762
rota_matrix   -0.5256    0.0000   -0.8508
rota_matrix    0.7238    0.5256   -0.4471
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.3618   -0.2630   -0.8944
rota_matrix   -0.5878    0.8090   -0.0001
rota_matrix    0.7236    0.5257   -0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.6709    0.6881   -0.2764
rota_matrix    0.1623    0.5000    0.8507
rota_matrix    0.7236    0.5259   -0.4471
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.0531    0.6881    0.7237
rota_matrix    0.6881   -0.5000    0.5259
rota_matrix    0.7237    0.5259   -0.4469
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.6379   -0.2630    0.7238
rota_matrix    0.2630   -0.8090   -0.5257
rota_matrix    0.7238    0.5257   -0.4469
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.0531   -0.6881   -0.7237
rota_matrix    0.6881   -0.5000    0.5259
rota_matrix   -0.7237   -0.5259    0.4469
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.6709   -0.6881    0.2764
rota_matrix    0.1623    0.5000    0.8507
rota_matrix   -0.7236   -0.5259    0.4471
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.3618    0.2630    0.8944
rota_matrix   -0.5878    0.8090   -0.0001
rota_matrix   -0.7236   -0.5257    0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.4471    0.8508    0.2762
rota_matrix   -0.5256    0.0000   -0.8508
rota_matrix   -0.7238   -0.5256    0.4471
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.6379    0.2630   -0.7238
rota_matrix    0.2630   -0.8090   -0.5257
rota_matrix   -0.7238   -0.5257    0.4469
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.3618    0.5878    0.7236
rota_matrix    0.2630    0.8090   -0.5257
rota_matrix   -0.8944    0.0001   -0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.1384    0.9511   -0.2763
rota_matrix    0.4255   -0.3090   -0.8506
rota_matrix   -0.8943    0.0002   -0.4474
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.4476   -0.0000   -0.8943
rota_matrix   -0.0000   -1.0000    0.0000
rota_matrix   -0.8943    0.0000   -0.4476
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.1384   -0.9511   -0.2763
rota_matrix   -0.4255   -0.3090    0.8506
rota_matrix   -0.8943   -0.0002   -0.4474
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.3618   -0.5878    0.7236
rota_matrix   -0.2630    0.8090    0.5257
rota_matrix   -0.8944   -0.0001   -0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.1384   -0.4255    0.8943
rota_matrix    0.9511   -0.3090    0.0002
rota_matrix    0.2763    0.8506    0.4474
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.8617   -0.4255    0.2765
rota_matrix    0.4255    0.3090   -0.8506
rota_matrix    0.2765    0.8506    0.4473
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.6709    0.1623   -0.7236
rota_matrix   -0.6881    0.5000   -0.5259
rota_matrix    0.2764    0.8507    0.4471
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.4471    0.5256   -0.7238
rota_matrix   -0.8508    0.0000    0.5256
rota_matrix    0.2762    0.8508    0.4471
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.9473    0.1623    0.2761
rota_matrix    0.1623   -0.5000    0.8507
rota_matrix    0.2761    0.8507    0.4473
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.0531    0.6881   -0.7237
rota_matrix   -0.6881   -0.5000   -0.5259
rota_matrix   -0.7237    0.5259    0.4469
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.6379   -0.2630   -0.7238
rota_matrix   -0.2630   -0.8090    0.5257
rota_matrix   -0.7238    0.5257    0.4469
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.4471   -0.8508    0.2762
rota_matrix    0.5256    0.0000    0.8508
rota_matrix   -0.7238    0.5256    0.4471
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.3618   -0.2630    0.8944
rota_matrix    0.5878    0.8090    0.0001
rota_matrix   -0.7236    0.5257    0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.6709    0.6881    0.2764
rota_matrix   -0.1623    0.5000   -0.8507
rota_matrix   -0.7236    0.5259    0.4471
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.3618    0.2630   -0.8944
rota_matrix    0.5878    0.8090    0.0001
rota_matrix    0.7236   -0.5257   -0.4472
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.4471    0.8508   -0.2762
rota_matrix    0.5256    0.0000    0.8508
rota_matrix    0.7238   -0.5256   -0.4471
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.6379    0.2630    0.7238
rota_matrix   -0.2630   -0.8090    0.5257
rota_matrix    0.7238   -0.5257   -0.4469
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.0531   -0.6881    0.7237
rota_matrix   -0.6881   -0.5000   -0.5259
rota_matrix    0.7237   -0.5259   -0.4469
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.6709   -0.6881   -0.2764
rota_matrix   -0.1623    0.5000   -0.8507
rota_matrix    0.7236   -0.5259   -0.4471
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.4471   -0.5256    0.7238
rota_matrix   -0.8508    0.0000    0.5256
rota_matrix   -0.2762   -0.8508   -0.4471
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.6709   -0.1623    0.7236
rota_matrix   -0.6881    0.5000   -0.5259
rota_matrix   -0.2764   -0.8507   -0.4471
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix   -0.8617    0.4255   -0.2765
rota_matrix    0.4255    0.3090   -0.8506
rota_matrix   -0.2765   -0.8506   -0.4473
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.1384    0.4255   -0.8943
rota_matrix    0.9511   -0.3090    0.0002
rota_matrix   -0.2763   -0.8506   -0.4474
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.9473   -0.1623   -0.2761
rota_matrix    0.1623   -0.5000    0.8507
rota_matrix   -0.2761   -0.8507   -0.4473
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000


"""


class ncs_group:
  """NCS group: one group of NCS operators and center and where it applies"""
  def __init__(self, ncs_rota_matr=None, center_orth=None, trans_orth=None,
      chain_residue_id=None,source_of_ncs_info=None,rmsd_list=None,
      ncs_domain_pdb=None,
      residues_in_common_list=None,cc=None,note=None,
       exclude_h=None,exclude_d=None):
    self._chain_residue_id=chain_residue_id  # just one of these
    self._rmsd_list=rmsd_list
    self._residues_in_common_list=residues_in_common_list
    self._centers=center_orth
    self._translations_orth=trans_orth
    self._rota_matrices=ncs_rota_matr
    if self._centers is not None:
      self._n_ncs_oper=len(self._centers)
    elif self._rmsd_list is not None:
      self._n_ncs_oper=len(self._rmsd_list)
    else:
      self._n_ncs_oper=0
    self._source_of_ncs_info=source_of_ncs_info
    self._ncs_domain_pdb=ncs_domain_pdb
    self._cc=cc
    self._note=note
    self._exclude_h=exclude_h
    self._exclude_d=exclude_d
    self._have_helical_symmetry=False
    self._have_point_group_symmetry=False

  def __repr__(self):
    """Return string representation of this NCS group"""
    return "NCS group with %s ops" %(self._n_ncs_oper)

  def apply_cob_to_vector(self,vector=None,
         change_of_basis_operator=None,
         coordinate_offset=None,
         unit_cell=None,new_unit_cell=None):
    """Apply change of basis to vector"""
    if coordinate_offset is not None:
      from scitbx.math import  matrix
      new_vector=matrix.col(vector)+matrix.col(coordinate_offset)
    elif change_of_basis_operator:
      frac=unit_cell.fractionalize(vector)
      new_frac = change_of_basis_operator.c() * frac
      new_vector=new_unit_cell.orthogonalize(new_frac)
    return new_vector


  def copy_rot_trans(self,list_of_matrices,list_of_translations,
      change_of_basis_operator=None,
      coordinate_offset=None,
      scale_factor=None,
      unit_cell=None,new_unit_cell=None):
    """Copy a list of matrices and translations, applying change of
    basis operator.  If change_of_basis_operator is None, then return copy
    of what we have"""
    from scitbx.math import  matrix
    new_list_of_matrices=[]
    new_list_of_translations=[]
    if change_of_basis_operator is not None:
      a=  matrix.sqr(new_unit_cell.orthogonalization_matrix()) \
        * change_of_basis_operator.c().as_rational().r \
        * matrix.sqr(unit_cell.fractionalization_matrix())
      a_inv=a.inverse()
    else:
      a=None
    for ncs_r,ncs_t in zip(list_of_matrices,list_of_translations):
      if scale_factor is not None:
        # expand by scale_factor about (0,0,0). Affects all translations
        new_list_of_matrices.append(deepcopy(ncs_r))
        new_list_of_translations.append(scale_factor*matrix.col(ncs_t))
      elif change_of_basis_operator is None and coordinate_offset is None:
        new_list_of_matrices.append(deepcopy(ncs_r))
        new_list_of_translations.append(deepcopy(ncs_t))
      elif coordinate_offset is not None:
        # TT 2015-11-02 special case of below where cob is just a translation
        # R' =  R
        # T' =  T + t - R t
        new_list_of_matrices.append(deepcopy(ncs_r))  # these are the same
        from scitbx import matrix
        delta = ncs_r * matrix.col(coordinate_offset)
        t_prime=matrix.col(ncs_t) + \
          matrix.col(coordinate_offset) - matrix.col(delta)
        new_list_of_translations.append(t_prime)

      else:
        # tt 2011-10-02
        # Formula for conversion of NCS rotation matrix and translation
        # relating two points in coordinate system to a second coordinate system
        # The change-of-basis operator is new_x = a x + t
        # The NCS operator is y = R x + T (in original coordinate system)
        # Then if NCS operator in new coordinate system is y' = R' x' + T':
        # R' = a R a_inv
        # T' = a T + t - a R a_inv t = transformed T minus R' * t
        #
        # Derivation:
        # x' and y' (values in new coordinate system) can be written:
        #   x'=a x + t
        #   y'=a y + t
        # Or rewriting:
        #   x= a_inv (x' - t)
        #   y= a_inv (y' - t)
        # Then as y = R x + T  (in original coordinate system), we can say:
        #   a_inv (y' - t) = R (a_inv (x' - t) ) + T
        # Or...
        #   y' = [a R a_inv] x' - [a R a_inv ] t + t + a t
        # So that:
        #   R' = a R a_inv
        #   T' = a T + t - a R a_inv t = transformed T minus R' * t

        # matrices are a ncs_r a_inv
        ncs_r_prime=a * ncs_r * a_inv
        new_list_of_matrices.append(ncs_r_prime)

        # translation vectors are partly a * ncs_t + t (the cob operator)
        frac=unit_cell.fractionalize(ncs_t)
        new_frac = change_of_basis_operator.c() * frac
        new_ncs_t=new_unit_cell.orthogonalize(new_frac)
        # translation vectors depend on the change of basis and the rotation
        # as well as the change-of-basis operator
        t_as_col=change_of_basis_operator.c().t().as_rational().as_float()
        # the basis translation in orig coordinate system
        cob_trans=unit_cell.orthogonalize(t_as_col)
        # correction for the basis translation in new coordinate system
        delta = ncs_r_prime * cob_trans
        t_prime=matrix.col(new_ncs_t) - matrix.col(delta)
        new_list_of_translations.append(t_prime)

    return new_list_of_matrices,new_list_of_translations

  def copy_vector_list(self,list_of_vectors,
      change_of_basis_operator=None,
      coordinate_offset=None,
      scale_factor=None,
         unit_cell=None,new_unit_cell=None):
    """Copy a list of vectors, optionally applying a change of basis,
       coordinate offset, scale factor"""
    from scitbx.math import  matrix
    new_vector_list=[]
    for vector in list_of_vectors:
      if scale_factor is not None:
        new_vector=scale_factor*matrix.col(vector)
      elif change_of_basis_operator is None and coordinate_offset is None:
        new_vector=deepcopy(vector)
      else:
        new_vector=self.apply_cob_to_vector(vector=vector,
          change_of_basis_operator=change_of_basis_operator,
          coordinate_offset=coordinate_offset,
           unit_cell=unit_cell,new_unit_cell=new_unit_cell)
      new_vector_list.append(new_vector)
    return new_vector_list

  def deep_copy(self,change_of_basis_operator=None,unit_cell=None,
      coordinate_offset=None,
      new_unit_cell=None,
      scale_factor=None,
      extract_point_group_symmetry=None,
      ops_to_keep=None,
      hierarchy_to_match_order=None):  # make full copy;
    """Make a deep copy of this NCS group"""
    # optionally apply change-of-basis operator (requires old, new unit cells)
    # optionally apply coordinate_offset (adding coordinate_offset to coords)
    # optionally sort operators to match order in hierarchy
    # optionally apply scale factor (magnification)
    # optionally keep only ops_to_keep operators
    # optionaly extract point group symmetry

    # Can do only one of the above five things at most
    assert [change_of_basis_operator,scale_factor,
      coordinate_offset,ops_to_keep,hierarchy_to_match_order,
      extract_point_group_symmetry].count(None)>=4

    if hierarchy_to_match_order and self._chain_residue_id is not None:
      return self.deep_copy_order(
         hierarchy_to_match_order=hierarchy_to_match_order)

    if ops_to_keep is not None:
      return self.deep_copy_ops_to_keep(ops_to_keep=ops_to_keep)

    if extract_point_group_symmetry:
      return self.extract_point_group_symmetry()

    from mmtbx.ncs.ncs import ncs
    new=ncs_group()
    new._chain_residue_id=self._chain_residue_id
    new._rmsd_list=deepcopy(self._rmsd_list)
    new._residues_in_common_list=deepcopy(self._residues_in_common_list)

    # centers simply get affected by the change of basis operator if present
    #   or scale_factor
    new._centers=self.copy_vector_list(self._centers,
      change_of_basis_operator=change_of_basis_operator,
      coordinate_offset=coordinate_offset,
      scale_factor=scale_factor,
         unit_cell=unit_cell,new_unit_cell=new_unit_cell)

    # matrices and translations may need to be adjusted if change of basis set
    #  Scale_factor applied to translations
    new._rota_matrices,new._translations_orth=self.copy_rot_trans(
       self._rota_matrices,self._translations_orth,
         scale_factor=scale_factor,
         coordinate_offset=coordinate_offset,
         unit_cell=unit_cell,new_unit_cell=new_unit_cell)

    new._n_ncs_oper=deepcopy(self._n_ncs_oper)
    new._source_of_ncs_info=self._source_of_ncs_info
    new._ncs_domain_pdb=deepcopy(self._ncs_domain_pdb)
    new._cc=deepcopy(self._cc)
    new._note=deepcopy(self._note)
    new._exclude_h=self._exclude_h
    new._exclude_d=self._exclude_d
    return new

  def deep_copy_ops_to_keep(self,ops_to_keep=None):
    """Deep copy, but keep only ops_to_keep operators"""
    assert ops_to_keep is not None

    new=self.deep_copy() # exact copy.  Now remove all except ops_to_keep

    new._n_ncs_oper=0
    new._rota_matrices=[]
    new._rota_matrices_inv=[]
    new._translations_orth=[]
    new._translations_orth_inv=[]
    new._residues_in_common_list=[]
    new._centers=[]

    new._rmsd_list=[]
    new._residues_in_common_list=[]

    if self._chain_residue_id:
      [group,residue_range_list]=self._chain_residue_id
    else:
      group=self._n_ncs_oper*[None]
      residue_range_list=self._n_ncs_oper*[None]
    new_residue_range_list=[]
    new_group=[]

    for i in range(self._n_ncs_oper):
      if i in ops_to_keep:
        new._n_ncs_oper+=1
        new._rota_matrices.append(self.rota_matrices()[i])
        new._rota_matrices_inv.append(self.rota_matrices_inv()[i])
        new._translations_orth.append(self.translations_orth()[i])
        new._translations_orth_inv.append(self.translations_orth_inv()[i])
        new._rmsd_list.append(self._rmsd_list[i])
        new._residues_in_common_list.append(self._residues_in_common_list[i])
        new._centers.append(self._centers[i])
        new_residue_range_list.append(residue_range_list[i])
        new_group.append(group[i])

    if self._chain_residue_id:
      new._chain_residue_id=[new_group,new_residue_range_list]
    else:
      new._chain_residue_id=None
    return new

  def get_order_dict(self,hierarchy_to_match_order=None):
    """Set up dicts for chain_id_from_index and index_from_chain_id"""
    self.chain_id_from_index={}
    self.index_from_chain_id={}
    i=0
    for m in hierarchy_to_match_order.models()[:1]:
      for chain in m.chains():
        id=remove_single_quotes(chain.id)
        self.chain_id_from_index[i]=id
        self.index_from_chain_id[id]=i
        i+=1

  def get_new_group(self,hierarchy_to_match_order=None):
    """Get a new group. Change the order of the operators to match hierarchy"""
    self.get_order_dict(hierarchy_to_match_order=hierarchy_to_match_order)

    # figure out what is the new order of groups
    sort_list=[]
    [group,residue_range_list] = self._chain_residue_id
    for x in group:
      id=self.index_from_chain_id.get(x)
      assert id is not None
      sort_list.append([id,x])
    sort_list.sort(key=itemgetter(0))
    new_group=[]
    for [id,x] in sort_list:
      new_group.append(x)

    # Identify which operator is the new first one
    first_chain_id=new_group[0]
    i=0
    first_op=None
    for chain_id in group:
      if chain_id==first_chain_id:
        first_op=i
      i+=1
    assert first_op is not None
    return new_group,first_op

  def deep_copy_order(self,hierarchy_to_match_order=None):
    """Make a full copy; reorder to match hierarchy"""
    assert self._chain_residue_id is not None

    # Get the new order of chain IDs
    new_group,first_op=self.get_new_group(
      hierarchy_to_match_order=hierarchy_to_match_order)

    from mmtbx.ncs.ncs import ncs
    new=ncs_group()
    new._chain_residue_id=self._chain_residue_id
    new._rmsd_list=deepcopy(self._rmsd_list)
    new._residues_in_common_list=deepcopy(self._residues_in_common_list)

    # now adjust all the operators so that the reference is first_op
    # operator j maps copy j on to copy 0.
    # we want operator j to map copy j onto copy first_op. First map to
    #  copy 0 then map copy 0 to copy first_op
    #  Rx + T =  Rinv_fo (Rj x + Tj) + Tinv_fo
    #  ==> R= Rinv_fo Rj   T = Rinv_fo Tj + Tinv_fo
    # check: if fo=0 then Rinv_fo=U Tinv_fo=0 -> R=Rj and T=Tj ok.

    # centers are the same
    new._centers=deepcopy(self._centers)
    new._rota_matrices=[]
    new._translations_orth=[]
    r_first_op_inv=self._rota_matrices[first_op].inverse()
    t_first_op_inv=-1.*r_first_op_inv*self._translations_orth[first_op]

    for r,t in zip(self._rota_matrices,self._translations_orth):
      new_r=r_first_op_inv*r
      new_t=r_first_op_inv*t + t_first_op_inv

      new._rota_matrices.append(new_r)
      new._translations_orth.append(new_t)
    new._n_ncs_oper=deepcopy(self._n_ncs_oper)
    new._source_of_ncs_info=self._source_of_ncs_info
    new._ncs_domain_pdb=deepcopy(self._ncs_domain_pdb)
    new._cc=deepcopy(self._cc)
    new._note=deepcopy(self._note)
    new._exclude_h=self._exclude_h
    new._exclude_d=self._exclude_d

    # Now just change the order of everything
    translations_orth=[]
    rota_matrices=[]
    centers=[]
    rmsd_list=[]
    residues_in_common_list=[]
    group=[]
    residue_range_list=[]
    # NOTE: [group,residue_range_list] = new._chain_residue_id

    for chain_id in new_group:
      for t,r,c,rmsd,res_in_common,g,res_range in zip(
          new._translations_orth,
          new._rota_matrices,
          new._centers,
          new._rmsd_list,
          new._residues_in_common_list,
          new._chain_residue_id[0],
          new._chain_residue_id[1]):
        if g==chain_id: # this is the one
          translations_orth.append(t)
          rota_matrices.append(r)
          centers.append(c)
          rmsd_list.append(rmsd)
          residues_in_common_list.append(res_in_common)
          group.append(g)
          residue_range_list.append(res_range)

    new._translations_orth=translations_orth
    new._rota_matrices=rota_matrices
    new._centers=centers
    new._rmsd_list=rmsd_list  # NOTE will not be correct perhaps leave out
    new._residues_in_common_list=residues_in_common_list
    new._chain_residue_id=[group,residue_range_list]

    assert group==new_group
    chain_residue_id=[group,residue_range_list]
    return new


  def display_summary(self,verbose=None):
    """Summarize this NCS group"""
    text=""
    text+="\nSummary of NCS group with "+str(self.n_ncs_oper())+" operators:"
    i=0
    if verbose:
      if self._chain_residue_id:
        text+="\nID of chain/residue where these apply: "+\
           str(self._chain_residue_id)
      if self._rmsd_list and self._chain_residue_id:
        text+="\nRMSD (A) from chain "+str(self._chain_residue_id[0][0])+':'+\
         self.print_list(self._rmsd_list)
      if self._residues_in_common_list and self._chain_residue_id:
        text+="\nNumber of residues matching chain "+\
            str(self._chain_residue_id[0][0])+':'+\
             str(self._residues_in_common_list)
      if self._source_of_ncs_info:
        text+="\nSource of NCS info: "+str(self._source_of_ncs_info)
      if self._ncs_domain_pdb:
        text+="\nNCS domains represented by: "+str(self._ncs_domain_pdb)
      if self._cc:
        text+="\nCorrelation of NCS: "+str(self._cc)
      if self._note:
        text+="\nNOTE: "+str(self._note)
      for center,trans_orth,ncs_rota_matr in zip (
         self._centers, self._translations_orth,self._rota_matrices):
        if center is None: continue
        i+=1
        text+="\n\nOPERATOR "+str(i)
        text+="\nCENTER: "+" %8.4f  %8.4f  %8.4f" %tuple(center)
        r = ncs_rota_matr.elems
        text+="\n\nROTA 1: "+" %8.4f  %8.4f  %8.4f" %tuple(r[0:3])
        text+="\nROTA 2: "+" %8.4f  %8.4f  %8.4f" %tuple(r[3:6])
        text+="\nROTA 3: "+" %8.4f  %8.4f  %8.4f" %tuple(r[6:9])
        text+="\nTRANS:  "+" %8.4f  %8.4f  %8.4f" %tuple(trans_orth)
      text+="\n"
    return text

  def format_group_specification(self):
    """Write out NCS group as text file with ncs_spec format"""
    if not self._chain_residue_id or len(self._chain_residue_id)<2:
      return ""

    # Need to test for existence because we might have operators or
    # chain specifications but not both

    if self._chain_residue_id is not None:
      [group,residue_range_list] = self._chain_residue_id
    else:
      group=self.n_ncs_oper()*[None]
      residue_range_list=self.n_ncs_oper()*[None]

    if self._centers is not None:
      [centers, translations_orth,rota_matrices]=\
         [self._centers, self._translations_orth,self._rota_matrices]
    else:
      centers=self.n_ncs_oper()*[None]
      translations_orth=self.n_ncs_oper()*[None]
      rota_matrices=self.n_ncs_oper()*[None]

    if self._rmsd_list is not None:
       rmsd_list=self._rmsd_list
    else:
      rmsd_list=self.n_ncs_oper()*[None]

    if self._residues_in_common_list is not None:
       residues_in_common_list=self._residues_in_common_list
    else:
      residues_in_common_list=self.n_ncs_oper()*[None]

    text="\nnew_ncs_group\n"
    if self._cc is not None: text+="NCS_CC "+str(self._cc)+"\n"
    if self._note is not None: text+="NOTE "+str(self._note)+"\n"
    if self._ncs_domain_pdb is not None:
      text+="  NCS_DOMAIN_PDB "+str(self._ncs_domain_pdb)+"\n"

    count=0
    for id,residue_ranges, center,trans_orth,ncs_rota_matr, \
        rmsd,common in zip (
        group,residue_range_list,
        centers, translations_orth,rota_matrices,
        rmsd_list,residues_in_common_list):
     count+=1
     text+='new_operator\n'
     if center is not None:
       for j in range(3):
         text+="\nrota_matrix "+" %8.4f  %8.4f  %8.4f" %tuple(
          ncs_rota_matr.elems[j*3:j*3+3])
       text+="\ntran_orth  "+" %8.4f  %8.4f  %8.4f" %tuple(trans_orth)
       text+="\n"
       text+="\ncenter_orth "+" %8.4f  %8.4f  %8.4f" %tuple(center)
       text+="\n"

     if id is not None: text+="CHAIN "+str(id)+ "\n"
     if rmsd is not None: text+="RMSD "+str(rmsd)+ "\n"
     if common is not None: text+="MATCHING "+str(common)+ "\n"

     if residue_ranges is not None and residue_ranges:
       for residue_range in residue_ranges:
         text+="  RESSEQ "
         text+=str(residue_range[0])+":"+str(residue_range[1])+"\n"
       text+="\n"
    return text

  def format_for_phenix_refine(self, prefix="pdb_interpretation.ncs_group"):
    """Write out NCS group formatted for phenix refine"""
    if not self._chain_residue_id or len(self._chain_residue_id)<2:
      return ""
    exclude=""
    if self._exclude_h: exclude+=" and (not element H) "
    if self._exclude_d: exclude+=" and (not element D) "
    [group,residue_range_list] = self._chain_residue_id
    count=0
    text = []
    for id,residue_ranges in zip (group,residue_range_list):
      count+=1
      if count==1:
        text.append("%s {"%prefix)
        l = "  reference = chain '{}'".format(id)
      else:
        l = "  selection = chain '{}'".format(id)
      if residue_ranges:
        first=True
        for residue_range in residue_ranges:
          if first:
            first=False
            l += " and (resseq "
          else:
            l += " or resseq  "
          l += str(residue_range[0])+":"+str(residue_range[1])
        l += " ) " + exclude
      text.append(l)
    text.append("}")
    text = '\n'.join(text)
    return text

  def format_for_biomt(self,crystal_number=None,skip_identity_if_first=False,
       ncs_domain_pdb=True):
    """Write out NCS group formatted for PDB BIOMT records"""

    serial_number=0
    from iotbx.mtrix_biomt import container

    result=container()
    for t,r in zip (
       self.translations_orth_inv(),self.rota_matrices_inv()):
      serial_number+=1
      result.add(r,t,serial_number, coordinates_present=False)

    return result.format_BIOMT_pdb_string()

  def format_for_resolve(self,crystal_number=None,skip_identity_if_first=False,
       ncs_domain_pdb=True):
    """Write out NCS group as text file for resolve"""
    text="new_ncs_group"
    if ncs_domain_pdb and self._ncs_domain_pdb is not None:
        text+="\nncs_domain_pdb "+str(self._ncs_domain_pdb)+"\n"
    i=0
    for center,trans_orth,ncs_rota_matr in zip (
       self._centers, self._translations_orth,self._rota_matrices):
      i+=1
      if i==1 and skip_identity_if_first and \
        is_identity(ncs_rota_matr,trans_orth): continue
      for j in range(3):
       text+="\nrota_matrix "+" %8.4f  %8.4f  %8.4f" %tuple(
          ncs_rota_matr[j*3:j*3+3])
      text+="\ntran_orth  "+" %8.4f  %8.4f  %8.4f" %tuple(trans_orth)
      text+="\n"
      text+="\ncenter_orth "+" %8.4f  %8.4f  %8.4f" %tuple(center)
      if crystal_number is not None:
        text+="\ncrystal_number "+str(crystal_number)
      text+="\n"
    return text

  def n_ncs_oper(self):
    """Return number of NCS operators"""
    return self._n_ncs_oper

  def chain_residue_id(self):
    """Return chain_residue_id dict"""
    return self._chain_residue_id

  def rmsd_list(self):
    """Return list of RMSD between NCS copies"""
    return self._rmsd_list

  def cc(self):
    """Return overall NCS CC if available"""
    return self._cc

  def note(self):
    """Return overall note if available"""
    return self._note

  def add_rmsd_list(self,rmsd_list):
    """Set the rmsd_list"""
    self._rmsd_list=rmsd_list

  def add_cc(self,cc):
    """Set the overall CC"""
    self._cc=cc

  def add_note(self,note):
    """Set the overall note"""
    self._note=note

  def residues_in_common_list(self):
    """Return the residues_in_common_list"""
    return self._residues_in_common_list

  def add_residues_in_common_list(self,residues_in_common_list):
    """Set the residues_in_common_list"""
    self._residues_in_common_list=residues_in_common_list

  def add_chain_residue_id(self,chain_residue_id):
    """Set the chain_residue_id"""
    self._chain_residue_id=chain_residue_id


  def centers(self):
    """Return the centers for the NCS groups"""
    return self._centers

  def translations_orth(self):
    """Return the orthogonal translations"""
    return self._translations_orth

  def rota_matrices(self):
    """Return the rotation matrices"""
    return self._rota_matrices

  def translations_orth_inv(self):
    """Return inverses of orthogonal translations"""
    if not hasattr(self,"_translations_orth_inv"):
      self.get_inverses()
    return self._translations_orth_inv

  def rota_matrices_inv(self):
    """Return inverses of the rotation matrices"""
    if not hasattr(self,"_rota_matrices_inv"):
      self.get_inverses()
    return self._rota_matrices_inv

  def delete_inv(self):
    """Remove inverses of rotation and translation matrices"""
    if hasattr(self,"_rota_matrices_inv"):
      del self._rota_matrices_inv
    if hasattr(self,"_translations_orth_inv"):
      del self._translations_orth_inv

  def adjust_magnification(self,magnification=None):
    """Set the magnification"""
    if not magnification or magnification==1:
      return # nothing to do

    self._translations_orth=self.copy_vector_list(self._translations_orth,
      scale_factor=magnification)
    self._centers=self.copy_vector_list(self._centers,
      scale_factor=magnification)

    self.get_inverses()

  def invert_matrices(self):
    """Not used"""
    self.get_inverses()
    # move the inverses to std
    self._translations_orth=deepcopy(self._translations_orth_inv)
    self._rota_matrices=deepcopy(self._rota_matrices_inv)
    self.get_inverses()

  def rotate_matrices(self,rot=None):
    """Rotate all matrices by rot"""
    translations_orth_rot=deepcopy(self._translations_orth)
    rota_matrices_rot=deepcopy(self._rota_matrices)

    self._translations_orth=[]
    self._rota_matrices=[]

    # create r_rot, t_rot such that r_rot x + t_rot is the same as
    #       rot * ( r [rot_inv x ] + t) : rotate x to orig, apply, r t, rotate back
    # r_rot(x)+t_rot ==  rot * ( r [rot_inv x ] + t)
    #  So:  t_rot=rot t
    #       r_rot=rot r rot_inv

    rot_inv=rot.inverse()
    for r,t in zip(rota_matrices_rot,translations_orth_rot):
      r_rot=rot *(r*rot_inv)
      t_rot=rot*t
      self._rota_matrices.append(r_rot)
      self._translations_orth.append(t_rot)
    self.get_inverses()

  def get_inverses(self):
    """Set up inverses of matrices and translations"""
    self._translations_orth_inv=[]
    self._rota_matrices_inv=[]
    for r,t in zip(self.rota_matrices(),self.translations_orth()):
      r_inv=r.inverse()
      t_inv=-1.*r_inv*t
      self._rota_matrices_inv.append(r_inv)
      self._translations_orth_inv.append(t_inv)

  def rotations_translations_forward_euler(self):
    """Get rotations_forward_euler,translations_forward_euler.
    Note usual rt is from molecule j to molecule 1. Here it is opposite."""
    from scitbx.math import euler_angles
    rotations_forward_euler=[]
    translations_forward_euler=[]
    for r,t in zip(self._rota_matrices,self._translations_orth):
      r_inv=r.inverse()
      t_inv=-1.*r_inv*t
      r_inv_euler=euler_angles.zyz_angles(r_inv)
      rotations_forward_euler.append(r_inv_euler)
      translations_forward_euler.append(t_inv)
    return rotations_forward_euler,translations_forward_euler

  def source_of_ncs_info(self):
    """Return the source of ncs information"""
    return self._source_of_ncs_info

  def ncs_domain_pdb(self):
    """Return name of file with NCS domains"""
    return self._ncs_domain_pdb

  def print_list(self,list_of_real):
    """Print list of real numbers"""
    text=""
    for number in list_of_real:
     text+="  "+str(self.round(number,2))
    return text

  def round(self,value,n_digit):  # round off value to n_digit digits
    """Round value to n_digit"""
    if type(value) == type(1):
       return self.round(float(value),n_digit)
    if type(value) != type(1.0):
       return self.round(0.0,1)

    if n_digit == 0:
      rounder=1
    else:
      rounder=10**n_digit
    if value >= 0:
      rounded=float(int(0.5+value*rounder))/rounder
    else:
      value1=-1.*value
      rounded=float(int(0.5+value1*rounder))/rounder
      rounded=-1.*rounded
    return rounded

  def extract_point_group_symmetry(self,
   tol_r=None,
   abs_tol_t=None,
   rel_tol_t=None):
   """Sequentially remove operators until pg symmetry is achieved or none
    are left"""
   ops_to_keep=[self.identity_op_id()]
   n_ops=len(self.rota_matrices_inv())
   for test_op in range(n_ops):
     if test_op in ops_to_keep: continue
     test_ops_to_keep=ops_to_keep+[test_op]
     new_group=self.deep_copy(
         ops_to_keep=test_ops_to_keep)
     if new_group.is_point_group_symmetry(
         tol_r=default_tol_r,
         abs_tol_t=default_abs_tol_t,
         rel_tol_t=default_rel_tol_t,
         symmetry_to_match=self):  # test_op keeps us in the original set
       ops_to_keep.append(test_op)
   new_group=self.deep_copy(
         ops_to_keep=ops_to_keep)
   if new_group.is_point_group_symmetry(
         tol_r=default_tol_r,
         abs_tol_t=default_abs_tol_t,
         rel_tol_t=default_rel_tol_t):
     return new_group
   else:
     return None

  def is_similar_ncs_group(self, other,
   tol_r=default_tol_r,
   abs_tol_t=default_abs_tol_t,
   rel_tol_t=default_rel_tol_t,
   allow_self_contained_in_other = True):
    '''
     Return True if all operations of self match one of other
    '''

    if (not allow_self_contained_in_other) and \
       len(self.rota_matrices_inv()) != len(other.rota_matrices_inv()):
      return False

    if len(self.rota_matrices_inv()) < 1:
      return False

    for r,t in zip(self.rota_matrices_inv(),self.translations_orth_inv()):
        is_similar=False
        for r2,t2 in zip(other.rota_matrices_inv(),
           other.translations_orth_inv()):
          if is_same_transform(r,t,r2,t2,tol_r=tol_r,
            abs_tol_t=abs_tol_t,rel_tol_t=rel_tol_t):
            is_similar=True
            break
        if not is_similar:
          return False
    return True

  def is_point_group_symmetry(self,
   tol_r=default_tol_r,
   abs_tol_t=default_abs_tol_t,
   rel_tol_t=default_rel_tol_t,
   symmetry_to_match=None):
    """Return True if any 2 sequential operations is a member of the
    set.  Test by sequentially applying all pairs of
    operators and verifying that the result is a member of the set"""

    # Allow checking self operators vs some other symmetry object if desired:

    if len(self.rota_matrices_inv()) < 2:
      return False

    if symmetry_to_match is None:
      symmetry_to_match=self

    for r,t in zip(self.rota_matrices_inv(),self.translations_orth_inv()):
      for r1,t1 in zip(self.rota_matrices_inv(),self.translations_orth_inv()):
        new_r = r1 * r
        new_t = (r1 * t) + t1
        is_similar=False
        for r2,t2 in zip(symmetry_to_match.rota_matrices_inv(),
           symmetry_to_match.translations_orth_inv()):
          if is_same_transform(new_r,new_t,r2,t2,tol_r=tol_r,
            abs_tol_t=abs_tol_t,rel_tol_t=rel_tol_t):
            is_similar=True
            break
        if not is_similar:
          return False
    self._have_point_group_symmetry=True
    self.tol_r=tol_r
    self.abs_tol_t=abs_tol_t
    self.rel_tol_t=rel_tol_t
    return True

  def sort_by_z_translation(self,tol_z=0.01,
       allow_negative_z_translation = False):
    """Sort rotation matrices and translations by z-translation"""
    n=len(self.rota_matrices_inv())
    z_translations=[]
    sort_z_translations=[]
    for i1 in range(n): # figure out if translation is along z
      z_translations.append(self.translations_orth_inv()[i1][2])
      sort_z_translations.append([self.translations_orth_inv()[i1][2],i1])
    rota_matrices_sav=deepcopy(self._rota_matrices)
    translations_orth_sav=deepcopy(self._translations_orth)

    # sort the z-translations to reorder the matrices. Could be backwards
    sort_z_translations = sorted(sort_z_translations, key = lambda s: s[0])
    #sort_z_translations.sort()
    sorted_indices=[]
    sorted_z=[]
    n_plus_one=0
    n_minus_one=0
    delta=None
    all_same_delta=True
    for z,i1 in sort_z_translations:
      if sorted_indices:
        if i1==sorted_indices[-1]+1: n_plus_one+=1
        if i1==sorted_indices[-1]-1: n_minus_one+=1
        delta_z=z-sorted_z[-1]
        if delta is None:
          delta=delta_z
        elif abs(delta-delta_z)>tol_z:
          is_helical=False # XX not used
      sorted_indices.append(i1)
      sorted_z.append(z)
    if allow_negative_z_translation and n_minus_one>n_plus_one:
      sorted_indices.reverse()
      self._helix_z_translation= -1 * delta
    else:
      self._helix_z_translation=delta

    # Reorder the operators:
    self._rota_matrices=len(rota_matrices_sav)*[None]
    self._translations_orth=len(rota_matrices_sav)*[None]
    for i1,i2 in zip(sorted_indices,range(n)):
      self._rota_matrices[i2]=rota_matrices_sav[i1]
      self._translations_orth[i2]=translations_orth_sav[i1]
    self.delete_inv() # remove the inv matrices/rotations so they regenerate
    if len(self._rota_matrices)<2:
      self._helix_theta=None
    else:
      self._helix_theta=self.get_theta_along_z(
        self._rota_matrices[0],self._rota_matrices[1])
    self.get_inverses()
    return sorted_indices

  def get_trans_along_z(self,t0,t1):
    """Get translation along z between t0 and t1"""
    return t1[2]-t0[2]

  def get_theta_along_z(self,m0,m1):
    """Get theta along z for m0, m1"""
    import math
    cost=m0[0]
    sint=m0[1]
    t0=180.*math.atan2(sint,cost)/3.14159
    cost=m1[0]
    sint=m1[1]
    t1=180.*math.atan2(sint,cost)/3.14159
    delta_rot = t1 - t0
    if delta_rot > 180:
      delta_rot = delta_rot - 360
    if delta_rot <= -180:
      delta_rot = delta_rot + 360
    return delta_rot

  def is_helical_along_z(self,tol_z=0.01,
   tol_r=default_tol_r,
   abs_tol_t=default_abs_tol_t,
   rel_tol_t=default_rel_tol_t):
    """Return True if this is helical symmetry along z.
    This assumes the operators are in order, but allow special case
    where the identity operator is placed at the beginning but belongs
    at the end.
    Also assumes the axis of helical symmetry is parallel to the Z-axis
    and returns False if not."""

    # For helical symmetry sequential application of operators moves up or
    #  down the list by an index depending on the indices of the operators.
    if len(self.rota_matrices_inv()) < 2:
      return False
    if self.is_point_group_symmetry(tol_r=tol_r,
            abs_tol_t=abs_tol_t,rel_tol_t=rel_tol_t):
      return False

    n=len(self.rota_matrices_inv())
    if n < 3:
      return False

    is_helical=True
    rota_matrices_sav=deepcopy(self._rota_matrices)
    translations_orth_sav=deepcopy(self._translations_orth)
    sorted_indices=self.sort_by_z_translation(tol_z=tol_z)

    offset_list=[]
    n_missing_list=[]
    self.helix_oper_forwards=None
    self.helix_oper_reverse=None
    self.helix_oper_identity=None
    for i1 in range(n): # figure out offset made by this self.helix_operator
      offset,n_missing=self.oper_adds_offset(i1,tol_r=tol_r,
          abs_tol_t=abs_tol_t,rel_tol_t=rel_tol_t)
      offset_list.append(offset)
      n_missing_list.append(n_missing)
      if offset==1:self.helix_oper_forwards=sorted_indices[i1]
      if offset==-1:self.helix_oper_reverse=sorted_indices[i1]
      if offset==0:self.helix_oper_identity=sorted_indices[i1]
    # offset_list should be one instance of each value and will be 0 at the
    #  operator that is unity
    if None in offset_list:
      is_helical=False

    if is_helical:
      ii=offset_list.index(0)
      if not is_identity(
          self.rota_matrices_inv()[ii],self.translations_orth_inv()[ii]):
        is_helical=False

    if is_helical:
      for i1 in range(n):
        if n_missing_list[i1] != abs(i1-ii):
          is_helical=False
          break
    if is_helical:
      offset_list.sort()
      start_value=offset_list[0]
      expected_list=list(range(start_value,start_value+n))
      if offset_list != expected_list:
        is_helical=False
    # restore
    sys.stdout.flush()
    self._rota_matrices=rota_matrices_sav
    self._translations_orth=translations_orth_sav
    self.delete_inv() # remove the inv matrices/rotations so they regenerate
    if is_helical:
      self._have_helical_symmetry=True
      self.tol_z=tol_z
      self.tol_r=tol_r
      self.abs_tol_t=abs_tol_t
      self.rel_tol_t=rel_tol_t
      return True
    else:
      return False

  def get_forwards_reverse_helix(self,r1=None,t1=None,r2=None,t2=None):
    """Get the forwards and reverse transforms, deciding which is which based
    on the order of operators supplied for r1 t1 and r2 t2"""
    assert self._have_helical_symmetry
    for dir in ['forwards','reverse']:
      if dir=='forwards':
        r_forwards,t_forwards=self.helix_rt_forwards()
        r_reverse,t_reverse=self.helix_rt_reverse()
      else:
        r_forwards,t_forwards=self.helix_rt_reverse()
        r_reverse,t_reverse=self.helix_rt_forwards()

      # apply to n-1 and see if we get n:
      new_r = r_forwards*r1
      new_t= (r_forwards* t1) + t_forwards
      if is_same_transform(new_r,new_t,r2,t2,
         tol_r=self.tol_r,
         abs_tol_t=self.abs_tol_t,
         rel_tol_t=self.rel_tol_t):
        return r_forwards,t_forwards,r_reverse,t_reverse
    from libtbx.utils import Sorry
    raise Sorry(
     "Unable to find forward and reverse operators for this helical symmetry")

  def get_helix_parameters(self,tol_z=default_tol_z):
    """Get helix parameters from z-translation and theta"""
    from libtbx import group_args
    helix_z_translation=self.get_helix_z_translation()
    helix_theta=self.get_helix_theta()
    if helix_z_translation is not None and helix_theta is not None:
      return group_args(helix_z_translation=helix_z_translation,
        helix_theta=helix_theta)

  def extend_helix_operators(self,z_range=None,tol_z=default_tol_z,
      max_operators=None):
    """Extend the operators to go from -z_range to z_range"""
    assert self._have_helical_symmetry
    rota_matrices_inv_sav=deepcopy(self.rota_matrices_inv())
    translations_orth_inv_sav=deepcopy(self.translations_orth_inv())
    # only apply centers if some existing ones are not zero
    have_non_zero_c=False
    from scitbx import matrix
    for c in self._centers:
      if list(c)!=[0,0,0]:
        have_non_zero_c=True
        break

    sort_list=[]
    for r,t,c in zip(rota_matrices_inv_sav,translations_orth_inv_sav,self._centers):
      z_value=t[2]
      sort_list.append([z_value,r,t,c])
    sort_list.sort(key=itemgetter(0))
    z_first,r_first,t_first,c_first=sort_list[0]
    z_last,r_last,t_last,c_last=sort_list[-1]
    z_next_to_last,r_next_to_last,t_next_to_last,c_next_to_last=sort_list[-2]
    z_translation=self.get_helix_z_translation()
    if not z_translation:
      from libtbx.utils import Sorry
      raise Sorry("Cannot apply extension of helical NCS operators with no"+
        " Z-translation")

    # Figure out which direction to add single shifts to each end
    r_forwards,t_forwards,r_reverse,t_reverse=self.get_forwards_reverse_helix(
      r1=r_next_to_last,t1=t_next_to_last,r2=r_last,t2=t_last)

    if max_operators:
      max_new_each_direction=max(0,
          (1+max_operators-len(rota_matrices_inv_sav))//2)
    else:
      max_new_each_direction=None

    # Add on at end until we get to z_max (or z_min if z_translation<0)
    r_list=[r_last]
    t_list=[t_last]
    c_list=[c_last]
    while 1:
      new_r = r_forwards*r_list[-1]
      new_t= (r_forwards * t_list[-1]) + t_forwards
      new_c= (r_forwards*c_list[-1])+t_forwards
      if max_new_each_direction and len(c_list)>=max_new_each_direction:
        break
      elif is_in_range(new_t[2],-z_range,z_range):
        r_list.append(new_r)
        t_list.append(new_t)
        c_list.append(new_c)
      else:
        break
    rota_matrices_inv=rota_matrices_inv_sav+r_list[1:]
    translations_orth_inv=translations_orth_inv_sav+t_list[1:]
    centers=self._centers+c_list[1:]

    # and for other end
    r_list=[r_first]
    t_list=[t_first]
    c_list=[c_first]
    while 1:
      new_r = r_reverse*r_list[-1]
      new_t= (r_reverse * t_list[-1]) + t_reverse
      new_c= (r_reverse*c_list[-1])+t_reverse
      if max_new_each_direction and len(c_list)>=max_new_each_direction:
        break
      elif is_in_range(new_t[2],-z_range,z_range):
        r_list.append(new_r)
        t_list.append(new_t)
        c_list.append(new_c)
      else:
        break
    rota_matrices_inv+=r_list[1:]
    translations_orth_inv+=t_list[1:]
    centers+=c_list[1:]
    # Now we have a new set...invert and save them
    self._n_ncs_oper=len(rota_matrices_inv)
    self._rota_matrices=[]
    self._translations_orth=[]
    self._centers=[]
    for r_inv,t_inv,c in zip(rota_matrices_inv,translations_orth_inv,centers):
      r_std=r_inv.inverse()
      t_std=-1.*r_std*t_inv
      self._rota_matrices.append(r_std)
      self._translations_orth.append(t_std)
      if have_non_zero_c:
        self._centers.append(c)
      else:
        self._centers.append(matrix.col((0,0,0)))

    if self._rmsd_list:
      self._rmsd_list=len(rota_matrices_inv)*[None]
    if self._residues_in_common_list:
       self._residues_in_common_list=len(rota_matrices_inv)*[None]

    if self._chain_residue_id:
       self._chain_residue_id= [
        len(rota_matrices_inv)*[None],
        len(rota_matrices_inv)*[[]]
         ]

    self.delete_inv()
    self.get_inverses()

    # And reorder them now...
    self.sort_by_z_translation(tol_z=tol_z)

  def get_helix_z_translation(self):
    """Return the helix z translation"""
    assert self._have_helical_symmetry
    if hasattr(self,'_helix_z_translation'):
      return self._helix_z_translation
    return None

  def get_helix_theta(self):
    """Return the helix theta"""
    assert self._have_helical_symmetry
    if hasattr(self,'_helix_theta'):
      return self._helix_theta
    return None

  def helix_rt_reverse(self):
    """Return r and t for moving one reverse in a helix"""
    assert self._have_helical_symmetry
    if not hasattr(self,'helix_oper_reverse') or not self.helix_oper_reverse:
      if not hasattr(self,'helix_oper_forwards') or \
         not self.helix_oper_forwards: # no info, quit
        return None,None
      else: # generate from forwards:
        r_forwards,t_forwards=self.helix_rt_forwards()
        r_reverse=r_forwards.inverse()
        t_reverse=-1.*r_reverse*t_forwards
        return r_reverse,t_reverse

    i1=self.helix_oper_reverse
    r1=self.rota_matrices_inv()[i1]
    t1=self.translations_orth_inv()[i1]
    return r1,t1

  def helix_rt_forwards(self):
    """Return r and t for moving one forwards in a helix"""
    assert self._have_helical_symmetry
    if not hasattr(self,'helix_oper_forwards') or not self.helix_oper_forwards:
      if not hasattr(self,'helix_oper_reverse') or \
         not self.helix_oper_reverse: # no info, quit
        return None,None
      else: # generate from reverse:
        r_reverse,t_reverse=self.helix_rt_reverse()
        r_forwards=r_reverse.inverse()
        t_forwards=-1.*r_forwards*t_reverse
        return r_forwards,t_forwards

    i1=self.helix_oper_forwards
    r1=self.rota_matrices_inv()[i1]
    t1=self.translations_orth_inv()[i1]
    return r1,t1

  def oper_adds_offset(self,i1,tol_r=None,abs_tol_t=None,rel_tol_t=None):
    """Figure out what operator is created from operator i1 + any other one"""
    n=len(self.rota_matrices_inv())
    r1=self.rota_matrices_inv()[i1]
    t1=self.translations_orth_inv()[i1]
    offset_list=[]
    missing_list=[]
    for i in range(n): # apply i1+i and see what k we get (k-i should be const)
      r=self.rota_matrices_inv()[i]
      t=self.translations_orth_inv()[i]
      new_r = r1 * r
      new_t = (r1 * t) + t1
      match_offset=None
      for j in range(n):
        r2=self.rota_matrices_inv()[j]
        t2=self.translations_orth_inv()[j]
        if is_same_transform(new_r,new_t,r2,t2,
           tol_r=tol_r,
           abs_tol_t=abs_tol_t,
           rel_tol_t=rel_tol_t):
          match_offset=j-i
      if match_offset is not None:
        if not match_offset in offset_list: offset_list.append(match_offset)
      else:
        missing_list.append(i)
    if len(offset_list)==1:
      return offset_list[0],len(missing_list)
    else:
      return None,None



  def identity_op_id(self):
    """Return id of identity operator"""
    id=0
    for center,trans_orth,ncs_rota_matr in zip (
       self._centers, self._translations_orth,self._rota_matrices):
      if is_identity(ncs_rota_matr,trans_orth):
         return id
      id+=1
    return None

  def map_inside_unit_cell(self,unit_cell=None):
    """Map all the operators inside the unit cell.  Must be supplied"""
    assert unit_cell is not None
    if len(self._centers)==0: return

    new_centers=[]
    new_translations_orth=[]
    from scitbx.math import  matrix
    #rotation matrices do not change, just translations
    # find the identity:
    first_coordinate_offset=None
    for center,trans_orth,ncs_rota_matr in zip (
       self._centers, self._translations_orth,self._rota_matrices):
      if is_identity(ncs_rota_matr,trans_orth):
        first_coordinate_offset=matrix.col(offset_inside_cell(
          center,unit_cell=unit_cell))
        break
    if first_coordinate_offset is None:
      raise Sorry("Identity not found in NCS matrices?")
    for center,trans_orth,ncs_rota_matr in zip (
       self._centers, self._translations_orth,self._rota_matrices):
      coordinate_offset=offset_inside_cell(center,unit_cell=unit_cell)

      new_centers.append(matrix.col(center)+matrix.col(coordinate_offset))
      #  T'=T - R x_i + x_1
      delta = matrix.col(ncs_rota_matr * matrix.col(coordinate_offset))
      t_prime=matrix.col(trans_orth) - delta + first_coordinate_offset
      new_translations_orth.append(t_prime)

    self._centers=new_centers
    self._translations_orth=new_translations_orth

  def add_identity_op(self):
    """Add identity op to group if not present"""
    if self.identity_op_id() is not None:
       return # nothing to do

    from scitbx import matrix

    self._rota_matrices.append(matrix.sqr(
      [1.,0.,0.]+[0.,1.,0.]+[0.,0.,1.]))
    self._translations_orth.append(matrix.col([0.,0.,0.]))
    self._centers.append(matrix.col([0.,0.,0.]))
    self._n_ncs_oper+=1

class ncs:
  """Class to hold any number of local symmetry (NCS) groups"""
  def __init__(self,exclude_h=None,exclude_d=None):
    self._ncs_groups=[]  # each group is an ncs_group object
    self.source_info=None
    self._ncs_read=False
    self._exclude_h=exclude_h
    self._exclude_d=exclude_d
    self._ncs_obj = None
    self._ncs_name = None # an optional name like "D3"
    self._shift_cart = (0,0,0)  # shift to place object in original location

  def __repr__(self):
    """Text representation of an NCS object containing some number of NCS groups"""
    text = "NCS object with %s groups: " %(len(self._ncs_groups))
    for g in self._ncs_groups:
      text+=str(g)
    return text

  def deep_copy(self,change_of_basis_operator=None,unit_cell=None,
      coordinate_offset=None,
      scale_factor=None,
      new_unit_cell=None,
      ops_to_keep=None,
      extract_point_group_symmetry=None,
      hierarchy_to_match_order=None):  # make a copy
    """Make new ncs object with same overall params as this one:"""
    from mmtbx.ncs.ncs import ncs

    new=ncs(exclude_h=self._exclude_h,exclude_d=self._exclude_d)
    new.source_info=self.source_info
    new._ncs_name=self._ncs_name
    new._ncs_read=self._ncs_read
    new._shift_cart=deepcopy(self._shift_cart)  # shift_cart is what was applied
    if coordinate_offset:
       new._shift_cart=tuple(
         [a+b for a,b in zip(new._shift_cart,coordinate_offset)])
    # deep_copy over all the ncs groups:
    for ncs_group in self._ncs_groups:
      new_group=ncs_group.deep_copy(
         change_of_basis_operator=change_of_basis_operator,
         coordinate_offset=coordinate_offset,
         scale_factor=scale_factor,
         unit_cell=unit_cell,new_unit_cell=new_unit_cell,
         ops_to_keep=ops_to_keep,
         extract_point_group_symmetry=extract_point_group_symmetry,
         hierarchy_to_match_order=hierarchy_to_match_order)

      new._ncs_groups.append(new_group)
    return new

  def change_of_basis(self,change_of_basis_operator=None,unit_cell=None,
      new_unit_cell=None):
    """Apply change of basis to all groups"""
    if change_of_basis_operator is None or unit_cell is None or\
        new_unit_cell is None:
       raise Sorry("For change of basis unit_cell, "+
           "new_unit_cell and operator are all required")
    return self.deep_copy(change_of_basis_operator=change_of_basis_operator,
      unit_cell=unit_cell,new_unit_cell=new_unit_cell)

  def magnification(self,scale_factor=None):
    """Apply magnification to all groups"""
    if scale_factor is None:
       raise Sorry("For magnification a scale factor is required.")
    return self.deep_copy(scale_factor=scale_factor)

  def set_shift_cart(self,shift_cart):
    """Set the shift_cart (origin offset)"""
    self._shift_cart=shift_cart

  def shifted(self, eps=1.e-3):
    ''' Return True if this ncs_object has been shifted from its original
     location (e.g., by boxing a map and model and this ncs object).
    '''

    r = self.shift_cart()
    if(r is None): return False
    from scitbx.array_family import flex
    if(flex.max(flex.abs(flex.double(r)))<=eps): return False
    return True

  def shift_cart(self):
    """Return the shirt_cart (origin offset)"""
    return self._shift_cart

  def coordinate_offset(self,coordinate_offset=None,unit_cell=None,
      new_unit_cell=None):
    """Offset coordinates by coordinate_offset.
    NOTE: Returns new object, self is unchanged."""
    if coordinate_offset is None:
       raise Sorry("For coordinate_offset an offset is required.")
    return self.deep_copy(coordinate_offset=coordinate_offset)

  def map_inside_unit_cell(self,unit_cell=None):
    """Map all the operators inside the unit cell.  Must be supplied and the
    centers for the operators must exist and not be zero"""
    for ncs_group in self._ncs_groups:
      ncs_group.map_inside_unit_cell(unit_cell=unit_cell)

  def ncs_read(self):
    """Return True if NCS has been read in"""
    return self._ncs_read

  def ncs_groups(self):
    """Return the NCS groups in this NCS object"""
    return self._ncs_groups

  def identity_op_id_in_first_group(self):
    """Return True if identity operator is in first (usually main) NCS group"""
    if self._ncs_groups:
      return self._ncs_groups[0].identity_op_id()
    else:
      return None

  def ncs_oper_in_first_group(self):
    """Return number of copies in first (usually main) NCS group"""
    if self._ncs_groups:
      return self._ncs_groups[0].n_ncs_oper()
    else:
      return None

  def rotate_about_z(self,rot_deg=None,invert_matrices=True):
    """Rotate all the ops by rot_deg about z"""
    if invert_matrices:
      rot_deg=-rot_deg
    oper=get_rot_z(rot_deg=rot_deg)
    self.rotate_matrices(rot=oper)

  def rotate_about_y(self,rot_deg=None,invert_matrices=True):
    """Rotate all the ops by rot_deg about y"""
    if invert_matrices:
      rot_deg=-rot_deg
    oper=get_rot_y(rot_deg=rot_deg)
    self.rotate_matrices(rot=oper)

  def ncs_from_pdb_input_BIOMT(self,pdb_inp=None,log=None,quiet=False,
     invert_matrices=True):
    """Obtain NCS object from BIOMT records in PDB file"""
    p=pdb_inp.process_BIOMT_records()
    if not p:
      print("No BIOMT records available", file=log)
      return

    self.ncs_from_import(rot_list=p.r,trans_list=p.t,invert_matrices=invert_matrices)

  def ncs_from_import(self,rot_list=None,trans_list=None,invert_matrices=True):
    """Obtain NCS object from a list of rotations, translations"""

    self.init_ncs_group()

    for r,t in zip(rot_list,trans_list):
      self._rota_matrix=r.as_list_of_lists()
      self._trans=list(t)
      self._center=[0.,0.,0.]
      self.save_oper()

    self.save_existing_group_info()
    # Invert them (we use mapping from operator i to 1 ; biomt is 1 to i
    if invert_matrices:
      self.invert_matrices()

    self._ncs_read=True


  def select_first_ncs_group(self):
    """Just keep the first ncs group and remove others"""
    self._ncs_groups=self._ncs_groups[:1]
    return self

  def select_first_ncs_operator(self):
    """Just keep the first ncs operator in the first group and remove others"""
    self.select_first_ncs_group()
    if self._ncs_groups:
      self._ncs_groups=[self._ncs_groups[0].deep_copy(ops_to_keep=[0])] #  keep first only
    return self

  def set_unit_ncs(self):
    """Just make a single ncs operator"""

    self.init_ncs_group()

    self._rota_matrix=[[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]
    self._trans=[0.,0.,0.]
    self._center=[0.,0.,0.]
    self.save_oper()

    self.save_existing_group_info()
    self._ncs_read=True

  def read_ncs(self,file_name=None,lines=[],source_info="",
       log=None,quiet=False):
    """Read NCS from a file"""
    if not log: log=sys.stdout
    if not quiet:
      if file_name:
        print("Reading NCS information from: ",file_name, file=log)
    if source_info:
       print(" based on ",source_info, file=log)
       self.source_info=source_info
    else:
       self.source_info=str(file_name)
    if file_name:
      if not os.path.isfile(file_name):
        raise Sorry("The file "+str(file_name)+" does not seem to exist?")
      else:
        with open(file_name) as f:
          lines=f.readlines()
    self.init_ncs_group()

    read_something=False

    for line in lines:
      if not line : continue
      spl=line.split()
      if len(spl)<1: continue
      key=spl[0].lower()
      if key=='transformations' and len(spl)>1 and \
         spl[1].lower()=='formatted':  # start all over!
         self._ncs_groups=[]
         self.init_ncs_group()
      elif key=='new_ncs_group': # start new group
        self.save_existing_group_info()
      # NOTE: new operator signified by rota_matrix or new_operator
      elif key=='new_operator':
        self.save_oper()
      elif key=='rota_matrix': # read set of rota
        if self._rota_matrix and \
            len(self._rota_matrix)==3 or len(self._rota_matrix)==0:
          self.save_oper()
        set=self.get_3_values_after_key(line)
        self._rota_matrix.append(set)

      elif key=='tran_orth': # read translation
        self._trans=self.get_3_values_after_key(line)
      elif key=='center_orth': # read translation
        self._center=self.get_3_values_after_key(line)
      elif key=='ncs_cc': # read  cc
        self._cc=self.get_1_value_after_key(line)
      elif key=='note' or key=='note:': # read anything
        self._note=" ".join(line.split()[1:])
      elif key=='chain':
        self._chain=self.get_1_char_after_key(line)
      elif key=='resseq':
        self._resseq_list.append(self.get_res_range_after_key(line))
      elif key=='rmsd':
        self._rmsd=self.get_1_value_after_key(line)
      elif key=='matching':
        self._residues_in_common=self.get_1_value_after_key(line)
      elif key=='source_info':
        self.source_info=self.get_1_char_after_key(line)
      elif key=='ncs_domain_pdb':
        self._ncs_domain_pdb=self.get_1_char_after_key(line)
      elif len(spl)==3 and spl[0]=='No' and spl[1]=='NCS' and spl[2]=='found':
        read_something=True
      else:
        pass
    self.save_existing_group_info()
    if read_something or len(self._ncs_groups) > 0:
      self._ncs_read=True
    else: # Try as biomtr
      import iotbx.pdb
      try:
        pdb_inp = iotbx.pdb.input(lines=lines,source_info=file_name)
        self.ncs_from_pdb_input_BIOMT(pdb_inp=pdb_inp,log=log,quiet=quiet)
      except Exception as e:
        pass

  def save_existing_group_info(self):
        """Save current group information"""

        self.save_oper()
        if self._n_ncs_oper > 0:  # save last-read ncs group.
          self.save_ncs_group()


  def get_res_range_after_key(self,line):
    """Get residue range after ':' in text string"""
    spl = line.replace(':', ' ').split()
    if  len(spl)<3:
      raise Sorry("Cannot interpret this NCS file"+"\n"+str(line))
    start,end=None,None
    try:
      start=int(spl[1])
      end=int(spl[2])
    except Exception:
      raise Sorry("Cannot interpret this NCS file"+"\n"+str(line))
    return [start,end]

  def get_1_char_after_key(self,line):
    """Get 1 (or more) characters in second word in line"""
    spl=line.split()
    if  len(spl)<2:
      raise Sorry("Cannot interpret this NCS file"+"\n"+str(line))
    char=None
    try:
      char=spl[1]
    except Exception:
      raise Sorry("Cannot interpret this NCS file"+"\n"+str(line))
    return char

  def get_1_value_after_key(self,line):
    """Get one value in second word in line"""
    spl=line.split()
    if  len(spl)<2:
      raise Sorry("Cannot interpret this NCS file"+"\n"+str(line))
    cc=None
    try:
      cc=float(spl[1])
    except Exception:
      raise Sorry("Cannot interpret this NCS file"+"\n"+str(line))
    return cc

  def get_3_values_after_key(self,line):
    """Get 3 values starting with 2nd word in line"""
    spl=line.split()
    if  len(spl)<4:
      raise Sorry("Cannot interpret this NCS file"+"\n"+str(line))
    set=[]
    try:
      for item in spl[1:4]:
        set.append(float(item))
    except Exception:
      raise Sorry("Cannot interpret this NCS file"+"\n"+str(line))
    return set

  def init_ncs_group(self):
     """Set up an NCS group"""
     self._n_ncs_oper=0
     self._ncs_trans_orth=[]
     self._ncs_rota_matr=[]
     self._ncs_center_orth=[]
     self.init_oper()
     self._rmsd_list=[]
     self._residues_in_common_list=[]
     self._cc=None
     self._note=None
     self._ncs_domain_pdb=None
     self._chain_residue_id=[]

     self._list_of_resseq_list=[]
     self._group=[]

  def init_oper(self):
     """Initialize operators"""
     self._rota_matrix=[]
     self._trans=None
     self._center=None
     self._rmsd=None
     self._residues_in_common=None
     self._resseq_list=[]
     self._chain=None

  def save_oper(self):
     """Save operators"""
     # decide if there is anything to save:

     have_oper=True
     for item in (self._trans,self._rota_matrix,self._center):
       if not item:
          have_oper=False
     if self._rota_matrix and len(self._rota_matrix)!=3:
          raise Sorry("Cannot interpret this NCS file (rotations not understood)")
     if self._trans and len(self._trans)!=3:
          raise Sorry("Cannot interpret this NCS file (translations not understood)")
     have_something=False
     if have_oper or self._rmsd or self._residues_in_common:
       have_something=True
     if not have_something: return
     self._n_ncs_oper+=1
     if have_oper:
       from scitbx import matrix
       self._ncs_trans_orth.append(matrix.col(self._trans))
       self._ncs_rota_matr.append(matrix.sqr(
         self._rota_matrix[0]+self._rota_matrix[1]+self._rota_matrix[2] ))
       self._ncs_center_orth.append(matrix.col(self._center))
     else:
       self._ncs_trans_orth.append(None)
       self._ncs_rota_matr.append(None)
       self._ncs_center_orth.append(None)
     self._rmsd_list.append(self._rmsd)
     self._residues_in_common_list.append(self._residues_in_common)
     self._list_of_resseq_list.append(self._resseq_list)
     self._group.append(self._chain)

     self.init_oper()

  def import_ncs_group(self,ncs_rota_matr=None,
       center_orth=None,
       trans_orth=None,
       chain_residue_id=None,
       residues_in_common_list=None,
       rmsd_list=None,
       ncs_domain_pdb=None,
       cc=None,
       source_of_ncs_info=None,
       ncs_group_object=None):

     """Import an NCS group"""
     if not ncs_group_object:
       list_length=None
       if center_orth is None and trans_orth:
         center_orth=len(trans_orth)*[(0,0,0)]
       for lst in [trans_orth,ncs_rota_matr,center_orth]:
         if not lst or len(lst)<1:
           print("Length too short:",type(lst),lst, end=' ')
           if lst is not None:
             print(len(lst))
           else:
             print("0")
           raise Sorry("The NCS operators in this file appear incomplete?")
         if not list_length: list_length=len(lst)
         if list_length!=len(lst):
           print("Length of list incorrect:",type(lst),lst,len(lst),list_length)
           raise Sorry("The NCS operators in this file appear incomplete?")
       ncs_group_object=ncs_group(
         ncs_rota_matr=ncs_rota_matr,
         center_orth=center_orth,
         trans_orth=trans_orth,
         chain_residue_id=remove_quotes_from_chain_id(chain_residue_id),
         residues_in_common_list=residues_in_common_list,
         rmsd_list=rmsd_list,
         source_of_ncs_info=source_of_ncs_info,
         ncs_domain_pdb=ncs_domain_pdb,
         cc=cc,
         exclude_h=self._exclude_h,exclude_d=self._exclude_d)
     self._ncs_groups.append(ncs_group_object)

  def save_ncs_group(self):
     """Save an NCS group"""
     # check that there is something  here:
     have_something=False
     for lst in [self._ncs_trans_orth,
         self._ncs_rota_matr,self._ncs_center_orth,
         self._residues_in_common_list,self._rmsd_list]:
        if lst and self._n_ncs_oper and \
           len(lst) != self._n_ncs_oper:
          print("Lengh of list does not match number of operators:",\
             type(lst),lst,self._n_ncs_oper)
          raise Sorry("The NCS operators in this file appear incomplete?")
        if lst is not None and len(lst)<1:
          print("Length of operators too short:",lst,len(lst))
          raise Sorry("The NCS operators in this file appear incomplete?")
        if lst is not None: have_something=True
     if not have_something: return
     self._chain_residue_id=[self._group,self._list_of_resseq_list]
     ncs_group_object=ncs_group(
       ncs_rota_matr=self._ncs_rota_matr,
       center_orth=self._ncs_center_orth,
       trans_orth=self._ncs_trans_orth,
       source_of_ncs_info=self.source_info,
       ncs_domain_pdb=self._ncs_domain_pdb, # 041309
       rmsd_list=self._rmsd_list,
       residues_in_common_list=self._residues_in_common_list,
       chain_residue_id=self._chain_residue_id,
       cc=self._cc,note=self._note)
     self._ncs_groups.append(ncs_group_object)
     self.init_ncs_group()

  def show_summary(self, verbose=True, log = None):
    """Summarize NCS groups"""
    return self.display_all(verbose=verbose, log = log)

  def display_all(self,verbose=True,log=None):
    """Long summary"""
    if log==None:
      log=sys.stdout
    count=0
    text=""
    if self._ncs_name:
      text+="NCS TYPE: %s " %(self._ncs_name)

    for ncs_group in self._ncs_groups:
      count+=1
      text+="\n\nGROUP "+str(count)
      text+=ncs_group.display_summary(verbose=verbose)
    text+="\n\n"
    log.write(text)
    return text

  def shift_cart(self):
    """Return the shift_cart (offset) value"""
    if self._shift_cart:
      return self._shift_cart
    else:
      return (0,0,0)

  def shift_back_cart(self):
    """Return inverse of shift_cart"""
    return tuple([-a for a in self.shift_cart()])

  def as_ncs_spec_string(self, format = 'ncs_spec'):
    '''
     Shifts to original location and returns text string
    '''
    assert format in ['ncs_spec','phil']
    shifted_ncs=self.coordinate_offset(coordinate_offset=self.shift_back_cart())
    if format == 'ncs_spec':
      return shifted_ncs.format_all_for_group_specification(
         log=null_out(),quiet=True,out=null_out())
    elif format == 'phil':
      return shifted_ncs.format_all_for_phenix_refine(
         quiet=True,out=null_out())


  def format_all_for_group_specification(self,log=None,quiet=True,out=None,
       file_name=None):
    """Return ncs_spec format of entire NCS object"""
    if file_name is not None:
       out=open(file_name,'w')
    if out==None:
       out=sys.stdout
    if log==None:
      log=sys.stdout
    elif hasattr(out,'name'):
      print("NCS written as ncs object information to:"\
        ,out.name, file=log)
    all_text=""
    text="Summary of NCS information\n"
    import time,os
    text+=time.ctime()+"\n"
    text+=os.getcwd()+"\n\n"
    if self.source_info is not None:
      text+="source_info "+str(self.source_info)+"\n"
    if not self._ncs_groups:
      text+="No NCS found\n"
    if out is not None or not quiet: out.write("\n"+text+"\n\n")
    for ncs_group in self._ncs_groups:
      text=ncs_group.format_group_specification()
      if out is not None or not quiet: out.write("\n"+text+"\n\n")
      all_text+="\n"+text
    all_text+="\n"
    if out is not None and file_name is not None:
      out.close()
    return all_text

  def format_all_for_biomt(self,log=None,quiet=False,out=None,):
    """Return BIOMT records for NCS operators in first NCS group"""
    if out==None:
       out=sys.stdout
    if log==None:
      log=sys.stdout
    else:
      print("\n\nNCS operators written in BIOMT format :",out.name, file=log)
    all_text=""
    if self._ncs_groups and len(self._ncs_groups)>1:
      print(
      "\nWARNING: BIOMT format cannot be used for more than one NCS group",
       "\nOnly writing out one NCS group",file=log)
    for ncs_group in self._ncs_groups[:1]:
      text=ncs_group.format_for_biomt()
      if not quiet: out.write("\n"+text+"\n\n")
      all_text+="\n"+text
    return all_text

  def format_all_for_resolve(self,log=None,quiet=False,out=None,
      crystal_number=None,skip_identity_if_first=False,ncs_domain_pdb=True):
    """Format NCS object for resolve"""
    if out==None:
       out=sys.stdout
    if log==None:
      log=sys.stdout
    else:
      print("\n\nNCS operators written in format for resolve to:",out.name, file=log)
    all_text=""
    for ncs_group in self._ncs_groups:
      text=ncs_group.format_for_resolve(crystal_number=crystal_number,
         skip_identity_if_first=skip_identity_if_first,
         ncs_domain_pdb=ncs_domain_pdb)
      if not quiet: out.write("\n"+text+"\n\n")
      all_text+="\n"+text
    return all_text

  def format_all_for_phenix_refine(self,quiet=False,out=None,
        prefix="refinement.pdb_interpretation.ncs_group"):
    '''
    This function is an older version of creating phil for phenix refine,
    it is modified to replicate a new phil parameters that can handle
    selection to the level of atoms, "format_phil_for_phenix_refine".

    When it will still can be used in the older form, which allows only
    residue level selection.
    '''
    if hasattr(self._ncs_obj,'show'):
      if prefix == 'refinement.pdb_interpretation.ncs_group':
        prefix="pdb_interpretation"
      if quiet:
        out = null_out()
      elif out is None:
        out=sys.stdout
      all_text = self._ncs_obj.show(format='phil',log=null_out(),header=False)
      # all_text = convert_phil_format(all_text,to_type=prefix)
      if all_text:
        if not quiet:
          print(all_text + '\n', file=out)
      return all_text
    else:
      # this is only being used when only a spec file is provided
      if out == None:
        out=sys.stdout
      all_text=""
      for ncs_group in self._ncs_groups:
        text= ncs_group.format_for_phenix_refine(prefix=prefix)
        if text:
          if not quiet: out.write(text+'\n')
          all_text+="\n"+text
      return all_text

  def format_phil_for_phenix_refine(self,log=None,quiet=False,out=None):
    """ Writes NCS phil selection in phenix_refine format """
    if out==None: out=sys.stdout
    phil_str = self._ncs_obj.show(format='phil',log=null_out(),header=False)
    # ncs_str = convert_phil_format(phil_str,to_type="ncs")
    ncs_str = phil_str
    if ncs_str:
      if not quiet: out.write(ncs_str + '\n')
    return ncs_str

  def format_phil_for_ncs(self,log=None,quiet=False,out=None):
    """ Writes NCS phil selection in NCS format """
    if out==None: out=sys.stdout
    if log==None:
      log=sys.stdout
    else:
      msg  = "NCS phil selection written in ncs selection format to:"
      print(msg,out.name, file=log)
    phil_str = self._ncs_obj.show(format='phil',log=null_out(),header=False)
    if phil_str:
      if not quiet: out.write(phil_str + '\n')
    return phil_str

  def set_ncs_name(self,ncs_name):
    """Set the ncs name"""
    self._ncs_name=ncs_name

  def get_ncs_name(self):
    """Return the ncs name"""
    return self._ncs_name

  def add_source_info(self,source_info):
    """Add source information"""
    if self.source_info is None:
       self.source_info=str(source_info)
    else:
       self.source_info+=str(source_info)

  def add_cc_list(self,cc_list):
    """Add list of CC values for each NCS group"""
    if len(self._ncs_groups) != len(cc_list):
      raise Sorry("Number of NCS groups does not match length of cc_list...")
    for ncs_group,cc in zip(self._ncs_groups,cc_list):
      ncs_group.add_cc(cc)

  def overall_note(self):
    """Add an overall note"""
    overall_note=""
    for ncs_group in self._ncs_groups:
      if ncs_group._note is not None:
        overall_note+=" "+ncs_group._note
    return overall_note

  def overall_cc(self):
    """Calculate overall cc from cc_list"""
    cc_all=0.
    n=0
    for ncs_group in self._ncs_groups:
      if ncs_group._cc is not None:
        cc_all+=ncs_group._cc
        n+=1
    if n>0:
      cc_all=cc_all/float(n)
    else:
      cc_all=None
    return cc_all

  def overall_rmsd(self):
    """Calculate overall rmsd from rmsd_list"""
    rmsd_all=0.
    n=0
    for ncs_group in self._ncs_groups:
      if ncs_group.rmsd_list() is not None:
        for rmsd in ncs_group.rmsd_list():
          if rmsd:
            rmsd_all+=rmsd
            n+=1
    if n>0:
      rmsd_all=rmsd_all/float(n)
    else:
      rmsd_all=None
    return rmsd_all

  def max_operators(self):
    """Return number of operators in NCS group with the most operators"""
    n_max=0
    for ncs_group in self._ncs_groups:
      if ncs_group and ncs_group.n_ncs_oper()>n_max:
        n_max=ncs_group.n_ncs_oper()
    return n_max

  def is_similar_ncs_object(self, other,
   tol_r=default_tol_r,
   abs_tol_t=default_abs_tol_t,
   rel_tol_t=default_rel_tol_t,
   allow_self_contained_in_other = True):
    '''
      Determine if self and other are similar ncs objects.
      ncs groups do not have to be in same order
    '''
    if not self._ncs_groups and not other._ncs_groups:
      return True  # nothing there for either one

    if self.shift_cart() != other.shift_cart():
      return False
    if not self._ncs_groups:
      return False
    if not other._ncs_groups:
      return False
    if len(self._ncs_groups) != len (other._ncs_groups):
      return False
    for ncs_group in self._ncs_groups:
      found=False
      for other_ncs_group in other._ncs_groups:
        if ncs_group.is_similar_ncs_group(other_ncs_group,tol_r=tol_r,
            abs_tol_t=abs_tol_t,rel_tol_t=rel_tol_t,
            allow_self_contained_in_other = allow_self_contained_in_other):
          found=True
          break
      if not found:
        return False
    return True

  def is_point_group_symmetry(self,
   tol_r=default_tol_r,
   abs_tol_t=default_abs_tol_t,
   rel_tol_t=default_rel_tol_t):
    """Return True if this NCS group has point symmetry"""
    if not self._ncs_groups:
      return False
    for ncs_group in self._ncs_groups[:1]:
      if not ncs_group.is_point_group_symmetry(tol_r=tol_r,
            abs_tol_t=abs_tol_t,rel_tol_t=rel_tol_t):
        return False
    return True

  def adjust_magnification(self,magnification=None):
    """Set magnification of each NCS group"""
    if not self._ncs_groups:
      return self
    for ncs_group in self._ncs_groups:
      ncs_group.adjust_magnification(magnification=magnification)
    return self

  def rotate_matrices(self,rot=None):
    """Rotate all matrices in all groups by rot"""
    if not self._ncs_groups:
      return
    for ncs_group in self._ncs_groups:
      ncs_group.rotate_matrices(rot=rot)

  def invert_matrices(self):
    """Invert matrices in all groups"""
    if not self._ncs_groups:
      return
    for ncs_group in self._ncs_groups:
      ncs_group.invert_matrices()

  def sort_by_z_translation(self,tol_z=default_tol_z):
    """Sort by z translation in all groups"""
    if not self._ncs_groups:
      return
    for ncs_group in self._ncs_groups:
      ncs_group.sort_by_z_translation(tol_z=tol_z)

  def extend_helix_operators(self,z_range=None,tol_z=default_tol_z,
      max_operators=None):
    """Extend helix operators in all groups"""
    if not self._ncs_groups:
      return
    for ncs_group in self._ncs_groups:
      ncs_group.extend_helix_operators(z_range=z_range,tol_z=tol_z,
        max_operators=max_operators)

  def get_helix_parameters(self,z_range=None,tol_z=default_tol_z,
      max_operators=None):
    """Get helix parameters in all groups"""
    if not self._ncs_groups:
      return
    # return values for group 0
    return self._ncs_groups[0].get_helix_parameters(tol_z=tol_z,)

  def is_helical_along_z(self,
   tol_r=default_tol_r,
   abs_tol_t=default_abs_tol_t,
   rel_tol_t=default_rel_tol_t):
    """Return True if all groups are helical symmetry along z"""
    if not self._ncs_groups:
      return False
    for ncs_group in self._ncs_groups:
      if not ncs_group.is_helical_along_z(tol_r=tol_r,
            abs_tol_t=abs_tol_t,rel_tol_t=rel_tol_t):
        return False
    return True

  def add_identity_op(self):
    """Add identity operator to all groups"""
    for ncs_group in self._ncs_groups:
      if ncs_group.identity_op_id() is None:
        ncs_group.add_identity_op()


  def apply_ncs_to_sites(self, sites_cart=None,ncs_obj=None,
      exclude_identity=False,ncs_id=None):
    """Apply NCS to sites_cart"""

    if type(sites_cart) == type((1,2,3)): # it was a single site
      from scitbx.array_family import flex
      sites_cart = flex.vec3_double((sites_cart,))

    if not ncs_obj:
      ncs_obj = self
    from scitbx.array_family import flex
    new_sites_cart=flex.vec3_double()
    if (not ncs_obj or ncs_obj.max_operators()<=1):
      if exclude_identity:
        return new_sites_cart  # nothing as we exclude identity
      else:
        return sites_cart

    ncs_group=ncs_obj.ncs_groups()[0]
    from scitbx.matrix import col

    if ncs_id is not None:
      t=ncs_group.translations_orth_inv()[ncs_id]
      r=ncs_group.rota_matrices_inv()[ncs_id]
      for site in sites_cart:
        new_sites_cart.append(r * col(site)  + t)
      return new_sites_cart

    if exclude_identity:
      identity_op=ncs_group.identity_op_id()
    else:
      identity_op=None
    ii=-1
    for t,r in zip(
       ncs_group.translations_orth_inv(),ncs_group.rota_matrices_inv()):
      ii+=1
      if exclude_identity and ii==identity_op: continue
      for site in sites_cart:
        new_sites_cart.append(r * col(site)  + t)
    return new_sites_cart


test_ncs_info="""

new_ncs_group
NCS_CC 0.92
new_operator

rota_matrix    1.0000    0.0000    0.0000
rota_matrix    0.0000    1.0000    0.0000
rota_matrix    0.0000    0.0000    1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth   30.2920   -2.8923   16.6160
CHAIN A
RMSD 0.2
MATCHING 12.0
  RESSEQ 1:26

new_operator

rota_matrix   -0.9971    0.0424   -0.0635
rota_matrix   -0.0297   -0.9816   -0.1889
rota_matrix   -0.0703   -0.1864    0.9800
tran_orth    70.9461    5.2622    3.7549

center_orth   39.8735    3.8824   16.7239
CHAIN B
RMSD 0.1
MATCHING 15.0
  RESSEQ 101:126



new_ncs_group
NCS_CC 0.95
new_operator

rota_matrix    1.0000    0.0000    0.0000
rota_matrix    0.0000    1.0000    0.0000
rota_matrix    0.0000    0.0000    1.0000
tran_orth     0.0000    0.0000    0.0001

center_orth   31.2920   -2.8923   16.6160
CHAIN A
RMSD 0.6
MATCHING 13.0
  RESSEQ 1:25

new_operator

rota_matrix   -0.9970    0.0424   -0.0635
rota_matrix   -0.0297   -0.9816   -0.1889
rota_matrix   -0.0703   -0.1864    0.9800
tran_orth    70.9461    5.2622    3.7549

center_orth   38.8735    3.8824   16.7239
CHAIN B
RMSD 0.5
MATCHING 11.0
  RESSEQ 101:124


"""

test1_ncs_info="""
new_ncs_group
new_operator

rota_matrix    1.0000    0.0000    0.0000
rota_matrix    0.0000    1.0000    0.0000
rota_matrix    0.0000    0.0000    1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth   0.0000    0.0000    0.0000
new_operator

rota_matrix    0.3090    0.9500   -0.0000
rota_matrix   -0.9500    0.3090    0.0000
rota_matrix    0.0000   -0.0000    1.0000
tran_orth   -23.0000  147.0000    0.0000

center_orth   0.0000    0.0000    0.0000
"""
def euler_frac_to_rot_trans(euler_values,frac,unit_cell):
    """Get RT in cctbx form from euler angles and fractional translation as
    used in phaser. Note: Specific for phaser EULER FRAC"""

    ncs_rota_matr_inv=scitbx.rigid_body.euler(
      euler_values[2],euler_values[1],euler_values[0],"zyz").rot_mat()
    ncs_rota_matr=ncs_rota_matr_inv.inverse()

    orth=unit_cell.orthogonalize(frac)
    trans_orth=-1.*ncs_rota_matr*orth
    return ncs_rota_matr,trans_orth
#####################################################



if __name__=="__main__":
  log=sys.stdout
  args=sys.argv[1:]
  if 'exercise' in args:
    file_name='TEST.NCS'
    f=open(file_name,'w')
    f.write(test_ncs_info)
    f.close()
    ncs_object=ncs()
    ncs_object.read_ncs(file_name,source_info=file_name)
    ncs_object.display_all()
    file2='TEST2.NCS'
    text=ncs_object.format_all_for_group_specification(file_name=file2)

    if not text or text != test_ncs_info:
     print("NOT OK ...please compare TEST.NCS (std) vs TEST2.NCS (output)")
     ff=open('txt.dat','w')
     ff.write(text)
     ff.close()
    else:

     print ("Running exercise_1")
     ncs_object=ncs()
     ncs_object.read_ncs(lines=test1_ncs_info.splitlines())
     ncs_lines=ncs_object.format_all_for_group_specification().splitlines()
     biomt_lines=ncs_object.format_all_for_biomt().splitlines()

     biomt_ncs_object=ncs()
     biomt_ncs_object.read_ncs(lines=biomt_lines)

     biomt_text=biomt_ncs_object.display_all()

     std_ncs_object=ncs()
     std_ncs_object.read_ncs(lines=ncs_lines)
     std_ncs_text=std_ncs_object.display_all()
     print ("STANDARD \n %s \n BIOMTR \n %s " %(
       std_ncs_text, biomt_text))
     for a,b in zip(std_ncs_text.splitlines(),biomt_text.splitlines()):
       assert a.strip()==b.strip()

     assert std_ncs_object.is_similar_ncs_object(biomt_ncs_object)
     new_ncs_obj=ncs()
     # Shift ncs coordinates, make ncs with both, make sure can find
     #   similarity even if order is different
     new_ncs_obj._ncs_groups=deepcopy(std_ncs_object._ncs_groups)
     offset_ncs=std_ncs_object.deep_copy(coordinate_offset=(10,10,10))
     assert not offset_ncs.is_similar_ncs_object(biomt_ncs_object)
     new_ncs_obj._ncs_groups.append(deepcopy(offset_ncs._ncs_groups[0]))
     second_ncs_obj=ncs()
     second_ncs_obj._ncs_groups=deepcopy(offset_ncs._ncs_groups)
     second_ncs_obj._ncs_groups.append(deepcopy(std_ncs_object._ncs_groups[0]))
     assert not new_ncs_obj.is_similar_ncs_object(biomt_ncs_object)
     assert new_ncs_obj.is_similar_ncs_object(second_ncs_obj)
     print("OK")

  elif len(args)>0 and args[0] and os.path.isfile(args[0]):
    ncs_object=ncs()
    ncs_object.read_ncs(args[0],source_info=args[0])
    ncs_object.display_all()
    if 1:
      file2='OUTPUT.NCS'
      text=ncs_object.format_all_for_group_specification(file_name=file2)
      print("IS point-group: ", end=' ')
      print(ncs_object.is_point_group_symmetry())
      print("IS helical:", end=' ')
      print(ncs_object.is_helical_along_z())
    if 1:
      file3='OUTPUT2.NCS'
      new_ncs_object=ncs_object.deep_copy(ops_to_keep=[0,1,6])
      text=new_ncs_object.format_all_for_group_specification(file_name=file3)
      print(text)
      print("IS point-group: ", end=' ')
      print(new_ncs_object.is_point_group_symmetry())
      print("IS helical:", end=' ')
      print(new_ncs_object.is_helical_along_z())
