import os
import time
import heapq
import math
import iotbx
import iotbx.pdb
import libtbx.load_env
from scitbx.array_family import flex
import operator
import numpy as np
from qrefine.super_cell import expand

import boost.python
ext = boost.python.import_ext("mmtbx_pair_interaction_ext")

dat_path = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests_data = os.path.join(dat_path,"tests","unit","data_files")
PI=3.141592920
A2B=1.8897161646320724
global results
results=[]
global mol_pair_has_interaction
mol_pair_has_interaction=[]

global wave_functions
wave_functions = []

#class wfc(object):
#  def __init__(self):
#    self.node_offsets=[
#                        [ 0, -2, -5 ],
#                        [ 1, -1, -4],
#                        [ 2, 0, -3 ],
#                        [ 3, 1, -2 ],
#                        [ 4, 2, -1],
#                        [5,3,0]]
#    self.coefficients_of_first_derivative = [
#                        [ -274, 6, -24 ],
#                        [ 600, -60, 150 ],
#                        [ -600, -40, -400 ],
#                        [ 400, 120, 600 ],
#                        [ -150, -30, -600 ],
#                        [24, 4, 274]]
#
#    self.coefficients_of_second_derivative = [
#                        [ 225, -5, -50 ],
#                        [ -770, 80, 305 ],
#                        [ 1070, -150, -780 ],
#                        [ -780, 80, 1070 ],
#                        [ 305, -5, -770 ],
#                        [ -50, 0, 225 ] ]
#    self.prefactor_of_first_derivative = 1.0 / 120.0
#    self.prefactor_of_second_derivative = 2.0 / 120.0
#    self.a=None
#    self.b=None
#    self.position_max=None
#    self.square_position_max=None
#    self.ngrid=None
#
#    self.grid_positions=None
#    self.grid_values=None
#    self.first_derivative_of_grid_values=None
#    self.second_derivative_of_grid_values=None
#
#    self.core_cutdens = 1E-12

def load_wfc(element):
  lines=open(os.path.join(dat_path,"./plugin/yoink/dat/"+element+"_lda.wfc")).readlines()
  num_orbitals=int(lines[0].strip().split()[0])
  occ_electrons_array = [[0]*num_orbitals]
  line_three = lines[2].strip().split()
  for i in range(len(line_three)):
    occ_electrons_array[0][i]=int(line_three[i])
  ##TODO check
  occ_electrons=np.array(occ_electrons_array[0])
  line_four=lines[3].strip().split()
  xmin=float(line_four[0])
  zz=float(line_four[1])
  dx=float(line_four[2])
  ngrid=int(line_four[3])
  ##read the grid and build the density
  temp_rr_array=[]
  for i in range(ngrid):
    temp_rr_array.append([0,0,0])
  temp_r_array=[0]*ngrid
  wfcin_array=[[0]*num_orbitals]

  node_offsets=[
                        [ 0, -2, -5 ],
                        [ 1, -1, -4],
                        [ 2, 0, -3 ],
                        [ 3, 1, -2 ],
                        [ 4, 2, -1],
                        [5,3,0]]
  coefficients_of_first_derivative = [
                        [ -274, 6, -24 ],
                        [ 600, -60, 150 ],
                        [ -600, -40, -400 ],
                        [ 400, 120, 600 ],
                        [ -150, -30, -600 ],
                        [24, 4, 274]]

  coefficients_of_second_derivative = [
                        [ 225, -5, -50 ],
                        [ -770, 80, 305 ],
                        [ 1070, -150, -780 ],
                        [ -780, 80, 1070 ],
                        [ 305, -5, -770 ],
                        [ -50, 0, 225 ] ]
  prefactor_of_first_derivative = 1.0 / 120.0
  prefactor_of_second_derivative = 2.0 / 120.0
  core_cutdens = 1E-12

  for j in range(ngrid):
    line=lines[4+j].strip().split()
    temp_r_array[j]=float(line[0])
    for m in range(1,len(line)):
      wfcin_array[0][m-1]=float(line[m])
    ##TODO check
    wfcin=np.array(wfcin_array[0])
    wfcin_square=np.multiply(wfcin,wfcin)
    temp_rr_array[j][0]=np.dot(occ_electrons,wfcin_square)
    if(temp_rr_array[j][0]/(4.0*PI*(temp_r_array[j]**2)) < core_cutdens):
      ngrid=j+1
      break
  rr_array=temp_rr_array[:ngrid]
  r_array =temp_r_array[:ngrid]

  #wfc_obj= wfc()

  # XXX TO MOVE TO C++ :

  wfc_obj = ext.wfc(
    node_offsets = node_offsets,
    coefficients_of_first_derivative = coefficients_of_first_derivative,
    coefficients_of_second_derivative = coefficients_of_second_derivative,
    prefactor_of_first_derivative = prefactor_of_first_derivative,
    prefactor_of_second_derivative = prefactor_of_second_derivative,
    core_cutdens = core_cutdens,
    rr_array = rr_array,
    r_array = r_array,
    ngrid = ngrid,
    zz = zz)
#  print "ngrid",ngrid
#
#  ##calculate derivatives
#  first_derivative_of_grid_values = [0]*ngrid
#  second_derivative_of_grid_values =  [0]*ngrid
#  grid_values = [0]*ngrid
#
#  noef=wfc_obj.node_offsets
#  coef1=wfc_obj.coefficients_of_first_derivative
#  coef2=wfc_obj.coefficients_of_second_derivative
#  fac1=wfc_obj.prefactor_of_first_derivative
#  fac2=wfc_obj.prefactor_of_second_derivative
#
#  for i in range(ngrid):
#    ic=1
#    if(i<=1): ic=0
#    elif(i>=ngrid-3): ic=2
#    for j in range(6):
#      rr_array[i][1] = rr_array[i][1] + coef1[j][ic]* rr_array[i + noef[j][ic]][0]
#      rr_array[i][2] = rr_array[i][2] + coef2[j][ic]* rr_array[i + noef[j][ic]][0]
#    rr_array[i][1] = rr_array[i][1] * fac1
#    rr_array[i][2] = rr_array[i][2] * fac2
#    r = r_array[i]
#    r1 = 1.0 / r
#    r2 = r1 * r1
#    r3 = r2 * r1
#    r4 = r3 * r1
#    delta = 1.0 / dx
#    delta2 = delta * delta
#    grid_values[i] = rr_array[i][0] * r2 / (4.0 * PI)
#    first_derivative_of_grid_values[i] = (rr_array[i][1] * delta - 2.0 * rr_array[i][0])* r3 / (4.0 * PI)
#    second_derivative_of_grid_values[i] = (rr_array[i][2] * delta2- 5.0 * rr_array[i][1] * delta + 6.0 * rr_array[i][0])* r4 / (4.0 * PI)
#
#  wfc_obj.grid_positions=r_array
#  wfc_obj.a=math.exp(xmin)/zz
#  wfc_obj.b=dx
#  wfc_obj.ngrid=ngrid
#  wfc_obj.position_max=r_array[ngrid-1]
#  wfc_obj.square_position_max=(r_array[ngrid-1]**2)
#  wfc_obj.first_derivative_of_grid_values=first_derivative_of_grid_values
#  wfc_obj.second_derivative_of_grid_values=second_derivative_of_grid_values
#  wfc_obj.grid_values=grid_values
#  return(wfc_obj)
  return wfc_obj

#def denisty_of_molecule():
#  xyz_type_list=[]
#  results[:]=[]
#  for a in atoms:
#    e=a.element.strip(" ").lower()
#    if(len(e)==1):e=e+"_"
#    xyz_type_list.append((a.xyz,e))

#def cal_hessian_CPP(distanceVector,distanceReciprocal,distanceUnitVector,fac1,fac2):
#  h = ext.hessian(
#    distanceVector     = distanceVector,
#    distanceReciprocal = distanceReciprocal,
#    distanceUnitVector = distanceUnitVector,
#    fac1               = fac1,
#    fac2               = fac2)
#  h = np.array(h).reshape(3,3)
#  return h

#def cal_hessian(distanceVector,distanceReciprocal,distanceUnitVector,fac1,fac2):
#  hessian=np.array([[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]])
#  #hessian[0,0]=9999
#  distanceUnitVector2 = np.multiply(distanceUnitVector,distanceUnitVector)
#  distanceReciprocal2 = distanceReciprocal * distanceReciprocal
#  fac11 = fac1 * distanceReciprocal  #first derivative of density divided by distance
#  fac3 = fac2 + fac11  # second derivative of density + fac11
#  #hessian is a 3x3 symmetrix matrix.get the upper right part at first.
#  for j in range(3):
#    hessian[j,j]= fac3 * distanceUnitVector2[j]- fac11
#    for k in range(j+1,3):
#       hessian[j,k]= distanceReciprocal2 * distanceVector[j]* distanceVector[k] * fac3
#  #get the lower left part of hessian
#  for j in range(3):
#     for  k in range(j+1,3):
#       hessian[k,j]=hessian[j,k]
#  return(hessian)

#class density_props_class(object):
#  def __init__(self,density=0,gradient_vector=np.array([0]*3),hessian=np.array([[0,0,0],[0,0,0],[0,0,0]])):
#    self.density=density
#    self.gradient_vector=gradient_vector
#    self.hessian=hessian
#
#  def add(self,density_props_obj):
#    self.density=self.density+density_props_obj.density
#    self.gradient_vector=self.gradient_vector+density_props_obj.gradient_vector
#    self.hessian=self.hessian+density_props_obj.hessian
#
#  def cal_silva(self):
#    silva=0
#    for i in range(3):
#      temp = 0.0
#      for j in range(3):
#        temp += self.gradient_vector[j] * self.hessian[i, j]
#      silva += (self.density * temp - self.gradient_vector[i]* self.gradient)**2
#    self.silva=silva
#
#  def has_silva_interaction(self,silva_type):
#    if(silva_type=='dori' and self.density<0.0001): return(False)
#    if(silva_type=='sedd' and self.density<0.1): return(False)
#    self.cal_silva()
#    silva_interaction=self.silva
#    if(silva_type=='dori'):
#      silva_interaction *=(4.0 / (self.gradient**3))
#      silva_interaction /=( 1.0 + silva_interaction)
#      if(silva_interaction>=0.9 ):
#        return(True)
#      else:
#        return(False)
#    elif(silva_type=='sedd'):
#      silva_interaction *= (4.0 / (self.density**8))
#      silva_interaction = math.log((1.0 + silva_interaction))
#      if(silva_interaction<=5 ):
#        return(True)
#      else: return(False)



#def atom_density_props(p,a_xyz,wfc_obj,density_only=False):
#  d_vector=np.array(list(a_xyz))-np.array(p)
#  d=max(np.linalg.norm(d_vector), 1.0E-10)
#  d_reciprocal=1.0/d
#  d_unit_vector=d_vector*d_reciprocal
#  f=0
#  fp=0
#  fpp=0
#  if(d<wfc_obj.position_max):
#    ir=0
#    r=0
#    grid_positions=wfc_obj.grid_positions
#    if(d<=grid_positions[0]):
#      ir=1
#      r=grid_positions[0]
#    else:
#      ir=int(1+math.floor(math.log(d/wfc_obj.a)/wfc_obj.b))
#      r=d
#    rr=[0]*4
#    dr1=[0]*4
#    #x1dr12=[[0,0,0,0]]*4
#    x1dr12=[]
#    for i in range(4):
#      x1dr12.append([0,0,0,0])
#    for i in range(4):
#      ii=min(max(ir, 2), wfc_obj.ngrid) - 3 + i
#      rr[i]=grid_positions[ii]
#      dr1[i] = r - rr[i]
#      for j in range(i):
#        x1dr12[i][j] = 1.0 / (rr[i] - rr[j])
#        x1dr12[j][i] = -x1dr12[i][j]
#    for i in range(4):
#      ii = min(max(ir, 2), wfc_obj.ngrid) - 3 + i
#      prod = 1.0
#      for j in range(4):
#        if (i == j):
#          continue
#        prod = prod * dr1[j] * x1dr12[i][j]
#      f = f + wfc_obj.grid_values[ii] * prod
#      fp=fp+wfc_obj.first_derivative_of_grid_values[ii]*prod
#      fpp=fpp+wfc_obj.second_derivative_of_grid_values[ii]*prod
#  density=f
#  if(density_only):return(density)
#  fac1=-fp
#  fac2=fpp
#  gradient_vector=np.array(d_unit_vector)*(-fac1)
#  hessian=cal_hessian_CPP(d_vector, d_reciprocal,d_unit_vector,fac1,fac2)
#  #hessian_old=cal_hessian(d_vector, d_reciprocal,d_unit_vector,fac1,fac2)
#
#  #print
#  #print "hessian c++:", hessian
#  #print "hessian py: ", hessian_old
#  #print np.allclose(hessian, hessian_old)
#
#  #return(density_props_class(density, gradient_vector, hessian))
#  return(ext.density_props(density=density, gradient_vector=gradient_vector, hessian=hessian.reshape(-1)))

def distance_between_atom_point(a,p):
  d_vector=np.array(list(a[0]))-np.array(p)
  return(np.linalg.norm(d_vector))

#def has_interaction_at_point(p, atoms_info, silva_type='dori',debug=False):
def has_interaction_at_point(p, atom_xyz, element_flags, atoms_info):

  #ext.has_interaction_at_point(
  #  p = p, atom_xyz=atom_xyz, element_flags=element_flags, wfc_obj=wave_functions)
  #print "HERE"

  density_props_obj=ext.density_props()
  for a in atoms_info:
    dist = np.linalg.norm( np.array(list(a[0]))-np.array(p) )
    if(dist<15):
      density_props_obj.add(ext.atom_density_props(p, a[0], element_wfc_dict[a[1]]))
      #print "density_props_obj.density:", density_props_obj.density


  #if(silva_type=='sedd'):
  #  atom_distance_dict={}
  #  for a in  atoms_included:
  #    d=distance_between_atom_point(a,p)
  #    if(d<10):
  #      atom_distance_dict[a]=d
  #  atom_distance_dict_sort=sorted(atom_distance_dict.items(), key=operator.itemgetter(1))
  #  if(len(atom_distance_dict_sort)<2):return(False)
  #  atoms_included=[atom_distance_dict_sort[0][0],atom_distance_dict_sort[1][0]]
  #density_props_obj=density_props_class()


  #[density_props_obj.add(atom_density_props(p,a[0],element_wfc_dict[a[1]])) for  a in atoms_included]

  #for a in atoms_included:
  #  density_props_obj.add(ext.atom_density_props(p, a[0], element_wfc_dict[a[1]]))


  density_props_obj.density=max(density_props_obj.density,1.0E-30)
  density_props_obj.gradient=np.dot(density_props_obj.gradient_vector,density_props_obj.gradient_vector)
  #return(density_props_obj.has_silva_interaction(silva_type='dori'))
  print "density_props_obj.gradient", density_props_obj.gradient
  return(density_props_obj.has_silva_interaction())


def point_in_pair(ix,iy,iz,step_size,xyz,xyz_min,atom_ids,atom_in_residue, atom_ids_unique):
   point=(xyz_min[0]+ix*step_size,xyz_min[1]+iy*step_size,xyz_min[2]+iz*step_size)
   ixyz=flex.vec3_double([point]*xyz.size())-xyz
   ixyz_dot=ixyz.dot()
   #print "-"*79
   # XXX Missing one residue ???
   shortest_distances_in_residues = flex.double()
   i = atom_in_residue[0]
   tmp = []
   tmp2 = []
   cntr = 0
   for i_atom in atom_in_residue:
     tmp.append(ixyz_dot[cntr])
     if(i!=i_atom):
       i=i_atom
       shortest_distances_in_residues.append(min(tmp))
       tmp2.append(tmp)
       tmp = []
     cntr+=1
   #print
   #print "tmp2",tmp2
   shortest_distances_in_residues.append(min(tmp2[-1]))
   #print
   #print "-"*80
   #print "shortest_distances_in_residues", shortest_distances_in_residues.size(), len(atom_ids_unique)
   #print list(ixyz_dot)
   #print shortest_distances_in_residues.size()
   #print "-"*80

   first, second = 999999, 999999
   resid_1, resid_2 = None, None
   for i,sdir in enumerate(shortest_distances_in_residues):
     #print "i",i
     sdir = shortest_distances_in_residues[i]
     if(sdir >= 200): continue
     if(sdir < first):
       second = first
       first = sdir
       resid_2 = resid_1
       resid_1 = i
       #print "resid_1", resid_1
     elif(sdir < second and sdir != first):
       if(resid_1 != resid_2):
         second = sdir
         resid_2 = i
         #print "resid_2", resid_2
   #print "first, second",first, second
   #print "resid_1, resid_2", resid_1, resid_2
   #print "atom_in_residue",atom_in_residue
   #print "atom_ids_unique", atom_ids_unique
   #print "ixyz_dot[61]", ixyz_dot[61]
   #print "atom_in_residue[61]",atom_in_residue[61]
   #print "xyz[61]",list(flex.double(xyz[61])/A2B)
   pair = [atom_ids_unique[resid_1],atom_ids_unique[resid_2]]
   pair.sort()
   result_ = (point, tuple(pair))


   result = (None,None)
   ixyz_dot_dict=dict(zip(list(range(xyz.size())),list(ixyz_dot)))
   ixyz_dot_dict_sort=heapq.nsmallest(2,ixyz_dot_dict.items(),key=lambda s:s[1])
   #print "ixyz_dot_dict_sort", ixyz_dot_dict_sort
   if(ixyz_dot_dict_sort[0][1]<200):
     atom_pair=(atom_ids[ixyz_dot_dict_sort[0][0]], atom_ids[ixyz_dot_dict_sort[1][0]])
     mol1=atom_pair[0]
     mol2=atom_pair[1]
     if(mol1!=mol2):
       mol_pair=[mol1,mol2]
       mol_pair.sort()
       if(mol_pair not in mol_pair_has_interaction):
         result = (point,tuple(mol_pair))
   #  else:
   #    return(None,None)
   #else:
   #  return(None,None)
   #print result
   #print result_
   #print
   return result_

def collect_result(result):
    #global results
    if(result is not None):
     results.append(result)

#def check_intermolecular_interaction(p,atoms_info,silva_type,mol_pair):
#  if(mol_pair not in results and has_interaction_at_point(p,atoms_info,silva_type)):
#    return(mol_pair)

def check_intermolecular_interaction(p, atom_xyz, element_flags, mol_pair, atoms_info):
  if(mol_pair not in results and has_interaction_at_point(p, atom_xyz, element_flags, atoms_info)):
    return(mol_pair)

def run(ph,core=None): #

  if(1):
    xyz=ph.atoms().extract_xyz()*A2B
    ph.atoms().set_xyz(xyz)
    atoms=ph.atoms()
    atoms_group_dict={}
    mols=[]
    element_types=set()
    for atom in atoms:
      e=atom.element.strip(" ")
      if(len(e)==1):e=e+"_"
      element_types.add(e.lower())
      if(atom.id_str()[10:] not in mols):
        atoms_group_dict[atom.id_str()[10:]]=[atom]
        mols.append(atom.id_str()[10:])
      else:
        atoms_group_dict[atom.id_str()[10:]].append(atom)
    mol_id_dict=dict(zip(mols,range(1,len(mols)+1)))
    id_mol_dict=dict(zip(range(1,len(mols)+1),mols))
    ##load wfc data for each element
    global element_wfc_dict
    element_wfc_dict={}
    for element in element_types:
      wfc_obj=load_wfc(element)
      element_wfc_dict[element]=wfc_obj
  # End of stage 1
  if(core is not None):
    core_atoms=[]
    for idx,item in enumerate(mols):
      if(idx+1 in core):
        core_atoms=core_atoms+atoms_group_dict[item]
    non_core_atoms=set(atoms)-set(core_atoms)
    non_core_atoms_filtered=[]
    for ai in core_atoms:
      for aj in non_core_atoms:
        if(ai.distance(aj)<15):
          non_core_atoms_filtered.append(aj)
    selection=flex.bool(len(atoms),False)
    atom_list=list(atoms)
    for atom in set(core_atoms+non_core_atoms_filtered):
      #selection[atom.serial_as_int()-1]=True
      selection[atom_list.index(atom)-1]=True
    sub_ph=ph.select(selection)
    del atom_list
    #
    interactions=get_interactions(sub_ph,silva_type='sedd',core=core)
    new_core=[]
    print("1. interactions got: ",len(interactions))
    for item in interactions:
      pair=[mol_id_dict[item[0]],mol_id_dict[item[1]]]
      if(len(set(pair)&set(core))>0):
        new_core+=pair
    new_core=list(set(new_core)|set(core))
    print("new core:",new_core)
    interaction_mols=new_core
    #
    interactions=get_interactions(sub_ph,silva_type='dori',core=new_core)
    print("2. interactions got: ",len(interactions))
    for item in interactions:
      pair=[mol_id_dict[item[0]],mol_id_dict[item[1]]]
      if(len(set(pair)&set(new_core))>0):
        interaction_mols+=pair
    interaction_mols=list(set(interaction_mols))
    interaction_atoms=[]
    for i in range(len(core_atoms)):
      core_atoms[i]=core_atoms[i].serial_as_int()
    for mol_id in interaction_mols:
        mol=id_mol_dict[mol_id]
        ams=[a.serial_as_int() for a in atoms_group_dict[mol]]
        interaction_atoms+=ams
    return(core_atoms,interaction_atoms,interaction_mols)
  else:
    interactions=get_interactions(ph)
    inter_pair=[]
    for item in interactions:
      tmp_pair=[mol_id_dict[item[0]],mol_id_dict[item[1]]]
      tmp_pair.sort()
      inter_pair.append(tmp_pair)
    return(inter_pair)

def get_interactions(ph,step_size=0.5*A2B,silva_type='dori',core=None):
  if(silva_type=='sedd'):step_size=0.3*A2B
  atoms=ph.atoms()
  print("num atoms:",len(atoms))
  atom_ids=[]
  for atom in atoms:
    atom_ids.append(atom.id_str()[10:])
  xyz=ph.atoms().extract_xyz()
  xyz_min=xyz.min()
  xyz_max=xyz.max()
  xyz_step=[ int(math.floor((xyz_max[i]-xyz_min[i])/step_size)+1) for i in range(3)]
  ##get points between each molecule pair
  import multiprocessing as mp
  results[:]=[]
  t0=time.time()
  print("points to process:",xyz_step[0]*xyz_step[1]*xyz_step[2])
  pool = mp.Pool(1)#mp.cpu_count())
  # XXX PVA: move to C++ start

  #print "atom_ids", atom_ids
  #mapping_dict = dict(enumerate(set(atom_ids)))
  #atom_in_residue = []
  #for atom in atom_ids:
  #  for k,v in zip(mapping_dict.keys(),mapping_dict.values()):
  #    if(v==atom):
  #      atom_in_residue.append(k)
  #print "atom_in_residue",atom_in_residue
  atom_in_residue = []
  cntr=0
  for rg in ph.residue_groups():
    for atom in rg.atoms():
      atom_in_residue.append(cntr)
    cntr+=1

  atom_ids_unique = []
  for atom_id in atom_ids:
    if(not atom_id in atom_ids_unique):
      atom_ids_unique.append(atom_id)

  for ix in range(xyz_step[0]):
    for iy in range(xyz_step[1]):
      for iz in range(xyz_step[2]):
        #point_in_pair(ix,iy,iz,step_size,xyz,xyz_min,mols)
        pool.apply_async(point_in_pair, args=(ix,iy,iz,step_size,xyz,xyz_min,atom_ids,atom_in_residue,atom_ids_unique), callback=collect_result)
  # XXX PVA: move to C++ end
  pool.close()
  pool.join()
  print("Time (3D loop): ",(time.time()-t0))
  mol_pair_point_dict={}
  for item in results:
    point,mol_pair=item[0],item[1]
    if(point is not None):
      if(mol_pair not in mol_pair_point_dict):
        mol_pair_point_dict[mol_pair]=[point]
      else:
        mol_pair_point_dict[mol_pair].append(point)
  xyz_type_list=[]
  results[:]=[]
  #
  #
  elements = atoms.extract_element()
  atom_xyz = atoms.extract_xyz()
  eldict = dict(enumerate(set(elements)))
  tmp = {}
  for k, v in zip(eldict.keys(), eldict.values()):
    tmp[v.strip()] = k
  eldict = tmp
#
#  print "eldict", eldict
  print "element_wfc_dict", element_wfc_dict
#
  element_flags = flex.int()
  for e in elements:
    element_flags.append( eldict[e.strip()] )
#  print "element_flags", element_flags

  #global wave_functions = []
  for e in eldict:
    k = e.lower()+"_"
    print "k", k
    wave_functions.append( element_wfc_dict[k] )
    #wave_functions.append( ext.wfc() )
#  print "wave_functions", wave_functions

  for a in atoms:
    e=a.element.strip(" ").lower()
    if(len(e)==1):e=e+"_"
    xyz_type_list.append((a.xyz,e))

 ###
  for i in range(3):
    for mol_pair, points in mol_pair_point_dict.items():
      for p in points:
        #check_intermolecular_interaction(p, atom_xyz, element_flags, mol_pair, wave_functions, xyz_type_list)

        print "p            ", p
        print "atom_xyz     ", atom_xyz
        print "element_flags", element_flags

        o = ext.has_interaction_at_point(
          p = p, a_xyz=atom_xyz, element_flags=element_flags, wfc_obj=wave_functions)
        print "HERE (start)", o
 ###


  print("mol pairs to process:",len(mol_pair_point_dict))
  t0=time.time()
  pool = mp.Pool(1)#mp.cpu_count())
  for mol_pair, points in mol_pair_point_dict.items():
    if(core is not None and len(set(mol_pair)&set(core))==2): # XXX & ?
      continue
    for p in points:
      pool.apply_async(check_intermolecular_interaction,
        args=(p, atom_xyz, element_flags, mol_pair, xyz_type_list), callback=collect_result)
      #pool.apply_async(check_intermolecular_interaction,args=(p,xyz_type_list,silva_type,mol_pair), callback=collect_result)
  pool.close()
  pool.join()
  print("Time ",(time.time()-t0))
  mol_pair_has_interaction=list(set(results))
  return(mol_pair_has_interaction)

if __name__=="__main__":
  f="b.pdb"
  print("1. clustering")
  pdb_inp = iotbx.pdb.input(f)
  ph = pdb_inp.construct_hierarchy()
  run(ph)
  #f=qr_unit_tests_data+"/1yjp.pdb"
  #print("1. clustering")
  #pdb_inp = iotbx.pdb.input(f)
  #ph = pdb_inp.construct_hierarchy()
  #run(ph)
  #print("2. fragment for residues 1 and 2")
  #pdb_inp = iotbx.pdb.input(f)
  #ph = pdb_inp.construct_hierarchy()
  #run(ph,core=[1,2])
  #print("3. clustering within expanded ph")
  #pdb_inp = iotbx.pdb.input(f)
  #cs=pdb_inp.crystal_symmetry()
  #ph = pdb_inp.construct_hierarchy()
  #expansion = expand(
  #    pdb_hierarchy        = ph,
  #    crystal_symmetry     = cs,
  #    select_within_radius = 10.0).ph_super_sphere
  #run(expansion)
  #print("4. fragment for residues 1 and 2 within expanded ph")
  #pdb_inp = iotbx.pdb.input(f)
  #ph = pdb_inp.construct_hierarchy()
  #cs=pdb_inp.crystal_symmetry()
  #expansion = expand(
  #    pdb_hierarchy        = ph,
  #    crystal_symmetry     = cs,
  #    select_within_radius = 10.0).ph_super_sphere
  #run(expansion,core=[1,2])
