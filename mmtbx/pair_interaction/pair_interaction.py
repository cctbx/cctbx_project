import os
import time
import math
import iotbx
import iotbx.pdb
import libtbx.load_env
from scitbx.array_family import flex
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

def load_wfc(element):
  lines=open(os.path.join(
    dat_path,"./plugin/yoink/dat/"+element+"_lda.wfc")).readlines()
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
  return wfc_obj

def run(ph, core=None):
  atom_in_residue = []
  atoms_group_dict={}
  mols=[]
  cntr=1
  element_types=set()
  for rg in ph.residue_groups():
    atoms_group_dict[cntr]=list(rg.atoms())
    for atom in rg.atoms():
      atom_in_residue.append(cntr)
      e=atom.element.strip(" ")
      if(len(e)==1):e=e+"_"
      element_types.add(e.lower())
    mols.append(cntr)
    cntr+=1
  xyz=ph.atoms().extract_xyz()*A2B
  ph.atoms().set_xyz(xyz)
  atoms=ph.atoms()
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
      selection[atom_list.index(atom)-1]=True
    sub_ph=ph.select(selection)
    del atom_list
    #
    interactions=get_interactions(sub_ph,atom_in_residue,silva_type='sedd',core=core)
    new_core=[]
    print("1. interactions got: ",len(interactions))
    for pair in interactions:
      if(len(set(pair).intersection(set(core)))>0):
        new_core+=pair
    new_core=list(set(new_core)|set(core))
    print("new core:",new_core)
    interaction_mols=new_core
    #
    interactions=get_interactions(sub_ph, atom_in_residue, silva_type='dori', core=new_core)
    print("2. interactions got: ",len(interactions))
    interaction_mols=[]
    for item in interactions:
      interaction_mols=interaction_mols+list(item)
    interaction_atoms=[]
    for i in range(len(core_atoms)):
      core_atoms[i]=core_atoms[i].serial_as_int()
    print(interaction_mols)
    for mol_id in interaction_mols:
      ams=[a.serial_as_int() for a in atoms_group_dict[mol_id]]
      interaction_atoms+=ams
    return(core_atoms, interaction_atoms, interaction_mols)
  else:
    return get_interactions(ph, atom_in_residue)

def get_interactions(ph, atom_in_residue, step_size=0.5*A2B, silva_type='dori',
      core=None):
  t0=time.time()
  atoms = ph.atoms()
  print("num atoms:",len(atoms))
  elements = atoms.extract_element()
  eldict = dict(enumerate(set(elements)))
  tmp = {}
  for k, v in zip(eldict.keys(), eldict.values()):
    tmp[v.strip()] = k
  eldict = tmp
  element_flags = flex.int()
  for e in elements:
    element_flags.append(eldict[e.strip()])
  for e in eldict:
    k = e.lower()
    if(len(k)==1): k+="_"
    wave_functions.append(element_wfc_dict[k])
  xyz = ph.atoms().extract_xyz()
  xyz_min=xyz.min()
  xyz_max=xyz.max()
  xyz_step = [
    int(math.floor((xyz_max[i]-xyz_min[i])/step_size)+1) for i in range(3)]
  print("points to process:",xyz_step[0]*xyz_step[1]*xyz_step[2])
  interacting_pairs = ext.points_and_pairs(
    ngrid           = xyz_step,
    step_size       = step_size,
    xyz             = xyz,
    xyz_min         = xyz_min,
    atom_in_residue = atom_in_residue,
    element_flags   = element_flags,
    wfc_obj         = wave_functions)
  tmp = []
  for it in interacting_pairs:
    pair = [int(it[0]),int(it[1])]
    tmp.append(tuple(pair))
  interacting_pairs = list(set(tmp))
  print("Time ",(time.time()-t0))
  return interacting_pairs

if __name__=="__main__":
  f=qr_unit_tests_data+"/1yjp.pdb"
  print("1. clustering")
  pdb_inp = iotbx.pdb.input(f)
  ph = pdb_inp.construct_hierarchy()
  run(ph)
  print("2. fragment for residues 1 and 2")
  pdb_inp = iotbx.pdb.input(f)
  ph = pdb_inp.construct_hierarchy()
  run(ph,core=[1,2])

  print("3. clustering within expanded ph")
  pdb_inp = iotbx.pdb.input(f)
  cs=pdb_inp.crystal_symmetry()
  ph = pdb_inp.construct_hierarchy()
  expansion = expand(
      pdb_hierarchy        = ph,
      crystal_symmetry     = cs,
      select_within_radius = 10.0).ph_super_sphere
  run(expansion)
  print("4. fragment for residues 1 and 2 within expanded ph")
  pdb_inp = iotbx.pdb.input(f)
  ph = pdb_inp.construct_hierarchy()
  cs=pdb_inp.crystal_symmetry()
  expansion = expand(
      pdb_hierarchy        = ph,
      crystal_symmetry     = cs,
      select_within_radius = 10.0).ph_super_sphere
  run(expansion,core=[1,2])
