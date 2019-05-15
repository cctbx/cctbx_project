from __future__ import absolute_import, division, print_function
import os
import time
import math
import iotbx
import iotbx.pdb
import libtbx.load_env
from scitbx.array_family import flex
import numpy as np
from qrefine.super_cell import expand
from libtbx import easy_pickle
from libtbx.utils import Sorry
import six
import boost.python
ext = boost.python.import_ext("mmtbx_pair_interaction_ext")

dat_path = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests_data = os.path.join(dat_path,"tests","unit","data_files")
A2B=1.8897161646320724
global results
results=[]
global mol_pair_has_interaction
mol_pair_has_interaction=[]

global wave_functions
wave_functions = []

def load_wfc(element):
  folder = libtbx.env.find_in_repositories("mmtbx/pair_interaction")
  for fn in os.listdir(folder):
    if(fn.startswith(element) and fn.endswith(".pkl")):
      fn = "/".join([folder,fn])
      wfc_obj = easy_pickle.load(fn)
      return wfc_obj
  #
  folder = libtbx.env.find_in_repositories("qrefine/plugin/yoink/dat")
  if(folder is None):
    raise Sorry("No _lda.wfc files found.")
  #
  lines=open(os.path.join(
    dat_path,"./plugin/yoink/dat/"+element+"_lda.wfc")).readlines()
  num_orbitals=int(lines[0].strip().split()[0])
  occ_electrons_array = [[0]*num_orbitals]
  line_three = lines[2].strip().split()
  for i in range(len(line_three)):
    occ_electrons_array[0][i]=int(line_three[i])
  ##TODO check
  occ_electrons=np.array(occ_electrons_array[0])

  line_four = lines[3].strip().split()
  xmin      = float(line_four[0])
  zz        = float(line_four[1])
  dx        = float(line_four[2])
  ngrid     = int(line_four[3])

  r_array = flex.double()
  wfcin_array  = []
  for line in lines[4:]:
    line = line.split()
    r_array.append(float(line[0]))
    tmp = []
    for i in range(1,len(line)):
      tmp.append(float(line[i]))
    wfcin_array.append(tmp)
  assert len(wfcin_array) == ngrid, [len(wfcin_array) , ngrid]

  wfc_obj = ext.wfc(
    ngrid         = ngrid,
    zz            = zz,
    r_array       = r_array,
    wfcin_array   = wfcin_array,
    occ_electrons = occ_electrons)

  easy_pickle.dump("%s_wfc_obj.pkl"%element, wfc_obj)
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
    #print element
    #print dir(wfc_obj)
    #STOP()


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
    print(("1. interactions got: ",len(interactions)))
    for pair in interactions:
      if(len(set(pair).intersection(set(core)))>0):
        new_core+=pair
    new_core=list(set(new_core)|set(core))
    print(("new core:",new_core))
    interaction_mols=new_core
    #
    interactions=get_interactions(sub_ph, atom_in_residue, silva_type='dori', core=new_core)
    print(("2. interactions got: ",len(interactions)))
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
  print(("num atoms:",len(atoms)))
  elements = atoms.extract_element()
  eldict = dict(enumerate(set(elements)))
  tmp = {}
  for k, v in six.iteritems( eldict):
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
  print(("points to process:",xyz_step[0]*xyz_step[1]*xyz_step[2]))
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
  print(("Time ",(time.time()-t0)))
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
