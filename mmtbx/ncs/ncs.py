#! /usr/bin/env python
import sys, os, string
from libtbx.utils import Sorry
 # hierarchy:  there can be any number of ncs groups.
 #   each group has a set of NCS operators and centers and may apply
 #      to a part of the structure.
 #  for group in ncs.ncs_groups(): returns list of groups
 #  id=group.chain_and_residue_id() returns id of where it applies
 #  for center in group.centers(): returns list of centers of ncs regions in group
 #  for rota_matr in group.rota_matrices(): returns rota matrices
 #  for trans_orth in group.translations_orth(): returns translation matrices


class ncs_group:  # one group of NCS operators and center and where it applies
  def __init__(self, ncs_rota_matr=None, center_orth=None, trans_orth=None,
      chain_residue_id=None,source_of_ncs_info=None,rmsd_list=None,
      ncs_domain_pdb=None,
      residues_in_common_list=None,cc=None,exclude_h=None,exclude_d=None):
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

    self._exclude_h=exclude_h
    self._exclude_d=exclude_d


  def apply_rt_mx_to_matrix(self,rt_mx=None,matrix=None):
    new_matrix=rt_mx.r*matrix
    return new_matrix

  def apply_cob_to_vector(self,vector=None,
         change_of_basis_operator=None,
         unit_cell=None,new_unit_cell=None):
    frac=unit_cell.fractionalize(vector)
    new_frac = change_of_basis_operator.c() * frac
    new_vector=new_unit_cell.orthogonalize(new_frac)
    return new_vector

  def copy_rot_trans(self,list_of_matrices,list_of_translations,
      change_of_basis_operator=None,
      unit_cell=None,new_unit_cell=None):
    # if change_of_basis_operator is None, then return copy of what we have
    from copy import deepcopy
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
    for r,t in zip(list_of_matrices,list_of_translations):
      if change_of_basis_operator is None:
        new_list_of_matrices.append(deepcopy(r))
        new_list_of_translations.append(deepcopy(t))
      else:  # translations get normal transformation, matrices get A R A_inv
        t_prime=self.apply_cob_to_vector(vector=t,
          change_of_basis_operator=change_of_basis_operator,
           unit_cell=unit_cell,new_unit_cell=new_unit_cell)
        new_list_of_translations.append(t_prime)
        r_prime=a * r * a_inv
        new_list_of_matrices.append(r_prime)
    return new_list_of_matrices,new_list_of_translations

  def copy_vector_list(self,list_of_vectors,
      change_of_basis_operator=None,
         unit_cell=None,new_unit_cell=None):
    from copy import deepcopy
    new_vector_list=[]
    for vector in list_of_vectors:
      if change_of_basis_operator is None:
        new_vector=deepcopy(vector)
      else:
        new_vector=self.apply_cob_to_vector(vector=vector,
          change_of_basis_operator=change_of_basis_operator,
           unit_cell=unit_cell,new_unit_cell=new_unit_cell)
      new_vector_list.append(new_vector)
    return new_vector_list

  def deep_copy(self,change_of_basis_operator=None,unit_cell=None,
      new_unit_cell=None):  # make full copy; 
    # optionally apply change-of-basis operator (requires old, new unit cells)
   
    from mmtbx.ncs.ncs import ncs
    from copy import deepcopy
    new=ncs_group()
    new._chain_residue_id=self._chain_residue_id
    new._rmsd_list=deepcopy(self._rmsd_list)
    new._residues_in_common_list=deepcopy(self._residues_in_common_list)

    # centers simply get affected by the change of basis operator if present
    new._centers=self.copy_vector_list(self._centers,
      change_of_basis_operator=change_of_basis_operator,
         unit_cell=unit_cell,new_unit_cell=new_unit_cell)

    # matrices and translations may need to be adjusted if change of basis set
    new._rota_matrices,new._translations_orth=self.copy_rot_trans(
       self._rota_matrices,self._translations_orth,
         change_of_basis_operator=change_of_basis_operator,
         unit_cell=unit_cell,new_unit_cell=new_unit_cell)

    new._n_ncs_oper=deepcopy(self._n_ncs_oper)
    new._source_of_ncs_info=self._source_of_ncs_info
    new._ncs_domain_pdb=deepcopy(self._ncs_domain_pdb)
    new._cc=deepcopy(self._cc)
    new._exclude_h=self._exclude_h
    new._exclude_d=self._exclude_d
    return new


  def display_summary(self):
    text=""
    text+="\nSummary of NCS group with "+str(self.n_ncs_oper())+" operators:"
    i=0
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
    if not self._chain_residue_id or len(self._chain_residue_id)<2:
      return ""

    # Need to test for existence because we might have operators or
    # chain specifications but not both

    if self._chain_residue_id is not None:
      [group,residue_range_list] = self._chain_residue_id
    else:
      group=self.n_ncs_oper*[None]
      residue_range_list=self.n_ncs_oper*[None]

    if self._centers is not None:
      [centers, translations_orth,rota_matrices]=\
         [self._centers, self._translations_orth,self._rota_matrices]
    else:
      centers=self.n_ncs_oper*[None]
      translations_orth=self.n_ncs_oper*[None]
      rota_matrices=self.n_ncs_oper*[None]

    if self._rmsd_list is not None:
       rmsd_list=self._rmsd_list
    else:
      rmsd_list=self.n_ncs_oper*[None]

    if self._residues_in_common_list is not None:
       residues_in_common_list=self._residues_in_common_list
    else:
      residues_in_common_list=self.n_ncs_oper*[None]

    text="\nnew_ncs_group\n"
    if self._cc is not None: text+="NCS_CC "+str(self._cc)+"\n"
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
       for j in xrange(3):
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

  def format_for_phenix_refine(self):
    if not self._chain_residue_id or len(self._chain_residue_id)<2:
      return ""
    exclude=""
    if self._exclude_h: exclude+=" and (not element H) "
    if self._exclude_d: exclude+=" and (not element D) "
    [group,residue_range_list] = self._chain_residue_id
    count=0
    for id,residue_ranges in zip (group,residue_range_list):
      count+=1
      if count==1:
        text="refinement.ncs.restraint_group { \n"
        text+="reference = chain '"+str(id)+"'"
      else:
        text+="selection = chain '"+str(id)+"'"
      if residue_ranges:
        first=True
        for residue_range in residue_ranges:
          if first:
            first=False
            text+=" and (resseq "
          else:
            text+=" or resseq  "
          text+=str(residue_range[0])+":"+str(residue_range[1])
        text+=" ) "+exclude+"\n"
      else :
        text += "\n"
    text+= "} \n"
    return text

  def format_for_resolve(self,crystal_number=None,skip_identity=False,
       ncs_domain_pdb=True):
    text="new_ncs_group"
    if ncs_domain_pdb and self._ncs_domain_pdb is not None:
        text+="\nncs_domain_pdb "+str(self._ncs_domain_pdb)+"\n"
    i=0
    for center,trans_orth,ncs_rota_matr in zip (
       self._centers, self._translations_orth,self._rota_matrices):
      i+=1
      if i==1 and skip_identity: continue
      for j in xrange(3):
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
    return self._n_ncs_oper

  def chain_residue_id(self):
    return self._chain_residue_id

  def rmsd_list(self):
    return self._rmsd_list

  def cc(self):
    return self._cc

  def add_rmsd_list(self,rmsd_list):
    self._rmsd_list=rmsd_list

  def add_cc(self,cc):
    self._cc=cc

  def residues_in_common_list(self):
    return self._residues_in_common_list

  def add_residues_in_common_list(self,residues_in_common_list):
    self._residues_in_common_list=residues_in_common_list

  def add_chain_residue_id(self,chain_residue_id):
    self._chain_residue_id=chain_residue_id


  def centers(self):
    return self._centers

  def translations_orth(self):
    return self._translations_orth

  def rotations_translations_forward_euler(self):
    # note usual rt is from molecule j to molecule 1. Here it is opposite.
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

  def rota_matrices(self):
    return self._rota_matrices

  def source_of_ncs_info(self):
    return self._source_of_ncs_info

  def ncs_domain_pdb(self):
    return self._ncs_domain_pdb

  def print_list(self,list_of_real):
    text=""
    for number in list_of_real:
     text+="  "+str(self.round(number,2))
    return text

  def round(self,value,n_digit):  # round off value to n_digit digits
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

class ncs:
  def __init__(self,exclude_h=None,exclude_d=None):
    self._ncs_groups=[]  # each group is an ncs_group object
    self.source_info=None
    self._ncs_read=False
    self._exclude_h=exclude_h
    self._exclude_d=exclude_d

  def deep_copy(self,change_of_basis_operator=None,unit_cell=None,
      new_unit_cell=None):  # make a copy
    from mmtbx.ncs.ncs import ncs

    # make new ncs object with same overall params as this one:
    new=ncs(exclude_h=self._exclude_h,exclude_d=self._exclude_d)
    new.source_info=self.source_info
    new._ncs_read=self._ncs_read

    # deep_copy over all the ncs groups:
    for ncs_group in self._ncs_groups:
      new._ncs_groups.append(ncs_group.deep_copy(
         change_of_basis_operator=change_of_basis_operator,
         unit_cell=unit_cell,new_unit_cell=new_unit_cell))
    return new

  def change_of_basis(self,change_of_basis_operator=None,unit_cell=None,
      new_unit_cell=None):
    if change_of_basis_operator is None or unit_cell is None or\
        new_unit_cell is None:
       raise Sorry("For change of basis unit_cell, "+
           "new_unit_cell and operator are all required")
    return self.deep_copy(change_of_basis_operator=change_of_basis_operator,
      unit_cell=unit_cell,new_unit_cell=new_unit_cell)

  def ncs_read(self):
    return self._ncs_read

  def ncs_groups(self):
    return self._ncs_groups

  def read_ncs(self,file_name=None,lines=[],source_info="",log=None,quiet=False):
    if not log: log=sys.stdout
    if not quiet:
      if file_name:
        print >>log,"Reading NCS information from: ",file_name
    if source_info:
       print >>log," based on ",source_info
       self.source_info=source_info
    else:
       self.source_info=str(file_name)
    if file_name:
      if not os.path.isfile(file_name):
        raise Sorry("The file "+str(file_name)+" does not seem to exist?")
      else:
        lines=open(file_name).readlines()
    self.init_ncs_group()

    read_something=False

    for line in lines:
      if not line : continue
      spl=string.split(line)
      if len(spl)<1: continue
      key=string.lower(spl[0])
      if key=='transformations' and len(spl)>1 and \
         string.lower(spl[1])=='formatted':  # start all over!
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

  def save_existing_group_info(self):

        self.save_oper()
        if self._n_ncs_oper > 0:  # save last-read ncs group.
          self.save_ncs_group()


  def get_res_range_after_key(self,line):
    spl=string.split(string.replace(line,':',' '))
    if  len(spl)<3:
      raise Sorry("Cannot interpret this NCS file"+"\n"+str(line))
    start,end=None,None
    try:
      start=string.atoi(spl[1])
      end=string.atoi(spl[2])
    except Exception:
      raise Sorry("Cannot interpret this NCS file"+"\n"+str(line))
    return [start,end]

  def get_1_char_after_key(self,line):
    spl=string.split(line)
    if  len(spl)<2:
      raise Sorry("Cannot interpret this NCS file"+"\n"+str(line))
    char=None
    try:
      char=spl[1]
    except Exception:
      raise Sorry("Cannot interpret this NCS file"+"\n"+str(line))
    return char

  def get_1_value_after_key(self,line):
    spl=string.split(line)
    if  len(spl)<2:
      raise Sorry("Cannot interpret this NCS file"+"\n"+str(line))
    cc=None
    try:
      cc=string.atof(spl[1])
    except Exception:
      raise Sorry("Cannot interpret this NCS file"+"\n"+str(line))
    return cc

  def get_3_values_after_key(self,line):
    spl=string.split(line)
    if  len(spl)<4:
      raise Sorry("Cannot interpret this NCS file"+"\n"+str(line))
    set=[]
    try:
      for item in spl[1:4]:
        set.append(string.atof(item))
    except Exception:
      raise Sorry("Cannot interpret this NCS file"+"\n"+str(line))
    return set

  def init_ncs_group(self):
     self._n_ncs_oper=0
     self._ncs_trans_orth=[]
     self._ncs_rota_matr=[]
     self._ncs_center_orth=[]
     self.init_oper()
     self._rmsd_list=[]
     self._residues_in_common_list=[]
     self._cc=None
     self._ncs_domain_pdb=None
     self._chain_residue_id=[]

     self._list_of_resseq_list=[]
     self._group=[]

  def init_oper(self):
     self._rota_matrix=[]
     self._trans=None
     self._center=None
     self._rmsd=None
     self._residues_in_common=None
     self._resseq_list=[]
     self._chain=None

  def save_oper(self):
     # decide if there is anything to save:

     have_oper=True
     for item in (self._trans,self._rota_matrix,self._center):
       if not item:
          have_oper=False
     if self._rota_matrix and len(self._rota_matrix)!=3:
          raise Sorry("Cannot interpret this NCS file")
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
       source_of_ncs_info=None):
     list_length=None
     for list in [trans_orth,ncs_rota_matr,center_orth]:
       if not list or len(list)<1:
          raise Sorry("The NCS operators in this file appear incomplete?")
       if not list_length: list_length=len(list)
       if list_length!=len(list):
          raise Sorry("The NCS operators in this file appear incomplete?")
     ncs_group_object=ncs_group(
       ncs_rota_matr=ncs_rota_matr,
       center_orth=center_orth,
       trans_orth=trans_orth,
       chain_residue_id=chain_residue_id,
       residues_in_common_list=residues_in_common_list,
       rmsd_list=rmsd_list,
       source_of_ncs_info=source_of_ncs_info,
       ncs_domain_pdb=ncs_domain_pdb,
       cc=cc,
       exclude_h=self._exclude_h,exclude_d=self._exclude_d)
     self._ncs_groups.append(ncs_group_object)

  def save_ncs_group(self):
     # check that there is something  here:
     have_something=False
     for list in [self._ncs_trans_orth,
         self._ncs_rota_matr,self._ncs_center_orth,
         self._residues_in_common_list,self._rmsd_list]:
        if list is not None and self._n_ncs_oper and \
           len(list) != self._n_ncs_oper:
          raise Sorry("The NCS operators in this file appear incomplete?")
        if list is not None and len(list)<2:
          raise Sorry("The NCS operators in this file appear incomplete?")
        if list is not None: have_something=True
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
       cc=self._cc)
     self._ncs_groups.append(ncs_group_object)
     self.init_ncs_group()


  def display_all (self,log=None):
    if log==None:
      log=sys.stdout
    count=0
    text=""
    for ncs_group in self._ncs_groups:
      count+=1
      text+="\n\nGROUP "+str(count)
      text+=ncs_group.display_summary()
    text+="\n\n"
    log.write(text)
    return text


  def format_all_for_group_specification(self,log=None,quiet=True,out=None,
       file=None):
    if file is not None:
       out=open(file,'w')
    if out==None:
       out=sys.stdout
    if log==None:
      log=sys.stdout
    else:
      print >>log,"NCS written as ncs object information to:"\
        ,out.name
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
    return all_text

  def format_all_for_resolve(self,log=None,quiet=False,out=None,
      crystal_number=None,skip_identity=False,ncs_domain_pdb=True):
    if out==None:
       out=sys.stdout
    if log==None:
      log=sys.stdout
    else:
      print >>log,"\n\nNCS operators written in format for resolve to:",out.name
    all_text=""
    for ncs_group in self._ncs_groups:
      text=ncs_group.format_for_resolve(crystal_number=crystal_number,
         skip_identity=skip_identity,ncs_domain_pdb=ncs_domain_pdb)
      if not quiet: out.write("\n"+text+"\n\n")
      all_text+="\n"+text
    return all_text

  def format_all_for_phenix_refine(self,log=None,quiet=False,out=None):
    if out==None:
      out=sys.stdout
    if log==None:
      log=sys.stdout
    else:
      print>>log,\
       "NCS operators written in format for phenix.refine to:",out.name
    all_text=""
    for ncs_group in self._ncs_groups:
      text= ncs_group.format_for_phenix_refine()
      if text:
        if not quiet: out.write("\n"+text+"\n\n")
        all_text+="\n"+text
    return all_text

  def add_source_info(self,source_info):
    if self.source_info is None:
       self.source_info=str(source_info)
    else:
       self.source_info+=str(source_info)

  def add_cc_list(self,cc_list):
   if len(self._ncs_groups) != len(cc_list):
     raise Sorry("Number of NCS groups does not match length of cc_list...")
   for ncs_group,cc in zip(self._ncs_groups,cc_list):
    ncs_group.add_cc(cc)

  def overall_cc(self):
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
    n_max=0
    for ncs_group in self._ncs_groups:
      if ncs_group.n_ncs_oper()>n_max:
        n_max=ncs_group.n_ncs_oper()
    return n_max

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
if __name__=="__main__":
  log=sys.stdout
  args=sys.argv[1:]
  if 'exercise' in args:
    file='TEST.NCS'
    f=open(file,'w')
    f.write(test_ncs_info)
    f.close()
    ncs_object=ncs()
    ncs_object.read_ncs(file,source_info=file)
    ncs_object.display_all()
    file2='TEST2.NCS'
    text=ncs_object.format_all_for_group_specification(file=file2)

    if not text or text != test_ncs_info:
     print "NOT OK ...please compare TEST.NCS (std) vs TEST2.NCS (output)"
     ff=open('txt.dat','w')
     ff.write(text)
     ff.close()
    else:
     print "OK"
  elif len(args)>0 and args[0] and os.path.isfile(args[0]):
    ncs_object=ncs()
    ncs_object.read_ncs(args[0],source_info=args[0])
    ncs_object.display_all()
    if 1:
      file2='OUTPUT.NCS'
      text=ncs_object.format_all_for_group_specification(file=file2)
