import sys, os, re, string
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
      residues_in_common_list=None,cc=None):
    self._chain_residue_id=chain_residue_id  # just one of these
    self._rmsd_list=rmsd_list
    self._residues_in_common_list=residues_in_common_list
    self._centers=center_orth
    self._translations_orth=trans_orth
    self._rota_matrices=ncs_rota_matr
    self._n_ncs_oper=len(self._centers)
    self._source_of_ncs_info=source_of_ncs_info
    self._cc=cc


  def display_summary(self):
    text=""
    text+="\nSummary of NCS group with "+str(self.n_ncs_oper())+" operators:"
    i=0
    if self._chain_residue_id:
      text+="\nID of chain/residue where these apply: "+\
         str(self._chain_residue_id)
    if self._rmsd_list:
      text+="\nRMSD (A) from chain "+str(self._chain_residue_id[0][0])+':'+\
       str(self._rmsd_list)
    if self._residues_in_common_list:
      text+="\nNumber of residues matching chain "+\
          str(self._chain_residue_id[0][0])+':'+\
           str(self._residues_in_common_list)
    if self._source_of_ncs_info:
      text+="\nSource of NCS info: "+str(self._source_of_ncs_info)
    if self._cc:
      text+="\nCorrelation of NCS: "+str(self._cc)

    for center,trans_orth,ncs_rota_matr in zip (
       self._centers, self._translations_orth,self._rota_matrices):
      i+=1
      text+="\n\nOPERATOR "+str(i)
      text+="\nCENTER: "+" %8.4f  %8.4f  %8.4f" %tuple(center)
      text+="\n\nROTA 1: "+" %8.4f  %8.4f  %8.4f" %tuple(ncs_rota_matr[0:3])
      text+="\nROTA 2: "+" %8.4f  %8.4f  %8.4f" %tuple(ncs_rota_matr[3:6])
      text+="\nROTA 3: "+" %8.4f  %8.4f  %8.4f" %tuple(ncs_rota_matr[6:9])
      text+="\nTRANS:  "+" %8.4f  %8.4f  %8.4f" %tuple(trans_orth)
    return text

  def format_for_ncs(self):
    if not self._chain_residue_id or len(self._chain_residue_id)<2:
      return ""
    [group,residue_range_list] = self._chain_residue_id
    count=0
    for id,residue_ranges in zip (group,residue_range_list):
     count+=1
     if count==1:
       text="refinement.ncs_restraint_group { \n"
       text+="reference = chain "+str(id)
     else:
       text+="selection = chain "+str(id)
     if residue_ranges:
       first=True
       for residue_range in residue_ranges:
         if first:
           first=False
           text+=" and (resseq "
         else:
           text+=" or resseq  "
         text+=str(residue_range[0])+":"+str(residue_range[1])
       text+=" ) \n"
    text+= "} \n"
    return text
  def format_for_resolve(self):
    text="new_ncs_group"
    for center,trans_orth,ncs_rota_matr in zip (
       self._centers, self._translations_orth,self._rota_matrices):
      for j in xrange(3):
       text+="\nrota_matrix "+" %8.4f  %8.4f  %8.4f" %tuple(
          ncs_rota_matr[j*3:j*3+3])
      text+="\ntran_orth  "+" %8.4f  %8.4f  %8.4f" %tuple(trans_orth)
      text+="\n"
      text+="\ncenter_orth "+" %8.4f  %8.4f  %8.4f" %tuple(center)
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

  def rota_matrices(self):
    return self._rota_matrices


class ncs:
  def __init__(self):
    self._ncs_groups=[]  # each group is a set of NCS operators and centers

  def ncs_groups(self):
    return self._ncs_groups

  def read_from_resolve(self,file_name,source_info="",log=None):
    if not log: log=sys.stdout
    print >>log,"Reading NCS matrices from resolve file: ",file_name
    if source_info: print >>log," based on ",source_info
    self.source_info=source_info
    if not os.path.isfile(file_name):
      raise Sorry("The file "+str(file_name)+" does not seem to exist?")
    self.init_ncs_group()
    for line in open(file_name).readlines():
      if not line: continue
      if re.search('new_ncs_group',string.lower(line)): # start new group
        if self._n_ncs_oper > 0:
          self.save_ncs_group()
      elif re.search('rota_matrix',string.lower(line)): # read set of rota
        if self._rota_matrix and len(self._rota_matrix)>2:
          self.save_oper()
        set=self.get_3_values_after_key(line)
        self._rota_matrix.append(set)
      elif re.search('tran_orth',string.lower(line)): # read translation
        set=self.get_3_values_after_key(line)
        self._trans=set
      elif re.search('center_orth',string.lower(line)): # read translation
        set=self.get_3_values_after_key(line)
        self._center=set
      else:
        pass

    if self._trans:
      self.save_oper()
    if self._n_ncs_oper > 0:
      self.save_ncs_group()

  def get_3_values_after_key(self,line):
    spl=string.split(line)
    if  len(spl)<4:
      raise Sorry("Cannot interpret this NCS file"+"\n"+str(line))
    set=[]
    try:
      for item in spl[1:4]:
        set.append(string.atof(item))
    except:
      raise Sorry("Cannot interpret this NCS file"+"\n"+str(line))
    return set

  def init_ncs_group(self):
     self._n_ncs_oper=0
     self._ncs_trans_orth=[]
     self._ncs_rota_matr=[]
     self._ncs_center_orth=[]
     self.init_oper()

  def init_oper(self):
     self._rota_matrix=[]
     self._trans=None
     self._center=None

  def save_oper(self):
     self._n_ncs_oper+=1
     for item in (self._trans,self._rota_matrix,self._center):
       if not item:
          raise Sorry("Cannot interpret this NCS file")
       if len(self._rota_matrix)!=3:
          raise Sorry("Cannot interpret this NCS file")
     from scitbx import matrix
     self._ncs_trans_orth.append(matrix.col(self._trans))
     self._ncs_rota_matr.append(matrix.sqr(
      self._rota_matrix[0]+self._rota_matrix[1]+self._rota_matrix[2] ))
     self._ncs_center_orth.append(matrix.col(self._center))
     self.init_oper()

  def import_ncs_group(self,ncs_rota_matr=None,
       center_orth=None,
       trans_orth=None,
       chain_residue_id=None,
       residues_in_common_list=None,
       rmsd_list=None,
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
       source_of_ncs_info=source_of_ncs_info)
     self._ncs_groups.append(ncs_group_object)

  def save_ncs_group(self):
     # check that everything is here:
     for list in [self._ncs_trans_orth,self._ncs_rota_matr,self._ncs_center_orth]:
        if len(list) != self._n_ncs_oper:
          raise Sorry("The NCS operators in this file appear incomplete?")

     ncs_group_object=ncs_group(
       ncs_rota_matr=self._ncs_rota_matr,
       center_orth=self._ncs_center_orth,
       trans_orth=self._ncs_trans_orth,
       source_of_ncs_info=self.source_info)
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
    log.write(text)
    return text


  def format_all_for_resolve(self,log=None,quiet=False,out=None):
    if out==None:
       out=sys.stdout
    if log==None:
      log=sys.stdout
    else:
      print >>log,"\n\nNCS operators written in format for resolve to:",out.name
    all_text=""
    for ncs_group in self._ncs_groups:
      text=ncs_group.format_for_resolve()
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
      text= ncs_group.format_for_ncs()
      if text:
        if not quiet: out.write("\n"+text+"\n\n")
        all_text+="\n"+text
    return all_text

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
    return cc_all
