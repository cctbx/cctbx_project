from __future__ import absolute_import, division, print_function
import os, sys
import copy
from cctbx.array_family import flex
from libtbx.utils import Sorry
from six.moves import cStringIO as StringIO
from libtbx import easy_run
from six.moves import zip
from six.moves import range

master_phil_str = """
  use_afitt = False
    .type = bool
  ligand_file_name = None
    .type = str
  ligand_names = None
    .type = str
  ff = 'mmff94s'
    .type = str
  scale = 10
    .type = str
"""

class covalent_object:
  def __init__(self):
    self.n_atoms = None
    self.resname = None
    self.res_id = None
    self.ligand_res_id = None
    self.charge = None
    self.partial_charges = []
    self.atom_elements = []
    self.bonds = []
    self.nbonds = None
    self.sites_cart_ptrs = []
    self.formal_charges = []
    self.res_bond = []


class afitt_object:
  def __init__(self,
               ligand_paths,  # ligand CIF restraints file
               ligand_names,  # ligand 3-codes
               pdb_hierarchy, #
               ff='mmff94s',     #
               scale=10, #
               ):
    self.n_atoms = []
    self.resname = ligand_names
    self.res_ids = [] #[chain, altloc,resseq]
    self.charge = []
    self.partial_charges = []
    self.atom_elements = []
    self.bonds = []
    self.nbonds = []
    self.sites_cart_ptrs = []
    self.formal_charges = []
    self.total_model_atoms = 0
    self.ff = ff
    self.scale = scale
    self.ligand_paths = ligand_paths
    self.pdb_hierarchy = pdb_hierarchy
    self.covalent_data = []
    self.occupancies = []

    cif_objects = []
    for ligand_path in ligand_paths:
      cif_objects.append( self.read_cif_file(ligand_path) )
    self.process_cif_objects( cif_objects, pdb_hierarchy )

  def __repr__(self):
    outl = "Afitt object"
    outl += "\n  %-15s : " %"ligand_files"
    for file in getattr(self, "ligand_paths"):
      outl += "%s " %file
    for attr in ["ff", "scale"]:
      outl += "\n  %-15s : %s" % (attr, getattr(self, attr))
    if self.sites_cart_ptrs:
      atoms = self.pdb_hierarchy.atoms()
      for resname, ptrs in zip(self.resname, self.sites_cart_ptrs):
        outl += "\n    %s" % (resname)
        for j, group in enumerate(ptrs):
          outl += "\n      Entity %s" % (j+1)
          for i in group:
            outl += "\n       %5d : %s" % (i, atoms[i].quote()[1:-1])
    else:
      for resname in self.resname:
        outl += "\n    %s" % resname
    return outl

  def read_cif_file(self, ligand_path):
    from iotbx import cif
    cif_object = cif.reader(file_path=ligand_path, strict=False).model()
    return cif_object

  def get_sites_cart_pointers(self,
                              atom_ids,
                              pdb_hierarchy,
                              chain_id,
                              altloc,
                              resseq,
                              ):
    sites_cart_ptrs=[-1]*len(atom_ids)
    #this should be simplified by using iotbx.pdb.atom_selection.cache
    for model in pdb_hierarchy.models():
      for chain in model.chains():
        if chain.id != chain_id: continue
        for conformer in chain.conformers():
          if conformer.altloc != altloc: continue
          for residue in conformer.residues():
            if residue.resseq != resseq: continue
            for atom in residue.atoms():
              for atom_id in atom_ids:
                if atom.name.strip() != atom_id.strip(): continue
                loc=atom_ids.index(atom_id)
                sites_cart_ptrs[loc] = atom.i_seq
    return sites_cart_ptrs

  def get_res_ids(self, pdb_hierarchy, resname):
    ids=[]
    atoms=[]
    #this should be simplified by using iotbx.pdb.atom_selection.cache
    for model in pdb_hierarchy.models():
      for chain in model.chains():
        for conformer in chain.conformers():
          for residue in conformer.residues():
            if residue.resname == resname:
              id_list=[chain.id,conformer.altloc,residue.resseq]
              ids.append(id_list)
              atoms.append([atom.i_seq for atom in residue.atoms()])

    #if different ligand residues have the same chain name and only one
    #of them has an altconf , multiple instances will be created of all
    #of them. This ugly piece of code removes the extra copies of residues
    #that have only one altconf. There must be a prettier way... Yes, I think
    #you need to use the 'pure main conf' and 'pure alt.conf' and 'proper
    #alt.conf' classification (see http://cci.lbl.gov/cctbx_docs/iotbx/iotbx.pdb.html#api-documentation)
    #but haven't done this yet.
    id_to_remove=[]
    for i in range(len(ids)-1):
      for j in range(i+1,len(ids)):
        atoms_i=atoms[i]
        atoms_j=atoms[j]
        if len(list(set(atoms_i) & set(atoms_j))) == len(atoms_i) and \
           len(list(set(atoms_i) & set(atoms_j))) == len(atoms_j):
          id_to_remove.append(i)
          break
    filtered_ids=[]
    for id in range(len(ids)):
      if id not in id_to_remove:
        filtered_ids.append(ids[id])
    return filtered_ids

  def get_occupancies(self, ptrs, pdb_hierarchy):
    for ptr in ptrs:
      for atom in pdb_hierarchy.atoms():
        if atom.i_seq == ptr:
          if 'occ' in locals():
            if occ > atom.occ:
              occ=atom.occ
          else:
            occ=atom.occ
    return occ

  def check_covalent(self, geometry):
    # for each ligand type
    for resname_i,resname in enumerate(self.resname):
      # covalent_data is a list of lists, has one list for each ligand type.
      # Each ligand_type list holds a covalent_object for each instance of that
      # ligand in the model.
      self.covalent_data.append([])
      # lig_atoms = []
      #for each instance of the ligand type
      for instance_i, instance in enumerate(self.res_ids[resname_i]):
        #search for covalent bonds between ligand and other residues
        # lig_atoms = lig_atoms + self.sites_cart_ptrs[resname_i][instance_i]
        # nonlig_atoms = [atom for atom in self.pdb_hierarchy.atoms()
        #                 if atom.i_seq not in lig_atoms ]
        nonlig_atoms = [atom for atom in self.pdb_hierarchy.atoms()
                        if atom.parent().parent().resseq != instance[2]]
        bond = []
        for lig_atm_iseq in self.sites_cart_ptrs[resname_i][instance_i]:
          for atom in nonlig_atoms:
            bond_t = geometry.bond_params_table.lookup(lig_atm_iseq, atom.i_seq)
            if bond_t != None:
              bond.append([atom, lig_atm_iseq])
        #if more than one covalent bond, exit
        assert len(bond) < 2, "Ligand %s has more than one covalent bond " \
                              "to the model. This is unsupported at " \
                              "present." %(resname)
        # if no covalent bonds, add None covalent_object and continue
        if len(bond) == 0:
          self.covalent_data[-1].append(None)
          continue

        #start populating the covalent object
        cov_obj = covalent_object()
        cov_res=bond[0][0].parent()
        cov_obj.resname = cov_res.resname
        cov_obj.res_id=[cov_res.parent().parent().id,
                       cov_res.altloc,
                       cov_res.parent().resseq]
        cov_obj.ligand_resname = resname
        cov_obj.ligand_res_id = instance
        cov_obj.sites_cart_ptrs = [atom.i_seq for atom in cov_res.atoms()]

        #get the cif file for the covalently bound residue (CBR)
        from mmtbx import monomer_library
        import mmtbx.monomer_library.server
        mon_lib_srv = monomer_library.server.server()
        get_func = getattr(mon_lib_srv, "get_comp_comp_id", None)
        if (get_func is not None):
          ml=get_func(comp_id=cov_obj.resname)
        else:
          ml=mon_lib_srv.get_comp_comp_id_direct(comp_id=cov_obj.resname)
        cif_object = ml.cif_object

        # atom id's in the CBR
        cif_atom_ids = [i for i in cif_object['_chem_comp_atom.atom_id']]
        model_atom_ids = [atom.name.strip() for atom in cov_res.atoms()]

        # because order of atoms in pdb can be different from cif, this
        # array points from the pdb atom to location of that atom in the cif array
        # this also takes care of the missings atoms because which, not
        # bein in the model_atom_ids, won't be in the cif_object_ptrs either
        cif_object_ptrs = [cif_atom_ids.index(atom) for atom in model_atom_ids]

        # get charges, elements, formal charges
        partial_charges = [float(i) for i in cif_object['_chem_comp_atom.partial_charge']]
        cov_obj.partial_charges = [partial_charges[ptr] for ptr in cif_object_ptrs]
        atom_elements = [i for i in cif_object['_chem_comp_atom.type_symbol']]
        cov_obj.atom_elements = [atom_elements[ptr] for ptr in cif_object_ptrs]
        if '_chem_comp_atom.charge' in cif_object:
          formal_charges = [float(i) for i in cif_object['_chem_comp_atom.charge']]
          cov_obj.formal_charges = [formal_charges[ptr] for ptr in cif_object_ptrs]

        # get bonds but exclude missing atom bonds and adjust atom indexes
        bond_atom_1 = [i for i in cif_object['_chem_comp_bond.atom_id_1']]
        bond_atom_2 = [i for i in cif_object['_chem_comp_bond.atom_id_2']]
        bond_dict={'single':1, 'double':2, 'triple':3, 'aromatic':4, 'coval':1,
                 'deloc':4}
        bond_type = [bond_dict[i] for i in cif_object['_chem_comp_bond.type']]
        bonds = list(zip(bond_atom_1, bond_atom_2, bond_type))
        for b in bonds:
          if b[0] in model_atom_ids and b[1] in model_atom_ids:
            cov_obj.bonds.append((model_atom_ids.index(b[0]), model_atom_ids.index(b[1]), b[2]))





        # import code; code.interact(local=dict(globals(), **locals())); sys.exit()
        # atom_ids_filt = [i for i in cif_atom_ids if i not in missing_atoms]

        # add the bond between ligand and CBR
        lig_atom_i = self.sites_cart_ptrs[resname_i][instance_i].index(bond[0][1])
        cov_atom_i = model_atom_ids.index(bond[0][0].name.strip())
        cov_obj.res_bond = [lig_atom_i, self.n_atoms[resname_i]+cov_atom_i, 1]

        #MISSING ATOMS
        # find atoms in the CBR cif that are misssing in the
        # model (covalent bond deletes a hydrogen)
        missing_atoms = [i for i in cif_atom_ids if i not in model_atom_ids]
        # assert that all missing atoms would've been bound to the CBR atom
        # linked to the ligand. Other atoms should not be missing!
        for missing_atom in missing_atoms:
          miss_atm_bound_to_cov_atm=False
          for b in bonds:
            if (b[0] == missing_atom and b[1] == bond[0][0].name.strip()) \
              or (b[1] == missing_atom and b[0] == bond[0][0].name.strip()):
                miss_atm_bound_to_cov_atm=True
          assert miss_atm_bound_to_cov_atm==True,\
            "Atom (%s) in residue (%s %s), which is covalently bound to " \
            "ligand (%s %s), is missing" \
            %(missing_atom, cov_res.resname,cov_res.parent().resseq.strip(),
              resname, instance[2].strip())

        cov_obj.charge = sum(cov_obj.partial_charges)
        cov_obj.nbonds = len(cov_obj.bonds)
        cov_obj.n_atoms = len(cov_obj.atom_elements)

        #sanity checks
        assert cov_obj.n_atoms == len(cov_res.atoms())
        assert cov_obj.n_atoms == len(cov_obj.sites_cart_ptrs)

        #add the covalent_object to covalent_data
        self.covalent_data[-1].append(cov_obj)

  def process_cif_objects(self, cif_objects, pdb_hierarchy):
    cif_object_d = {}
    for res in self.resname:
      for cif_object in cif_objects:
        for i, id in enumerate(cif_object['comp_list']['_chem_comp.id']):
          if res == id:
            assert res not in cif_object_d, "More than one cif file containing residue %s!" %res
            cif_object_d[res] = cif_object
      if res not in cif_object_d:
        raise Sorry("No restraints file containing residue %s!" %res)

    for res in self.resname:
      cif_object = cif_object_d[res]
      for i, id in enumerate(cif_object['comp_list']['_chem_comp.id']):
        if res == id:
          self.n_atoms.append(
            int(cif_object['comp_list']['_chem_comp.number_atoms_all'][i]) )
      comp_rname='comp_%s' %res
      assert comp_rname in cif_object, "Residue %s not in cif file!" %res
      try:
        self.partial_charges.append(
          [float(i) for i in cif_object[comp_rname]['_chem_comp_atom.partial_charge']]
          )
      except Exception:
        self.partial_charges.append( [0]*self.n_atoms[-1] )
      self.atom_elements.append(
        [i for i in cif_object[comp_rname]['_chem_comp_atom.type_symbol']]
        )
      atom_ids = \
        [i for i in cif_object[comp_rname]['_chem_comp_atom.atom_id']]
      bond_atom_1 = \
        [atom_ids.index(i) for i in cif_object[comp_rname]['_chem_comp_bond.atom_id_1']]
      bond_atom_2 = \
        [atom_ids.index(i) for i in cif_object[comp_rname]['_chem_comp_bond.atom_id_2']]
      bond_dict={'single':1, 'double':2, 'triple':3, 'aromatic':4, 'coval':1,
                 'deloc':4}
      bond_type = \
        [bond_dict[i] for i in cif_object[comp_rname]['_chem_comp_bond.type']]
      self.bonds.append( list(zip(bond_atom_1, bond_atom_2, bond_type)) )
      self.charge.append( sum(self.partial_charges[-1]) )
      self.nbonds.append ( len(self.bonds[-1]) )
      res_ids = self.get_res_ids(pdb_hierarchy, res)
      self.res_ids.append(res_ids)
      this_res_sites_cart_ptrs=[]
      for residue_instance in self.res_ids[-1]:
        this_res_sites_cart_ptrs.append( self.get_sites_cart_pointers(
                                          atom_ids,
                                          pdb_hierarchy,
                                          chain_id=residue_instance[0],
                                          altloc=residue_instance[1],
                                          resseq=residue_instance[2])
                                        )
        for ptrs in this_res_sites_cart_ptrs:
          if ptrs.count(-1)>2:
            raise Sorry("Atoms missing from %s. Likely need to add hydrogens." % res)
      self.sites_cart_ptrs.append( this_res_sites_cart_ptrs )
      this_occupancies=[]
      for ptrs in this_res_sites_cart_ptrs:
        this_occupancies.append( self.get_occupancies(ptrs, pdb_hierarchy) )
      self.occupancies.append( this_occupancies )
      if '_chem_comp_atom.charge' in cif_object[comp_rname]:
        self.formal_charges.append(
          [float(i) for i in cif_object[comp_rname]['_chem_comp_atom.charge']]
          )
      else:
        self.formal_charges.append([])

    self.total_model_atoms=pdb_hierarchy.atoms_size()


  def make_afitt_input(self, sites_cart, resname_i, instance_i):
    r_i=resname_i
    i_i=instance_i
    sites_cart_ptrs=self.sites_cart_ptrs[r_i][i_i]
    elements=self.atom_elements[r_i]
    assert len(elements) ==  len(sites_cart_ptrs), \
           "No. of atoms in residue %s, instance %d does not equal to \
           number of atom seq pointers." %(self.resname[resname_i], instance_i)
    f = StringIO()
    if len(self.covalent_data) == 0 or self.covalent_data[r_i][i_i] == None:
    # if True:
      f.write(  '%d\n' %self.n_atoms[r_i])
      f.write('residue_type %s chain %s number %d total_charge %d\n'
              %(self.resname[r_i], self.res_ids[r_i][i_i][0],1,self.charge[r_i] ))
      #~ import code; code.interact(local=dict(globals(), **locals()))
      for atom,ptr in zip(elements, sites_cart_ptrs):

        f.write('%s   %20.16f   %20.16f   %20.16f\n' %(atom,
              sites_cart[ptr][0], sites_cart[ptr][1], sites_cart[ptr][2]) )
      f.write('bond_table_nbonds %d\n' %self.nbonds[r_i])
      for bond in self.bonds[r_i]:
        f.write('%d %d %d\n' %(bond[0], bond[1], bond[2]))
      if self.formal_charges[r_i]:
        n_non_zero_charges = len([ch for ch in self.formal_charges[r_i] if ch != 0])
        f.write("formal_charges %d\n" %n_non_zero_charges)
        if self.formal_charges[r_i]:
          for i,fcharge in enumerate(self.formal_charges[r_i]):
            if fcharge != 0:
              f.write ('%d %d\n' %(i,fcharge))
      f.write('fixed_atoms 0\n')
    else:
      cov_obj =  self.covalent_data[r_i][i_i]
      # print '%d %d\n ' %(self.n_atoms[r_i] , cov_obj.n_atoms)
      f.write('%d\n' %(self.n_atoms[r_i] + cov_obj.n_atoms) )
      f.write('residue_type %s chain %s number %d total_charge %d\n'
              %(self.resname[r_i],
                self.res_ids[r_i][i_i][0],
                1,
                self.charge[r_i] + cov_obj.charge ))
      for atom,ptr in zip(elements, sites_cart_ptrs):
        f.write('%s   %20.16f   %20.16f   %20.16f\n' %(atom,
              sites_cart[ptr][0], sites_cart[ptr][1], sites_cart[ptr][2]) )
      for atom,ptr in zip(cov_obj.atom_elements, cov_obj.sites_cart_ptrs):
        f.write('%s   %20.16f   %20.16f   %20.16f\n' %(atom,
              sites_cart[ptr][0], sites_cart[ptr][1], sites_cart[ptr][2]) )
      # import code; code.interact(local=dict(globals(), **locals()))
      f.write('bond_table_nbonds %d\n'
              %(self.nbonds[r_i]+cov_obj.nbonds+1) )
      for bond in self.bonds[r_i]:
        f.write('%d %d %d\n' %(bond[0], bond[1], bond[2]))
      for bond in cov_obj.bonds:
        f.write('%d %d %d\n' %(
          bond[0] + self.n_atoms[r_i],
          bond[1] + self.n_atoms[r_i],
          bond[2]))
      f.write('%d %d %d\n' %(
          cov_obj.res_bond[0],
          cov_obj.res_bond[1],
          cov_obj.res_bond[2]))
      if self.formal_charges[r_i] or cov_obj.formal_charges:
        n_non_zero_charges = 0
        if self.formal_charges[r_i]:
          n_non_zero_charges += len([ch for ch in self.formal_charges[r_i] if ch != 0])
        if cov_obj.formal_charges:
          n_non_zero_charges += len([ch for ch in cov_obj.formal_charges if ch != 0])
        f.write("formal_charges %d\n" %n_non_zero_charges)
        if self.formal_charges[r_i]:
          for i,fcharge in enumerate(self.formal_charges[r_i]):
            if fcharge != 0:
              f.write ('%d %d\n' %(i,fcharge))
        if cov_obj.formal_charges:
          for i, fcharge in enumerate(cov_obj.formal_charges):
            if fcharge != 0:
              f.write('%d %d\n' %(i+self.n_atoms[r_i], fcharge))
      f.write('fixed_atoms %d\n' %cov_obj.n_atoms)
      for i in range(cov_obj.n_atoms):
        f.write('%d\n' %(i+self.n_atoms[r_i]))
    ##
    # ofile=open('tmpfile','w')
    # ofile.write(f.getvalue())
    # ofile.close()
    # sys.exit()
    # print f.getvalue()
    return f.getvalue()

def get_afitt_command():
  cmd = "afitt_helper_mmff help" # used because afitt_helper_mmff hangs on no input
  ero = easy_run.fully_buffered(command=cmd,
                               )
  out = StringIO()
  ero.show_stderr(out=out)
  exe = "afitt_helper_mmff"
  if out.getvalue().find("OpenEye")>-1:
    return exe
  if os.environ.get("OE_EXE", False):
    oe_dir = os.environ.get("OE_EXE")
    exe = os.path.join(oe_dir, exe)
    if os.path.exists(exe):
      return exe
  return None

def call_afitt(afitt_input, ff):
  exe = get_afitt_command()
  if exe is None:
    raise Sorry("AFITT command not found. Add to path or correctly set OE_EXE")
  cmd = '%s -ff %s' % (exe, ff)
  #print cmd
  #print afitt_input
  ero = easy_run.fully_buffered(command=cmd,
                                stdin_lines=afitt_input,
                               )
  out = StringIO()
  ero.show_stdout(out=out)
  if 'ENERGYTAG' not in out.getvalue().split():
    ero.show_stderr()
    print("AFITT energy call exited with errors printed above.")
    sys.exit()
  return out

def process_afitt_output(afitt_output,
                         geometry,
                         afitt_o,
                         resname_i,
                         instance_i,
                         afitt_allgradients,
                         afitt_alltargets,
                         verbose=False,
                         phenix_gnorms=None):
  r_i=resname_i
  i_i=instance_i
  ptrs = afitt_o.sites_cart_ptrs[r_i][i_i]
  afitt_gradients = flex.vec3_double()
  for line in afitt_output.getvalue().splitlines():
    if line.startswith('ENERGYTAG'):
       afitt_energy=float(line.split()[1])
    elif line.startswith('GRADIENTTAG'):
       afitt_gradients.append (
          (float(line.split()[1]),
           float(line.split()[2]),
           float(line.split()[3]) ) )
  ### debug_stuff
  if verbose:
    print ("AFITT_ENERGY %s_%d_%s: %10.4f\n"
                  %(afitt_o.resname[r_i],
                    int(afitt_o.res_ids[r_i][i_i][2]),
                    afitt_o.res_ids[r_i][i_i][1],
                    afitt_energy ))
  ### end_debug
  #geometry.residual_sum += afitt_energy
  #~ import inspect
  #~ for i in inspect.stack():
    #~ print i[1], i[2], i[4]
  #~ print "\n\n\n\n"
  cov_ptrs=[]
  if afitt_o.covalent_data[r_i][i_i] is not None:
    cov_ptrs= afitt_o.covalent_data[r_i][i_i].sites_cart_ptrs
  if (geometry.gradients is not None):
    # AFITT prints gradient lines for fixed atoms too so I need to check
    # that no. of gradient == ligand atoms + fixed atoms. But since the
    # fixed atom gradient's are always zero and we add gradients to the
    # Phenix gradients, I don't add code to actually do anything with the
    # fixed atom gradients. NOTE: if one day for some reason we decide to
    # replace Phenix gradients with AFITT gradients, this would need to be added.
    assert afitt_gradients.size() == len(ptrs)  +  len(cov_ptrs)
    if afitt_o.scale == 'gnorm':
      from math import sqrt
      # phenix_norm=phenix_gnorms[r_i][i_i]
      phenix_norm=0
      afitt_norm=0
      for afitt_gradient, ptr in zip(afitt_gradients, ptrs):
        phenix_norm += geometry.gradients[ptr][0]**2+geometry.gradients[ptr][1]**2+geometry.gradients[ptr][2]**2
        afitt_norm += afitt_gradient[0]**2+afitt_gradient[1]**2+afitt_gradient[2]**2
      phenix_norm = sqrt(phenix_norm)
      afitt_norm = sqrt(afitt_norm)
      gr_scale = phenix_norm/afitt_norm
      ### debug_stuff
      if verbose:
        print ("GRNORM_RATIO %s_%d_%s: %10.4f\n"
                    %(afitt_o.resname[r_i],
                      int(afitt_o.res_ids[r_i][i_i][2]),
                      afitt_o.res_ids[r_i][i_i][1],
                      gr_scale ))
        # print phenix_norm, afitt_norm
      ### end_debug
    elif afitt_o.scale == 'noafitt':
      gr_scale = None
    else:
      gr_scale = float(afitt_o.scale)

    ### debug_stuff
    print_gradients = False
    if print_gradients:
      print("\n\nGRADIENTS BEFORE AFTER AFITT\n")
      # print "NORMS: %10.4f         %10.4f\n" %(phenix_norm, afitt_norm)
      for afitt_gradient, ptr in zip(afitt_gradients, ptrs):
        print("(%10.4f %10.4f %10.4f) (%4.4f %4.4f %4.4f)" \
            %(geometry.gradients[ptr][0], geometry.gradients[ptr][1], geometry.gradients[ptr][2],
            afitt_gradient[0], afitt_gradient[1], afitt_gradient[2]))
      sys.exit()
    ### end_debug
    if gr_scale:
      scaled_gradients = []
      # occupancy = afitt_o.occupancies[r_i][i_i]
      for afitt_gradient in afitt_gradients:
        scaled_gradient = (afitt_gradient[0]*gr_scale,
                           afitt_gradient[1]*gr_scale,
                           afitt_gradient[2]*gr_scale)
        scaled_gradients.append(scaled_gradient)
      afitt_allgradients[(r_i,i_i)] = scaled_gradients
      afitt_alltargets[(r_i,i_i)] = gr_scale*afitt_energy

def apply_target_gradients(afitt_o, geometry, afitt_allgradients, afitt_alltargets):
  # import code; code.interact(local=dict(globals(), **locals()))
  # sys.exit()
  if (geometry.gradients is not None):
    for key in afitt_allgradients:
      r_i = key[0]
      i_i = key[1]
      gradients = afitt_allgradients[key]
      target = afitt_alltargets[key]
      ptrs=afitt_o.sites_cart_ptrs[r_i][i_i]
      cov_ptrs=[]
      if afitt_o.covalent_data[r_i][i_i] is not None:
        cov_ptrs= afitt_o.covalent_data[r_i][i_i].sites_cart_ptrs
      for i_seq, gradient in zip(ptrs+cov_ptrs,gradients):
        gx = gradient[0] + geometry.gradients[i_seq][0]
        gy = gradient[1] + geometry.gradients[i_seq][1]
        gz = gradient[2] + geometry.gradients[i_seq][2]
        geometry.gradients[i_seq] = (gx,gy,gz)
      geometry.residual_sum += target

  #   for i_seq in afitt_allgradients:
  #     gradient=[0,0,0]
  #     for loc in afitt_allgradients[i_seq]:
  #       gradient[0] += loc[0]
  #       gradient[1] += loc[1]
  #       gradient[2] += loc[2]
  #     if len(afitt_allgradients[i_seq]) >1:
  #       for r in range(3):
  #         gradient[r] /= len(afitt_allgradients[i_seq])
  #     gx = gradient[0] + geometry.gradients[i_seq][0]
  #     gy = gradient[1] + geometry.gradients[i_seq][1]
  #     gz = gradient[2] + geometry.gradients[i_seq][2]
  #     geometry.gradients[i_seq] = (gx,gy,gz)
  # for target in afitt_alltargets:
  #   geometry.residual_sum += afitt_alltargets[target]

  return geometry

def get_afitt_energy(cif_file,
                     ligand_names,
                     pdb_hierarchy,
                     ff,
                     sites_cart,
                     geometry=None):
  afitt_o = afitt_object(
                cif_file,
                ligand_names,
                pdb_hierarchy,
                ff)

  if geometry is not None:
    afitt_o.check_covalent(geometry)

  energies=[]
  for resname_i,resname in enumerate(afitt_o.resname):
    for instance_i, instance in enumerate(afitt_o.res_ids[resname_i]):
      #~ import code; code.interact(local=dict(globals(), **locals()))
      afitt_input = afitt_o.make_afitt_input(sites_cart,
                                             resname_i,
                                             instance_i,
                                             )
      lines = call_afitt(afitt_input, ff)
      for line in lines.getvalue().splitlines():
        if line.startswith('ENERGYTAG'):
          energy=float(line.split()[1])
      energies.append([resname, int(instance[2]), instance[1].strip(), energy] )
  return energies

def validate_afitt_params(params):
  if params.ligand_names is None:
    raise Sorry("Ligand name(s) not specified\n\t afitt.ligand_names=%s" %
                params.ligand_names)
  if params.ff not in ["mmff94", "mmff94s", "pm3", "am1"]:
    raise Sorry("Invalid force field\n\t afitt.ff=%s" % params.ff)

def get_non_afitt_selection(restraints_manager,
                            sites_cart,
                            hd_selection,
                            ignore_hd,
                            verbose=False,
                            ):
  if ignore_hd:
    general_selection = ~hd_selection
  else:
    general_selection = hd_selection|~hd_selection
  ligand_i_seqs = []
  for ligand in restraints_manager.afitt_object.sites_cart_ptrs:
    for group in ligand:
      ligand_i_seqs += group
  for i_seq in ligand_i_seqs:
    general_selection[i_seq] = False
  if verbose:
    print(restraints_manager.afitt_object)
    print("\nNumber of atoms in selection : %d" % len(list(filter(None, general_selection))))
    print(list(general_selection))
  return general_selection

def get_afitt_selection(restraints_manager,
                        sites_cart,
                        hd_selection,
                        ignore_hd,
                        verbose=False,
                        ):
  if ignore_hd:
    hd_selection = ~hd_selection
  else:
    hd_selection = hd_selection|~hd_selection
  general_selection = hd_selection&~hd_selection
  ligand_i_seqs = []
  for ligand in restraints_manager.afitt_object.sites_cart_ptrs:
    for group in ligand:
      ligand_i_seqs += group
  for i_seq in ligand_i_seqs:
    general_selection[i_seq] = True
  rc = general_selection&hd_selection
  if verbose:
    print(restraints_manager.afitt_object)
    print("\nNumber of atoms in selection : %d" % len(list(filter(None, general_selection))))
  return rc

def write_pdb_header(params, out=sys.stdout, remark="REMARK   3  "):
  print("%sAFITT PARAMETERS" % (remark), file=out)
  for attr in params.__dict__:
    if attr.find("__")==0: continue
    print("%s  %s: %s" % (remark,
                                  attr.upper(),
                                  str(getattr(params, attr)).upper(),
                                 ), file=out)
  print("%s" % remark, file=out)

def _show_gradient(g):
  return "(%9.3f %9.3f %9.3f)" % (g)

def adjust_energy_and_gradients(result,
                                restraints_manager,
                                sites_cart,
                                hd_selection,
                                afitt_o,
                                verbose=False):
  # import code; code.interact(local=dict(globals(), **locals()))
  # sys.exit()
  if result.afitt_residual_sum<1e-6:
    if verbose: 'returning without adjusting energy and gradients'
    return result
  general_selection = get_non_afitt_selection(restraints_manager, sites_cart, hd_selection, False)
  rm = restraints_manager.select(general_selection)
  old_normalisation = getattr(rm, "normalization", None)
  if old_normalisation is None:
    es = rm.energies_sites(
      sites_cart = sites_cart.select(general_selection),
      compute_gradients = True,
      normalization = False,
    )
  else:
    rm.normalization=False
    es = rm.energies_sites(
      sites_cart = sites_cart.select(general_selection),
      compute_gradients = True,
    )
    rm.normalization = old_normalisation
  protein_residual_sum = es.residual_sum
  protein_gradients = es.gradients
  #
  general_selection = get_afitt_selection(restraints_manager,
                                          sites_cart,
                                          hd_selection,
                                          False,
                                         )
  rm = restraints_manager.select(general_selection)
  old_normalisation = getattr(rm, "normalization", None)
  if old_normalisation is None:
    es = rm.energies_sites(
      sites_cart = sites_cart.select(general_selection),
      compute_gradients = True,
      normalization = False,
    )
  else:
    rm.normalization=False
    es = rm.energies_sites(
      sites_cart = sites_cart.select(general_selection),
      compute_gradients = True,
    )
    rm.normalization = old_normalisation
  ligand_residual_sum = es.residual_sum
  ligand_gradients = es.gradients
  #
  if verbose:
    print('gradients')
    print('phenix + afitt')
    for i, s in enumerate(general_selection):
      ls = ""
      if s: ls = "*"
      print("%3d %s %s" % (i+1,_show_gradient(result.gradients[i]), ls))
    print('protein-ligand complex')
    for i, s in enumerate(general_selection):
      ls = ""
      if s: ls = "*"
      print("%3d %s %s" % (i+1,_show_gradient(result.complex_gradients[i]), ls))
    print('protein only')
    for i, s in enumerate(protein_gradients):
      print("%3d %s" % (i+1,_show_gradient(s)))
    print('ligand only')
    for i, s in enumerate(ligand_gradients):
      print("%3d %s" % (i+1,_show_gradient(s)))

  result.residual_sum -= ligand_residual_sum

  ligand_i = 0
  # protein_i = 0
  if verbose:
    print("%-40s %-40s %-40s %-40s" % ("phenix protein+ligand",
                                       "phenix+afitt",
                                       "phenix ligand only",
                                       "phenix+afitt final",
                                       ))
  for i, g in enumerate(result.complex_gradients):
    if verbose:
      outl = "%5d %s %s" % (i,_show_gradient(g),str(general_selection[i])[0])
      outl += _show_gradient(result.gradients[i])
    if general_selection[i]:
      # ligand
      result.gradients[i] = (
        result.gradients[i][0] - ligand_gradients[ligand_i][0],
        result.gradients[i][1] - ligand_gradients[ligand_i][1],
        result.gradients[i][2] - ligand_gradients[ligand_i][2],
      )
      if verbose:
        outl += " %3d %s %s" % ( ligand_i,
                                 _show_gradient(ligand_gradients[ligand_i]),
                                 _show_gradient(result.gradients[i]),
                                 )
      ligand_i+=1
    if verbose: print(outl)

  if verbose:
    print('total (phenix+afitt) residual_sum',result.residual_sum)
    print(result.complex_residual_sum)
    print('complex_residual_sum', result.complex_residual_sum)
    print('afitt_residual_sum', result.afitt_residual_sum)
    print('ligand_residual_sum',ligand_residual_sum)
    print('protein_residual_sum',protein_residual_sum)
    #print 'nonbonded_residual_sum',nonbonded_residual_sum
    print('\n\n')
    print('gradients')
    print('protein only')
    for i, s in enumerate(protein_gradients):
      print("%3d %s" % (i,_show_gradient(s)))
    print('ligand only')
    for i, s in enumerate(ligand_gradients):
      print("%3d %s" % (i,_show_gradient(s)))
    print('protein-ligand complex')
    for i, s in enumerate(result.complex_gradients):
      print("%3d %s" % (i,_show_gradient(s)))
    print('unadjusted')
    for i, s in enumerate(result.gradients):
      print("%3d %s" % (i,_show_gradient(s)))

  return result

def adjusted_phenix_g_norm(geometry,
                            restraints_manager,
                            sites_cart,
                            hd_selection,
                            afitt_o,
                            verbose=False):
  from math import sqrt
  general_selection = get_afitt_selection(restraints_manager,
                                          sites_cart,
                                          hd_selection,
                                          False,
                                         )
  rm = restraints_manager.select(general_selection)
  old_normalisation = getattr(rm, "normalization", None)
  if old_normalisation is None:
    es = rm.energies_sites(
      sites_cart = sites_cart.select(general_selection),
      compute_gradients = True,
      normalization = False,
    )
  else:
    rm.normalization=False
    es = rm.energies_sites(
      sites_cart = sites_cart.select(general_selection),
      compute_gradients = True,
    )
    rm.normalization = old_normalisation
  ligand_gradients = es.gradients
  ligand_gradients_and_sites={}
  i=0
  for site in range(len(general_selection)):
    if general_selection[site]:
      ligand_gradients_and_sites[site]=ligand_gradients[i]
      i+=1



  gnorms=[]
  for resname_i,resname in enumerate(afitt_o.resname):
    gnorms.append([])
    for instance_i, instance in enumerate(afitt_o.res_ids[resname_i]):
      gnorm = 0
      ptrs = afitt_o.sites_cart_ptrs[resname_i][instance_i]
      # import code; code.interact(local=dict(globals(), **locals()))
      for ptr in ptrs:
        l1 = ligand_gradients_and_sites[ptr][0]
        l2 = ligand_gradients_and_sites[ptr][1]
        l3 = ligand_gradients_and_sites[ptr][2]
        gnorm += l1**2+l2**2+l3**2
        i+=1
      gnorm = sqrt(gnorm)
      gnorms[resname_i].append(gnorm)
  return gnorms


def finite_difference_test(pdb_file,
                           cif_file,
                           ligand_names,
                           atom,
                           scale=1,
                           verbose=False):
  from mmtbx import monomer_library
  import mmtbx.monomer_library.server
  import mmtbx.monomer_library.pdb_interpretation
  import iotbx.pdb

  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv    = mon_lib_srv,
    ener_lib       = ener_lib,
    file_name      = pdb_file,
    raw_records    = None,
    force_symmetry = True)
  pdb_inp = iotbx.pdb.input(file_name=pdb_file)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  pdb_hierarchy.atoms().reset_i_seq()
  xrs = pdb_hierarchy.extract_xray_structure()
  sites_cart=xrs.sites_cart()

  grm = processed_pdb_file.geometry_restraints_manager(
    show_energies = False,
    plain_pairs_radius = 5.0,
    )
  afitt_o = afitt_object(
              [cif_file],
              ligand_names,
              pdb_hierarchy,
              scale=scale)
  afitt_o.check_covalent(grm)

  if verbose: print("Analytical Gradient")

  geometry = grm.energies_sites(
    sites_cart        = sites_cart,
    compute_gradients = True)
  if verbose: print("  phenix target:   %10.16f" %geometry.target)
  if verbose: print("  phenix gradient: %10.16f" %geometry.gradients[atom][0])

  geometry.complex_residual_sum = geometry.residual_sum
  geometry.complex_gradients = copy.deepcopy(geometry.gradients)
  afitt_allgradients = {}
  afitt_alltargets = {}
  for resname_i,resname in enumerate(afitt_o.resname):
    for instance_i, instance in enumerate(afitt_o.res_ids[resname_i]):
      afitt_input = afitt_o.make_afitt_input(sites_cart,
                                             resname_i,
                                             instance_i)

      lines = call_afitt(afitt_input,
                         afitt_o.ff)

      # print "pawel"
      # print afitt_o.covalent_data

      process_afitt_output(
          lines, geometry, afitt_o,
          resname_i, instance_i, afitt_allgradients, afitt_alltargets)
  if verbose: print("  afitt target:    %10.16f" %afitt_alltargets[(0,0)])
  if verbose:
    if atom in afitt_o.sites_cart_ptrs[0][0]:
      i = afitt_o.sites_cart_ptrs[0][0].index(atom)
      print("  afitt gradients: %10.16f" %afitt_allgradients[(0,0)][i][0])

  geometry = apply_target_gradients(afitt_o,
                                    geometry,
                                    afitt_allgradients,
                                    afitt_alltargets)
  geometry.afitt_residual_sum = geometry.residual_sum -\
                                geometry.complex_residual_sum
  grm.afitt_object = afitt_o
  geometry = adjust_energy_and_gradients(geometry,
                                         grm,
                                         xrs.sites_cart(),
                                         xrs.hd_selection(),
                                         afitt_o,
                                       )
  geometry.target = geometry.residual_sum


  if verbose: print("  final target:    %10.16f" %geometry.target)
  if verbose: print("  final gradient:  %10.16f" %geometry.gradients[atom][0])
  ana_gradient = geometry.gradients[atom][0]
  print("-> %10.9f"%(ana_gradient))

  if verbose: print("\nFinite Diff. Gradient")
  # finite differences
  e = 1.e-5
  site_cart_o = sites_cart[atom]
  ts = []
  phts = []
  afts = []
  for e_ in [e, -1*e]:
    if verbose: print("e = %f" %e_)
    afitt_allgradients = {}
    afitt_alltargets = {}
    site_cart = [site_cart_o[0]+e_,site_cart_o[1],site_cart_o[2]]
    sites_cart[atom] = site_cart
    geometry = grm.energies_sites(
      sites_cart        = sites_cart,
      compute_gradients = True)
    if verbose: print("  phenix target:   %10.16f" %geometry.target)
    phts.append(geometry.target)
    geometry.complex_residual_sum = geometry.residual_sum
    geometry.complex_gradients = copy.deepcopy(geometry.gradients)
    for resname_i,resname in enumerate(afitt_o.resname):
      for instance_i, instance in enumerate(afitt_o.res_ids[resname_i]):
        afitt_input = afitt_o.make_afitt_input(sites_cart,
                                               resname_i,
                                               instance_i)
        lines = call_afitt(afitt_input,
                           afitt_o.ff)
        process_afitt_output(
            lines, geometry, afitt_o,
            resname_i, instance_i, afitt_allgradients, afitt_alltargets)
    if verbose: print("  afitt target:    %10.16f" %afitt_alltargets[(0,0)])
    afts.append(afitt_alltargets[(0,0)])
    geometry = apply_target_gradients(
        afitt_o, geometry, afitt_allgradients, afitt_alltargets)
    geometry.afitt_residual_sum = geometry.residual_sum -\
                                geometry.complex_residual_sum
    grm.afitt_object = afitt_o
    geometry = adjust_energy_and_gradients(
      geometry,
      grm,
      sites_cart,
      xrs.hd_selection(),
      afitt_o,
      )
    geometry.target = geometry.residual_sum

    if verbose: print("  final target:    %10.16f" %geometry.target)
    t=geometry.target
    ts.append(t)
  if verbose: print("  phenix finite diff.: %10.16f" %((phts[0]-phts[1])/(2*e)))
  if verbose: print("  afitt finite diff.: %10.16f" %((afts[0]-afts[1])/(2*e)))
  num_gradient = (ts[0]-ts[1])/(2*e)
  print("-> %10.9f" %(num_gradient))
  gradient_diff = num_gradient - ana_gradient
  assert abs(gradient_diff) <= 1e-4, \
    "TEST FAILS: (analytical - numerical)= %10.9f" %(gradient_diff)
  print("TEST PASSES: (analytical - numerical)= %10.9f" %(gradient_diff))
  return 0

def apply(result, afitt_o, sites_cart,phenix_gnorms=None):
  result.complex_residual_sum = result.geometry.residual_sum
  # needs to be more selective!!!
  result.complex_gradients = copy.deepcopy(result.geometry.gradients)
  afitt_allgradients = {}
  afitt_alltargets = {}
  for resname_i,resname in enumerate(afitt_o.resname):
    for instance_i, instance in enumerate(afitt_o.res_ids[resname_i]):
      afitt_input = afitt_o.make_afitt_input(sites_cart,
                                                  resname_i,
                                                  instance_i,
      )
      lines = call_afitt(afitt_input, afitt_o.ff)
      process_afitt_output(lines,
                           result.geometry,
                           afitt_o,
                           resname_i,
                           instance_i,
                           afitt_allgradients,
                           afitt_alltargets,
                           verbose=False,
                           phenix_gnorms=phenix_gnorms)
  result.geometry = apply_target_gradients(afitt_o, result.geometry,
                                           afitt_allgradients,
                                           afitt_alltargets)

  # used as a trigger for adjust the energy and gradients
  result.afitt_residual_sum = result.geometry.residual_sum -\
                              result.complex_residual_sum
  #print 'Afitt'
  #print afitt_o
  return result

def bond_test(model):
  rm = model.restraints_manager
  bond_params_table = rm.geometry.bond_params_table
  bond = bond_params_table.lookup(0,1)
  print(bond.distance_ideal,bond.weight)
  bond = bond_params_table.lookup(0,10)
  print(bond)
  assert 0

def run(pdb_file, cif_file, ligand_names, ff='mmff94s',covalent=False):
  import iotbx.pdb
  assert os.path.isfile(pdb_file), "File %s does not exist." %pdb_file
  assert os.path.isfile(cif_file), "File %s does not exist." %cif_file
  pdb_inp = iotbx.pdb.input(file_name=pdb_file)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  pdb_hierarchy.atoms().reset_i_seq()
  xrs = pdb_hierarchy.extract_xray_structure()
  sites_cart=xrs.sites_cart()
  grm=None
  if covalent:
    from mmtbx import monomer_library
    import mmtbx.monomer_library.server
    import mmtbx.monomer_library.pdb_interpretation
    mon_lib_srv = monomer_library.server.server()
    ener_lib = monomer_library.server.ener_lib()
    processed_pdb_file = monomer_library.pdb_interpretation.process(
      mon_lib_srv    = mon_lib_srv,
      ener_lib       = ener_lib,
      file_name      = pdb_file,
      raw_records    = None,
      force_symmetry = True)
    grm = processed_pdb_file.geometry_restraints_manager(
      show_energies = False,
      plain_pairs_radius = 5.0,
      )
  energies = get_afitt_energy([cif_file],
                              ligand_names,
                              pdb_hierarchy,
                              ff,
                              sites_cart,
                              grm)

  for energy in energies:
    print("%s_%d_%s AFITT_ENERGY: %10.4f" %(energy[0], energy[1], energy[2], energy[3]))

def run2():
  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument("pdb_file", help="pdb file")
  parser.add_argument("cif_file", help="cif file", default=0)
  parser.add_argument("ligand_names", help="3-letter ligand names separated by commas")
  parser.add_argument("-ff", help="afitt theory: mmff94, mmff94s pm3 or am1", default='mmff94s')
  parser.add_argument('-covalent', dest='covalent', action='store_true', help="calculate covalent energy (only for debugging)")
  args = parser.parse_args()
  ligand_names=args.ligand_names.split(',')
  run(args.pdb_file, args.cif_file, ligand_names, args.ff, args.covalent)

if (__name__ == "__main__"):
  run2()
