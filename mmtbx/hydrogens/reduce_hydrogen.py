from __future__ import absolute_import, division, print_function
import six
import sys, time
from libtbx.utils import Sorry
import mmtbx.model
import iotbx.pdb
import boost_adaptbx.boost.python as bp
from libtbx.utils import null_out
from libtbx import group_args
from scitbx import matrix
from cctbx.array_family import flex
from mmtbx.ligands.ready_set_utils import add_n_terminal_hydrogens_to_residue_group
from cctbx.geometry_restraints.linking_class import linking_class
#
from cctbx.maptbx.box import shift_and_box_model

ext = bp.import_ext("cctbx_geometry_restraints_ext")
get_class = iotbx.pdb.common_residue_names_get_class

# ==============================================================================

def get_h_restraints(resname, strict=True):
  from mmtbx.monomer_library import cif_types
  from mmtbx.chemical_components import get_cif_dictionary
  from mmtbx.ligands.rdkit_utils import get_molecule_from_resname
  molecule = get_molecule_from_resname(resname)
  if molecule is None: return None
  cc_cif = get_cif_dictionary(resname)
  cc = cc_cif['_chem_comp'][0]
  hs = []
  hsi = []
  chem_comp = cif_types.chem_comp(
    id=cc.id,
    three_letter_code=cc.three_letter_code,
    name=cc.name,
    group=cc.type,
    number_atoms_all=0, #cc.number_atoms_all,
    number_atoms_nh=0, #cc.number_atoms_nh,
    desc_level=".")
  comp_comp_id = cif_types.comp_comp_id(source_info=None, chem_comp=chem_comp)
  lookup = {}
  for i, a in enumerate(cc_cif.get('_chem_comp_atom',[])):
    lookup[a.atom_id]=i
    lookup[i]=a.atom_id
    if a.type_symbol in ['H', 'D']:
      hs.append(a.atom_id)
      hsi.append(i)
    comp_comp_id.atom_list.append(cif_types.chem_comp_atom(
      atom_id=a.atom_id,
      type_symbol=a.type_symbol,
      # type_energy=a.type_energy,
      # partial_charge=a.partial_charge,
      ))
  conf =  molecule.GetConformer()
  from rdkit import Chem # needed import
  from rdkit.Chem import rdMolTransforms
  for b in cc_cif.get('_chem_comp_bond',[]):
    if strict:
      if (b.atom_id_1 not in hs and
          b.atom_id_2 not in hs): continue
    if ( b.atom_id_1 not in lookup or
         b.atom_id_2 not in lookup): continue
    atom_idx1=lookup[b.atom_id_1]
    atom_idx2=lookup[b.atom_id_2]
    bl = rdMolTransforms.GetBondLength(conf, atom_idx1, atom_idx2)
    comp_comp_id.bond_list.append(cif_types.chem_comp_bond(
      atom_id_1=b.atom_id_1,
      atom_id_2=b.atom_id_2,
      type=b.value_order,
      value_dist='%0.3f' % (bl*.9),
      value_dist_esd=".1"))

  from mmtbx.ligands.rdkit_utils import enumerate_angles
  for angle in enumerate_angles(molecule):
    if strict:
      if angle[0] in hsi or angle[2] in hsi:
        av = rdMolTransforms.GetAngleDeg(conf, angle[0], angle[1], angle[2])
      else: continue
    else:
      av = rdMolTransforms.GetAngleDeg(conf, angle[0], angle[1], angle[2])
    if ( angle[0] not in lookup or
         angle[2] not in lookup): continue
    comp_comp_id.angle_list.append(cif_types.chem_comp_angle(
      atom_id_1=lookup[angle[0]],
      atom_id_2=lookup[angle[1]],
      atom_id_3=lookup[angle[2]],
      value_angle='%0.1f' % av,
      value_angle_esd="1"))

  from mmtbx.ligands.rdkit_utils import enumerate_torsions
  for i, angle in enumerate(enumerate_torsions(molecule)):
    if strict:
      if angle[0] in hsi or angle[3] in hsi:
        av = rdMolTransforms.GetDihedralDeg(conf, angle[0], angle[1], angle[2], angle[3])
      else: continue
    else:
      av = rdMolTransforms.GetDihedralDeg(conf, angle[0], angle[1], angle[2], angle[3])
    if ( angle[0] not in lookup or
         angle[3] not in lookup): continue
    comp_comp_id.tor_list.append(cif_types.chem_comp_tor(
      id='Var_%03d' % i,
      atom_id_1=lookup[angle[0]],
      atom_id_2=lookup[angle[1]],
      atom_id_3=lookup[angle[2]],
      atom_id_4=lookup[angle[3]],
      value_angle='%0.1f' % av,
      value_angle_esd='1',
      period='1'))
  return comp_comp_id

def atom_in_restraints(name, cc_cif):
  for a in cc_cif.get('_chem_comp_atom',[]):
    if a.atom_id==name:
      return a
  return None

def bonds_in_restraints(atom, exclude_hydrogens=False):
  from mmtbx.chemical_components import get_cif_dictionary
  cc_cif = get_cif_dictionary(atom.parent().resname)
  rc=[]
  for b in cc_cif.get('_chem_comp_bond',[]):
    if b.atom_id_1.strip()==atom.name.strip():
      if exclude_hydrogens:
        a=atom_in_restraints(b.atom_id_2, cc_cif)
        if a.type_symbol in ['H', 'D']: continue
      rc.append(b.atom_id_2)
    if b.atom_id_2.strip()==atom.name.strip():
      if exclude_hydrogens:
        a=atom_in_restraints(b.atom_id_1, cc_cif)
        if a.type_symbol in ['H', 'D']: continue
      rc.append(b.atom_id_1)
  return rc

# ==============================================================================

def mon_lib_query(residue, mon_lib_srv, construct_h_restraints=True, raise_sorry=True):
  # if get_class(residue.resname) in ['common_rna_dna']:
  #   md = get_h_restraints(residue.resname)
  #   return md
  # if print_time: print(residue.resname, get_class(residue.resname))
  if residue.resname == 'UNL':
    return None, None
  md, ani = mon_lib_srv.get_comp_comp_id_and_atom_name_interpretation(
    residue_name=residue.resname,
    atom_names=residue.atoms().extract_name())
  cif_object=None
  # if md is None:
  #   md, ani = mon_lib_srv.get_comp_comp_id_and_atom_name_interpretation(
  #     residue_name='%s_EL' % residue.resname,
  #     atom_names=residue.atoms().extract_name(),
  #     ad_hoc_single_atom_residues=True)
  if md is None:
    md = get_h_restraints(residue.resname, strict=False)
    if md is None:
      if raise_sorry:
        raise Sorry('Entity "%s" not found in CCD (or GeoStd). Please supply restraints.' % residue.resname)
      else:
        return None, None
    from six.moves import cStringIO as StringIO
    input_string='data_comp_list\n'
    input_string+=str(md.chem_comp.as_cif_loop())
    f=StringIO()
    md.show(f=f)
    # use strip in case 3-letter code has only 2 letters (e.g. DI)
    input_string += '\ndata_comp_%s\n' % residue.resname.strip()
    input_string += '\n%s' % f.getvalue()
    cif_object = iotbx.cif.reader(input_string=input_string).model()
  return md, cif_object

# ==============================================================================

def get_reduce_pdb_interpretation_params(use_neutron_distances):
  '''
  Create pdb_interpretation parameter scope.
  Do this in a function so other programs (reduce2) can use the same parameters
  '''
  p = mmtbx.model.manager.get_default_pdb_interpretation_params()
  p.pdb_interpretation.restraints_library.cdl=False # XXX this triggers a bug !=360
  p.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
  p.pdb_interpretation.disable_uc_volume_vs_n_atoms_check=True
  p.pdb_interpretation.use_neutron_distances = use_neutron_distances
  p.pdb_interpretation.proceed_with_excessive_length_bonds=True
  p.pdb_interpretation.allow_polymer_cross_special_position=True
  #p.pdb_interpretation.automatic_linking.link_metals = True
  p.pdb_interpretation.automatic_linking.link_residues = True
  p.pdb_interpretation.automatic_linking.exclude_hydrogens_from_bonding_decisions = True
  #
  return p

# ==============================================================================

class place_hydrogens():
  '''
  Add H atoms to a model

  Parameters
  ----------
  use_neutron_distances : bool
    use neutron distances instead of X-ray

  adp_scale : float
    scale factor for isotropic B of H atoms.
    B(H-atom) = adp_scale * B(parent non-H atom)

  keep_existing_H : bool
    keep existing H atoms in model, only place missing H
  '''

# ------------------------------------------------------------------------------

  def __init__(self,
               model,
               use_neutron_distances = False,
               n_terminal_charge     = 'residue_one',
               adp_scale             = 1,
               exclude_water         = True,
               stop_for_unknowns     = False,
               keep_existing_H       = False,
               validate_e            = False,
               print_time            = False):

    self.model                 = model
    self.use_neutron_distances = use_neutron_distances
    self.n_terminal_charge     = n_terminal_charge
    self.adp_scale             = adp_scale
    self.exclude_water         = exclude_water
    self.stop_for_unknowns     = stop_for_unknowns
    self.keep_existing_H       = keep_existing_H
    self.validate_e            = validate_e
    self.print_time            = print_time
    #
    self.no_H_placed_mlq        = list()
    self.site_labels_disulfides = list()
    self.site_labels_no_para    = list()
    #self.charged_atoms          = list()
    self.sl_removed             = list()
    self.n_H_initial            = 0
    self.n_H_final              = 0

    if self.print_time:
      self.time_rebox_model        = None
      self.time_remove_element_X   = None
      self.time_add_missing_H      = None
      self.time_terminal_propeller = None
      self.time_make_grm           = None
      self.time_remove_isolated    = None
      self.time_riding_manager     = None
      self.time_remove_H_nopara    = None
      self.time_reset              = None
      self.time_idealize           = None
      self.time_remove_H_on_links  = None

# ------------------------------------------------------------------------------

  def run(self):
    '''
    Function that places H atoms
    '''

    # Create symmetry if necessary
    # ------------------------------
    model_has_bogus_cs = False
    t0 = time.time()
    cs = self.model.crystal_symmetry()
    if (cs is None) or (cs.unit_cell() is None):
      self.model = shift_and_box_model(model = self.model)
      model_has_bogus_cs = True
      #self.model.add_crystal_symmetry_if_necessary() # this is slower than shift_and_box_model!!!!
    self.time_rebox_model = round(time.time()-t0, 2)


    # Don't stop if model contains element X atoms
    # This needs more discussion before being made final
    # ---------------------------------------------
    t0 = time.time()
    if ' X' in self.model.get_hierarchy().atoms().extract_element():
      self.model = self.model.select(~self.model.selection('element X'))
    self.time_remove_element_X = round(time.time()-t0, 2)


    # Remove existing H if requested
    # ------------------------------
    self.model.get_xray_structure()
    self.n_H_initial = self.model.get_hd_selection().count(True)
    if not self.keep_existing_H:
      self.model = self.model.select(~self.model.get_hd_selection())

    # Add missing H atoms and place them at bogus position
    # ----------------------------------------------------
    t0 = time.time()
    pdb_hierarchy = self.add_missing_H_atoms_at_bogus_position()
    self.time_add_missing_H = round(time.time()-t0, 2)
    # DEBUG
    #print(pdb_hierarchy.composition().n_hd)
    #f = open("intermediate1.pdb","w")
    #f.write(self.model.model_as_pdb())

    # Place N-terminal propeller hydrogens
    # TODO double check N-terminal position for PRO residues
    # ------------------------------------
    if self.n_terminal_charge in ['residue_one', 'first_in_chain']:
      t0 = time.time()
      self.place_n_terminal_propeller(pdb_hierarchy = pdb_hierarchy)
      self.time_terminal_propeller = round(time.time()-t0, 2)

    pdb_hierarchy.sort_atoms_in_place()
    pdb_hierarchy.atoms().reset_serial()

    # DEBUG
    #f = open("intermediate2.pdb","w")
    #f.write(self.model.model_as_pdb())

    # Make new model obj and get restraints manager
    # ---------------------------------------------
    p = get_reduce_pdb_interpretation_params(self.use_neutron_distances)
    ro = self.model.get_restraint_objects()
    t0 = time.time()
    self.model = mmtbx.model.manager(
      model_input       = None,
      pdb_hierarchy     = pdb_hierarchy,
      stop_for_unknowns = self.stop_for_unknowns,
      crystal_symmetry  = self.model.crystal_symmetry(),
      restraint_objects = ro,
      log               = null_out())
    self.model.process(pdb_interpretation_params=p,
                       make_restraints=True,
                       # retain_zero_dihedrals=True,
                       )
    #self.model.idealize_h_minimization()
    #STOP()
    self.time_make_grm = round(time.time()-t0, 2)

    #f = open("intermediate3.pdb","w")
    #f.write(self.model.model_as_pdb())

    # Return if no H have been placed
    sel_h = self.model.get_hd_selection()
    if sel_h.count(True) == 0: return

    # Remove isolated H atoms
    # -----------------------
    # (when heavy atom is missing, H needs not to be placed)
    t0 = time.time()
    sel_isolated = self.model.isolated_atoms_selection()
    sel_lone_H = sel_h & sel_isolated
    # As h_parameterization will not include these, they can be removed in the
    # next step; for book-keeping it is useful to keep track of lone H as a
    # selection
    #if not sel_lone_H.all_eq(False):
    #  self.model = self.model.select(~sel_lone_H)
    self.time_remove_isolated = round(time.time()-t0, 2)

    sel_h = self.model.get_hd_selection()

    #f = open("intermediate3a.pdb","w")
    #f.write(self.model.model_as_pdb())

    # Setup riding H manager
    # ----------------------
    t0 = time.time()
    self.model.setup_riding_h_manager(use_ideal_dihedral = True)
    riding_h_manager = self.model.riding_h_manager
    if riding_h_manager is None:
      return
    self.time_riding_manager = round(time.time()-t0, 2)

    # Remove H that could not be parameterized
    # ----------------------------------------
    t0 = time.time()
    sel_h_in_para = flex.bool(
      [bool(x) for x in riding_h_manager.h_parameterization])
    sel_h_not_in_para = sel_h_in_para.exclusive_or(sel_h)
    # no need to display lone H atoms in the log, so remove from labels
    sel_h_not_in_para_but_not_lone = sel_h_not_in_para.exclusive_or(sel_lone_H)
    self.site_labels_no_para = [atom.id_str().replace('pdb=','').replace('"','')
      for atom in self.model.get_hierarchy().atoms().select(sel_h_not_in_para_but_not_lone)]
    if not sel_h_not_in_para.all_eq(False):
      self.model = self.model.select(~sel_h_not_in_para)
    self.time_remove_H_nopara = round(time.time()-t0, 2)

    #f = open("intermediate4.pdb","w")
    #f.write(self.model.model_as_pdb())

# to be removed; was for curiosity only
#    if self.validate_e:
#      t0 = time.time()
#      self.validate_electrons()
#      if self.print_time:
#        print("validate electrons:", round(time.time()-t0, 2))

    # Reset occupancies, ADPs and idealize H atom positions
    # -----------------------------------------------------
    t0 = time.time()
    self.model.reset_adp_for_hydrogens(scale = self.adp_scale, keep_aniso=True)
    self.model.reset_occupancy_for_hydrogens_simple()
    self.time_reset = round(time.time()-t0, 2)
    t0 = time.time()
    self.model.idealize_h_riding()
    self.time_idealize = round(time.time()-t0, 2)

    # Remove H atoms that are involved in links (bonds, metal coordination, etc)
    # --------------------------------------------------------------------------
    t0 = time.time()
    self.exclude_H_on_links()
    self.time_remove_H_on_links = round(time.time()-t0, 2)


    # TODO: this should be ideally done *after* reduce optimization
    #if not self.exclude_water:
    #  self.model.add_hydrogens(1., occupancy=0.)

    self.n_H_final = self.model.get_hd_selection().count(True)

    if self.print_time:
      self.print_times()

  # ----------------------------------------------------------------------------

  def place_n_terminal_propeller(self, pdb_hierarchy):
    '''
    Place NH3 at residue #1 or at first residue in chain
    Changes hierarchy in place
    '''
    for m in pdb_hierarchy.models():
      for chain in m.chains():
        rgs = chain.residue_groups()[0]
        # by default, place NH3 only at residue with resseq 1
        if (self.n_terminal_charge == 'residue_one' and rgs.resseq_as_int() != 1):
          continue
        elif (self.n_terminal_charge == 'first_in_chain'):
          pass
        for ag in rgs.atom_groups():
          # SAC in 5xdq, 5zcp. Never needs propeller. Also AYA
          for ag in rgs.atom_groups():
            n=ag.get_atom('N') # assumes atom name "N"
            if n: break
          if not n: continue
          bonds=bonds_in_restraints(n, exclude_hydrogens=True)
          heavies=2
          if ag.resname in ['PRO']: # needs a PRO child lookup
            heavies=3
          if len(bonds)>=heavies: continue
          if (get_class(name=ag.resname) in
              ['common_amino_acid', 'modified_amino_acid', 'd_amino_acid']):
            if ag.get_atom('H'):
              ag.remove_atom(ag.get_atom('H'))
          # TODO make the function below smart, so it
          # 1) knows when to add H1H2H3 or not
          # 2) renames H to H1 (so no need to remove it beforehand)
          rc = add_n_terminal_hydrogens_to_residue_group(rgs) # rc is always empty list?

  # ----------------------------------------------------------------------------

  def add_missing_H_atoms_at_bogus_position(self):
    '''
    Add missing H atoms at bogus positions to the pdb_hierarchy

    This procedure changes the hierarchy in place.
    All H atoms are placed at center of coordinates + (0.5, 0.5, 0.5)
    The translation is necessary because sometimes the center of coordinates
    coincides with the position of a heavy atom.

    In one residue/entity, all newly placed H are superposed, they will be
    moved to their expected position later.
    '''
    # TODO temporary fix until v3 names are in mon lib
    alternative_names = [
      ('HA1', 'HA2', 'HA3'),
      ('HB1', 'HB2', 'HB3'),
      ('HG1', 'HG2', 'HG3'),
      ('HD1', 'HD2', 'HD3'),
      ('HE1', 'HE2', 'HE3'),
      ('HG11', 'HG12', 'HG13')
      ]
    # end TODO
    pdb_hierarchy = self.model.get_hierarchy()
    mon_lib_srv = self.model.get_mon_lib_srv()
    #XXX This breaks for 1jxt, residue 2, TYR
    for m in pdb_hierarchy.models():
      for chain in m.chains():
        for rg in chain.residue_groups():
          n_atom_groups = len(rg.atom_groups())
          for ag in rg.atom_groups():
            if n_atom_groups > 2 and ag.altloc == '':
              continue
            #print list(ag.atoms().extract_name())
            if(get_class(name=ag.resname) == "common_water"): continue
            actual = [a.name.strip().upper() for a in ag.atoms()]
            #
            mlq, cif_object = mon_lib_query(residue=ag, mon_lib_srv=mon_lib_srv, raise_sorry=False)
            if mlq is None:
              self.no_H_placed_mlq.append(ag.resname)
              continue

            if cif_object:
              ro = self.model.get_restraint_objects()
              if ro is None: ro=[]
              ro.append(('auto_%s' % ag.resname, cif_object))
              self.model.set_restraint_objects(ro)

            expected_h = []
            #expected_ha = []
            atom_dict = mlq.atom_dict()
            def _remove_atoms(atom_dict, names):
              remove=[]
              for k,v in atom_dict.items():
                if k in names:
                  remove.append(k)
              if remove:
                for r in remove:
                  del atom_dict[r]
              return atom_dict
            #
            # don't add polymer H atoms. Terminal H atoms added elsewhere
            #
            if mlq.test_for_peptide(atom_dict):
              atom_dict = _remove_atoms(atom_dict, ['H2', 'HXT'])
            elif mlq.test_for_rna_dna(atom_dict):
              atom_dict = _remove_atoms(atom_dict, ["HO3'", 'HO3*'])
            for k, v in six.iteritems(atom_dict):
              if(v.type_symbol=="H"):
                expected_h.append(k)
              #else:
              #  expected_ha.append(k)
            #print('expected H', expected_h)
            #
            # TODO start
            # temporary fix until v3 names are in mon lib
            if (get_class(name=ag.resname) in
                ['common_amino_acid', 'modified_amino_acid', 'd_amino_acid']):
              for altname in alternative_names:
                if (altname[0] in expected_h and altname[1] in expected_h):
                  if (atom_dict[altname[0]].type_energy == 'HCH2' and
                      atom_dict[altname[1]].type_energy == 'HCH2'):
                    expected_h.append(altname[2])
                    expected_h.remove(altname[0])
                    #print('renamed %s to %s' % (altname[0], altname[2]))
            # TODO end
            missing_h = list(set(expected_h).difference(set(actual)))
            if 0: print(ag.resname, missing_h)
            #new_xyz = ag.atoms().extract_xyz().mean()
            new_xyz = flex.double(ag.atoms().extract_xyz().mean()) + \
              flex.double([0.5,0.5,0.5])
            new_xyz = tuple(new_xyz)

            hetero = ag.atoms()[0].hetero
            segid = ag.atoms()[0].segid

            for mh in missing_h:
              if len(mh) < 4: mh = (' ' + mh).ljust(4)
              a = (iotbx.pdb.hierarchy.atom()
                .set_name(new_name=mh)
                .set_element(new_element="H")
                .set_xyz(new_xyz=new_xyz)
                .set_hetero(new_hetero=hetero)
                .set_segid(new_segid=segid))

              ag.append_atom(a)

    return pdb_hierarchy

# ------------------------------------------------------------------------------
#
#  def validate_electrons(self):
#    from elbow.quantum import electrons
#    atom_valences = electrons.electron_distribution(
#      self.model.get_hierarchy(), # needs to be altloc free
#      self.model.get_restraints_manager().geometry,
#      verbose=False,
#    )
#    atom_valences.validate(ignore_water=True, raise_if_error=False)
#    self.charged_atoms = atom_valences.get_charged_atoms()

# ------------------------------------------------------------------------------

  def exclude_H_on_links(self, verbose=False):
    """Remove H atoms bound to heavy atoms that form a link

    An exception are HD1 and HE2 of HIS. The mover functionality in reduce will
    take care of those.

    TODO: Could restraints manager have a list of links with relevant information?
          Then we don't have to loop through all proxies here.
    """
    from mmtbx.ligands.chemistry import get_valences
    origin_ids = linking_class()
    grm = self.model.get_restraints_manager()
    bond_proxies_simple, asu = grm.geometry.get_all_bond_proxies(
      sites_cart = self.model.get_sites_cart())
    atoms = self.model.get_atoms()
    elements = self.model.get_hierarchy().atoms().extract_element()
    exclusion_iseqs = list()
    exclusion_dict = dict()
    all_proxies = [p for p in bond_proxies_simple]
    for proxy in asu:
      all_proxies.append(proxy)
    # Loop through bond proxies to find links (i.e. proxies with origin_id != 0)
    for proxy in all_proxies:
      if(  isinstance(proxy, ext.bond_simple_proxy)): i,j=proxy.i_seqs
      elif(isinstance(proxy, ext.bond_asu_proxy)):    i,j=proxy.i_seq,proxy.j_seq
      else: assert 0 # never goes here
      if proxy.origin_id != 0:
        exclusion_iseqs.extend([i,j])
        exclusion_dict[i] = proxy.origin_id
        exclusion_dict[j] = proxy.origin_id
    sel_remove = flex.size_t()

    # Find H atoms bound to linked atoms
    removed_dict = {}
    parent_dict = {}
    bonds = {}
    bond_lengths = {}
    for proxy in all_proxies:
      if(  isinstance(proxy, ext.bond_simple_proxy)): i,j=proxy.i_seqs
      elif(isinstance(proxy, ext.bond_asu_proxy)):    i,j=proxy.i_seq,proxy.j_seq
      else: assert 0 # never goes here
      bonds.setdefault(i,[])
      bonds[i].append(j)
      bonds.setdefault(j,[])
      bonds[j].append(i)
      # Exception for HIS HD1 and HE2
      if (atoms[i].parent().resname == 'HIS' and
        atoms[i].name.strip() in ['HD1','DD1', 'HE2', 'DE2']): continue
      if (atoms[j].parent().resname == 'HIS' and
        atoms[j].name.strip() in ['HD1','DD1', 'HE2', 'DE2']): continue
      if(elements[i] in ["H","D"] and j in exclusion_iseqs):
        if i not in sel_remove:
          sel_remove.append(i)
          bond_lengths[i] = proxy.distance_ideal
          removed_dict[i] = exclusion_dict[j]
          parent_dict[i]=j
      if(elements[j] in ["H","D"] and i in exclusion_iseqs):
        if j not in sel_remove:
          sel_remove.append(j)
          bond_lengths[j] = proxy.distance_ideal
          removed_dict[j] = exclusion_dict[i]
          parent_dict[j]=i
    # remove H atoms NOT to remove - double negative!
    #verbose=True
    if verbose:
      print('removed_dict',removed_dict)
      for i_seq in sel_remove:
        print('remove?',atoms[i_seq].quote())
    remove_from_sel_remove=[]
    for ii, i_seq in reversed(list(enumerate(sel_remove))):
      j_seq=parent_dict[i_seq]
      # need to add the use of atomic charge
      valences=get_valences(elements[j_seq])
      number_of_bonds=len(bonds[j_seq])
      if number_of_bonds in valences:
        # remove this H from delection
        remove_from_sel_remove.append(i_seq) # ??
        del removed_dict[i_seq]
      else:
        bonds[j_seq].remove(i_seq)
        bonds[i_seq].remove(j_seq)

    #print(remove_from_sel_remove)
    fsc0=grm.geometry.shell_sym_tables[0].full_simple_connectivity()
    fsc1=grm.geometry.shell_sym_tables[1].full_simple_connectivity()
    #fsc2=grm.geometry.shell_sym_tables[2].full_simple_connectivity()
    for _i in remove_from_sel_remove:
      #print(list(fsc0[_i]))
      parent = fsc0[_i][0]
      #print('origin_id', exclusion_dict[parent], origin_ids.get_origin_key(exclusion_dict[parent]))
      first_neighbors = fsc1[_i]
      fn_filtered = [item for item in first_neighbors if item not in sel_remove]
      #print(list(first_neighbors), list(sel_remove), fn_filtered )
      # now improve geometry of the H being kept
      # TODO make sure all atoms are in same conformer
      # TODO check that neighbor atoms are all non H
      # TODO what about two tetrahedral geometry?
      if len(fn_filtered) == 3:
        #print('tetrahedral geometry')
        coord1 = matrix.col(atoms[fn_filtered[0]].xyz)
        coord2 = matrix.col(atoms[fn_filtered[1]].xyz)
        coord3 = matrix.col(atoms[fn_filtered[2]].xyz)
        coordp = matrix.col(atoms[parent].xyz)
        orth = (coord2-coord1).cross(coord3-coord1).normalize()
        vol = orth.dot(coordp-coord1)
        if vol > 0: orth = -orth
        atoms[_i].xyz = coordp - orth*bond_lengths[_i]
      if len(fn_filtered) == 2:
        #print('flat geometry')
        coord1 = matrix.col(atoms[fn_filtered[0]].xyz)
        coord2 = matrix.col(atoms[fn_filtered[1]].xyz)
        coordp = matrix.col(atoms[parent].xyz)
        half = ((coord1 - coordp).normalize() + (coord2 - coordp).normalize())
        atoms[_i].xyz = coordp-half.normalize()*bond_lengths[_i]

    if remove_from_sel_remove:
      sel_remove=list(sel_remove)
      for r in remove_from_sel_remove:
        sel_remove.remove(r)
        if verbose: print('keep',atoms[r].quote())
      sel_remove=flex.size_t(sel_remove)
    #
    sl_removed = [(atom.id_str().replace('pdb=','').replace('"',''),
                   origin_ids.get_origin_key(removed_dict[atom.i_seq]))
        for atom in self.model.get_hierarchy().atoms().select(sel_remove)]
    #
    if sel_remove:
      self.model = self.model.select(~flex.bool(self.model.size(), sel_remove))
    self.sl_removed = sl_removed
    self.exclusion_iseqs = exclusion_iseqs


# ------------------------------------------------------------------------------

  def show(self, log):
    '''
    Informative output
    '''
    if log is None: log = sys.stdout
    #
    if (not self.keep_existing_H and self.n_H_initial):
      msg = 'Number of hydrogen atoms trimmed from input model: %s \n'
      print(msg % self.n_H_initial, file=log)
    #
    msg = 'Number of hydrogen atoms added to the input model: %s \n'
    print(msg % self.n_H_final, file=log)
    #
    if self.no_H_placed_mlq:
      msg = '''
No H atoms were placed on the following residues because no restraints
were found:'''
      print(msg, file=log)
      for resname in self.no_H_placed_mlq:
        print(resname, file=log)
    #
    if self.site_labels_disulfides:
      msg = '''
The following cysteine HG atoms were not placed because the sulfur atom is
involved in a disulfide bond'''
      print(msg, file=log)
      for label in self.site_labels_disulfides:
        print(label, file=log)
    #
    if self.site_labels_no_para:
      msg = '''
The following H atoms were not placed because they could not be parameterized
(not enough restraints information)'''
      print(msg, file=log)
      for label in self.site_labels_no_para:
        print(label, file=log)
#    if self.charged_atoms:
#      msg = '''
#The following heavy atom have an unusual electron count. This could be because
#heavy atoms or H atoms are missing.'''
#      print(msg, file=log)
#      for item in self.charged_atoms:
#        idstr = item[0].id_str().replace('pdb=','').replace('"','')
#        if 'HOH' in idstr: continue
#        print(idstr, item[1])

    if self.sl_removed:
      print()
      msg = '''Atom %s was not placed because it is involved in %s'''
      for item in self.sl_removed:
        print(msg % (item[0], item[1]), file=log)

# ------------------------------------------------------------------------------

  def get_model(self):
    return self.model

# ------------------------------------------------------------------------------

  def get_counts(self):
    return group_args(
      number_h_final  = self.n_H_final,
      no_H_placed_mlq = self.no_H_placed_mlq,
      site_labels_disulfides = self.site_labels_disulfides,
      site_labels_no_para = self.site_labels_no_para)

# ------------------------------------------------------------------------------

  def get_times(self):
    return group_args(
      time_rebox_model        = self.time_rebox_model,
      time_remove_element_X   = self.time_remove_element_X,
      time_add_missing_H      = self.time_add_missing_H,
      time_terminal_propeller = self.time_terminal_propeller,
      time_make_grm           = self.time_make_grm,
      time_remove_isolated    = self.time_remove_isolated,
      time_riding_manager     = self.time_riding_manager,
      time_remove_H_nopara    = self.time_remove_H_nopara,
      time_reset              = self.time_reset,
      time_idealize           = self.time_idealize,
      time_remove_H_on_links  = self.time_remove_H_on_links)

# ------------------------------------------------------------------------------

  def print_times(self):
    print('Detailed timings:')
    print("Rebox model:", self.time_rebox_model)
    print('Remove element X:', self.time_remove_element_X)
    print("Add missing H at bogus position:", self.time_add_missing_H)
    print('Add N-terminal propeller:', self.time_terminal_propeller)
    print("Get new model obj and grm:", self.time_make_grm )
    print("Remove isolated H:", self.time_remove_isolated)
    print("Setup Riding manager:", self.time_riding_manager)
    print("Remove H that were not parameterized:", self.time_remove_H_nopara)
    print("Reset adp, occ:", self.time_reset)
    print("idealize H positions:", self.time_idealize)
    print("Remove H on links:", self.time_remove_H_on_links)
    print()

# ==============================================================================

## stub for reduce parameters
## TODO can be parameters or phil, depending on how many options are really needed
#reduce_master_params_str = """
#flip_NQH = True
#  .type = bool
#  .help = add H and rotate and flip NQH groups
#search_time_limit = 600
#  .type = int
#  .help = max seconds to spend in exhaustive search (default=600)
#"""
#
#def optimize(model):
#  """
#  Carry out reduce optimization
#
#  Parameters
#  ----------
#  model
#      mmtbx model object that contains H atoms
#      H atoms should be at approprite distances
#
#  Returns
#  -------
#  model
#      mmtbx model object with optimized H atoms
#  """
#  # hierarchy object --> has hierarchy of structure
#  pdb_hierarchy = model.get_hierarchy()
#  # geometry restraints manager --> info about ideal bonds, angles; what atoms are bonded, etc.
#  grm = model.get_restraints_manager()
#
#  print("Reduce optimization happens here")
#
#  return model
