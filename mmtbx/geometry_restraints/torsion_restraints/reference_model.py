from __future__ import absolute_import, division, print_function
import cctbx.geometry_restraints
from mmtbx.validation import rotalyze
from mmtbx.utils import rotatable_bonds
from mmtbx.rotamer.sidechain_angles import SidechainAngles
import mmtbx.monomer_library
from cctbx.array_family import flex
import iotbx.phil
import libtbx.load_env
from libtbx.utils import Sorry
from mmtbx import secondary_structure
from scitbx.matrix import rotate_point_around_axis
from libtbx.str_utils import make_sub_header
from mmtbx.geometry_restraints.torsion_restraints import utils
import sys
import time
import iotbx
from six.moves import zip
from six.moves import range

TOP_OUT_FLAG = True

renamed_ncs_search_str = iotbx.ncs.ncs_search_options.replace(
    "ncs_search", "search_options")
renamed_ncs_search_str = renamed_ncs_search_str.replace(
    "chain_max_rmsd = 2.", "chain_max_rmsd = 100.")
renamed_ncs_search_str = renamed_ncs_search_str.replace(
    "residue_match_radius = 4.0", "residue_match_radius = 1000")
# removing 'enabled' parameter
i1 = renamed_ncs_search_str.find("enabled")
i2 = renamed_ncs_search_str.find("exclude_selection")
renamed_ncs_search_str = renamed_ncs_search_str[:i1]+renamed_ncs_search_str[i2:]

reference_model_str = """
reference_model
    .caption = The reference torsion restraints are used to steer refinement \
      of the  working model.  This technique is advantageous in cases where \
      the working data set is low resolution, but there is a known related \
      structure solved at higher resolution.  The higher resolution \
      reference model is used to generate a set of dihedral restraints \
      that are applied to each matching dihedral in the working model. \
      To specify a PDB file as the reference model, add it to the list of \
      input files in the main window, then change the data type from \
      "Input model" to "Reference model".
    .style = box auto_align menu_item
    .short_caption = Reference model restraints
{
  enabled = False
    .short_caption = Reference model restraints
    .type = bool
    .help = Restrains the dihedral angles to a high-resolution reference \
      structure to reduce overfitting at low resolution.  You will need to \
      specify a reference PDB file (in the input list in the main window) \
      to use this option.
    .style = bold noauto
  file = None
    .type = path
    .short_caption = Reference model
    .style = bold file_type:pdb hidden
    .multiple = True
  use_starting_model_as_reference = False
    .type = bool
    .short_caption = use starting model as reference
  sigma = 1.0
    .type = float(value_min=0.001)
  limit = 15.0
    .type = float
  hydrogens = False
    .type = bool
    .help = Include dihedrals with hydrogen atoms
  main_chain = True
    .type = bool
    .help = Include dihedrals formed by main chain atoms
  side_chain = True
    .type = bool
    .help = Include dihedrals formed by side chain atoms
  fix_outliers = True
    .type = bool
    .help = Try to fix rotamer outliers in refined model
  strict_rotamer_matching = False
    .type = bool
    .help = Make sure that rotamers in refinement model matches those in \
      reference model even when they are not outliers
  auto_shutoff_for_ncs = False
    .type = bool
    .help = Do not apply to parts of structure covered by NCS restraints
  secondary_structure_only = False
    .type = bool
    .help = Only apply reference model restraints to secondary structure \
      elements (helices and sheets)
  reference_group
    .multiple=True
    .optional=True
    .short_caption=Reference group
    .style = noauto auto_align menu_item parent_submenu:reference_model
  {
    reference=None
      .type=atom_selection
      .short_caption=Selection in the reference model
    selection=None
      .type=atom_selection
      .short_caption=Selection in the refined model
    file_name=None
      .type=path
      .optional=True
      .short_caption = Reference model for this restraint group
      .style = bold hidden
      .help = this is to used internally to disambiguate cases where multiple \
              reference models contain the same chain ID. This normally does \
              not need to be set by the user
  }
  %s
}
""" % renamed_ncs_search_str

reference_model_params = iotbx.phil.parse(
    reference_model_str)

def add_reference_dihedral_restraints_if_requested(
    model,
    geometry,
    params=None,
    selection=None,
    log=None):
  if not params.enabled:
    return 0
  if (params.use_starting_model_as_reference and
    (len(params.file) > 0) and params.file[0] is not None):
    raise Sorry("Cannot not restrain working model to self and a "+
                    "reference model simultaneously")
  reference_file_list = set()
  reference_hierarchy_list = None
  if params.use_starting_model_as_reference:
    reference_hierarchy_list = [model.get_hierarchy()]
    reference_file_list = None
    print("*** Restraining model using starting model ***", file=log)
  else:
    for file_name in params.file:
      reference_file_list.add(file_name)
    for rg in params.reference_group:
      if rg.file_name is not None:
        reference_file_list.add(rg.file_name)
  print("*** Adding Reference Model Restraints (torsion) ***", file=log)
  #test for inserted TER cards in working model
  ter_indices = model._ter_indices
  if ter_indices is not None:
    utils.check_for_internal_chain_ter_records(
      pdb_hierarchy=model.get_hierarchy(),
      ter_indices=ter_indices)
  if reference_file_list is not None:
    reference_file_list = list(reference_file_list)
  rm = reference_model(
    model,
    reference_file_list=reference_file_list,
    reference_hierarchy_list=reference_hierarchy_list,
    params=params,
    selection=selection,
    log=log)
  rm.show_reference_summary(log=log)
  geometry.adopt_reference_dihedral_manager(rm)

class reference_model(object):

  def __init__(self,
               model,
               reference_hierarchy_list=None,
               reference_file_list=None,
               params=None,
               selection=None,
               log=None):
    assert [reference_hierarchy_list,
            reference_file_list].count(None) == 1
    if(log is None):
      log = sys.stdout
    self.log=log
    # self.model = model
    self.params = params
    self.selection = selection
    self.mon_lib_srv = model.get_mon_lib_srv()
    self.ener_lib = model.get_ener_lib()
    self.pdb_hierarchy = model.get_hierarchy()
    self.pdb_hierarchy.reset_i_seq_if_necessary()
    sites_cart = self.pdb_hierarchy.atoms().extract_xyz()
    if self.selection is None:
      self.selection = flex.bool(len(sites_cart), True)
    if reference_hierarchy_list is None:
      reference_hierarchy_list = \
        utils.process_reference_files(
          reference_file_list=reference_file_list,
          log=log)
    if reference_file_list is None:
      reference_file_list = \
          ["ref%d" % x for x in range(len(reference_hierarchy_list))]
    #
    # this takes 20% of constructor time.
    self.dihedral_proxies_ref = utils.get_reference_dihedral_proxies(
        reference_hierarchy_list=reference_hierarchy_list,
        reference_file_list=reference_file_list,
        mon_lib_srv=self.mon_lib_srv,
        ener_lib=self.ener_lib,
        restraint_objects=model.get_restraint_objects(),
        monomer_parameters=model.get_monomer_parameters(),
        log=log)
    self.i_seq_name_hash = utils.build_name_hash(
                             pdb_hierarchy=self.pdb_hierarchy)
    #reference model components
    self.sites_cart_ref = {}
    self.pdb_hierarchy_ref = {}
    self.i_seq_name_hash_ref = {}
    self.reference_dihedral_hash = {}
    self.reference_file_list = reference_file_list
    #triage reference model files
    for file, hierarchy in zip(reference_file_list,
                               reference_hierarchy_list):
      self.sites_cart_ref[file] = hierarchy.atoms().extract_xyz()
      self.pdb_hierarchy_ref[file] = hierarchy
      self.i_seq_name_hash_ref[file] = \
        utils.build_name_hash(
          pdb_hierarchy=hierarchy)
      self.reference_dihedral_hash[file] = \
        self.build_dihedral_hash(
          dihedral_proxies=self.dihedral_proxies_ref[file],
          sites_cart=self.sites_cart_ref[file],
          pdb_hierarchy=hierarchy,
          include_hydrogens=self.params.hydrogens,
          include_main_chain=self.params.main_chain,
          include_side_chain=self.params.side_chain)
    self.match_map = None
    self.proxy_map = None
    self.build_reference_dihedral_proxy_hash()
    #
    # This takes 80% of constructor time!!!
    self.residue_match_hash = {} # {key_model: ('file_name', key_ref)}
    self.match_map = {} # {'file_name':{i_seq_model:i_seq_ref}}
    if params.use_starting_model_as_reference:
      self.get_matching_from_self()
    else:
      self.get_matching_from_ncs(log=log)
    if self.match_map == {}:
      # making empty container
      new_ref_dih_proxies = self.reference_dihedral_proxies = \
          cctbx.geometry_restraints.shared_dihedral_proxy()
    else:
      new_ref_dih_proxies = self.get_reference_dihedral_proxies(model)

  def get_matching_from_self(self):
    """ Shortcut for the case when restraining on starting model """
    if self.reference_file_list[0] not in self.match_map.keys():
      self.match_map[self.reference_file_list[0]] = {}
    for chain in self.pdb_hierarchy.only_model().chains():
      for rg in chain.residue_groups():
        # Filling out self.residue_match_hash
        key_model = "%s %s" % (rg.unique_resnames()[0], rg.id_str().strip())
        if key_model not in self.residue_match_hash:
          self.residue_match_hash[key_model] = (self.reference_file_list[0], key_model)
        # Filling out self.match_map
        for atom in rg.atoms():
          self.match_map[self.reference_file_list[0]][atom.i_seq] = atom.i_seq



  def _make_matching_and_fill_dictionaries(self, model_h, ref_h, fn,
      m_cache, model_selection_str="all", ref_selection_str="all"):
    ref_cache = self.pdb_hierarchy_ref[fn].atom_selection_cache()
    model_selection_str += " and not element H and not element D and not water"
    ref_selection_str += " and not element H and not element D and not water"
    m_sel = m_cache.selection(model_selection_str)
    ref_sel = ref_cache.selection(ref_selection_str)
    combined_h = model_h.select(m_sel).deep_copy()
    if combined_h.atoms_size() == 0:
      msg = "Selection '%s' selected 0 atoms in refined model.\n" % (model_selection_str) +\
          "Please check if the selection provided is correct."
      raise Sorry(msg)
    ref_h = ref_h.select(ref_sel).deep_copy()
    if ref_h.atoms_size() == 0:
      msg = "Reference selection '%s' selected 0 atoms in %s.\n" % (ref_selection_str, fn) +\
          "Please check if the selection provided is correct."
      raise Sorry(msg)
    for chain in ref_h.only_model().chains():
      chain.id +="ref"
    combined_h.transfer_chains_from_other(ref_h)
    combined_h.reset_atom_i_seqs()
    temp_h = combined_h.deep_copy()
    temp_h.atoms().reset_i_seq()

    # combined_h.write_pdb_file(fn+"_combined.pdb")
    ncs_obj = iotbx.ncs.input(
        hierarchy=temp_h,
        params = self.params.search_options,
        log = self.log)
    # For each found NCS group we going to do matching procedure between
    # copies
    for group_list in ncs_obj.get_ncs_restraints_group_list():
      # combine selections from master and copies into one list...
      n_total_selections = len(group_list.copies) + 1
      ncs_iselections = [group_list.master_iselection]
      ncs_residue_groups = []
      ncs_hierarchys = []
      for i in range(len(group_list.copies)):
        ncs_iselections.append(group_list.copies[i].iselection)
      # and loop over it making new hierarchies.
      for i, isel in enumerate(ncs_iselections):
        ncs_h = combined_h.select(isel)
        ncs_hierarchys.append(ncs_h)
        rgs = list(ncs_h.residue_groups())
        if len(rgs)>0:
          ncs_residue_groups.append(rgs)

      n_total_ncs_residue_groups = len(ncs_residue_groups)
      if n_total_ncs_residue_groups == 0:
        continue
      len_ncs_rg = len(ncs_residue_groups[0])
      for ncs_rg in ncs_residue_groups:
        assert len(ncs_rg) == len_ncs_rg
      if fn not in self.match_map.keys():
        self.match_map[fn] = {}
      ref_indeces = []
      for i in range(n_total_ncs_residue_groups):
        if (len(ncs_residue_groups[i][0].parent().id) > 2 and
            ncs_residue_groups[i][0].parent().id[-3:] == 'ref'):
          ref_indeces.append(i)
      if len(ref_indeces) == 0:
        # Reference model does not participate in particular found NCS copy.
        # If it is in another copy, that's fine.
        continue
      for i in range(n_total_ncs_residue_groups):
        # Figuring out what is reference and what is model
        if i in ref_indeces:
          continue
        a = ncs_residue_groups[i]
        for j in range(len_ncs_rg):
          model_rg = a[j]
          # Here we want to be smarter and find more appropriate reference,
          # in case A, Aref, B, Bref we want to match A-Aref, B-Bref
          current_ref_index = ref_indeces[0]
          if len(ref_indeces) > 1:
            model_chain_id = model_rg.parent().id
            for ri in ref_indeces:
              ref_chain_id = ncs_residue_groups[ri][j].parent().id
              if (ref_chain_id[:-3] == model_chain_id
                  or ref_chain_id == model_chain_id):
                current_ref_index = ri
          reference_rg = ncs_residue_groups[current_ref_index][j]
          # Filling out self.residue_match_hash
          if len(reference_rg.parent().id) > 2 and reference_rg.parent().id[-3:] == 'ref':
            reference_rg.parent().id = reference_rg.parent().id[:-3]
          key_model = "%s %s" % (model_rg.unique_resnames()[0],
              model_rg.id_str().strip())
          key_ref = "%s %s" % (reference_rg.unique_resnames()[0],
              reference_rg.id_str().strip())
          if key_model not in self.residue_match_hash:
            self.residue_match_hash[key_model] = (self.reference_file_list[0], key_ref)
          # Filling out self.match_map
          info_rgs = [[],[]]
          assert reference_rg.atoms_size() == model_rg.atoms_size(), "%s, %s" % (
              model_rg.id_str(), reference_rg.id_str())
          for rg, info in [(model_rg, info_rgs[0]), (reference_rg, info_rgs[1])]:
            info.append(rg.parent().id)
            info.append(rg.unique_resnames()[0])
            info.append(rg.resseq)
            info.append(rg.icode)
          m_str = "chain '%s' and resseq '%s' and icode '%s'" % (
              info_rgs[0][0], info_rgs[0][2], info_rgs[0][3])
          ref_str = "chain '%s' and resseq '%s' and icode '%s'" % (
              info_rgs[1][0], info_rgs[1][2], info_rgs[1][3])
          m_sel = m_cache.selection(m_str)
          ref_sel = ref_cache.selection(ref_str)
          for m_atom, ref_atom in zip(self.pdb_hierarchy.select(m_sel).atoms(),
              self.pdb_hierarchy_ref[fn].select(ref_sel).atoms()):
            self.match_map[fn][m_atom.i_seq] = ref_atom.i_seq

  def is_reference_groups_provided(self):
    if hasattr(self.params, "reference_group"):
      if (len(self.params.reference_group) == 0 or
          (len(self.params.reference_group) == 1 and
          self.params.reference_group[0].reference is None and
          self.params.reference_group[0].selection is None and
          self.params.reference_group[0].file_name is None)):
        return False
      else:
        for i, rg in enumerate(self.params.reference_group):
          if rg.file_name is None:
            if len(self.reference_file_list) > i:
              rg.file_name = self.reference_file_list[i]
            elif len(self.params.file) != 1:
              raise Sorry("Ambigous definition of groups in reference model")
            else:
              rg.file_name = self.params.file[0]
        return True
    else:
      return False

  def get_matching_from_ncs(self, log):
    import iotbx.ncs
    m_cache = self.pdb_hierarchy.atom_selection_cache()
    if not self.is_reference_groups_provided():
      # only files are specified
      for fn in self.reference_file_list:
        print("\nreference file: %s" % fn, file=log)
        print("Model:              Reference:", file=log)
        self._make_matching_and_fill_dictionaries(
            model_h=self.pdb_hierarchy,
            ref_h=self.pdb_hierarchy_ref[fn],
            fn=fn,
            m_cache=m_cache)
    else:
      # We got reference_group section
      for rg in self.params.reference_group:
        file_name = rg.file_name
        print("\nreference file: %s" % file_name, file=log)
        print("Model:              Reference:", file=log)
        self._make_matching_and_fill_dictionaries(
            model_h=self.pdb_hierarchy,
            ref_h=self.pdb_hierarchy_ref[file_name],
            fn=file_name,
            m_cache=m_cache,
            model_selection_str=rg.selection,
            ref_selection_str=rg.reference)

  def proxy_remove(self, selection=None):
    if self.reference_dihedral_proxies is not None:
      self.reference_dihedral_proxies = \
          self.reference_dihedral_proxies.proxy_remove(selection=selection)

  def proxy_select(self, n_seq, iselection):
    import copy
    new_proxies = None
    if self.reference_dihedral_proxies is not None:
      new_proxies = self.reference_dihedral_proxies.proxy_select(
          n_seq, iselection)
    new_manager = copy.copy(self)
    new_manager.reference_dihedral_proxies = new_proxies
    return new_manager

  def show_sorted(self, by_value, sites_cart, site_labels, proxy_label, f):
    if self.reference_dihedral_proxies is not None:
      self.reference_dihedral_proxies.show_sorted(
          by_value=by_value,
          sites_cart=sites_cart,
          site_labels=site_labels,
          proxy_label=proxy_label,
          f=f)

  def target_and_gradients(self, unit_cell, sites_cart, gradient_array):
    res_sum = 0
    if self.reference_dihedral_proxies is not None:
      if unit_cell is None:
        res_sum = cctbx.geometry_restraints.dihedral_residual_sum(
            sites_cart=sites_cart,
            proxies=self.reference_dihedral_proxies,
            gradient_array=gradient_array)
      else:
        res_sum = cctbx.geometry_restraints.dihedral_residual_sum(
            unit_cell=unit_cell,
            sites_cart=sites_cart,
            proxies=self.reference_dihedral_proxies,
            gradient_array=gradient_array)
    return res_sum

  def get_n_proxies(self):
    if self.reference_dihedral_proxies is not None:
      return self.reference_dihedral_proxies.size()
    return 0

  def top_out_function(self, x, weight, top):
    return top*(1-exp(-weight*x**2/top))

  def top_out_gradient(self, x, weight, top):
    return (2*weight*x)*exp(-(weight*x**2)/top)

  def top_out_curvature(self, x, weight, top):
    return (2*weight*(top - 2*weight*x**2))/top**2*exp(-(weight*x**2)/top)

  def build_reference_dihedral_proxy_hash(self):
    self.reference_dihedral_proxy_hash = {}
    for ref in self.dihedral_proxies_ref.keys():
      self.reference_dihedral_proxy_hash[ref] = {}
      proxies = self.dihedral_proxies_ref[ref]
      for dp in proxies:
        key = ""
        for i_seq in dp.i_seqs:
          key += self.i_seq_name_hash_ref[ref][i_seq]
        self.reference_dihedral_proxy_hash[ref][key] = dp

  def build_dihedral_hash(self,
                          dihedral_proxies=None,
                          sites_cart=None,
                          pdb_hierarchy=None,
                          include_hydrogens=False,
                          include_main_chain=True,
                          include_side_chain=True):
    if not include_hydrogens:
      i_seq_element_hash = \
        utils.build_element_hash(pdb_hierarchy=pdb_hierarchy)
    i_seq_name_hash = \
      utils.build_name_hash(pdb_hierarchy=pdb_hierarchy)
    dihedral_hash = dict()

    for dp in dihedral_proxies:
      try:
        #check for H atoms if required
        if not include_hydrogens:
          for i_seq in dp.i_seqs:
            if i_seq_element_hash[i_seq] == " H":
              raise StopIteration()
        #ignore backbone dihedrals
        if not include_main_chain:
          sc_atoms = False
          for i_seq in dp.i_seqs:
            if i_seq_name_hash[i_seq][0:4] not in [' CA ',' N  ',' C  ',' O  ']:
              sc_atoms = True
              break
          if not sc_atoms:
            raise StopIteration()
        if not include_side_chain:
          sc_atoms = False
          for i_seq in dp.i_seqs:
            if i_seq_name_hash[i_seq][0:4] \
              not in [' CA ', ' N  ', ' C  ', ' O  ']:
              sc_atoms = True
              break
          if sc_atoms:
            raise StopIteration()
        key = ""
        for i_seq in dp.i_seqs:
          key = key+i_seq_name_hash[i_seq]
        di = \
          cctbx.geometry_restraints.dihedral(
            sites_cart=sites_cart,
            proxy=dp)
        dihedral_hash[key] = di.angle_model
      except StopIteration:
        pass
    return dihedral_hash

  def _is_proxy_already_present(self, proxy_array, proxy):
    for p in proxy_array:
      if proxy.i_seqs == p.i_seqs:
        return True
    return False

  def get_reference_dihedral_proxies(self, model):
    complete_dihedral_proxies = utils.get_dihedrals_and_phi_psi(
        model=model)
    generated_reference_dihedral_proxies = \
      cctbx.geometry_restraints.shared_dihedral_proxy()
    sigma = self.params.sigma
    limit = self.params.limit
    ref_ss_m = None
    ss_selection = None
    if self.params.secondary_structure_only:
      if (not libtbx.env.has_module(name="ksdssp")):
        raise RuntimeError(
          "ksdssp module is not configured, "+\
          "cannot generate secondary structure reference")
      ref_ss_m = {}
      ss_selection = {}
      for file in self.reference_file_list:
        ref_ss_m[file] = secondary_structure.manager(
          pdb_hierarchy=self.pdb_hierarchy_ref[file],
          sec_str_from_pdb_file=None)
        sec_str_from_pdb_file = ref_ss_m[file].actual_sec_str
        if sec_str_from_pdb_file != None:
          overall_selection = sec_str_from_pdb_file.overall_selection()
          sel_cache_ref = self.pdb_hierarchy_ref[file].atom_selection_cache()
          bsel = sel_cache_ref.selection(string=overall_selection)
          if bsel.all_eq(False):
            raise Sorry("No atom selected")
          ss_selection[file] = bsel

    t_sum = 0
    for dp in complete_dihedral_proxies:
      key_work = ""
      complete = True
      for i_seq in dp.i_seqs:
        if not self.selection[i_seq]:
          complete = False
      if not complete:
        continue

      for i_seq in dp.i_seqs:
        key_work = key_work + self.i_seq_name_hash[i_seq]
      #find matching key
      key = None
      file_match = None
      for file in self.reference_file_list:
        if key is not None or file not in self.match_map:
          continue
        else:
          key = ""
          ref_match = True
          for i_seq in dp.i_seqs:
            if ref_match:
              map_part = self.match_map[file].get(i_seq)
              if map_part is not None:
                key_part = self.i_seq_name_hash_ref[file].get(map_part)
                if key_part is None:
                  ref_match = False
                  key = None
                  file_match = None
                else:
                  key = key+key_part
              else:
                ref_match = False
                key = None
                file_match = None
            if key is not None:
              file_match = file
      try:
        reference_angle = self.reference_dihedral_hash[file_match][key]
      except Exception:
        continue
      w_limit = limit
      w_weight = 1/sigma**2
      if (self.params.secondary_structure_only and ss_selection is not None
          and ss_selection[file_match] is not None):
        limit2 = 15.0
        w_weight = 0.04
        if (ss_selection[file_match][self.match_map[file_match][dp.i_seqs[0]]] and
            ss_selection[file_match][self.match_map[file_match][dp.i_seqs[1]]] and
            ss_selection[file_match][self.match_map[file_match][dp.i_seqs[2]]] and
            ss_selection[file_match][self.match_map[file_match][dp.i_seqs[3]]]):
          limit2 = 30.0
          w_weight = 1
        dp_add = cctbx.geometry_restraints.dihedral_proxy(
            i_seqs=dp.i_seqs,
            angle_ideal=reference_angle,
            weight=w_weight,
            limit=w_limit,
            top_out=TOP_OUT_FLAG)
        generated_reference_dihedral_proxies.append(dp_add)
      else:
        dp_add = cctbx.geometry_restraints.dihedral_proxy(
            i_seqs=dp.i_seqs,
            angle_ideal=reference_angle,
            weight=1/sigma**2,
            limit=limit,
            top_out=TOP_OUT_FLAG)
        # print "Already_there:", self._is_proxy_already_present(generated_reference_dihedral_proxies,dp_add)
        # if not self._is_proxy_already_present(generated_reference_dihedral_proxies,dp_add):
        generated_reference_dihedral_proxies.append(dp_add)
    self.reference_dihedral_proxies = generated_reference_dihedral_proxies

  def show_reference_summary(self, log=None):
    if log is None:
      log = sys.stdout
    print("--------------------------------------------------------", file=log)
    print("Reference Model Matching Summary:", file=log)
    keys = list(self.residue_match_hash.keys())
    def get_key_chain_num(res):
      return res[4:]
    keys.sort(key=get_key_chain_num)
    for file in self.reference_file_list:
      print("\nreference file: %s\n" % file, file=log)
      print("Model:              Reference:", file=log)
      for key in keys:
        if self.residue_match_hash[key][0] == file:
          print("%s  <=====>  %s" % \
            (key, self.residue_match_hash[key][1]), file=log)
    print("\nTotal # of matched residue pairs: %d" % len(keys), file=log)
    print("Total # of reference model restraints: %d" % \
      len(self.reference_dihedral_proxies), file=log)
    print("--------------------------------------------------------", file=log)

  def set_rotamer_to_reference(self,
                               xray_structure,
                               mon_lib_srv=None,
                               log=None,
                               quiet=False):
    if self.mon_lib_srv is None:
      self.mon_lib_srv = mon_lib_srv
    assert isinstance(self.mon_lib_srv, mmtbx.monomer_library.server.server)
    if(log is None): log = sys.stdout
    make_sub_header(
      "Correcting rotamer outliers to match reference model",
      out=log)
    sa = SidechainAngles(False)
    r = rotalyze.rotalyze(pdb_hierarchy=self.pdb_hierarchy)
    rot_list_reference = {}
    coot_reference = {}
    for key in self.pdb_hierarchy_ref.keys():
      hierarchy = self.pdb_hierarchy_ref[key]
      rot_list_reference[key] = \
        rotalyze.rotalyze(pdb_hierarchy=hierarchy)
    model_hash = {}
    model_chis = {}
    reference_hash = {}
    reference_chis = {}
    model_outliers = 0
    for rot in r.results:
      model_hash[rot.id_str()] = rot.rotamer_name
      if rot.rotamer_name == "OUTLIER":
        model_outliers += 1

    for key in rot_list_reference.keys():
      reference_hash[key] = {}
      for rot in rot_list_reference[key].results:
        reference_hash[key][rot.id_str()] = rot.rotamer_name

    print("** evaluating rotamers for working model **", file=log)
    for model in self.pdb_hierarchy.models():
      for chain in model.chains():
        for residue_group in chain.residue_groups():
            all_dict = rotalyze.construct_complete_sidechain(residue_group)
            for atom_group in residue_group.atom_groups():
              try:
                atom_dict = all_dict.get(atom_group.altloc)
                chis = sa.measureChiAngles(atom_group, atom_dict)
                if chis is not None:
                  key = utils.id_str(
                          chain_id=chain.id,
                          resseq=residue_group.resseq,
                          resname=atom_group.resname,
                          icode=residue_group.icode,
                          altloc=atom_group.altloc)
                  model_chis[key] = chis
              except Exception:
                print('  %s%5s %s is missing some sidechain atoms, **skipping**' % (
                      chain.id, residue_group.resid(),
                      atom_group.altloc+atom_group.resname), file=log)
    if model_outliers == 0:
      print("No rotamer outliers detected in working model", file=log)
      return
    else:
      print("Number of rotamer outliers: %d" % model_outliers, file=log)

    print("\n** evaluating rotamers for reference model **", file=log)
    for file in self.pdb_hierarchy_ref.keys():
      hierarchy = self.pdb_hierarchy_ref[file]
      reference_chis[file] = {}
      for model in hierarchy.models():
        for chain in model.chains():
          for residue_group in chain.residue_groups():
              all_dict = rotalyze.construct_complete_sidechain(residue_group)
              for atom_group in residue_group.atom_groups():
                try:
                  atom_dict = all_dict.get(atom_group.altloc)
                  chis = sa.measureChiAngles(atom_group, atom_dict)
                  if chis is not None:
                    key = utils.id_str(
                            chain_id=chain.id,
                            resseq=residue_group.resseq,
                            resname=atom_group.resname,
                            icode=residue_group.icode,
                            altloc=atom_group.altloc)
                    reference_chis[file][key] = chis
                except Exception:
                  print('  %s%5s %s is missing some sidechain atoms, **skipping**' % (
                        chain.id, residue_group.resid(),
                        atom_group.altloc+atom_group.resname), file=log)

    print("\n** fixing outliers **", file=log)
    sites_cart_start = xray_structure.sites_cart()
    for model in self.pdb_hierarchy.models():
      for chain in model.chains():
        for residue_group in chain.residue_groups():
          if len(residue_group.conformers()) > 1:
            print("  %s%5s %s has multiple conformations, **skipping**" % (
              chain.id, residue_group.resid(),
              " "+residue_group.atom_groups()[0].resname), file=log)
            continue
          for conformer in residue_group.conformers():
            for residue in conformer.residues():
              if residue.resname == "PRO":
                continue
              key = utils.id_str(
                      chain_id=chain.id,
                      resseq=residue_group.resseq,
                      resname=residue_group.atom_groups()[0].resname,
                      icode=residue_group.icode,
                      altloc=conformer.altloc)
              if len(chain.id) == 1:
                chain_id = " "+chain.id
              else:
                chain_id = chain.id
              file_key = '%s%s%s' %(residue.resname,
                                    chain_id,
                                    residue_group.resid())
              file_key = file_key.strip()
              file_match = self.residue_match_hash.get(file_key)
              if file_match is not None:
                file = file_match[0]
              else:
                continue
              model_rot = model_hash.get(key)
              reference_rot = reference_hash[file].get(self.one_key_to_another(file_match[1]))
              m_chis = model_chis.get(key)
              r_chis = reference_chis[file].get(self.one_key_to_another(file_match[1]))
              if model_rot is not None and reference_rot is not None and \
                  m_chis is not None and r_chis is not None:
                if (model_rot == 'OUTLIER' and \
                    reference_rot != 'OUTLIER'): # or \
                    #atom_group.resname in ["LEU", "VAL", "THR"]:
                  self.change_residue_rotamer_in_place(
                      sites_cart_start,residue, m_chis,r_chis,self.mon_lib_srv)
                  xray_structure.set_sites_cart(sites_cart_start)

                elif self.params.strict_rotamer_matching and \
                  (model_rot != 'OUTLIER' and reference_rot != 'OUTLIER'):
                  if model_rot != reference_rot:
                    self.change_residue_rotamer_in_place(
                        sites_cart_start,residue, m_chis,r_chis,self.mon_lib_srv)
                    xray_structure.set_sites_cart(sites_cart_start)

  def one_key_to_another(self,key):
    # Work-around function, don't have time to dig into 10 different cache
    # types... Probably better data structure could be suggested to handle
    # all needed information.

    # sp = key.split()
    # assert len(sp) == 3
    var1 = " %s  %s" % (key[4:], key[:3])
    return var1

  def change_residue_rotamer_in_place(self,sites_cart, residue,
      m_chis, r_chis, mon_lib_srv):
    assert m_chis.count(None) == 0
    assert r_chis.count(None) == 0
    axis_and_atoms_to_rotate= \
      rotatable_bonds.axes_and_atoms_aa_specific(
          residue=residue,
          mon_lib_srv=mon_lib_srv,
          remove_clusters_with_all_h=True,
          log=None)
    if axis_and_atoms_to_rotate is None:
      return
    assert len(m_chis) == len(axis_and_atoms_to_rotate)
    assert len(r_chis) >= len(m_chis)
    counter = 0
    residue_iselection = residue.atoms().extract_i_seq()
    sites_cart_residue = sites_cart.select(residue_iselection)
    for aa in axis_and_atoms_to_rotate:
      axis = aa[0]
      atoms = aa[1]
      residue.atoms().set_xyz(new_xyz=sites_cart_residue)
      new_xyz = flex.vec3_double()
      angle_deg = r_chis[counter] - m_chis[counter]
      if angle_deg < 0:
        angle_deg += 360.0
      for atom in atoms:
        new_xyz = rotate_point_around_axis(
                    axis_point_1=sites_cart_residue[axis[0]],
                    axis_point_2=sites_cart_residue[axis[1]],
                    point=sites_cart_residue[atom],
                    angle=angle_deg, deg=True)
        sites_cart_residue[atom] = new_xyz
      sites_cart = sites_cart.set_selected(
            residue_iselection, sites_cart_residue)
      counter += 1

  def remove_restraints_with_ncs_matches(self,
                                         ncs_dihedral_proxies,
                                         ncs_match_hash,
                                         log=None):
    if not self.params.auto_shutoff_for_ncs:
      return
    if log is None:
      log = sys.stdout
    proxy_list = []
    remaining_proxies = cctbx.geometry_restraints.shared_dihedral_proxy()
    remaining_match_hash = {}
    for dp in ncs_dihedral_proxies:
      proxy_list.append(dp.i_seqs)
    for dp in self.reference_dihedral_proxies:
      if dp.i_seqs not in proxy_list:
        remaining_proxies.append(dp)
    for key in self.residue_match_hash:
      found_match = False
      for key2 in ncs_match_hash:
        if key == key2:
          found_match = True
        else:
          for match in ncs_match_hash[key2]:
            if key == match:
              found_match = True
      if not found_match:
        remaining_match_hash[key] = self.residue_match_hash[key]
    self.reference_dihedral_proxies = remaining_proxies
    self.residue_match_hash = remaining_match_hash
    print("\n**Removed reference restraints that overlap "+ \
                       "with torsion NCS restraints**\n", file=log)
    print("Updated Reference Model Restraints:", file=log)
    self.show_reference_summary(log=log)
