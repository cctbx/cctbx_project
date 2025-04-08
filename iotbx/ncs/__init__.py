"""Routines for reading and writing and manipulating symmetry files (NCS
for crystallography, symmetry for cryo-EM)
"""
from __future__ import absolute_import, division, print_function
from iotbx.pdb.atom_selection import selection_string_from_selection
from scitbx.array_family import flex
from mmtbx.ncs import ncs_search
from libtbx.utils import Sorry, null_out
import libtbx.phil
import iotbx.pdb.hierarchy
from mmtbx.ncs import ncs
from mmtbx.ncs.ncs_restraints_group_list import class_ncs_restraints_group_list, \
    NCS_restraint_group, NCS_copy
from scitbx import matrix
import sys
from iotbx.pdb.utils import all_chain_ids
from time import time
from six.moves import cStringIO as StringIO
from six.moves import zip

ncs_search_options = """\
ncs_search
  .short_caption = Search options
  .style = box
  .help = Set of parameters for NCS search procedure. Some of them also used \
      for filtering user-supplied ncs_group.
{
  enabled = False
    .type = bool
    .help = Enable NCS restraints or constraints in refinement (in some cases \
        may be switched on inside refinement program).
    .short_caption = Use NCS
    .style = noauto bold
  exclude_selection = "element H or element D or water"
    .type = atom_selection
    .help = Atoms selected by this selection will be excluded from the model \
      before any NCS search and/or filtering procedures. There is no way \
      atoms defined by this selection will be in NCS.
    .expert_level = 2
  chain_similarity_threshold = 0.85
    .type=float
    .short_caption = Sequence alignment threshold
    .help='''Threshold for sequence similarity between matching chains.
      A smaller value may cause more chains to be grouped together and can lower
      the number of common residues'''
    .expert_level = 0
  chain_max_rmsd = 2.
    .type = float
    .short_caption = Max RMSD between matching chains
    .help = '''limit of rms difference between chains to be considered
       as copies'''
    .expert_level = 0
  residue_match_radius = 4.0
    .type = float
    .help = Maximum allowed distance difference between pairs of matching \
        atoms of two residues
    .expert_level = 0
  try_shortcuts = False
    .type = bool
    .help = Try very quick check to speed up the search when chains are identical. \
        If failed, regular search will be performed automatically.
    .expert_level = 0
  minimum_number_of_atoms_in_copy = 3
    .type = int
    .help = Do not create ncs groups where master and copies would contain \
        less than specified amount of atoms
    .expert_level = 3
  validate_user_supplied_groups = True
    .type = bool
    .help = Enable validation of user-supplied ncs_group. Need to exercise \
        a lot of caution turning this off. This option is for \
        developers only.
    .expert_level = 3
}
"""

# parameters for manual specification of NCS - ASU mapping
ncs_group_phil_str = '''\
ncs_group
  .multiple = True
  .short_caption = NCS group definition
  .help = The definition of one NCS group. Note, that almost always in \
    refinement programs they will be checked and filtered if needed.
  .style = auto_align
  .expert_level=0
  {
  reference = None
    .type = str
    .short_caption = Reference selection
    .help = 'Residue selection string for the complete master NCS copy'
    .expert_level=0
  selection = None
    .type = str
    .short_caption = NCS related selection
    .help = 'Residue selection string for each NCS copy location in ASU'
    .multiple = True
    .expert_level=0
  }
'''

ncs_group_master_phil = libtbx.phil.parse(ncs_group_phil_str)

class input(object):
  def __init__(self,
          hierarchy=None,
          # XXX warning, ncs_phil_groups can be changed inside...
          ncs_phil_groups = None,
          params = None,
          log=None,
          ):
    """
    TODO:
    1. Transfer get_ncs_info_as_spec() to ncs/ncs.py:ncs

    Select method to build ncs_group_object

    order of implementation:
    1) ncs_phil_groups - user-supplied definitions are filtered
    2) hierarchy only - Performing NCS search

    Args:
    -----
      ncs_phil_groups:  iotbx.phil.parse(ncs_group_phil_str).extract().ncs_group
      chain_max_rmsd (float): limit of rms difference between chains to be considered
        as copies
      min_percent (float): Threshold for similarity between chains
        similarity define as:
        (number of matching res) / (number of res in longer chain)
      chain_similarity_threshold (float): min similarity between matching chains
      residue_match_radius (float): max allow distance difference between pairs of matching
        atoms of two residues
    """

    self.number_of_ncs_groups = 0 # consider removing/replacing with function

    self.ncs_restraints_group_list = class_ncs_restraints_group_list()
    # keep hierarchy for writing (To have a source of atoms labels)
    self.hierarchy = hierarchy
    # residues common to NCS copies. Used for .spec representation
    self.common_res_dict = {}
    # Collect messages, recommendation and errors
    self.messages = '' # Not used outside...
    self.old_i_seqs = None
    self.original_hierarchy = None
    self.truncated_hierarchy = None
    self.truncated_h_asc = None
    self.chains_info = None
    self.phil_groups_modified = False # Switched if modifying input phil groups

    extension = ''
    # set search parameters
    self.params = params
    if self.params is None:
      self.params = input.get_default_params().ncs_search
    #
    if log is None:
      self.log = sys.stdout
    else:
      self.log = log

    if hierarchy:
      # for a in hierarchy.atoms():
      #   print "oo", a.i_seq, a.id_str()
      # print "====="
      hierarchy.atoms().reset_i_seq()
      self.original_hierarchy = hierarchy.deep_copy()
      self.original_hierarchy.reset_atom_i_seqs()
      if self.params.exclude_selection is not None:
        # pdb_hierarchy_inp.hierarchy.write_pdb_file("in_ncs_pre_before.pdb")
        cache = hierarchy.atom_selection_cache()
        sel = cache.selection("not (%s)" % self.params.exclude_selection)
        self.truncated_hierarchy = hierarchy.select(sel)
      else:
        # this could be to save iseqs but I'm not sure
        self.truncated_hierarchy = hierarchy.select(flex.size_t_range(hierarchy.atoms_size()))
      self.old_i_seqs = self.truncated_hierarchy.atoms().extract_i_seq()
      # print "self.old_i_seqs", list(self.old_i_seqs)
      # self.truncated_hierarchy.atoms().reset_i_seq()
      self.truncated_hierarchy.reset_atom_i_seqs()
      self.truncated_h_asc = self.truncated_hierarchy.atom_selection_cache()
      # self.truncated_hierarchy.write_pdb_file("in_ncs_pre_after.pdb")
      self.chains_info = ncs_search.get_chains_info(self.truncated_hierarchy)


      if self.truncated_hierarchy.atoms_size() == 0:
        return

    #
    # print "ncs_groups before validation", ncs_phil_groups
    # validated_ncs_phil_groups is for correct handling of validation results:
    # it ends up as None in many corner cases where we actually want to do
    # automatic search. Return value is never None once the actual validation
    # takes place.
    validated_ncs_phil_groups = None
    validated_ncs_phil_groups = self.validate_ncs_phil_groups(
      pdb_h = self.truncated_hierarchy,
      ncs_phil_groups   = ncs_phil_groups,
      asc = self.truncated_h_asc,
      validate_user_supplied_groups = self.params.validate_user_supplied_groups)
    if validated_ncs_phil_groups is None:
      # print "Last chance, building from hierarchy"
      self.build_ncs_obj_from_pdb_asu(
          pdb_h=self.truncated_hierarchy,
          asc=self.truncated_h_asc)

    # error handling
    if self.ncs_restraints_group_list.get_n_groups() == 0:
      print('========== WARNING! ============\n', file=self.log)
      print('  No NCS relation were found !!!\n', file=self.log)
      print('================================\n', file=self.log)
    if self.messages != '':
      print(self.messages, file=self.log)

  @staticmethod
  def get_default_params():
    """
    Get parsed parameters (in form of Python objects). Use this function to
    avoid importing ncs_search phil strings and think about how to
    parse it. Does not need the instance of class (staticmethod).
    Then modify what needed to be modified and init this class normally.
    """
    import iotbx.phil
    return iotbx.phil.parse(
          input_string=ncs_search_options,
          process_includes=True).extract()

  def pdb_h_into_chain(self, pdb_h, ch_id="A"):
    new_chain = iotbx.pdb.hierarchy.chain(id=ch_id)
    n_res_groups = 0
    for chain in pdb_h.only_model().chains():
      n_res_groups += chain.residue_groups_size()
    new_chain.pre_allocate_residue_groups(
        number_of_additional_residue_groups=n_res_groups)
    new_resseq = 1
    for chain in pdb_h.only_model().chains():
      for rg in chain.residue_groups():
        new_rg = rg.detached_copy()
        new_rg.resseq = new_resseq
        original_iseqs = rg.atoms().extract_i_seq()
        for atom, orig_iseq in zip(new_rg.atoms(), original_iseqs):
          atom.tmp = orig_iseq
        new_resseq += 1
        new_chain.append_residue_group(residue_group=new_rg)
    return new_chain

  def validate_ncs_phil_groups(
      self,
      pdb_h,
      ncs_phil_groups,
      asc,
      validate_user_supplied_groups=True):
    """
    Note that the result of this procedure is corrected ncs_phil_groups.
    These groups will be later submitted to build_ncs_obj_from_phil
    procedure. This is sub-optimal and should be changed because
    everything is already processed here and ready to build proper
    NCS_restraint_group object.
    add filtered groups in self.ncs_restraints_group_list
    """
    def show_particular_ncs_group(ncs_gr):
      p_obj = ncs_group_master_phil.extract()
      p_obj.ncs_group[0].reference = ncs_gr.reference
      p_obj.ncs_group[0].selection = ncs_gr.selection
      to_show = ncs_group_master_phil.format(python_object=p_obj)
      to_show.show(out=self.log)

    def show_empty_selection_error_message(ng, where="reference"):
      print("  Missing or corrupted %s field:" % where, file=self.log)
      print("  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", file=self.log)
      print("      _ALL_ user-supplied groups will be ignored", file=self.log)
      print("  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", file=self.log)
      show_particular_ncs_group(ng)

    # Massage NCS groups
    # return ncs_phil_groups
    validated_ncs_groups = []
    if ncs_phil_groups is None:
      return None
    if(ncs_phil_groups is not None and len(ncs_phil_groups)==0):
      # print "exiting here"
      ncs_phil_groups=None
      return None
    if (ncs_phil_groups is not None and
        len(ncs_phil_groups)==1 and
        ncs_phil_groups[0].reference is None and
        len(ncs_phil_groups[0].selection) == 1 and
        ncs_phil_groups[0].selection[0] is None):
      # This is empty ncs_group definition somehow creeped into here.
      # Not a big deal.
      return None
    if(ncs_phil_groups is not None):
      print("Validating user-supplied NCS groups...", file=self.log)
      empty_cntr = 0
      for ng in ncs_phil_groups:
        if ng.reference is None or len(ng.reference.strip())==0:
          show_empty_selection_error_message(ng, where="reference")
          empty_cntr += 1
        for s in ng.selection:
          if s is None or len(s.strip())==0:
            show_empty_selection_error_message(ng, where="selection")
            empty_cntr += 1
      if(empty_cntr>0):
        print("  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", file=self.log)
        print("      _ALL_ user-supplied groups are ignored.", file=self.log)
        print("  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", file=self.log)
        ncs_phil_groups=None
        return None
    # Verify NCS selections
    msg="Empty selection in NCS group definition: %s"
    if not validate_user_supplied_groups:
      for ncs_group in ncs_phil_groups:
        print("  Copying user-supplied groups without validation:", file=self.log)
        show_particular_ncs_group(ncs_group)
        m_isel = asc.iselection(ncs_group.reference)
        ng = NCS_restraint_group(
            master_iselection = m_isel,
            str_selection = ncs_group.reference)
        for s_string in ncs_group.selection:
          c_isel = asc.iselection(s_string)
          c = NCS_copy(
              copy_iselection=c_isel,
              rot = None,
              tran = None,
              str_selection=s_string,
              rmsd=999)
          ng.append_copy(c)
        self.ncs_restraints_group_list.append(ng)
        validated_ncs_groups.append(ng)
      master_sel = flex.bool(pdb_h.atoms_size(), False)
      for gr in self.ncs_restraints_group_list:
        master_sel.set_selected(gr.master_iselection, True)
      # asu_sites = pdb_h.atoms().extract_xyz().select(master_sel)
      # self.ncs_restraints_group_list._show(hierarchy=pdb_h,brief=False)
      # STOP()
      self.ncs_restraints_group_list.recalculate_ncs_transforms(asu_site_cart=pdb_h.atoms().extract_xyz())
      # validated_ncs_groups.recalculate_ncs_transforms(asu_site_cart=pdb_h.atoms().extract_xyz())

    else:
      for ncs_group in ncs_phil_groups:
        print("  Validating:", file=self.log)
        show_particular_ncs_group(ncs_group)
        selection_list = []
        # first, check for selections producing 0 atoms
        user_original_reference_iselection = None
        user_original_copies_iselections = []
        n_atoms_in_user_ncs = 0
        s_string = ncs_group.reference
        if s_string is not None:
          sel = asc.iselection(s_string)
          selection_list.append(s_string)
          n_atoms_in_user_ncs = sel.size()
          if(n_atoms_in_user_ncs==0):
            raise Sorry(msg%s_string)
          user_original_reference_iselection = sel
        for s_string in ncs_group.selection:
          if(s_string is not None):
            sel = asc.iselection(s_string)
            selection_list.append(s_string)
            n_copy = sel.size()
            if(n_copy==0):
              raise Sorry(msg%s_string)
            user_original_copies_iselections.append(sel)
        #
        # The idea for user's groups is to pick them one by one,
        # select only reference and selections from the model,
        # If there are multiple chains in ref or selection -
        # combine them in one chain,
        # save atom original i_seq in atom.tmp
        # run searching procedure for the resulting hierarchy
        # if the user's selections were more or less OK - there should be
        # one group, get atom.tmp values for the selected atoms and using
        # original hierarchy convert them into string selections when needed.
        # If multiple groups produced - use them, most likely the user
        # provided something really wrong.
        # Need to pay some attention to what came out as master and what order
        # of references.
        #
        combined_h = iotbx.pdb.hierarchy.root()
        combined_h.append_model(iotbx.pdb.hierarchy.model())
        all_c_ids = all_chain_ids()
        cur_ch_id_n = 0
        master_chain = self.pdb_h_into_chain(pdb_h.select(
            user_original_reference_iselection),ch_id=all_c_ids[cur_ch_id_n])
        # print "tmp in master chain:", list(master_chain.atoms().extract_tmp_as_size_t())
        cur_ch_id_n += 1
        combined_h.only_model().append_chain(master_chain)

        # combined_h = iotbx.pdb.hierarchy.new_hierarchy_from_chain(master_chain)
        # print "tmp combined_h1:", list(combined_h.atoms().extract_tmp_as_size_t())
        for uocis in user_original_copies_iselections:
          # print "adding selection to combined:", s_string
          sel_chain = self.pdb_h_into_chain(pdb_h.select(
            uocis),ch_id=all_c_ids[cur_ch_id_n])
          combined_h.only_model().append_chain(sel_chain)
          cur_ch_id_n += 1

        combined_h.reset_atom_i_seqs()
        # combined_h.write_pdb_file("combined_in_validation.pdb")
        # print "tmp:", list(combined_h.atoms().extract_tmp_as_size_t())


        # XXX Here we will regenerate phil selections using the mechanism
        # for finding NCS in this module. Afterwards we should have perfectly
        # good phil selections, and later the object will be created from
        # them.
        # Most likely this is not the best way to validate user selections.

        # selection_list
        nrgl_fake_iseqs = ncs_search.find_ncs_in_hierarchy(
            ph=combined_h,
            chains_info=None,
            chain_max_rmsd=max(self.params.chain_max_rmsd, 10.0),
            log=None,
            chain_similarity_threshold=min(self.params.chain_similarity_threshold, 0.5),
            residue_match_radius=max(self.params.residue_match_radius, 1000.0))
        # hopefully, we will get only 1 ncs group
        # ncs_group.selection = []
        if nrgl_fake_iseqs.get_n_groups() == 0:
          # this means that user's selection doesn't match
          # print "ZERO NCS groups found"
          rejected_msg = "  REJECTED because copies don't match good enough.\n" + \
          "Try to revise selections or adjust chain_similarity_threshold or \n" + \
          "chain_max_rmsd parameters."
          print(rejected_msg, file=self.log)
          continue
        # User triggered the fail of this assert!
        selections_were_modified = False
        #
        for ncs_gr in nrgl_fake_iseqs:
          new_gr = ncs_gr.deep_copy()
          new_ncs_group = ncs_group_master_phil.extract().ncs_group[0]
          for i, isel in enumerate(ncs_gr.get_iselections_list()):
            m_all_isel = isel.deep_copy()
            original_m_all_isel = combined_h.atoms().\
                select(m_all_isel).extract_tmp_as_size_t()
            if n_atoms_in_user_ncs > original_m_all_isel.size():
              selections_were_modified = True
            # print "new isels", list(m_all_isel)
            # print "old isels", list(original_m_all_isel)
            all_m_select_str = selection_string_from_selection(
                pdb_h=pdb_h,
                selection=original_m_all_isel,
                chains_info=self.chains_info,
                atom_selection_cache=asc)
            # print "all_m_select_str", all_m_select_str
            if i == 0:
              new_gr.master_iselection = original_m_all_isel
              new_gr.master_str_selection = all_m_select_str
              new_ncs_group.reference=all_m_select_str
            else:
              new_gr.copies[i-1].iselection = original_m_all_isel
              new_gr.copies[i-1].str_selection = all_m_select_str
              new_ncs_group.selection.append(all_m_select_str)
          self.ncs_restraints_group_list.append(new_gr)
          new_ncs_group.selection = new_ncs_group.selection[1:]
          validated_ncs_groups.append(new_ncs_group)
        # Finally, we may check the number of atoms in selections that will
        # go further.
        # XXX Deleted, because this is taken care of previously
        ok_msg = "  OK. All atoms were included in" +\
        " validated selection.\n"
        modified_msg = "  MODIFIED. Some of the atoms were excluded from" + \
        " your selection.\n  The most common reasons are:\n" + \
        "    1. Missing residues in one or several copies in NCS group.\n" + \
        "    2. Presence of alternative conformations (they are excluded).\n" + \
        "    3. Residue mismatch in requested copies.\n" + \
        "  Please check the validated selection further down.\n"
        if selections_were_modified:
          print(modified_msg, file=self.log)
          self.phil_groups_modified = True
        else:
          print(ok_msg, file=self.log)
    # print "len(validated_ncs_groups)", len(validated_ncs_groups)
    # for ncs_gr in validated_ncs_groups:
    #   print "  reference:", ncs_gr.reference
    #   print "  selection:", ncs_gr.selection
    self.finalize_nrgl()
    return validated_ncs_groups

  def finalize_nrgl(self):
    self.ncs_restraints_group_list = self.ncs_restraints_group_list.\
        filter_out_small_groups(min_n_atoms=self.params.minimum_number_of_atoms_in_copy)
    self.number_of_ncs_groups = self.ncs_restraints_group_list.get_n_groups()
    #
    # Warning! str_selections updating with truncated hierarchy because
    # there are tests
    # phenix_regression/refinement/torsion_ncs/tst_refinement_torsion_ncs.py
    # phenix_regression/refinement/ncs/tst_ncs_5.py
    # that fail otherwise. It was ok for couple of years, so not rushing to fix.
    #
    self.ncs_restraints_group_list.update_str_selections_if_needed(
        hierarchy=self.truncated_hierarchy,
        asc=self.truncated_h_asc,
        chains_info=self.chains_info)
    self.ncs_restraints_group_list.update_i_seqs(self.old_i_seqs)


  def build_ncs_obj_from_pdb_asu(self,pdb_h, asc):
    """
    Build transforms objects and NCS <-> ASU mapping from a complete ASU
    Note that the MTRIX record are ignored, they are produced in the
    process of identifying the master NCS

    Args::
      pdb_h : pdb_hierarchy
    """
    if len(pdb_h.models()) > 1:
      raise Sorry('Multi-model PDB (with MODEL-ENDMDL) is not supported.')
    chain_ids = {x.id for x in pdb_h.models()[0].chains()}
    if len(chain_ids) > 1:
      t0 = time()
      if self.params.try_shortcuts:
        # probably the most parameters are not necessary, since
        # we are going after cases where molecule was multiplied using
        # BIOMT or MTRIX, so we are expecting identical chains with 0 rmsd.
        self.ncs_restraints_group_list = ncs_search.shortcut_1(
            hierarchy=pdb_h,
            chains_info=self.chains_info,
            chain_similarity_threshold=self.params.chain_similarity_threshold,
            chain_max_rmsd=self.params.chain_max_rmsd,
            log=null_out(),
            residue_match_radius=self.params.residue_match_radius)
      print ("Time spend for trying shortcut: %.2f" % (time()-t0), file=self.log)
      if self.ncs_restraints_group_list.get_n_groups() == 0:
        # shortcuts failed
        self.ncs_restraints_group_list = ncs_search.find_ncs_in_hierarchy(
          ph=pdb_h,
          chains_info=self.chains_info,
          chain_similarity_threshold=self.params.chain_similarity_threshold,
          chain_max_rmsd=self.params.chain_max_rmsd,
          log=self.log,
          residue_match_radius=self.params.residue_match_radius)
      self.finalize_nrgl()

  def get_ncs_restraints_group_list(self):
    return self.ncs_restraints_group_list

  def set_ncs_restraints_group_list(self, nrgl):
    self.ncs_restraints_group_list = nrgl

  def get_ncs_info_as_spec(
          self,
          exclude_h=None,
          exclude_d=None,
          stem=None,
          write_ncs_domain_pdb=False,
          log = None):
    """
    XXX This function should be transfered to mmtbx/ncs/ncs.py:ncs class as
    its classmethod, because it creates an object and this is the task of
    a constructor. And it definetely should be decoupled from file creation!


    Returns ncs spec object and can prints ncs info in a ncs_spec,
    format_all_for_resolve or format_all_for_phenix_refine format

    Note that while ncs_groups can master ncs can be comprised from several
    chains, the spec groups can not. So groups with multiple chains in the
    master selection are splitted

    Note that spec format does not support insertions notation
    for example "resseq 49" will include "resid 49" and "resid 49A"

    Args:
      write: (bool) when False, will not write to file or print
      exclude_h,exclude_d : parameters of the ncs object
    Return:
      spec_object
    """
    log = log or self.log
    if not stem : stem =''
    else: stem += '_'
    spec_object = ncs.ncs(exclude_h=exclude_h,exclude_d=exclude_d)
    xyz = self.truncated_hierarchy.atoms().extract_xyz()
    #===============================================================
    # New implementation
    # Here we have original i_seqs already in ncs_restraints_group_list,
    # so we should use self.hierarchy for everything
    xyz = self.hierarchy.atoms().extract_xyz()
    assert self.ncs_restraints_group_list is not None
    splitted_nrgl = self.ncs_restraints_group_list.split_by_chains(
        hierarchy=self.hierarchy)
    for i_group, group in enumerate(splitted_nrgl):
      center_orth = []
      rotations = []
      translations = []
      # chain id
      chain_id_list = []
      # this is [ [[1, 2] [5, 6]] ]
      residue_range_list = []
      rmsd_list = []
      # number of residues
      residues_count = []

      # Putting master in:
      center_orth.append(get_center_orth(xyz,group.master_iselection))
      rotations.append(matrix.sqr([1,0,0,0,1,0,0,0,1]))
      translations.append(matrix.col([0,0,0]))
      chain_id, ranges, count = get_chain_and_ranges(
          self.hierarchy.select(group.master_iselection))
      chain_id_list.append(chain_id)
      residue_range_list.append(ranges)
      residues_count.append(count)
      rmsd_list.append(0)
      for c in group.copies:
        center_orth.append(get_center_orth(xyz,c.iselection))
        # in spec files transform is copy -> master, not master -> copy
        r,t = inverse_transform(c.r,c.t)
        rotations.append(r)
        translations.append(t)

        chain_id, ranges, count = get_chain_and_ranges(
            self.hierarchy.select(c.iselection))
        chain_id_list.append(chain_id)
        residue_range_list.append(ranges)
        residues_count.append(count)
        rmsd_list.append(c.rmsd)
      # XXX This should be consistent with full_file_name parameter in
      # simple_ncs_from_pdb.py: create_ncs_domain_pdb_files()
      # This is here just because we need to output filename of the domain
      # into the spec file if pdb file is going to be created...
      ncs_domain_pdb = None
      if write_ncs_domain_pdb:
        ncs_domain_pdb = stem+'group_'+str(i_group+1)+'.pdb'
      spec_object.import_ncs_group(
        center_orth = center_orth,
        ncs_rota_matr = rotations,
        trans_orth = translations,
        rmsd_list = rmsd_list,
        chain_residue_id = [chain_id_list,residue_range_list],
        residues_in_common_list = residues_count,
        ncs_domain_pdb = ncs_domain_pdb)
    spec_object._ncs_obj = self
    return spec_object

  def print_ncs_phil_param(self,write=False,log=None):
    """
    Prints NCS information in the phil parameters format
    lines longer that 80 characters are folded

    Phil structure example:
      ncs_group {
        reference = 'chain A'
        selection = 'chain C'
        selection = 'chain E'
      }
      ncs_group {
        reference = 'chain B'
        selection = 'chain D'
        selection = 'chain F'
      }

    Args:
      write (bool): when true, print to log
      log : location of output, an open file or sys.stdout

    Returns:
      (str): NCS phil parameter string
    """
    if not log: log = sys.stdout
    groups = []
    for gr in self.ncs_restraints_group_list:
      master = format_80(gr.master_str_selection)
      groups.append('ncs_group {')
      groups.append("  reference = {}".format(master))
      for c in gr.copies:
        cp = format_80(c.str_selection)
        groups.append("  selection = {}".format(cp))
      groups.append('}')
    gr = '\n'.join(groups)
    gr += '\n'
    if write:
      print(gr, file=log)
    return gr

  def show(self,
           format=None,
           verbose=False,
           prefix='',
           header=True,
           log=None):

    """
    Display NCS object

    Args:
      format (str): "phil" : phil file representation
                    "spec" : spec representation out of NCS groups
                    "cctbx": cctbx representation out of NCS groups
                    "restraints"  : .ncs (phenix refine) format
                    "constraints" : .ncs (phenix refine) format
      verbose (bool): when True, will print selection strings, rotation and
        translation info
      prefix (str): a string to be added, padding the output, at the left of
        each line
      header (bool): When True, include header
      log: where to log the output, by default set to sys.stdout
    """
    if not log: log = self.log
    out_str = ''
    if (not format) or (format.lower() == 'cctbx'):
      out_str = self.__repr__(prefix)
      print(out_str, file=log)
      if verbose:
        print(self.show_ncs_selections(prefix), file=log)
    elif format.lower() == 'phil':
      out_str = self.show_phil_format(prefix=prefix,header=header)
      print(out_str, file=log)
    elif format.lower() == 'spec':
      # Does not add prefix in SPEC format
      out_str = self.show_search_parameters_values(prefix) + '\n'
      out_str += self.show_chains_info(prefix) + '\n'
      out_str += '\n' + prefix + 'NCS object "display_all"'
      spec_obj = self.get_ncs_info_as_spec()
      strio = StringIO()
      spec_obj.display_all(log=strio)
      out_str += strio.getvalue()
      out_str += '\n\n' + prefix + 'NCS object "display_inverse"'
      inv_spec_obj=spec_obj.deep_copy()
      inv_spec_obj.invert_matrices()
      strio = StringIO()
      inv_spec_obj.display_all(log=strio)
      out_str += strio.getvalue()
      print(out_str, file=log)
    return out_str

  def show_phil_format(self,prefix='',header=True,group_prefix=''):
    """
    Returns a string of NCS groups phil parameters

    Args:
      prefix (str): a string to be added, padding the output, at the left of
        each line
      header (bool): When True, include header
      group_prefix (str): prefix for the group only
    """
    str_out = []
    if header:
      msg = '\n{}NCS phil parameters:'
      str_out = [msg.format(prefix),'-'*len(msg)]
    str_line = prefix + '  {:s} = {}'
    str_ncs_group =  prefix + group_prefix + 'ncs_group {\n%s' + prefix + '\n}'
    for gr in self.ncs_restraints_group_list:
      str_gr = [str_line.format('reference',gr.master_str_selection)]
      for c in gr.copies:
        str_gr.append(str_line.format('selection',c.str_selection))
      str_gr = '\n'.join(str_gr)
      str_out.append(str_ncs_group%str_gr)
    str_out = '\n'.join(str_out)
    return str_out

  def show_search_parameters_values(self,prefix=''):
    """
    Returns a string of search parameters values

    Args:
      prefix (str): a string to be added, padding the output, at the left of
        each line
    """
    strio = StringIO()
    ncs_group_master_phil.format(self.params).show(out=strio)
    return strio.getvalue()

  def show_chains_info(self,prefix=''):
    """
    Returns formatted string for print out, string containing chains IDs in a
    table format, padded from the left with "prefix"

    Args:
      prefix (str): a string to be added, padding the output, at the left of
        each line
    """
    model = self.truncated_hierarchy.models()[0]
    chain_ids = {x.id for x in model.chains()}
    model_unique_chains_ids = tuple(sorted(chain_ids))
    ids = sorted(model_unique_chains_ids)
    str_out = ['\n{}Chains in model:'.format(prefix),'-'*51]
    n = len(ids)
    item_in_row = 10
    n_rows = n // item_in_row
    last_row = n % item_in_row
    str_ids = [prefix + '{:5s}' * item_in_row] * n_rows
    str_ids_last = prefix + '{:5s}' * last_row
    # connect all output stings
    str_out.extend(str_ids)
    str_out.append(str_ids_last)
    str_out.append('. '*26)
    str_out = '\n'.join(str_out)
    str_out = str_out.format(*ids)
    return str_out

  def show_transform_info(self,prefix=''):
    """
    Returns formatted string for print out, string containing chains IDs in a
    table format, padded from the left with "prefix"

    Args:
      prefix (str): a string to be added, padding the output, at the left of
        each line
    """
    str_out = ['\n{}Transforms:'.format(prefix),'-'*51]
    str_line = prefix + '{:<25s}:   {}'
    str_r = prefix + 'ROTA  {:2}{:10.4f}{:10.4f}{:10.4f}'
    str_t = prefix + 'TRANS   {:10.4f}{:10.4f}{:10.4f}'
    for i_gr, gr in enumerate(self.ncs_restraints_group_list):
      str_out.append(str_line.format('Group #',i))
      for j, c in enumerate(gr.copies):
        str_out.append(str_line.format('Transform #',j + 1))
        str_out.append(str_line.format('RMSD',c.rmsd))
        rot = [str_r.format(k,*x) for k,x in enumerate(c.r.as_list_of_lists())]
        str_out.extend(rot)
        tran = str_t.format(*[x for xi in c.t.as_list_of_lists() for x in xi])
        str_out.append(tran)
        str_out.append('~ '*20)
    str_out.pop()
    str_out = '\n'.join(str_out)
    return str_out

  def show_ncs_selections(self,prefix=''):
    """
    Return NCS selection strings as a string, for printing

    Args:
     prefix (str): a string to be added, padding the output, at the left of
       each line
    """
    str_out = ['\n{}NCS selections:'.format(prefix),'-'*51]
    str_line = prefix + '{:<25s}:   {}'
    for i, gr in enumerate(self.ncs_restraints_group_list):
      str_out.append(str_line.format('Group #',i))
      str_out.append(str_line.format('Master selection string',gr.master_str_selection))
      for c in gr.copies:
        str_out.append(str_line.format('Copy selection string',c_str))
    transforms_info = self.show_transform_info(prefix)
    str_out.append(transforms_info)
    str_out.append('. '*26)
    str_out = '\n'.join(str_out)
    return str_out

  def __repr__(self,prefix=''):
    """ print NCS object info, padded with "prefix" on the left """
    str_out = [self.show_search_parameters_values(prefix)]
    str_out.append(self.show_chains_info(prefix))
    # print transforms
    str_out = '\n'.join(str_out)
    return str_out

def get_chain_and_ranges(hierarchy):
  """
  Helper function to get info from hierarchy to create spec object.
  Used in get_ncs_info_as_spec
  hierarchy is already selected
  """
  c_ids = [c.id for c in hierarchy.only_model().chains()]
  assert len(set(c_ids)) == 1
  c_id = c_ids[0]
  c = hierarchy.only_model().chains()[0]
  rgs_all = []
  for chain in hierarchy.chains():
    rgs_all +=list(chain.residue_groups())
  ranges = []
  in_range = False
  first_id = None
  last_id = None
  for rg in rgs_all:
    resseq = int(rg.resseq)
    if in_range:
      if resseq - last_id > 1:
        # terminate range, start new one
        ranges.append([first_id, last_id])
        first_id = resseq
        last_id = resseq
      else:
        # continue range
        last_id = resseq
    else:
      # start new range
      first_id = resseq
      last_id = resseq
      in_range = True
  # dumping rest:
  ranges.append([first_id, last_id])
  return c_id, ranges, len(c.residue_groups())

def format_80(s):
  """
  Split string that is longer than 80 characters to several lines

  Args:
    s (str)

  Returns:
    ss (str): formatted string
  """
  i = 0
  ss = ''
  for x in s:
    ss += x
    i += 1
    if i == 80:
      i = 0
      ss += ' \\ \n'
  return ss

def inverse_transform(r,t):
  r_inv = r.transpose()
  t_inv = - r_inv*t
  return r_inv,t_inv

def get_center_orth(xyz,selection):
  """
  Compute the center of coordinates of selected coordinates

  Args:
    xyz (flex.vec3_double): Atoms coordinates
    selection (flex.bool): Atom selection array

  Returns:
    (tuple) center of coordinates for the selected coordinates
    Returns (-100,-100,-100) when selection is bad
  """
  try:
    new_xyz = xyz.select(selection)
    mean = new_xyz.mean()
  except RuntimeError:
    mean = (-100,-100,-100)
  return mean
