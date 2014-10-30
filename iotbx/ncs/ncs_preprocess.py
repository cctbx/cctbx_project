from __future__ import division
from mmtbx.ncs.ncs_from_pdb import get_ncs_object_from_pdb
from mmtbx.ncs.ncs_utils import apply_transforms
from scitbx.array_family import flex
from scitbx.math import superpose
from mmtbx.ncs import ncs_search
from libtbx.utils import Sorry
from libtbx.phil import parse
import iotbx.pdb.hierarchy
from mmtbx.ncs import ncs
from scitbx import matrix
from iotbx import pdb
import string
import math
import re
import os
import sys


# parameters for manual specification of NCS - ASU mapping
master_phil = parse("""
  ncs_group
    .multiple = True
    {
    master_selection = ''
      .type = str
      .help = 'Residue selection string for the complete master NCS copy'
    copy_selection = ''
      .type = str
      .help = 'Residue selection string for each NCS copy location in ASU'
      .multiple = True
  }
  """)


class ncs_group_object(object):

  def __init__(self):
    """
    process MTRIX, BOIMT and PHIL parameters and produce an object
    with information for NCS refinement

    Argument:
    transform_info:  an object produced from the PDB MTRIX or BIOMT records
    ncs_refinement_params: an object produced by the PHIL parameters
    """
    self.total_asu_length = None
    # iselection maps, each master ncs to its copies position in the asu
    self.ncs_to_asu_map = {}
    # iselection maps of each ncs copy to its master ncs
    self.asu_to_ncs_map = {}
    # keys are items in ncs_chain_selection, values are lists of selection str
    self.ncs_to_asu_selection = {}
    self.ncs_copies_chains_names = {}
    self.tr_id_to_selection = {}
    # dictionary of transform names, same keys as ncs_to_asu_map
    self.number_of_ncs_groups = 0
    self.ncs_group_map = {}
    # map transform name (s1,s2,...) to transform object
    self.ncs_transform = {}
    # list of transform to chain assignment
    self.transform_chain_assignment = []
    # map transform to list of master ncs parts in its ncs groups
    self.transform_to_ncs = {}
    # master ncs and non-ncs selection in a string and a flex.bool types
    self.ncs_atom_selection = None
    self.ncs_selection_str = ''
    # iselection of all part in ASU that are not related via NCS operators
    self.non_ncs_region_selection = flex.size_t([])
    self.all_master_ncs_selections = flex.size_t([])
    # list of selection strings of master NCS
    self.ncs_chain_selection = []
    # unique chains or selection identifiers
    self.model_unique_chains_ids = tuple()
    self.selection_ids = set()
    # transform application order
    self.model_order_chain_ids = []
    self.transform_to_be_used = set()
    # order of transforms  - used when concatenating or separating them
    self.transform_order = []
    # keep hierarchy for writing (To have a source of atoms labels)
    self.hierarchy = None
    # residues common to NCS copies. Used for .spec representation
    self.common_res_dict = {}
    # flag indicating if ncs operation found
    self.found_ncs_transforms = False
    # Collect messages, recommendation and errors
    self.messages = ''
    self.write_messages = False
    self.log = sys.stdout

  def preprocess_ncs_obj(self,
                         pdb_hierarchy_inp=None,
                         pdb_inp=None,
                         transform_info=None,
                         rotations = None,
                         translations = None,
                         ncs_selection_params = None,
                         ncs_phil_groups = None,
                         file_name=None,
                         file_path='',
                         spec_file_str='',
                         spec_source_info='',
                         cif_string = '',
                         quiet=True,
                         spec_ncs_groups=None,
                         pdb_string=None,
                         use_cctbx_find_ncs_tools=True,
                         use_simple_ncs_from_pdb=True,
                         use_minimal_master_ncs=True,
                         max_rmsd=0.02,
                         write_messages=False,
                         process_similar_chains=True,
                         min_percent=0.95,
                         min_contig_length=10,
                         log=sys.stdout,
                         check_atom_order=False,
                         exclude_misaligned_residues=False,
                         max_dist_diff=4.0):
    """
    Select method to build ncs_group_object

    order of implementation:
    1) rotations,translations
    2) transform_info
    3) ncs_selection_params
    4) ncs_phil_groups
    5) spec file
    6) iotbx.pdb.hierarchy.input object

    Args:
    -----
      pdb_hierarchy_inp: iotbx.pdb.hierarchy.input
      transform_info: object containing MTRIX or BIOMT transformation info
      rotations: matrix.sqr 3x3 object
      translations: matrix.col 3x1 object
      ncs_selection_params: Phil parameters
        Phil structure
           ncs_group (multiple)
           {
             master_selection = ''
             copy_selection = ''   (multiple)
           }
      ncs_phil_groups: a list of ncs_groups_container object, containing
        master NCS selection and a list of NCS copies selection
      file_name: (str) .ncs_spec or .mmcif  or .pdb file name
      file_path: (str)
      spec_file_str: (str) spec format data
      spec_source_info:
      quiet: (bool) When True -> quiet output when processing files
      spec_ncs_groups: ncs_groups object as produced by simple_ncs_from_pdb
      cif_string: (str) string of cif type data
      use_cctbx_find_ncs_tools: (bool) Enable using of chain base NCS search
      use_simple_ncs_from_pdb: (bool) Enable using use_simple_ncs_from_pdb
      use_minimal_master_ncs: (bool) use maximal or minimal common chains
        in master ncs groups
      max_rmsd (float): limit of rms difference between chains to be considered
        as copies
      write_messages (bool): When True, right messages to log
        nearly the same length (but not exactly the same) and are NCS related.
        Raise error if NCS relations are not found
      process_similar_chains (bool): When True, process chains that are close
        in length without raising errors
      min_percent (float): Threshold for similarity between chains
        similarity define as:
        (number of matching res) / (number of res in longer chain)
      min_contig_length (int): minimum length of matching chain segments
      check_atom_order (bool): check atom order in matching residues.
        When False, matching residues with different number of atoms will be
        excluded from matching set
      exclude_misaligned_residues (bool): check and exclude individual residues
        alignment quality
      max_dist_diff (float): max allow distance difference between pairs of matching
        atoms of two residues
    """
    extension = ''
    self.write_messages = write_messages
    self.check_atom_order = check_atom_order
    self.log = log
    if file_name: extension = os.path.splitext(file_name)[1]
    if pdb_hierarchy_inp:
      msg = 'pdb_hierarchy_inp is not iotbx.pdb.hierarchy.input object\n'
      assert isinstance(pdb_hierarchy_inp,iotbx.pdb.hierarchy.input),msg
    elif pdb_inp:
      ph = pdb_inp.construct_hierarchy()
      pdb_hierarchy_inp = iotbx.pdb.hierarchy.input_hierarchy_pair(pdb_inp,ph)
    if extension.lower() in ['.pdb','.cif', '.mmcif', '.gz']:
      pdb_hierarchy_inp = pdb.hierarchy.input(file_name=file_name)
    elif (not pdb_hierarchy_inp) and (pdb_string or cif_string):
      if pdb_string: input_str = pdb_string
      else: input_str = cif_string
      pdb_hierarchy_inp = pdb.hierarchy.input(pdb_string=input_str)
    if pdb_hierarchy_inp and (not transform_info):
      transform_info = pdb_hierarchy_inp.input.process_mtrix_records(eps=0.01)
      if transform_info.as_pdb_string() == '': transform_info = None
      else:
        transform_info = insure_identity_is_in_transform_info(transform_info)
    if transform_info or rotations:
      if ncs_only(transform_info) or rotations:
        self.build_ncs_obj_from_pdb_ncs(
          pdb_hierarchy_inp = pdb_hierarchy_inp,
          rotations=rotations,
          translations=translations,
          transform_info=transform_info)
      else:
        # in the case that all ncs copies are in pdb
        assert [use_cctbx_find_ncs_tools,use_simple_ncs_from_pdb].count(True)>0
        self.build_ncs_obj_from_pdb_asu(
            pdb_hierarchy_inp=pdb_hierarchy_inp,
            use_cctbx_find_ncs_tools=use_cctbx_find_ncs_tools,
            use_simple_ncs_from_pdb=use_simple_ncs_from_pdb,
            use_minimal_master_ncs=use_minimal_master_ncs,
            process_similar_chains=process_similar_chains,
            min_percent=min_percent,
            min_contig_length=min_contig_length,
            max_rmsd=max_rmsd,
            exclude_misaligned_residues=exclude_misaligned_residues,
            max_dist_diff=max_dist_diff)
    elif ncs_selection_params or ncs_phil_groups:
      self.build_ncs_obj_from_phil(
        ncs_selection_params=ncs_selection_params,
        ncs_phil_groups=ncs_phil_groups,
        pdb_hierarchy_inp=pdb_hierarchy_inp,
        min_percent=min_percent)
    elif extension.lower() == '.ncs_spec' or \
            spec_file_str or spec_source_info or spec_ncs_groups:
      self.build_ncs_obj_from_spec_file(
        file_name=file_name,
        file_path=file_path,
        file_str=spec_file_str,
        source_info=spec_source_info,
        pdb_hierarchy_inp=pdb_hierarchy_inp,
        spec_ncs_groups=spec_ncs_groups,
        quiet=quiet)
    elif pdb_hierarchy_inp:
      self.build_ncs_obj_from_pdb_asu(
        pdb_hierarchy_inp=pdb_hierarchy_inp,
        use_cctbx_find_ncs_tools=use_cctbx_find_ncs_tools,
        use_simple_ncs_from_pdb=use_simple_ncs_from_pdb,
        use_minimal_master_ncs=use_minimal_master_ncs,
        min_percent=min_percent,
        min_contig_length=min_contig_length,
        process_similar_chains=process_similar_chains,
        max_rmsd=max_rmsd,
        exclude_misaligned_residues=exclude_misaligned_residues,
        max_dist_diff=max_dist_diff)
    else:
      raise Sorry('Please provide one of the supported input')
    # error handling
    self.found_ncs_transforms = (len(self.transform_to_be_used) > 0)
    if write_messages:
      if self.found_ncs_transforms == 0:
        print >> log,'No NCS relation were found !!!\n'
      if self.messages != '':
        print >> log,self.messages

  def build_ncs_obj_from_pdb_ncs(self,
                                 pdb_hierarchy_inp,
                                 transform_info=None,
                                 rotations = None,
                                 translations = None):
    """
    Build transforms objects and NCS <-> ASU mapping using PDB file containing
    a single NCS copy and MTRIX  or BIOMT records

    Args:
      pdb_hierarchy_inp : iotbx.pdb.hierarchy.input
      transform_info : an object containing MTRIX or BIOMT transformation info
      rotations : matrix.sqr 3x3 object
      translations : matrix.col 3x1 object
    """
    self.collect_basic_info_from_pdb(pdb_hierarchy_inp=pdb_hierarchy_inp)
    msg = 'No NCS transform information\n'
    assert bool(transform_info or (rotations and translations)), msg
    if rotations:
      # add rotations,translations to ncs_refinement_groups
      self.add_transforms_to_ncs_refinement_groups(
        rotations=rotations,
        translations=translations)
    else:
      # use only MTRIX/BIOMT records from PDB
      self.process_pdb(transform_info=transform_info)
    self.transform_chain_assignment = get_transform_order(self.transform_to_ncs)
    self.ncs_copies_chains_names = self.make_chains_names(
      transform_assignment=self.transform_chain_assignment,
      unique_chain_names = self.model_unique_chains_ids)
    # build self.ncs_to_asu_selection
    for k in self.transform_chain_assignment:
      selection_str = 'chain ' + self.ncs_copies_chains_names[k]
      key =  k.split('_')[0]
      self.tr_id_to_selection[k] = (key,selection_str)
      self.ncs_to_asu_selection = add_to_dict(
        d=self.ncs_to_asu_selection,k=key,v=selection_str)
    # add the identity case to tr_id_to_selection
    for key in  self.ncs_copies_chains_names.iterkeys():
      if not self.tr_id_to_selection.has_key(key):
        sel = key.split('_')[0]
        self.tr_id_to_selection[key] = (sel,sel)
    self.finalize_pre_process(pdb_hierarchy_inp=pdb_hierarchy_inp)

  def build_ncs_obj_from_phil(self,
                              ncs_selection_params = None,
                              ncs_phil_groups = None,
                              pdb_hierarchy_inp = None,
                              process_similar_chains=True,
                              min_percent=0.95):
    """
    Build transforms objects and NCS <-> ASU mapping using phil selection
    strings and complete ASU

    Args:
      ncs_selection_params : Phil parameters
      pdb_hierarchy_inp : iotbx.pdb.hierarchy.input
      process_similar_chains (bool): When True, process chains that are close
        in length without raising errors
      min_percent (float): similarity between chains when looking for
        ncs relations

    Phil structure
    ncs_group (multiple)
    {
      master_selection = ''
      copy_selection = ''   (multiple)
    }
    """

    if not process_similar_chains:
      min_percent = 1.0
    # process params
    if ncs_selection_params:
      if isinstance(ncs_selection_params,str):
        ncs_selection_params = parse(ncs_selection_params)
      phil_param =  master_phil.fetch(
        source=ncs_selection_params,track_unused_definitions=True)
      working_phil = phil_param[0].extract()
      assert  phil_param[1] == [],'Check phil parameters...\n'
      ncs_phil_groups = working_phil.ncs_group
    else:
      assert not (ncs_phil_groups is None)
    assert self.ncs_selection_str == ''
    unique_selections = set()
    transform_sn = 0
    ncs_group_id = 0
    # populate ncs selection and ncs to copies location
    for group in ncs_phil_groups:
      gns = group.master_selection
      self.ncs_chain_selection.append(gns)
      unique_selections = uniqueness_test(unique_selections,gns)
      ncs_group_id += 1
      transform_sn += 1
      self.add_identity_transform(
        ncs_selection=gns,
        ncs_group_id=ncs_group_id,
        transform_sn=transform_sn)
      key = format_num_as_str(transform_sn)
      # update with identity transform
      self.update_ncs_copies_chains_names(
            masters = gns,copies = gns, tr_id = key)
      self.update_tr_id_to_selection(gns,gns,key)
      asu_locations = []
      for asu_select in group.copy_selection:
        unique_selections = uniqueness_test(unique_selections,asu_select)
        r, t, rmsd, msg = get_rot_trans(
          ph=pdb_hierarchy_inp,
          master_selection=gns,
          copy_selection=asu_select,
          max_rmsd=100,
          min_percent=min_percent)
        self.messages += msg
        if r.is_zero():
          msg = 'Master NCS and Copy are very poorly related, check selection.'
          self.messages += msg + '\n'
        asu_locations.append(asu_select)
        transform_sn += 1
        key = format_num_as_str(transform_sn)
        self.update_tr_id_to_selection(gns,asu_select,key)
        tr = Transform(
          rotation = r,
          translation = t,
          serial_num = transform_sn,
          coordinates_present = True,
          ncs_group_id = ncs_group_id,
          rmsd=rmsd)
        # Update ncs_group dictionary and transform_to_ncs list
        self.build_transform_dict(
          transform_id = key,
          transform = tr,
          selection_id = gns)
        self.ncs_group_map = update_ncs_group_map(
          ncs_group_map=self.ncs_group_map,
          ncs_group_id = ncs_group_id,
          selection_ids = gns,
          transform_id = key)
        assert not self.ncs_transform.has_key(key)
        self.ncs_transform[key] = tr
        self.selection_ids.add(gns)
        self.update_ncs_copies_chains_names(
            masters = gns,copies = asu_select, tr_id = key)
      self.ncs_to_asu_selection[gns] = asu_locations
      self.number_of_ncs_groups = ncs_group_id

    self.ncs_selection_str = '('+ self.ncs_chain_selection[0] +')'
    for i in range(1,len(self.ncs_chain_selection)):
      self.ncs_selection_str += ' or (' + self.ncs_chain_selection[i] + ')'

    self.transform_chain_assignment = get_transform_order(self.transform_to_ncs)
    self.finalize_pre_process(pdb_hierarchy_inp=pdb_hierarchy_inp)

  def build_ncs_obj_from_pdb_asu(self,pdb_hierarchy_inp,
                                 use_cctbx_find_ncs_tools=True,
                                 use_simple_ncs_from_pdb=True,
                                 use_minimal_master_ncs=True,
                                 max_rmsd=0.5,
                                 process_similar_chains=True,
                                 min_contig_length=10,
                                 min_percent=0.95,
                                 exclude_misaligned_residues=False,
                                 max_dist_diff=4.0):
    """
    Build transforms objects and NCS <-> ASU mapping from a complete ASU
    Note that the MTRIX record are ignored, they are produced in the
    process of identifying the master NCS

    Args::
      pdb_hierarchy_inp : iotbx.pdb.hierarchy.input
      use_cctbx_find_ncs_tools (bool): Enable using of chain base NCS search
      use_simple_ncs_from_pdb (bool): Enable using use_simple_ncs_from_pdb
      use_minimal_master_ncs (bool): indicate if using maximal or minimal NCS
        groups
      max_rmsd (float): limit of rms difference between chains to be considered
        as copies
      process_similar_chains (bool): When True, process chains that are close
      in length without raising errors
      min_contig_length (int): minimum length of matching chain segments
      min_percent (float): Threshold for similarity between chains
        similarity define as:
        (number of matching res) / (number of res in longer chain)
      exclude_misaligned_residues (bool): check and exclude individual residues
      alignment quality
      max_dist_diff (float): max allow distance difference between pairs of matching
        atoms of two residues
    """
    if not process_similar_chains:
      min_percent = 1.0
      min_contig_length = 100000
      self.check_atom_order = True
    if use_cctbx_find_ncs_tools:
      group_dict = ncs_search.find_ncs_in_hierarchy(
        ph=pdb_hierarchy_inp.hierarchy,
        min_contig_length=min_contig_length,
        min_percent=min_percent,
        use_minimal_master_ncs=use_minimal_master_ncs,
        rmsd_eps=max_rmsd,
        write=self.write_messages,
        log=self.log,
        check_atom_order=self.check_atom_order,
        exclude_misaligned_residues=exclude_misaligned_residues,
        max_dist_diff=max_dist_diff)
      # process atom selections
      self.total_asu_length = pdb_hierarchy_inp.hierarchy.atoms().size()
      self.build_ncs_obj_from_group_dict(group_dict,pdb_hierarchy_inp)
    elif use_simple_ncs_from_pdb:
      # Fixme: Remove this after testing simple ncs from pdb
      ncs_object = get_ncs_object_from_pdb(
        pdb_inp=pdb_hierarchy_inp.input,
        hierarchy=pdb_hierarchy_inp.hierarchy)
      if ncs_object:
        spec_ncs_groups = ncs_object.ncs_groups()
        self.build_ncs_obj_from_spec_file(
          pdb_hierarchy_inp=pdb_hierarchy_inp,
          spec_ncs_groups=spec_ncs_groups)

  def build_ncs_obj_from_group_dict(self,group_dict,pdb_hierarchy_inp):
    """
    Use group_dict to build ncs object

    Args:
      pdb_hierarchy_inp : iotbx.pdb.hierarchy.input
      group_dict (dict):
        keys: tuple of master chain IDs
        values: NCS_groups_container objects with Attributes:
          iselections (list of lists of flex.size_t):
            selection array for the complete ASU
          residue_index_list (list): list of list of matching residues indices
          copies (list of lists):List of lists of the ncs copies chain IDs
          transforms (list of transform objects):
            object with attributes:
              rotation : Rotation matrix object
              translation: Translation matrix object
              serial_num : (int) Transform serial number
              coordinates_present (bool): True when coordinates are presents
              ncs_group_id (int): group ID of the group containing this transform
              rmsd (float): RMS distance between ncs copies
    """
    if hasattr(pdb_hierarchy_inp,"hierarchy"):
      ph = pdb_hierarchy_inp.hierarchy
    else:
      ph = pdb_hierarchy_inp
    self.ncs_atom_selection = flex.bool([False]*self.total_asu_length)
    chains_info = ncs_search.get_chains_info(ph)
    ncs_related_atoms = flex.bool([False]*self.total_asu_length)
    sorted_group_keys = sorted(group_dict)
    for gr_n,key in enumerate(sorted_group_keys):
      ncs_gr = group_dict[key]
      transform_id = set()
      # get all master selection string
      m_all_list = [x for ix in ncs_gr.iselections[0] for x in list(ix)]
      m_all_list.sort()
      m_all_isel = flex.size_t(m_all_list)
      all_m_select_str = selection_string_from_selection(
        ph,m_all_isel,chains_info=chains_info)
      self.ncs_to_asu_selection[all_m_select_str] = []
      #
      for i in xrange(len(ncs_gr.copies)):
        # iterate of ncs copies in group
        tr = ncs_gr.transforms[i]
        tr_id = format_num_as_str(tr.serial_num)
        self.ncs_transform[tr_id] = tr
        for j in xrange(len(ncs_gr.copies[i])):
          # iterate over chains in ncs copy
          m_isel = ncs_gr.iselections[0][j]
          m_ch_id = ncs_gr.copies[0][j]
          m_select_str = selection_string_from_selection(
            ph,m_isel,chains_info=chains_info)
          c_isel = ncs_gr.iselections[i][j]
          c_select_str = selection_string_from_selection(
            ph,c_isel,chains_info=chains_info)
          transform_id.add(tr_id)
          key0 = 'chain {}_{}'.format(m_ch_id,tr_id)
          key1 = m_select_str
          key2 = key1 + '_' + tr_id
          self.asu_to_ncs_map[key1] = m_isel
          self.ncs_to_asu_map[key2] = c_isel
          self.tr_id_to_selection[key0] = (m_select_str,c_select_str)
          self.selection_ids.add(m_select_str)
          self.update_ncs_copies_chains_names(
            masters = m_select_str,
            copies = c_select_str,
            tr_id = tr_id)
          self.ncs_group_map = update_ncs_group_map(
            ncs_group_map=self.ncs_group_map,
            ncs_group_id = gr_n + 1,
            selection_ids = m_select_str,
            transform_id = tr_id)
          if i == 0:
            # master copy
            master_sel = flex.bool(self.total_asu_length,m_isel)
            self.ncs_atom_selection |=  master_sel
            ncs_related_atoms |=  master_sel
          else:
            # non-identity transforms
            self.transform_to_be_used.add(tr.serial_num)
            # example key: "chain A_s002"
            self.transform_to_ncs = add_to_dict(
              d=self.transform_to_ncs,k=tr_id,v=key2)
            copy_sel = flex.bool(self.total_asu_length,c_isel)
            ncs_related_atoms |=  copy_sel
        # Get complete master and copies selections
        if i != 0:
          c_all_list = [x for ix in ncs_gr.iselections[i] for x in list(ix)]
          c_all_list.sort()
          c_all_isel = flex.size_t(c_all_list)
          c_select_str = selection_string_from_selection(
            ph,c_all_isel,chains_info=chains_info)
          self.ncs_to_asu_selection[all_m_select_str].append(c_select_str)
    #
    self.number_of_ncs_groups = len(group_dict)
    ncs_selection_str_list = []
    selection_ids = sorted(self.selection_ids)
    for sel in selection_ids:
      ncs_selection_str_list.append('(' + sel + ')')
    self.ncs_selection_str = ' or '.join(ncs_selection_str_list)
    self.transform_chain_assignment = get_transform_order(self.transform_to_ncs)

    self.all_master_ncs_selections = self.ncs_atom_selection.iselection(True)
    # add the ncs_atom_selection all the regions that are not NCS related
    self.ncs_atom_selection = self.ncs_atom_selection | ~ncs_related_atoms
    self.finalize_pre_process(pdb_hierarchy_inp=pdb_hierarchy_inp)

  def build_ncs_obj_from_spec_file(self,file_name=None,file_path='',
                                   file_str='',source_info="",quiet=True,
                                   pdb_hierarchy_inp=None,
                                   spec_ncs_groups=None,
                                   join_same_spec_groups = True):
    """
    read .spec files and build transforms object and NCS <-> ASU mapping

    Arguments:
    spec file information using a file, file_str (str) of a source_info (str)
    quite: (bool) quite output when True
    pdb_hierarchy_inp: pdb hierarchy input
    spec_ncs_groups: ncs_groups object
    join_same_spec_groups: (bool) True: combine groups with similar transforms
    """
    if not spec_ncs_groups: spec_ncs_groups = []
    if (not bool(spec_ncs_groups)) and (file_name or file_str):
      fn = ""
      if file_name:
        fn = os.path.join(file_path,file_name)
      lines = file_str.splitlines()
      assert lines != [] or os.path.isfile(fn)
      spec_object = ncs.ncs()
      spec_object.read_ncs(
        file_name=fn,lines=lines,source_info=source_info,quiet=quiet)
      spec_ncs_groups = spec_object.ncs_groups()

    if isinstance(spec_ncs_groups, ncs.ncs):
      spec_ncs_groups = spec_ncs_groups.ncs_groups()
    if spec_ncs_groups:
      transform_sn = 0
      ncs_group_id = 0
      for gr in spec_ncs_groups:
        # create selection
        spec_group_list =get_ncs_group_selection(gr.chain_residue_id())
        gs = spec_group_list[0]
        if join_same_spec_groups:
          # leave groups with the same transforms separate
          group_exist =self.look_and_combine_groups(gr,spec_group_list)
          if group_exist: continue
        ncs_group_id += 1
        self.ncs_chain_selection.append(gs)
        asu_locations = []
        for i,ncs_copy_select in enumerate(spec_group_list):
          # invert transform - the rotation in gr is from the copy to the master
          r = gr.rota_matrices()[i]
          t = gr.translations_orth()[i]
          r,t = inverse_transform(r,t)
          rmsd = round(gr.rmsd_list()[i],4)
          transform_sn += 1
          key = format_num_as_str(transform_sn)
          self.update_tr_id_to_selection(gs,ncs_copy_select,key)
          if not is_identity(r,t):
            asu_locations.append(ncs_copy_select)
          tr = Transform(
            rotation = r,
            translation = t,
            serial_num = transform_sn,
            coordinates_present = True,
            ncs_group_id = ncs_group_id,
            rmsd=rmsd)
          # Update ncs_group dictionary and transform_to_ncs list
          self.build_transform_dict(
            transform_id = key,
            transform = tr,
            selection_id = gs)
          self.ncs_group_map = update_ncs_group_map(
            ncs_group_map=self.ncs_group_map,
            ncs_group_id = ncs_group_id,
            selection_ids = gs,
            transform_id = key)
          assert not self.ncs_transform.has_key(key)
          self.ncs_transform[key] = tr
          self.selection_ids.add(gs)
          self.update_ncs_copies_chains_names(
            masters = gs,copies = ncs_copy_select,tr_id = key)
        self.ncs_to_asu_selection[gs] = asu_locations
      self.number_of_ncs_groups = ncs_group_id

      self.ncs_selection_str = '('+ self.ncs_chain_selection[0] +')'
      for i in range(1,len(self.ncs_chain_selection)):
        self.ncs_selection_str += ' or (' + self.ncs_chain_selection[i] + ')'
      self.transform_chain_assignment=get_transform_order(self.transform_to_ncs)
      self.finalize_pre_process(pdb_hierarchy_inp=pdb_hierarchy_inp)

  def look_and_combine_groups(self,gr_new,spec_group_list):
    """
    In spec files groups of different masters and copies listed separately,
    even if they have the same rotation/translation and can be combined.
    This function combines them, updates the relevant object attributes and
    returns True/False to indicate if group found

    Args:
      spec_group_list (list): selection string for each ncs copy in the group
      gr_new (object): ncs group object

    Returns:
      found_same_group (bool): indicate if group found
    """
    gs_new = spec_group_list[0]
    found_same_group = False
    gr_r_list = gr_new.rota_matrices()
    gr_t_list = gr_new.translations_orth()
    # in spec files transforms are inverted
    gr_new_list = [inverse_transform(r,t) for (r,t) in zip(gr_r_list,gr_t_list)]
    for k,[_, tr_set] in self.ncs_group_map.iteritems():
      # all transforms need to be the same to join
      if len(gr_r_list) != len(tr_set): continue
      same_transforms = [-1,]*len(tr_set)
      tr_list = list(tr_set)
      tr_list.sort()
      for tr_key1 in tr_list:
        r1 = self.ncs_transform[tr_key1].r
        t1 = self.ncs_transform[tr_key1].t
        for i,(r2,t2) in enumerate(gr_new_list):
          if same_transforms[i] != -1: continue
          same,transpose = ncs_search.is_same_transform(r1,t1,r2,t2)
          test = (same and (not transpose))
          if test:
            same_transforms[i] = i
            break
      found_same_group = (same_transforms.count(-1) == 0)
      # update dictionaries only if same group was found and then break
      if found_same_group:
        self.selection_ids.add(gs_new)
        asu_locations = []
        for i in same_transforms:
          transform_id = tr_list[i]
          ncs_copy_select = spec_group_list[i]
          key = gs_new + '_' + transform_id
          # look at self.ncs_copies_chains_names
          self.update_ncs_copies_chains_names(
            masters=gs_new, copies=ncs_copy_select, tr_id=transform_id)
          r,t = gr_new_list[i]
          self.update_tr_id_to_selection(gs_new,ncs_copy_select,transform_id)
          if not is_identity(r,t):
            self.transform_to_ncs = add_to_dict(
              d=self.transform_to_ncs,k=transform_id,v=key)
            asu_locations.append(ncs_copy_select)

        self.ncs_to_asu_selection[gs_new] = asu_locations
        self.ncs_group_map[k][0].add(gs_new)
        break
    return found_same_group

  def update_ncs_copies_chains_names(self,masters, copies, tr_id):
    masters = get_list_of_chains_selection(masters)
    copies = get_list_of_chains_selection(copies)
    for m,c in zip(masters,copies):
      self.ncs_copies_chains_names[m +'_' + tr_id] = c.replace('chain ','')

  def update_tr_id_to_selection(self,masters, copies,tr_id):
    """
    Args:
      masters: (str) selection of master ncs
      copies: (str) selection of copy
      tr_id: (str) string like "chain A_001" where 001 is the transform number
    """
    tr_keys = get_list_of_chains_selection(masters)
    master_selection_list = separate_selection_string(masters)
    copies_selection_list = separate_selection_string(copies)
    for k,c,m in zip(tr_keys,copies_selection_list,master_selection_list):
      key = k + '_' + tr_id
      self.tr_id_to_selection[key] = (m,c)

  def add_transforms_to_ncs_refinement_groups(self,rotations,translations):
    """
    Add rotation matrices and translations vectors
    to ncs_refinement_groups
    """
    assert len(rotations) == len(translations)
    assert not self.ncs_transform, 'ncs_transform should be empty\n'
    sn = {1}
    self.add_identity_transform(ncs_selection=self.ncs_selection_str)
    n = 1
    for (r,t) in zip(rotations,translations):
      # check if transforms are the identity transform
      if not is_identity(r,t):
        n += 1
        sn.add(n)
        tr_sn = n
      else:
        tr_sn = 1
      key = format_num_as_str(tr_sn)
      tr = Transform(
        rotation = r,
        translation = t,
        serial_num = tr_sn,
        coordinates_present = False,
        ncs_group_id = 1)
      self.ncs_transform[key] = tr
      for select in self.ncs_chain_selection:
        self.build_transform_dict(
          transform_id = key,
          transform = tr,
          selection_id = select)
        self.selection_ids.add(select)
        chain_id = get_list_of_chains_selection(select)
        assert len(chain_id) == 1
        chain_id = chain_id[0]
        self.tr_id_to_selection[chain_id + '_' + key] = (select,select)
      self.ncs_group_map = update_ncs_group_map(
        ncs_group_map=self.ncs_group_map,
        ncs_group_id = 1,
        selection_ids = self.ncs_chain_selection,
        transform_id = key)

  def collect_basic_info_from_pdb(self,pdb_hierarchy_inp):
    """
    Build chain selection string and collect chains IDs from pdb
    Consider that chains can be not continuous
    """
    if pdb_hierarchy_inp:
      model  = pdb_hierarchy_inp.hierarchy.models()[0]
      chain_ids = {x.id for x in model.chains()}
      # Collect order if chains IDs and unique IDs
      self.model_unique_chains_ids = tuple(sorted(chain_ids))
      model_order_ch_ids = [(x.id,x.atoms().size()) for x in model.chains()]
      ch_n_atoms = {x:None for x in self.model_unique_chains_ids}
      for (ch,n) in model_order_ch_ids:
        if ch_n_atoms[ch] is None:
          ch_n_atoms[ch] = [(0,n)]
        else:
          _,last_n = ch_n_atoms[ch][-1]
          ch_n_atoms[ch].append((last_n, last_n + n))
      for ch,n in model_order_ch_ids:
        selection_range = ch_n_atoms[ch].pop(0)
        self.model_order_chain_ids.append((ch,selection_range))
      s = ' or chain '.join(self.model_unique_chains_ids)
      self.ncs_selection_str = 'chain ' + s
      assert self.ncs_chain_selection == []
      self.ncs_chain_selection =\
        ['chain ' + s for s in self.model_unique_chains_ids]
      self.ncs_chain_selection.sort()

  def compute_ncs_asu_coordinates_map(self,pdb_hierarchy_inp):
    """ Calculates coordinates maps from ncs to asu and from asu to ncs """
    # check is coordinates maps already calculated
    t1 = not bool(self.ncs_atom_selection)
    t2 = not bool(self.asu_to_ncs_map)
    t3 = not bool(self.ncs_to_asu_map)
    if t1 and t2 and t3:
      temp = pdb_hierarchy_inp.hierarchy.atom_selection_cache()
      # check if pdb_hierarchy_inp contain only the master NCS copy
      pdb_length = len(pdb_hierarchy_inp.hierarchy.atoms())
      self.ncs_atom_selection = temp.selection(self.ncs_selection_str)
      ncs_length = self.ncs_atom_selection.count(True)
      # keep track on the asu copy number
      copy_count = {}
      if pdb_length > ncs_length:
        self.total_asu_length = pdb_length
        selection_ref = flex.bool([False]*pdb_length)
        for k in self.transform_chain_assignment:
          key =  k.split('_')[0]
          ncs_selection = temp.selection(key)
          if not self.asu_to_ncs_map.has_key(key):
            copy_count[key] = 0
            selection_ref = (selection_ref | ncs_selection)
            self.asu_to_ncs_map[key] = ncs_selection.iselection(True)
          else:
            copy_count[key] += 1
          # ncs_to_asu_selection is a list of all the copies of a master
          asu_copy_ref = self.ncs_to_asu_selection[key][copy_count[key]]
          asu_selection = temp.selection(asu_copy_ref)
          selection_ref = update_selection_ref(selection_ref,asu_selection)
          self.ncs_to_asu_map[k] = asu_selection.iselection(True)
        self.non_ncs_region_selection = (~selection_ref).iselection(True)
        # add the non ncs regions to the master ncs copy
        self.all_master_ncs_selections = self.ncs_atom_selection.iselection(True)
        self.ncs_atom_selection |= ~selection_ref
        assert set(self.non_ncs_region_selection).intersection(
          set(self.all_master_ncs_selections)) == set()
        # Check that all the NCS related copies have the same atom numbers

      elif pdb_length == ncs_length:
        # this case is when the pdb hierarchy contain only the master NCS copy
        self.total_asu_length = self.get_asu_length(temp)
        ns = [True]*pdb_length + [False]*(self.total_asu_length - pdb_length)
        self.ncs_atom_selection = flex.bool(ns)
        self.all_master_ncs_selections=self.ncs_atom_selection.iselection(True)
        sorted_keys = sorted(self.transform_to_ncs)
        for i,k in enumerate(sorted_keys):
          v = self.transform_to_ncs[k]
          for transform_key in v:
            key =  transform_key.split('_')[0]
            ncs_selection =flex.bool(self.total_asu_length,temp.iselection(key))
            if not self.asu_to_ncs_map.has_key(key):
              self.asu_to_ncs_map[key] = ncs_selection.iselection(True)
            # make the selection at the proper location at the ASU
            temp_iselection = self.asu_to_ncs_map[key] + ((i + 1) * ncs_length)
            asu_selection = flex.bool(self.total_asu_length,temp_iselection)
            self.ncs_to_asu_map[transform_key] = asu_selection.iselection(True)

  def add_identity_transform(self,ncs_selection,ncs_group_id=1,transform_sn=1):
    """    Add identity transform
    Argument:

    ncs_selection: (str) selection string for the NCS master copy
    ncs_group_id: (int) the NCS group ID
    transform_sn: (int) Over all transform serial number
    """
    transform_obj = Transform(
      rotation = matrix.sqr([1,0,0,0,1,0,0,0,1]),
      translation = matrix.col([0,0,0]),
      serial_num = transform_sn,
      coordinates_present = True,
      ncs_group_id = ncs_group_id)
    id_str = format_num_as_str(transform_sn)
    self.ncs_transform[id_str] = transform_obj
    self.build_transform_dict(
      transform_id = id_str,
      transform = transform_obj,
      selection_id = ncs_selection)
    self.selection_ids.add(ncs_selection)
    self.ncs_group_map = update_ncs_group_map(
      ncs_group_map=self.ncs_group_map,
      ncs_group_id = ncs_group_id,
      selection_ids = ncs_selection,
      transform_id = id_str)

  def process_pdb(self,transform_info):
    """
    Process PDB Hierarchy object
    Args:
      transform_info

    Remarks:
    The information on a chain in a PDB file does not have to be continuous.
    Every time the chain name changes in the pdb file, a new chain is added
    to the model, even if the chain ID already exist. so there model.
    chains() might contain several chains that have the same chain ID
    """
    transform_info_available = bool(transform_info) and bool(transform_info.r)
    if transform_info_available:
      ti = transform_info
      for (r,t,n,cp) in zip(ti.r,ti.t,ti.serial_number,ti.coordinates_present):
        n = int(n)
        key = format_num_as_str(n)
        tr = Transform(
          rotation = r,
          translation = t,
          serial_num = n,
          coordinates_present = cp,
          ncs_group_id = 1)
        for select in self.ncs_chain_selection:
          self.build_transform_dict(
            transform_id = key,
            transform = tr,
            selection_id = select)
          self.selection_ids.add(select)
        self.ncs_group_map = update_ncs_group_map(
          ncs_group_map=self.ncs_group_map,
          ncs_group_id = 1,
          selection_ids = self.ncs_chain_selection,
          transform_id = key)
        # if ncs selection was not provided in phil parameter
        assert not self.ncs_transform.has_key(key)
        self.ncs_transform[key] = tr

  def build_transform_dict(self,
                           transform_id,
                           transform,
                           selection_id):
    """
    Apply all non-identity transforms
    Build transform_to_ncs dictionary, which provides the location of the
    particular chains or selections in the NCS (example key: "s002")
    and updates transform_to_be_used set

    Args:
      transform_id (str): s001,s002...
      transform : transform object, containing information on transformation
      selection_id (str): NCS selection string
    """
    if not is_identity(transform.r,transform.t):
      self.transform_to_be_used.add(transform.serial_num)
      key = selection_id + '_' + format_num_as_str(transform.serial_num)
      # example key: "chain A_s002"
      self.transform_to_ncs = add_to_dict(
        d=self.transform_to_ncs,k=transform_id,v=key)

  def get_asu_length(self,atom_selection_cache):
    """" Collect the length of all ncs copies """
    asu_total_length = self.ncs_atom_selection.count(True)
    for k in self.transform_chain_assignment:
      key =  k.split('_')[0]
      ncs_selection = atom_selection_cache.selection(key)
      asu_total_length += ncs_selection.count(True)
    return asu_total_length

  def build_MTRIX_object(self,ncs_only=True):
    """
    Build a MTRIX object from ncs_group_object
    Used for testing
    """
    assert  self.number_of_ncs_groups < 2
    result = iotbx.pdb.mtrix_and_biomt_records_container()
    tr_dict = self.ncs_transform
    tr_sorted = sorted(tr_dict,key=lambda k:tr_dict[k].serial_num)
    for key in tr_sorted:
      transform = self.ncs_transform[key]
      r = transform.r
      t = transform.t
      identity_test = is_identity(r,t)
      cp = (not ncs_only) or identity_test
      result.add(
        r=r,
        t=t,
        coordinates_present=cp,
        serial_number=transform.serial_num)
    return result

  def make_chains_names(self,
                        transform_assignment,
                        unique_chain_names):
    """
    Create a dictionary names for the new NCS copies
    keys: (str) chain_name + '_' + serial_num
    values: (str) (one or two chr long)

    Chain names might repeat themselves several times in a pdb file
    We want copies of chains with the same name to still have the
    same name after similar BIOMT/MTRIX transformation

    Arguments:
    transform_assignment : (list) transformation assignment order
    unique_chain_names : (tuple) a set of unique chain names

    Returns:
    new_names : a dictionary. {'A_1': 'G', 'A_2': 'H',....} map a chain
    name and a transform number to a new chain name
    """
    if not transform_assignment or not unique_chain_names: return {}
    # create list of character from which to assemble the list of names
    # total_chains_number = len(i_transforms)*len(unique_chain_names)
    total_chains_number = len(transform_assignment)
    # start naming chains with a single letter
    chr_list1 = list(set(string.ascii_uppercase) - set(unique_chain_names))
    chr_list2 = list(set(string.ascii_lowercase) - set(unique_chain_names))
    chr_list1.sort()
    chr_list2.sort()
    new_names_list = chr_list1 + chr_list2
    # check if we need more chain names
    if len(new_names_list) < total_chains_number:
      n_names =  total_chains_number - len(new_names_list)
      # the number of character needed to produce new names
      chr_number = int(math.sqrt(n_names)) + 1
      # build character list
      chr_list = list(string.ascii_uppercase) + \
                 list(string.ascii_lowercase) + \
                 list(string.digits)
      # take only as many characters as needed
      chr_list = chr_list[:chr_number]
      extra_names = set([ x+y for x in chr_list for y in chr_list])
      # make sure not using existing names
      extra_names = list(extra_names - set(unique_chain_names))
      extra_names.sort()
      new_names_list.extend(extra_names)
    assert len(new_names_list) >= total_chains_number
    dictionary_values = new_names_list[:total_chains_number]
    # create the dictionary
    zippedlists = zip(transform_assignment,dictionary_values)
    new_names_dictionary = {x:y for (x,y) in zippedlists}
    # add the master NCS to dictionary
    tr_set  = {format_num_as_str(x) for x in self.transform_to_be_used}
    for k,v in self.ncs_group_map.iteritems():
      tr_str = (v[1] - tr_set)
      assert len(tr_str) == 1
      tr_str = tr_str.pop()
      for ch_sel in v[0]:
        if not ' or ' in ch_sel:
          new_names_dictionary[ch_sel+'_'+tr_str] = ch_sel.replace('chain ','')
    return new_names_dictionary

  def finalize_pre_process(self,pdb_hierarchy_inp=None):
    """
    Steps that are common to most method of transform info
    """
    if pdb_hierarchy_inp:
      self.compute_ncs_asu_coordinates_map(pdb_hierarchy_inp=pdb_hierarchy_inp)
      # keep hierarchy for writing
      self.hierarchy = pdb_hierarchy_inp.hierarchy
      self.set_common_res_dict()
    self.transform_order = sort_dict_keys(self.transform_to_ncs)

  def set_common_res_dict(self):
    """
    Build common residues list and related RMSD
    for use when writing spec files list of common residues
    for each chain - transform pair
    """
    # collect all master chain IDs
    sorted_keys = sort_dict_keys(self.ncs_copies_chains_names)
    only_master_ncs_in_hierarchy = False
    if self.ncs_atom_selection.count(True) == self.hierarchy.atoms().size():
      only_master_ncs_in_hierarchy = True
    sc = self.hierarchy.atom_selection_cache()
    #
    for key in sorted_keys:
      master_sel_str, ncs_sel_str = self.tr_id_to_selection[key]
      if only_master_ncs_in_hierarchy:
        # use master ncs for ncs copy residues indices
        copy_selection_indices = sc.selection(master_sel_str).iselection(True)
        rmsd = 0
      else:
        copy_selection_indices = sc.selection(ncs_sel_str).iselection(True)
        tr_num = key.split('_')[1]
        tr = self.ncs_transform[tr_num]
        rmsd = tr.rmsd
      # get continuous res ids
      range_list = []
      t_ph = self.hierarchy.select(copy_selection_indices).models()[0].chains()
      for chain in t_ph:
        res_id = []
        # for rs in chain.residues():
        for rs in chain.residue_groups():
          resid = rs.resid().strip()
          j = rs.resseq_as_int()
          # check if we have insertions
          if str(j) == resid:
            if res_id:
              # step larger than one residue -> close segment
              if (res_id[1] + 1) < j:
                range_list.append(res_id)
                res_id = [j,j]
              else:
                # increase segment range by one
                res_id[1] += 1
            else:
              # start new segment
              res_id = [j,j]
          else:
            # This representation does not handle insertions !!!
            msg = "Sequence may contain insertions and can't be "
            msg += "represented using only residue ID. \n"
            self.messages += msg
        if res_id and (res_id[1] == j):
          # close the last segment
          range_list.append(res_id)
      range_list.sort()
      self.common_res_dict[key] = ([range_list,copy_selection_indices],rmsd)

  def get_ncs_restraints_group_list(self):
    """    a list of ncs_restraint_group objects    """
    ncs_restraints_group_list = []
    group_id_list = sort_dict_keys(self.ncs_group_map)
    for k in group_id_list:
      v = self.ncs_group_map[k]
      master_isel = flex.size_t([])
      for key in sorted(list(v[0])):
        if self.asu_to_ncs_map.has_key(key):
          master_isel.extend(self.asu_to_ncs_map[key])
      new_nrg = NCS_restraint_group(master_isel)

      for tr in sorted(list(v[1])):
        if self.transform_to_ncs.has_key(tr):
          r = self.ncs_transform[tr].r
          t = self.ncs_transform[tr].t
          ncs_isel = flex.size_t([])
          for sel in self.transform_to_ncs[tr]:
            ncs_isel.extend(self.ncs_to_asu_map[sel])
          new_ncs_copy = NCS_copy(copy_iselection=ncs_isel, rot=r, tran=t)
          new_nrg.copies.append(new_ncs_copy)
      # compare master_isel_test and master_isel
      ncs_restraints_group_list.append(new_nrg)
    return ncs_restraints_group_list

  def update_using_ncs_restraints_group_list(self,ncs_restraints_group_list):
    """
    Update ncs_group_object rotations and transformations.

    Note that to insure proper assignment the ncs_restraints_group_list
    should be produced using the get_ncs_restraints_group_list method

    Args:
      ncs_restraints_group_list: a list of ncs_restraint_group objects
    """
    assert len(ncs_restraints_group_list) == len(self.ncs_group_map)
    group_id_list = sort_dict_keys(self.ncs_group_map)
    for k in group_id_list:
      v = self.ncs_group_map[k]
      nrg = ncs_restraints_group_list.pop(0)
      for tr in sorted(list(v[1])):
        if self.transform_to_ncs.has_key(tr):
          ncs_copy = nrg.copies.pop(0)
          self.ncs_transform[tr].r = ncs_copy.r
          self.ncs_transform[tr].t = ncs_copy.t
          # Test that the correct transforms are updated
          ncs_isel = flex.size_t([])
          for sel in self.transform_to_ncs[tr]:
            ncs_isel.extend(self.ncs_to_asu_map[sel])
          assert ncs_copy.iselection == ncs_isel

  def get_transform_records(self, file_name=None,
                          ncs_only=True,
                          pdb_hierarchy=None,
                          xrs=None,
                          fmodel=None,
                          crystal_symmetry=None,
                          mtrix=None,
                          biomt=None,
                          write=False,
                          log = sys.stdout):
    """
    Write to a file or prints transformation records.
    with or without PDB atoms and Cryst records.
    If no pdb_hierarchy, xray structure or fmodel are provided, the function
    will return only the MTRIX/BIOMT records

    Args:
      file_name: (str) output file name
      ncs_only: (bool) When False, the comple ASU will be printed (applicable
                only with MTRIX records)
      pdb_hierarchy: (pdb_hierarchy object)
      xrs: (xray structure) for crystal symmetry
      fmodel: (fmodel object)
      crystal_symmetry: crystal symmetry records
      mtrix: (bool) When True -> write MTRIX records
      biomt: (bool) When True -> write BIOMT records
      write: (bool) when False, will will not write to file or print

    Return:
      PDB string
    """
    if (not mtrix) and (not biomt):
      mtrix = True
      biomt = False
    assert bool(mtrix) == (not bool(biomt))
    if biomt: ncs_only = True
    mtrix_object = self.build_MTRIX_object(ncs_only=ncs_only)
    pdb_header_str = ''
    new_ph_str = ''
    transform_rec = ''
    #
    if fmodel:
      xrs = fmodel.xray_structure
    if xrs and self.hierarchy and (not pdb_hierarchy):
      pdb_hierarchy = self.hierarchy
      pdb_str = xrs.as_pdb_file()
      pdb_header_str = get_pdb_header(pdb_str)
      xyz = pdb_hierarchy.atoms().extract_xyz()
      new_xyz = xrs.sites_cart()
      if new_xyz.size() > xyz.size():
        ncs_only = True
        xrs = xrs.select(self.ncs_atom_selection.iselection(True))
        new_xyz = xrs.sites_cart()
      assert new_xyz.size() == xyz.size()
      pdb_hierarchy.atoms().set_xyz(new_xyz)
    if pdb_hierarchy:
      ph = pdb_hierarchy
      if xrs:
        crystal_symmetry = xrs.crystal_symmetry()
      pdb_str = ph.as_pdb_string(crystal_symmetry=crystal_symmetry)
      if not pdb_header_str:
       pdb_header_str = get_pdb_header(pdb_str)
      if ncs_only:
        new_ph = ph.select(self.ncs_atom_selection.iselection(True))
      else:
        msg = 'The complete ASU hierarchy need to be provided !!!\n'
        assert len(self.ncs_atom_selection) == len(ph.atoms()),msg
        new_ph = ph
      new_ph_str = new_ph.as_pdb_string(crystal_symmetry=None)
    #
    if mtrix:
      transform_rec = mtrix_object.as_pdb_string()
    elif biomt:
      transform_rec = mtrix_object.format_BOIMT_pdb_string()
    #
    if write:
      if file_name:
        f = open(file_name,'w')
        print >> f, pdb_header_str
        print >> f, transform_rec
        print >> f, new_ph_str
        f.close()
      else:
        print >> log,pdb_header_str
        print >> log,transform_rec
        print >> log,new_ph_str
    return '\n'.join([pdb_header_str,transform_rec,new_ph_str])

  def get_ncs_info_as_spec(
          self,
          pdb_hierarchy_asu=None,
          xrs=None,
          fmodel=None,
          write=False):
    """
    Returns ncs spec object and can prints ncs info in a ncs_spec,
    format_all_for_resolve or format_all_for_phenix_refine format

    Note that while ncs_groups can master ncs can be comprised from several
    chains, the spec groups can not. So groups with multiple chains in the
    master selection are splitted

    Note that spec format does not support insertions notation
    for example "resseq 49" will include "resid 49" and "resid 49A"

    Args:
      file_name: (str) output file name
      pdb_hierarchy: (pdb_hierarchy object)
      xrs: (xray structure) for crystal symmetry
      fmodel: (fmodel object)
      write: (bool) when False, will not write to file or print

    Return:
      spec_object
    """
    spec_object = ncs.ncs()

    if [bool(xrs),bool(pdb_hierarchy_asu)].count(True) == 0:
      if self.hierarchy and \
              (self.hierarchy.atoms().size() == self.total_asu_length):
        pdb_hierarchy_asu = self.hierarchy
    if fmodel:
      xrs = fmodel.xray_structure
    if xrs and (not pdb_hierarchy_asu):
      xyz = xrs.sites_cart()
    else:
      xyz = pdb_hierarchy_asu.atoms().extract_xyz()

    for gn,gr in self.ncs_group_map.iteritems():
      for gr_chains in gr[0]:
        gr_dict = {}
        # gr -> [master selection str, set of transforms]
        # Process one chain, in the master ncs, at a time
        for gr_chain in get_list_of_chains_selection(gr_chains):
          for s_str in gr[1]:
            # get ncs copy chain name, that corresponds to master and transform
            gr_key = gr_chain + '_' + s_str
            gr_dict[gr_key] = self.ncs_copies_chains_names[gr_key]
        sorted_keys = sort_dict_keys(gr_dict)
        center_orth = []
        rotations = []
        translations = []
        chain_residue_id_list = []
        chain_residue_range_list = []
        rmsd_list = []
        residues_in_common_list = []
        for k in sorted_keys:
          chain_id = gr_dict[k]
          chain_residue_id_list.append(chain_id)
          [range_list,ncs_sel],rmsd = self.common_res_dict[k]
          chain_residue_range_list.append(range_list)
          center_orth.append(get_center_orth(xyz,ncs_sel))
          transform_key = k.split('_')[1]
          tr = self.ncs_transform[transform_key]
          # use the rmsd of the ncs related atoms rather than the transform
          # rmsd_list.append(tr.rmsd)
          rmsd_list.append(rmsd)
          # in spec files transform is copy -> master, not master -> copy
          r,t = inverse_transform(tr.r,tr.t)
          rotations.append(r)
          translations.append(t)
          res_num = sum([y-x+1 for [x,y] in range_list])
          residues_in_common_list.append(res_num)
        # build group
        spec_object.import_ncs_group(
          center_orth = center_orth,
          ncs_rota_matr = rotations,
          trans_orth = translations,
          rmsd_list = rmsd_list,
          chain_residue_id = [chain_residue_id_list,chain_residue_range_list],
          residues_in_common_list = residues_in_common_list,
          ncs_domain_pdb = None)

    if write:
      spec_object.display_all(log=self.log)
      f=open("simple_ncs_from_pdb.resolve",'w')
      spec_object.format_all_for_resolve(log=self.log,out=f)
      f.close()
      f=open("simple_ncs_from_pdb.ncs",'w')
      spec_object.format_all_for_phenix_refine(log=self.log,out=f)
      f.close()
      f=open("simple_ncs_from_pdb.ncs_spec",'w')
      spec_object.format_all_for_group_specification(log=self.log,out=f)
      f.close()
    return spec_object

  def print_ncs_phil_param(self,write=False,log=sys.stdout):
    """
    Prints NCS information in the phil parameters format
    lines longer that 80 characters are folded

    Phil structure example:
      ncs_group {
        master_selection = 'chain A'
        copy_selection = 'chain C'
        copy_selection = 'chain E'
      }
      ncs_group {
        master_selection = 'chain B'
        copy_selection = 'chain D'
        copy_selection = 'chain F'
      }

    Args:
      write (bool): when true, print to log
      log : location of output, an open file or sys.stdout

    Returns:
      (str): NCS phil parameter string
    """
    groups = []
    for master, copies in self.ncs_to_asu_selection.iteritems():
      master = format_80(master)
      groups.append('ncs_group {')
      groups.append("  master_selection = {}".format(master))
      for cp in copies:
        cp = format_80(cp)
        groups.append("  copy_selection = {}".format(cp))
      groups.append('}')
    gr = '\n'.join(groups)
    gr += '\n'
    if write:
      print >> log,gr
    return gr

  def build_asu_hierarchy(self,
                          pdb_hierarchy,
                          round_coordinates=True):
    """
    Build ASU hierarchy

    Arguments:
    pdb_hierarchy: pdb hierarchy of the master NCS
    round_coordinates: (bool) round coordinates of new NCS copies,
                        for sites_cart constancy
    Return:
    ASU hierarchy
    """
    # Build only for PDB when there is a single NCS group
    assert self.number_of_ncs_groups < 2
    new_ph = pdb_hierarchy.deep_copy()
    ncs_restraints_group_list = self.get_ncs_restraints_group_list()
    new_sites = apply_transforms(
      ncs_coordinates = pdb_hierarchy.atoms().extract_xyz(),
      ncs_restraints_group_list = ncs_restraints_group_list,
      total_asu_length =  self.total_asu_length,
      extended_ncs_selection = flex.size_t_range(pdb_hierarchy.atoms().size()),
      round_coordinates = round_coordinates)
    model = new_ph.models()[0]
    tr_assignment_order = []
    for tr in self.transform_order:
      for (ch_id, (sel_start,sel_end)) in self.model_order_chain_ids:
        key = 'chain ' + ch_id
        tr_key  =  key + '_' + tr

        ncs_selection = self.asu_to_ncs_map[key][sel_start:sel_end]
        tr_assignment_order.append([tr_key,ncs_selection])

    for tr,ncs_selection in tr_assignment_order:
      new_part = pdb_hierarchy.select(ncs_selection).deep_copy()
      new_chain = iotbx.pdb.hierarchy.ext.chain()
      new_chain.id = self.ncs_copies_chains_names[tr]
      for res in new_part.residue_groups():
        new_chain.append_residue_group(res.detached_copy())
      model.append_chain(new_chain)
    new_ph.atoms().set_xyz(new_sites)
    return new_ph

def add_to_dict(d,k,v):
  if d.has_key(k):
    d[k].append(v)
  else:
    d[k] = [v]
  return d

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
      ss += ' \ \n'
  return ss

def get_residue_sequence(chain_ph):
  """
  Get a list of pairs (residue name, residue number) and the chain ID
  Exclude water

  Args:
    chain_ph (iotbx_pdb_hierarchy_ext): hierarchy of a single chain

  Returns:
    chain_id (str): Chain ID
    res_list_new (list of str): list of residues names
    resid_list_new (list of str): list of residues number
  """
  res_list_new = []
  resid_list_new = []
  for res_info in chain_ph.atom_groups():
    x = res_info.resname
    if x.lower() != 'hoh':
      # Exclude water
      res_list_new.append(x)
      resid_list_new.append(res_info.parent().resseq)
  #
  atoms = chain_ph.atoms()
  chain_id = atoms[0].chain().id
  return chain_id, res_list_new,resid_list_new

def selection_string_from_selection(pdb_hierarchy_inp,selection,
                                    chains_info=None):
  """
  Convert a selection array to a selection string

  Args:
    pdb_hierarchy_inp : iotbx.pdb.hierarchy.input (or iotbx.pdb.hierarchy)
    selection (flex.bool or flex.size_t)
    chains_info : object containing
      chains (str): chain IDs OR selections string
      res_name (list of str): list of residues names
      resid (list of str): list of residues sequence number, resid
      atom_names (list of list of str): list of atoms in residues
      atom_selection (list of list of list of int): the location of atoms in ph
      chains_atom_number (list of int): list of number of atoms in each chain

  Returns:
    sel_str (str): atom selection string
  """
  # Todo: consider moving to another location
  # create a hierarchy from the selection
  if hasattr(pdb_hierarchy_inp,"hierarchy"):
    pdb_hierarchy_inp = pdb_hierarchy_inp.hierarchy
  # pdb_hierarchy_inp is a hierarchy
  if isinstance(selection,flex.bool): selection = selection.iselection(True)
  selection_set = set(selection)
  sel_list = []
  if not chains_info:
    chains_info = ncs_search.get_chains_info(pdb_hierarchy_inp)
  chain_ids = sorted(chains_info)
  for ch_id in chain_ids:

    a_sel = {x for xi in chains_info[ch_id].atom_selection for x in xi}
    if len(a_sel) != chains_info[ch_id].chains_atom_number:
      s = ' and (not resname hoh)'
    else: s = ''
    ch_sel = 'chain {}{}'.format(ch_id,s)
    test_set = a_sel.intersection(selection_set)
    if not test_set: continue
    complete_ch_not_present = (test_set != a_sel)
    res_sel = []
    first_n = None
    pre_res_n = -10000
    if complete_ch_not_present:
      # collect continuous ranges of residues when possible
      res_len = len(chains_info[ch_id].resid)
      for i in xrange(res_len):
        # test that all atoms in residue are included in selection
        a_sel = set(chains_info[ch_id].atom_selection[i])
        test_set = a_sel.intersection(selection_set)
        if not test_set: continue
        all_atoms_present = (test_set == a_sel)
        res_id = chains_info[ch_id].resid[i]
        # ensure that insertion are not included if shouldn't
        next_res_id = '0'
        if i < (res_len - 1):
          next_res_id = chains_info[ch_id].resid[i+1]
        try:
          # check that res_id is a number and not insertion residue
          res_num = int(res_id)
          # int will fail if the next residue is insertion
          int(next_res_id)
        except ValueError:
          # res_id is an insertion type residue
          res_num = -10000
        if all_atoms_present:
          if res_num != -10000:
            # normal case
            if pre_res_n == -10000:
              # start new range
              first_n = res_num
              pre_res_n = res_num
            elif res_num == (pre_res_n + 1):
              pre_res_n += 1
            else:
              res_seq = resseq_string(first_n,pre_res_n)
              res_sel.append(res_seq)
              first_n = res_num
              pre_res_n = res_num
          else:
            # insertion in sequence
            res_sel,first_n,pre_res_n =update_res_sel(res_sel,first_n,pre_res_n)
            res_sel.append('resid ' + res_id )
        else:
          # not all residue's atoms are in selection
          res_sel,first_n,pre_res_n = update_res_sel(res_sel,first_n,pre_res_n)
          s = '(resid ' + res_id + ' and (name '
          # get present atoms
          atom_name = chains_info[ch_id].atom_names[i]
          dx = min(test_set)
          atom_name = [atom_name[x-dx] for x in test_set]
          atom_str = ' or name '.join(atom_name)
          res_sel.append(s + atom_str + '))')
      res_sel,first_n,pre_res_n = update_res_sel(res_sel,first_n,pre_res_n)
    s = get_clean_selection_string(ch_sel,res_sel)
    sel_list.append(s)
  # add parenthesis what selection is more than just a chain
  s_l = []
  sel_list.sort()
  for s in sel_list:
    if len(s) > 8:
      s = '(' + s + ')'
    s_l.append(s)
  sel_str = ' or '.join(s_l)
  return sel_str

def update_res_sel(res_sel,first_res_n,pre_res_n):
  """ update the residue selection list and markers of continuous section """
  if pre_res_n != -10000:
    res_seq = resseq_string(first_res_n,pre_res_n)
    first_res_n = None
    pre_res_n = -10000
    res_sel.append(res_seq)
  return res_sel,first_res_n,pre_res_n

def inverse_transform(r,t):
  r = r.transpose()
  t = - r*t
  return r,t

def get_list_of_chains_selection(selection_str):
  """
  Args:
    selection_str: (str) selection string

  Returns:
    (list of str) of the format ['chain X', 'chain Y',...]
  """
  sstr = selection_str.replace(')',' ')
  sstr = sstr.replace('CHAIN','chain')
  sstr = sstr.replace('Chain','chain') + ' '
  pos_list = [x.end() for x in re.finditer('chain ',sstr)]
  ch_id_list = [sstr[i:sstr.find(' ',i)] for i in pos_list]
  chain_list = ['chain ' + x for x in ch_id_list]
  return chain_list

def separate_selection_string(s):
  s.replace('Chain ','chain ')
  s.replace('CHAIN ','chain ')
  s.replace(' or chain ','chain ')
  s.replace(') or (chain ','chain ')
  s.replace(')or(chain ','chain ')
  if (s[-1] == ')') and (s[0] == '('): s = s[:-1]
  s = s.split('chain ')
  s_list = []
  for sel_str in s:
    sel_str = sel_str.strip()
    if sel_str.endswith(' or'): sel_str = sel_str[:-3]
    if not sel_str in ['','(']:
      s_list.append('chain ' + sel_str)
  return s_list

def get_pdb_header(pdb_str):
  """
  Collect CRYST and SCALE records

  Args:
    pdb_str: (str) pdb type string

  Returns:
    the portion of the pdb_str till the first ATOM
  """
  pdb_str = pdb_str.splitlines()
  pdb_header_lines = []
  for x in pdb_str:
    if x.startswith('ATOM'): break
    else: pdb_header_lines.append(x)
  return '\n'.join(pdb_header_lines)

def get_center_orth(xyz,selection):
  """
  Args:
    xyz:(flex.vec3_double) the complete asu sites cart (atoms coordinates)
    selection:

  Returns:
    (flex.vec3_double) center of coordinates for the selected coordinates
  """
  new_xyz = xyz.select(selection)
  return new_xyz.mean()

def format_num_as_str(n):
  """  return a 3 digit string of n  """
  if n > 999 or n < 0:
    raise IOError('Input out of the range 0 - 999.')
  else:
    s1 = n//100
    s2 = (n%100)//10
    s3 = n%10
    return str(s1) + str(s2) + str(s3)

def ncs_only(transform_info):
  """
  Verify that all transforms are not present
  (excluding the identity transform)

  Args:
    transform_info: (transformation object)

  Returns:
    (bool): True if all transforms are not present
  """
  present = False
  if transform_info:
    ti = transform_info
    for (r,t,n,cp) in zip(ti.r,ti.t,ti.serial_number,ti.coordinates_present):
      if not is_identity(r,t):
        present = present or cp
  return not present

def is_identity(r,t):
  """ test if r, rotation matrix is identity, and t, translation is zero """
  return r.is_r3_identity_matrix() and t.is_col_zero()

def all_ncs_copies_present(transform_info):
  """
  Check if all transforms coordinates are present,
  if the complete ASU is present

  Args:
    transform_info: (transformation object)

  Returns:
    (bool)
  """
  test = True
  for cp in transform_info.coordinates_present:
    test = test and cp
  return test

def order_transforms(transform_to_ncs):
  """ set transformation application order  """
  return sorted(transform_to_ncs,key=lambda key: int(key[1:]))

def uniqueness_test(unique_selection_set,new_item):
  """
  When processing phil parameters. Insert new item to set, if not there,
  raise an error if already in set

  Args:
    unique_selection_set: (set)
    new_item: (str)

  Returns:
    unique_selection_set: updated set
  """
  if new_item in unique_selection_set:
    raise IOError,'Phil selection strings are not unique !!!'
  else:
    unique_selection_set.add(new_item)
    return unique_selection_set

def get_rot_trans(ph=None,
                  master_selection=None,
                  copy_selection=None,
                  max_rmsd=0.02,
                  min_percent = 0.95):
  """
  Get rotation and translation using superpose.
  Can raise "Sorry" if chains are not exactly the same length, but a large
  subset of them is NCS related. Return the subset of atoms indices used in
  such a case.
  Note that if the master and copy have the same length it will not check
  for similarity.

  Args:
    ph : pdb hierarchy input object
    master/copy_selection (str): master and copy selection strings
    max_rmsd (float): limit of rms difference between chains to be considered
      as copies
    min_percent (float): Threshold for similarity between chains

  Returns:
    r: rotation matrix
    t: translation vector
    rmsd (float): RMSD between master and copy
    msg (str): error messages
  """
  msg = ''
  r_zero = matrix.sqr([0]*9)
  t_zero = matrix.col([0,0,0])
  #
  if hasattr(ph,'hierarchy'):
    ph = ph.hierarchy
  if ph:
    cache = ph.atom_selection_cache().selection
    master_ncs_ph = ph.select(cache(master_selection))
    ncs_copy_ph = ph.select(cache(copy_selection))
    chain_id_m,seq_m,res_ids_m  = get_residue_sequence(master_ncs_ph)
    chain_id_c,seq_c,res_ids_c = get_residue_sequence(ncs_copy_ph)
    res_sel_m, res_sel_c, similarity = ncs_search.res_alignment(
      seq_a=seq_m,seq_b=seq_c,
      min_contig_length=0,min_percent=0)
    if similarity == 0 :
      # similarity between chains is small, do not consider as same chains
      return r_zero,t_zero,0,''
    #
    sel_m, sel_c,res_sel_m,res_sel_c,msg1 = get_matching_atoms(
      ph,res_sel_m,res_sel_c,res_ids_m,res_ids_c,
      master_selection,copy_selection)
    msg += msg1
    sel_m = sel_m.iselection(True)
    sel_c = sel_c.iselection(True)
    # master
    other_sites = ph.select(sel_m).atoms().extract_xyz()
    # copy
    ref_sites = ph.select(sel_c).atoms().extract_xyz()
    if ref_sites.size() > 0:
      lsq_fit_obj = superpose.least_squares_fit(
          reference_sites = ref_sites,
          other_sites     = other_sites)
      r = lsq_fit_obj.r
      t = lsq_fit_obj.t
      rmsd = ref_sites.rms_difference(lsq_fit_obj.other_sites_best_fit())
      if rmsd > max_rmsd:
        return r_zero,t_zero,0,msg
    else:
      return r_zero,t_zero,0,'No sites to compare.\n'
    #
    if similarity < min_percent:
      msg='Chains {0} and {1} appear to be NCS related but are dis-similar..\n'
      msg = msg.format('master','copy')
    return r,t,round(rmsd,4),msg
  else:
    return r_zero,t_zero,0,msg

def update_selection_ref(selection_ref,new_selection):
  """
  Test for overlapping selection and then updates the selection_ref
  with the new_selection

  Both received and return arguments are flex.bool
  """
  test = (selection_ref & new_selection).count(True) == 0
  assert test,'Overlapping atom selection. Check phil parameters...\n'
  return selection_ref | new_selection

def get_ncs_group_selection(chain_residue_id):
  """
  Args:
    chain_residue_id: [[chain id's],[[[residues range]],[[...]]]

  Returns:
    selection lists, with selection string for each ncs copy in the group
  """
  chains = chain_residue_id[0]
  res_ranges = chain_residue_id[1]
  assert len(chains) == len(res_ranges)
  ncs_selection = []
  for c,rr in zip(chains, res_ranges):
    c = c.strip()
    assert c.find(' ') < 0,'Multiple chains in a single spec ncs group\n'
    ch_selection = 'chain ' + c
    res_range = ['resseq {0}:{1}'.format(s,e) for s,e in rr]
    res_range = '(' + ' or '.join(res_range) + ')'
    ncs_selection.append(ch_selection + ' and ' + res_range)
  return ncs_selection

def get_transform_order(transform_to_ncs):
  """ order transforms mainly for proper chain naming """
  transform_order = order_transforms(transform_to_ncs)
  transform_chain_assignment = []
  for tr_id in transform_order:
    for tr_selection in transform_to_ncs[tr_id]:
      transform_chain_assignment.append(tr_selection)
  return transform_chain_assignment

def update_ncs_group_map(
        ncs_group_map, ncs_group_id, selection_ids, transform_id):
  """
  Update ncs_group_map that maps a group ID to a list:
   [master selection,list of objects containing transforms and copies
   selections]
  """
  if isinstance(selection_ids, str): selection_ids = [selection_ids]
  if ncs_group_map.has_key(ncs_group_id):
    ncs_group_map[ncs_group_id][0].update(set(selection_ids))
    ncs_group_map[ncs_group_id][1].add(transform_id)
  else:
    ncs_group_map[ncs_group_id] = [set(selection_ids),{transform_id}]
  return ncs_group_map

def sort_dict_keys(d):
  """ sort dictionary d by values """
  return sorted(d,key=lambda k:d[k])

def get_matching_atoms(ph,res_num_a,res_num_b,res_ids_a,res_ids_b,
                       master_selection,copy_selection):
  """
  Get selection of matching chains, match residues atoms
  We keep only residues with continuous matching atoms

  Args:
    ph (hierarchy): all chains
    res_num_a/b (list of int): indices of matching residues position
    res_ids_a/b (list if str): residues IDs
    master/copy_selection (str): master and copy selection strings

  Returns:
    sel_a, sel_b (flex.bool): matching atoms selection
    res_num_a/b (list of int): updated res_num_a/b
    msg (str): message regarding matching residues with different atom number
  """
  atom_cache = ph.atom_selection_cache().selection
  l = ph.atoms().size()
  sel_a = flex.bool([False,] * l)
  sel_b = flex.bool([False,] * l)
  #
  res_num_a_updated = []
  res_num_b_updated = []
  residues_with_different_n_atoms = []
  for i in xrange(len(res_num_a)):
    na = res_num_a[i]
    nb = res_num_b[i]
    resid_a = res_ids_a[na]
    resid_b = res_ids_b[nb]
    sa = atom_cache('({}) and resid {}'.format(master_selection,resid_a))
    sb = atom_cache('({}) and resid {}'.format(copy_selection,resid_b))
    sa = sa.iselection(True)
    sb = sb.iselection(True)
    res_a = ph.select(sa)
    res_b = ph.select(sb)
    if sa.size() != sb.size():
      # collect residue IDs of matching residues with different number of atoms
      residues_with_different_n_atoms.append(resid_a)
    # select only atoms that exist in both residues
    atoms_names_a = list(res_a.atoms().extract_name())
    atoms_names_b = list(res_b.atoms().extract_name())
    atoms_a,atoms_b,similarity = ncs_search.res_alignment(
      seq_a=atoms_names_a, seq_b=atoms_names_b,
      min_contig_length=100)
    # get the number of the atom in the chain
    sa = flex.size_t(atoms_a) + min(sa)
    sb = flex.size_t(atoms_b) + min(sb)
    # keep only residues with continuous matching atoms
    if sa.size() != 0:
      res_num_a_updated.append(na)
      res_num_b_updated.append(nb)
    sa = flex.bool(l,sa)
    sb = flex.bool(l,sb)
    sel_a |= sa
    sel_b |=  sb
  if residues_with_different_n_atoms:
    problem_res_nums = {x.strip() for x in residues_with_different_n_atoms}
    msg = "NCS related residues with different number of atoms, selection"
    msg += master_selection + ':' + copy_selection + '\n['
    msg += ','.join(problem_res_nums) + ']\n'
  else:
    msg = ''
  return sel_a,sel_b,res_num_a_updated,res_num_b_updated,msg

def resseq_string(first_res_num,previous_res_num):
  if previous_res_num > first_res_num:
    res_seq = 'resseq {}:{}'.format(first_res_num,previous_res_num)
  else:
    res_seq = 'resseq {}'.format(first_res_num)
  return res_seq

def get_clean_selection_string(ch_sel,res_selection):
  """
  Args:
    ch_sel (str): such as 'chain A'
    res_selection (list of str): such as ['resseq 1:10','resid 27c and name CA']

  Returns:
    s (str): such as 'chain A and (resseq 1:10 or (resid 27c and name CA))'
  """
  if res_selection:
    if len(res_selection) > 1:
      s = ch_sel + ' and (' + ' or '.join(res_selection) + ')'
    else:
      s = ch_sel + ' and ' + res_selection[0]
  else:
    s = ch_sel
  # remove extra spaces
  s = s.replace('  ',' ')
  s = s.replace('  ',' ')
  return s

def insure_identity_is_in_transform_info(transform_info):
  """
  Add the identity matrix in cases where the pdb or mmcif files do not
  contain it

  Args:
    transform_info (transformation object): contain rotation, translation,
      serial number and indication if present

  Return:
    transform_info : Add or reorder the transform_info so that the
      identity matrix has serial number 1
  """
  ti = transform_info
  ti_zip =  zip(ti.r,ti.t,ti.serial_number,ti.coordinates_present)
  identity_sn = []
  t_r = []
  t_t = []
  t_cp = []
  for i,(r,t,sn,cp) in enumerate(ti_zip):
    if is_identity(r=r,t=t):
      identity_sn.append([i,sn])
      if (i == 0) and (sn == 1):
        t_r.append(r)
        t_t.append(t)
        t_cp.append(cp)
    else:
      t_r.append(r)
      t_t.append(t)
      t_cp.append(cp)
  if identity_sn == [[0,1]]: return transform_info
  # identity transform is missing or not in the first location
  # add identity transform as the first transform
  t_r.insert(0,matrix.sqr([1,0,0,0,1,0,0,0,1]))
  t_t.insert(0,matrix.col([0,0,0]))
  t_cp.insert(0,True)
  # re-assign serial numbers
  s = '{0:3d}'
  t_sn = [s.format(i+1) for i in range(len(t_cp))]
  ti.r = t_r
  ti.t = t_t
  ti.serial_number = t_sn
  ti.coordinates_present = t_cp
  return ti

class NCS_copy():

  def __init__(self,copy_iselection, rot, tran):
    """
    used for NCS groups list copies

    Attributes:
      iselection (flex.size_t): NCS copy selection
      r (matrix obj): rotation matrix from master to this copy
      t (matrix obj): translation vector from master to this copy
    """
    self.iselection = copy_iselection
    self.r = rot
    self.t = tran

class NCS_restraint_group(object):

  def __init__(self,master_iselection):
    """
    used for NCS groups list

    Attributes:
      master_iselection (flex.size_t): NCS group master copy selection
      copies (list): list of NCS_copy objects
    """
    self.master_iselection = master_iselection
    self.copies = []

class Transform(object):
  """ Transformation object """
  def __init__(self,
               rotation = None,
               translation = None,
               serial_num = None,
               coordinates_present = None,
               ncs_group_id = None,
               rmsd = 0.0):
    """
    Basic transformation properties

    Argument:
    rotation : Rotation matrix object
    translation: Translation matrix object
    serial_num : (int) Transform serial number
    coordinates_present: equals 1 when coordinates are presents in PDB file
    ncs_group_id : (int) ncs groups in, all the selections that the same
                   transforms are being applied to
    """
    self.r = rotation
    self.t = translation
    self.serial_num = serial_num
    self.coordinates_present = bool(coordinates_present)
    self.ncs_group_id = ncs_group_id
    self.rmsd = rmsd

class selections(object):
  """
  object for grouping multiple selection information for master ncs and copy
  """

  def __init__(self,
               master_chain_id='',
               copy_chain_id='',
               master_sel=flex.size_t([]),
               copy_sel=flex.size_t([]),
               not_in_master=flex.size_t([]),
               not_in_copy=flex.size_t([])):
    """
    Arg:
      master_chain_id (str)
      opy_chain_id (str)
      master_sel (flex.size_t): master selection
      copy_sel (flex.size_t): copy selection
      not_in_master (flex.size_t): atoms in master to exclude
      not_in_copy (flex.size_t): atoms in copy to exclude
    """
    self.master_chain_id  = master_chain_id
    self.copy_chain_id    = copy_chain_id
    self.master_sel       = master_sel
    self.copy_sel         = copy_sel
    self.not_in_master    = not_in_master
    self.not_in_copy      = not_in_copy
