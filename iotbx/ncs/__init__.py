from __future__ import division
import ncs_preprocess
import sys
import iotbx.pdb
from scitbx.array_family import flex
from libtbx import group_args
import scitbx.matrix

ncs_search_options = """\
ncs_search
  .short_caption = NCS search options
  .style = box noauto
{
  enabled = False
    .type = bool
    .help = Use NCS restraints or constraints in refinement (can be \
              determined automatically)
    .short_caption = Use NCS
    .style = noauto bold
  check_atom_order = True
    .type = bool
    .help = '''check atom order in matching residues'''
    .style = noauto
  exclude_misaligned_residues = True
    .type = bool
    .help = '''check and exclude individual residues
        alignment quality'''
    .style = noauto
  allow_different_size_res = True
    .type = bool
    .help = '''keep matching residue with different
        number of atoms'''
    .style = noauto
  process_similar_chains = True
    .type = bool
    .help = '''When True, process chains that are close
       in length without raising errors'''
    .style = noauto
  match_radius = 4.0
    .type = float
    .help = '''max allow distance difference between pairs of matching
        atoms of two residues'''
  similarity_threshold = 0.95
    .type=float
    .help='''Threshold for similarity between matching chains.
      A smaller value cause more chains to be grouped together and can lower
      the number of common residues'''
  min_contig_length = 10
    .type=int
    .help = "segments < min_contig_length rejected"
  minimize_param = *chains transforms
    .type = choice(multi=False)
    .help = '''chains : minimize the number of NCS related chains in each NCS
      groups. Transforms : minimize the number of NCS operations'''
    .style = noauto
  min_percent = 85.
    .help = '''Threshold for similarity between chains, where similarity
    define as: (number of matching res) / (number of res in longer chain)'''
    .type = float
    .short_caption = Min. percent identity
    .expert_level = 0
    .style = bold noauto
  max_rmsd = 2.
    .type = float
    .short_caption = Max. RMSD
    .help = '''limit of rms difference between chains to be considered
       as copies'''
    .expert_level = 0
    .style = bold noauto
}
"""

# This is not used anywhere
ncs_groups = '''\
refinement
{
  ncs
  .short_caption = NCS restraints and constraints groups selections
  {
    restraint_group
    .multiple = True
    .optional = True
    .short_caption = Restraint group
    {
      reference = None
        .type = str
        .short_caption = Master NCS selection string
      selection = None
        .type = str
        .multiple = True
        .short_caption = NCS copy selection string
    }
    constraint_group
    .multiple = True
    .optional = True
    .short_caption = Constraint group
    {
      reference = None
        .type = str
        .short_caption = Master NCS selection string
      selection = None
        .type = str
        .multiple = True
        .short_caption = NCS copy selection string
    }
  }
}
'''

def input(pdb_hierarchy_inp=None,
          pdb_inp=None,
          hierarchy=None,
          transform_info=None,
          rotations = None,
          translations = None,
          ncs_phil_string = None,
          ncs_phil_groups = None,
          file_name=None,
          file_path='',
          spec_file_str='',
          spec_source_info='',
          cif_string = '',
          quiet=True,
          spec_ncs_groups=None,
          pdb_string=None,
          use_minimal_master_ncs=True,
          max_rmsd=2.0,
          write_messages=False,
          log=None,
          process_similar_chains=True,
          min_percent=0.85,
          chain_similarity_limit=0.95,
          min_contig_length=10,
          check_atom_order=True,
          allow_different_size_res=True,
          exclude_misaligned_residues=True,
          max_dist_diff=4.0,
          ignore_chains=None):
    """
    Select method to build ncs_group_object

    order of implementation:
    1) rotations,translations
    2) transform_info
    3) ncs_phil_string
    4) ncs_phil_groups
    5) spec file
    6) mmcif file
    7) iotbx.pdb.hierarchy.input object

    Args:
    -----
      pdb_hierarchy_inp: iotbx.pdb.hierarchy.input
      transform_info: object containing MTRIX or BIOMT transformation info
      rotations: matrix.sqr 3x3 object
      translations: matrix.col 3x1 object
      ncs_phil_string: Phil parameters
        Phil structure
           ncs_group (multiple)
           {
             reference = ''
             selection = ''   (multiple)
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
      use_minimal_master_ncs (bool): use maximal or minimal common chains
        in master ncs groups
      max_rmsd (float): limit of rms difference between chains to be considered
        as copies
      write_messages (bool): When True, write messages to log
        nearly the same length (but not exactly the same) and are NCS related.
        Raise error if NCS relations are not found
      process_similar_chains (bool): When True, process chains that are close
       in length without raising errors
      min_percent (float): Threshold for similarity between chains
        similarity define as:
        (number of matching res) / (number of res in longer chain)
      chain_similarity_limit (float): min similarity between matching chains
      min_contig_length (int): minimum length of matching chain segments
      check_atom_order (bool): check atom order in matching residues.
        When False, matching residues with different number of atoms will be
        excluded from matching set
      allow_different_size_res (bool): keep matching residue with different
        number of atoms
      exclude_misaligned_residues (bool): check and exclude individual residues
        alignment quality
      max_dist_diff (float): max allow distance difference between pairs of matching
        atoms of two residues
      ignore_chains (set of str): set of chain IDs to exclude
    """
    if not log: log = sys.stdout
    ncs_group_obj = ncs_preprocess.ncs_group_object()
    ncs_group_obj.preprocess_ncs_obj(
      pdb_hierarchy_inp=pdb_hierarchy_inp,
      pdb_inp=pdb_inp,
      hierarchy=hierarchy,
      transform_info=transform_info,
      rotations=rotations,
      translations=translations,
      ncs_phil_string=ncs_phil_string,
      ncs_phil_groups=ncs_phil_groups,
      file_name=file_name,
      file_path=file_path,
      spec_file_str=spec_file_str,
      spec_source_info=spec_source_info,
      cif_string=cif_string,
      quiet=quiet,
      spec_ncs_groups=spec_ncs_groups,
      pdb_string=pdb_string,
      use_minimal_master_ncs=use_minimal_master_ncs,
      max_rmsd=max_rmsd,
      write_messages=write_messages,
      log=log,
      process_similar_chains=process_similar_chains,
      min_percent=min_percent,
      min_contig_length=min_contig_length,
      check_atom_order=check_atom_order,
      allow_different_size_res=allow_different_size_res,
      exclude_misaligned_residues=exclude_misaligned_residues,
      max_dist_diff=max_dist_diff,
      chain_similarity_limit=chain_similarity_limit,
      ignore_chains=ignore_chains)
    return ncs_group_obj

def extract_tncs_groups(distance_threshold, file_name=None, pdb_inp=None,
                        show=True):
  assert [file_name, pdb_inp].count(None) == 1
  # Compute mean distance between master and copy, after moving to origin
  def mean_dist_due_to_pure_rotation(sites_cart_1, sites_cart_2):
    s1 = sites_cart_1 - sites_cart_1.mean()
    s2 = sites_cart_2 - sites_cart_2.mean()
    distances = flex.sqrt((s1 - s2).dot())
    return flex.mean(distances)
  #
  if(file_name is not None):
    pdb_inp = iotbx.pdb.input(file_name=file_name)
    ncs_inp = iotbx.ncs.input(pdb_inp = pdb_inp,
      exclude_misaligned_residues=False)
  else:
    ncs_inp = iotbx.ncs.input(pdb_inp = pdb_inp,
      exclude_misaligned_residues=False)
  sites_cart = pdb_inp.atoms().extract_xyz()
  if(show):
    ncs_groups_str = ncs_inp.get_ncs_info_as_spec().\
      format_all_for_phenix_refine(quiet=True, prefix="ncs_group")
    print ncs_groups_str
  ncs_groups = ncs_inp.get_ncs_restraints_group_list()
  # Find ptNCS
  result = []
  for g in ncs_groups:
    m_sel = g.master_iselection
    for c in g.copies:
      c_sel = c.iselection
      dist_m_c = mean_dist_due_to_pure_rotation(
        sites_cart_1 = sites_cart.select(m_sel),
        sites_cart_2 = sites_cart.select(c_sel))
      if(dist_m_c < distance_threshold):
        pair = group_args(
          master_selection = m_sel,
          copy_selection   = c_sel,
          r                = c.r,
          t                = c.t,
          dist_m_c         = dist_m_c)
        result.append(pair)
  # Find tNCS
  cs = pdb_inp.crystal_symmetry()
  o = cs.unit_cell().orthogonalize
  f = cs.unit_cell().fractionalize
  # Loop over NCS groups
  for g in ncs_groups:
    all_selections = []
    all_selections.append(g.master_iselection)
    for c in g.copies:
      all_selections.append(c.iselection)
    # For each NCS copy loop over symmetry operations (except identity)
    for i, sel_i in enumerate(all_selections):
      master_sites_cart = sites_cart.select(sel_i)
      for smx in cs.space_group().smx():
        r = scitbx.matrix.sqr(smx.r().as_double())
        t = scitbx.matrix.col(smx.t().as_double())
        if(r.is_r3_identity_matrix()): # skip identity operator
          assert t.is_col_zero()
          continue
        copy_sites_cart = o(r.elems*f(master_sites_cart)) # no t required!
        # Loop over NCS copies again to see if symmetry related copy matches
        # any of remaining copies
        for j, sel_j in enumerate(all_selections):
          if(i != j):
            sites_cart_j = sites_cart.select(sel_j)
            dist_m_c = mean_dist_due_to_pure_rotation(
              sites_cart_1 = copy_sites_cart,
              sites_cart_2 = sites_cart_j)
            if(dist_m_c < distance_threshold):
              pair = group_args(
                master_selection = sel_i,
                copy_selection   = None, # means symmetry related
                r                = r.transpose(), # since need to apply to copy
                t                = t,
                dist_m_c         = dist_m_c)
              result.append(pair)
  return result
