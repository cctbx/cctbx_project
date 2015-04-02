from __future__ import division
import ncs_preprocess
import sys

simple_ncs_phil_params = """\
simple_ncs_from_pdb
  .short_caption = Simple NCS from PDB file
{

 pdb_in = None
   .type = path
   .help = 'Input PDB file to be used to identify ncs'
   .short_caption = PDB file
   .style = bold noauto OnChange:load_ncs_pdb_file file_type:pdb no_map
 temp_dir = ""
   .type = path
   .help = "temporary directory (ncs_domain_pdb will be written there)"
   .style = noauto
 min_percent = 85.
   .help = '''Threshold for similarity between chains, where similarity
   define as: (number of matching res) / (number of res in longer chain)'''
   .type = float
   .short_caption = Min. percent identity
   .expert_level = 0
   .style = bold noauto
 max_rmsd = 2.
   .help = '''limit of rms difference between chains to be considered
      as copies'''
   .type = float
    .short_caption = Max. RMSD
   .expert_level = 0
   .style = bold noauto
 max_rmsd_user = None
   .help = '''When given, overrides the max_rmsd'''
   .type = float
   .short_caption = User-specified Max. RMSD
   .expert_level = 0
   .style = bold noauto
 ncs_domain_pdb_stem  = None
   .type = str
   .help = '''NCS domains will be written to ncs_domain_pdb_stem+"group_"+nn'''
   .style = noauto
 write_ncs_domain_pdb = False
   .type = bool
   .help = '''You can write out PDB files representing NCS domains for
        density modification if you want'''
   .style = noauto
 domain_finding_parameters
    .style = box auto_align noauto
  {
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
   .help = "segments < min_contig_length rejected"
   .type=int
 }
 verbose = False
   .type = bool
   .help = "Verbose output"
   .short_caption = Debugging output
   .style = noauto
 show_ncs_phil = False
   .type = bool
   .help = "Show phil parameters for the NCS groups"
   .short_caption = Show phil parameters
   .style = noauto
 show_summary = False
   .type = bool
   .help = "Show NCS information summary"
   .short_caption = Show NCS information summary
   .style = noauto
 raise_sorry = False
   .type = bool
   .help = "Raise sorry if problems"
   .short_caption = raise sorry
 debug = False
   .type = bool
   .help = "Debugging output"
   .short_caption = Debugging output
   .style = noauto
 dry_run = False
   .type = bool
   .help = '''Just read in and check parameter names'''
   .style = noauto
 ncs_file_format = *restraints constraints
   .type = choice(multi=False)
   .help = '''control the .ncs file format
                restraints:  "refinement.ncs.restraint_group"
                constraints: "refinement.ncs.constraint_group"'''
   .style = noauto
 minimize_param = *chains transforms
   .type = choice(multi=False)
   .help = '''chains : minimize the number of NCS related chains in each NCS
   groups.
   transforms : minimize the number of NCS operations'''
   .style = noauto
 exclude_misaligned_residues = True
   .type = bool
   .help = "check and exclude individual residues alignment quality"
   .style = noauto
 check_atom_order = False
   .type = bool
   .help = "check atom order in matching residues"
   .style = noauto
    }
"""

ncs_search_options = """\
  check_atom_order = False
   .type = bool
   .help = '''check atom order in matching residues'''
   .style = noauto
  use_minimal_master_ncs = True
   .type = bool
   .help = '''Minimize number of chains in master ncs groups'''
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
  quiet = True
   .type = bool
   .help = '''quiet output when processing files'''
   .style = noauto
"""

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
          check_atom_order=False,
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

