"""
  simple_ncs_from_pdb.py
  tct 2006-12-12

**************************************
To use as a method, specify input PDB and any other commands as text arguments:

ncs_search=simple_ncs_from_pdb(args=['input_pdb.pdb','max_rmsd=4.'])

Alternatively, you can specify inputs with a params phil object:
ncs_search=simple_ncs_from_pdb(params=params)

Now ncs_search.ncs_object is an object from "sources/mmtbx/mmtbx/ncs.py"
You can get a summary or text for resolve or phenix.refine with:

ncs_search.ncs_object.display_all()
text=ncs_search.ncs_object.format_all_for_resolve()
text=ncs_search.ncs_object.format_all_for_phenix_refine()
**************************************

Purpose: Figure out the ncs from the chains in this pdb object
and return an "ncs" object from "ncs.py" that represents it. This ncs object
can write out a little file for resolve and a file for phenix.refine
specifying this ncs.

Major assumption: residue numbers are consistent among chains

Approach: Use residue numbers to align the residue names, identify
pairs of chains that can match. Choose groupings of chains that maximize the
smallest number of matching residues between each member of a group and the
first (reference) member of the group. Within a pair of chains, allow some
segments to match and others not. Each pair of segments must have a
length >= min_length and percent identity >=min_percent.  A pair of segments
may not end in a mismatch. An overall pair of chains must have an rmsd
of CA atoms of <= rmsd_max.

If find_invariant_domain is specified then once all chains that can be matched
with the above algorithm are identified, all remaining chains are matched,
allowing the break-up of chains into invariant domains. The invariant
domains each get a separate NCS group.

If residue numbers are not the same for corresponding chains, but
they are simply offset by a constant for each chain, this will be
recognized and the chains will be aligned.

"""
import iotbx
import iotbx.phil
from iotbx import pdb
from iotbx.pdb import resseq_decode
import libtbx.phil
import libtbx.phil.command_line
from libtbx.utils import Sorry,null_out
from libtbx import adopt_init_args
from mmtbx.ncs.ncs import ncs
import sys, os, re, string, time
from cctbx.array_family import flex
from scitbx.math import superpose
from scitbx.math import matrix
from mmtbx.ncs.ncs import ncs_group
import mmtbx.monomer_library.pdb_interpretation
from mmtbx.invariant_domain import find_domain
import math
from phenix.utilities.arg_display_methods import arg_display_methods
from phenix.utilities.composite_params import get_composite_master_params
from phenix.utilities.list_methods import list_methods
from phenix.utilities.headers import headers
from phenix.utilities.is_debug import is_debug
from iotbx.pdb import resseq_encode, resseq_decode
from phenix.utilities.catenate_equals import catenate_equals

ncs_master_params = iotbx.phil.parse("""
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
 min_length= 10
   .help = "minimum number of matching residues in a segment"
   .type = int
    .short_caption = Min. number of matching residues per segment
   .expert_level = 2
 min_fraction_represented = 0.10
   .help = "Minimum fraction of residues represented by NCS to keep."
           "If less...skip ncs entirely"
   .type = float
    .short_caption = Min. fraction represented
   .expert_level = 2
 njump   = 1
   .help = "Take every njumpth residue instead of each 1"
   .type = int
   .short_caption = Number of residues per jump
    .expert_level = 2
 njump_recursion   = 10
   .help = "Take every njump_recursion residue instead of each 1 on recursive call"
   .type = int
   .expert_level = 2
 min_length_recursion = 50
   .help = "minimum number of matching residues in a segment for recursive call"
   .type = int
   .short_caption = Min. length recursion
   .expert_level = 2
 min_percent= 95.
   .help = "min percent identity of matching residues"
   .type = float
   .short_caption = Min. percent identity
   .expert_level = 0
   .style = bold noauto
 max_rmsd = 2.
   .help = "max rmsd of 2 chains. If 0, then only search for domains"
   .type = float
    .short_caption = Max. RMSD
   .expert_level = 0
   .style = bold noauto
 quick = True
   .type = bool
   .help = "If quick is set and all chains match, just look for 1 NCS group"
   .short_caption = Quick search
   .style = noauto
 max_rmsd_user = 3.
   .help=max rmsd of chains suggested by user (i.e., if called from phenix.refine \
         with suggested ncs groups)

   .type = float
   .short_caption = Max. RMSD for user-specified chains
   .expert_level = 2
 maximize_size_of_groups = True
   .type = bool
   .help = '''You can request that the scoring be set up to maximize
       the number of members in NCS groups (maximize_size_of_groups=True)
       or that scoring is set up to maximize the length of the matching
       segments in the NCS group (maximize_size_of_groups=False)'''
 require_equal_start_match = True
   .type = bool
   .help = '''You can require that all matching segments start at the same
           relative residue number for all members of an NCS group,
           trimming the matching region as necessary. This
           is required if residue numbers in different chains are not the
           same, but not otherwise'''
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
   find_invariant_domains = True
   .type = bool
   .help = "Find the parts of a set of chains that follow NCS"
   initial_rms = 0.5
   .type=float
   .help="Guess of RMS among chains"
   match_radius = 2.0
   .type = float
   .help = "Keep atoms that are within match_radius of NCS-related atoms"
   similarity_threshold = 0.75
   .type=float
   .help="Threshold for similarity between segments"
   smooth_length = 0
   .help = "two segments separated by smooth_length or less get connected"
   .type=int
   min_contig_length = 3
   .help = "segments < min_contig_length rejected"
   .type=int
   min_fraction_domain = 0.2
     .help = "domain must be this fraction of a chain"
     .type = float
   max_rmsd_domain = 2.
     .help = "max rmsd of domains"
     .type = float
 }
 verbose = False
   .type = bool
   .help = "Verbose output"
   .short_caption = Debugging output
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
    }
""")
master_params = ncs_master_params

restraint_group_params = iotbx.phil.parse("""
 restraint_group
  .multiple=True
  .optional=True
  .short_caption=Restraint group
  .style = noauto auto_align
{
  reference=None
    .type=atom_selection
    .optional=True
    .short_caption=Reference selection
    .input_size=400
    .style = bold
  selection=None
    .type=atom_selection
    .multiple=True
    .short_caption=Restrained selection
    .input_size=400
    .style = bold
  coordinate_sigma=0.05
    .type=float
  b_factor_weight=10
    .type=float
}
""")

class simple_ncs_from_pdb(arg_display_methods,list_methods,headers):

  def __init__(self, args   = None,
                     params = None,
                     ignore_chains   = [],
                     required_chains = [],
                     exclude_chains = [],
                     ncs_master_params = ncs_master_params,
                     command_name   = "simple_ncs_from_pdb",
                     all_chain_proxies = None,
                     pdb_inp = None,
                     hierarchy = None,
                     suppress_print = False,
                     source_info    = None,
                     pdb_file       = None,
                     groups_only = None,
                     suggested_ncs_groups = None,
                     log = sys.stdout,
                     quiet=False,
                     exclude_h=False,
                     exclude_d=False,
                     write_ncs_domain_pdb=None,
                     ncs_domain_pdb_stem=None,
                     temp_dir=None
               ):
    self.log=log
    self.quiet=quiet
    self.exclude_h=exclude_h
    self.exclude_d=exclude_d
    self.exclude_chains=exclude_chains
    self.required_chains=required_chains
    args=catenate_equals(args).new_args()
    self.process_inputs(args)
    if self.args==[]:
       self.args=None
    args=self.args
    if suggested_ncs_groups is None: # take it from args if not directly given
      suggested_ncs_groups=self.suggested_ncs_groups
    else:
      self.suggested_ncs_groups=suggested_ncs_groups
    allow_recursion=self.allow_recursion
    exact_match_only=self.exact_match_only

    master_params_name='simple_ncs_from_pdb'
    self.Name='simple_ncs_from_pdb'

    from phenix.utilities import citations
    if not suppress_print:
      citations.add_citation('phenix','simple_ncs_from_pdb')

    # run a test case
    if args is not None and 'exercise' in args:
      self.exercise()
      return

    args=self.special_cases(args)

    master_params=get_composite_master_params(
         method_list=['simple_ncs_from_pdb'],
         location_list=['phenix.command_line'])

    args=self.get_keyword_table(args,out=self.log)       # set self.keyword_table

    summary,header=self.get_summary_and_header(command_name)

    done,master_params,new_params,changed_params,help=self.get_params(
        command_name,master_params,args,out=self.log)
    if params is None :
      params = new_params
    if done: return
    if not quiet: print >>self.log, header

    if help or (params and params.simple_ncs_from_pdb.verbose):
      print >>self.log, "Values of all params:"
      master_params.format(python_object=params).show(out=log)

    if help or params is None: return

    # Done with standard processing of inputs
    # overwrite with direct inputs, if any:
    if write_ncs_domain_pdb is not None:
      params.simple_ncs_from_pdb.write_ncs_domain_pdb=write_ncs_domain_pdb
    if ncs_domain_pdb_stem is not None:
      params.simple_ncs_from_pdb.ncs_domain_pdb_stem=ncs_domain_pdb_stem
    if temp_dir is not None:
      params.simple_ncs_from_pdb.temp_dir=temp_dir

    # Things that must be defined...

    self.params=params
    if not suppress_print:
      print >>self.log,"Parameters used for simple_ncs_from_pdb:"
      master_params.format(python_object=params).show(out=self.log)
      print >>self.log

    if params.simple_ncs_from_pdb.dry_run:
      print "ARGS: ",args
      return

    # read in the PDB file if needed
    if((all_chain_proxies is None) and (pdb_inp is None and hierarchy is None)):
      if(pdb_file is None):
        if args is not None and args and args[0] and os.path.isfile(args[0]):
          pdb_file=args[0]
        elif params.simple_ncs_from_pdb.pdb_in is not None:
          pdb_file=params.simple_ncs_from_pdb.pdb_in
        else:
          raise Sorry("\nNeed PDB file for simple_ncs_from_pdb"+
             "\n\nPlease make the PDB file the first argument like this: \n"+
             "phenix.simple_ncs_from_pdb mypdb.pdb ...\n")
      if not os.path.isfile(pdb_file):
         raise Sorry("The file "+str(pdb_file)+" is missing?")

      raw_records = flex.std_string()
      raw_records.extend(flex.split_lines(open(pdb_file).read()))
      if pdb_inp is None:
        pdb_inp = iotbx.pdb.input(source_info=None, lines=raw_records)
      mon_lib_srv = mmtbx.monomer_library.server.server()
      ener_lib = mmtbx.monomer_library.server.ener_lib()
      processed_pdb= mmtbx.monomer_library.pdb_interpretation.process(
         mon_lib_srv=mon_lib_srv,
         ener_lib=ener_lib,
         params=None,
         raw_records=raw_records,
         strict_conflict_handling=False,
         max_atoms=None,
         log=null_out())

      all_chain_proxies=processed_pdb.all_chain_proxies
      self.source_info=pdb_file
    if hierarchy is None:
      hierarchy = pdb_inp.construct_hierarchy()

    # set input params
    self.verbose=params.simple_ncs_from_pdb.verbose
    if self.verbose: hierarchy.show(out=self.log)
    if source_info:
      self.source_info=source_info
    if not hasattr(self,'source_info') or not self.source_info:
      self.source_info="None"
    self.find_invariant_domains=\
       params.simple_ncs_from_pdb.domain_finding_parameters.find_invariant_domains
    self.min_fraction_domain= \
       params.simple_ncs_from_pdb.domain_finding_parameters.min_fraction_domain
    self.max_rmsd_domain=params.simple_ncs_from_pdb.domain_finding_parameters.max_rmsd_domain
    self.min_percent=params.simple_ncs_from_pdb.min_percent
    self.quick=params.simple_ncs_from_pdb.quick
    self.maximize_size_of_groups=params.simple_ncs_from_pdb.maximize_size_of_groups
    self.require_equal_start_match=params.simple_ncs_from_pdb.require_equal_start_match
    self.write_ncs_domain_pdb=params.simple_ncs_from_pdb.write_ncs_domain_pdb
    self.ncs_domain_pdb_stem=params.simple_ncs_from_pdb.ncs_domain_pdb_stem
    self.all_chain_proxies=all_chain_proxies
    self.hierarchy=hierarchy
    self.pdb_inp=pdb_inp
    self.njump=params.simple_ncs_from_pdb.njump
    if self.njump<1:
      raise Sorry("njump must be >=1")
    self.min_length=params.simple_ncs_from_pdb.min_length
    self.min_fraction_represented=params.simple_ncs_from_pdb.min_fraction_represented


    # identify chains in the PDB file
    chains,chain_ids,starting_residue_numbers,offset_dict=\
        self.get_chain_list(hierarchy)
    self.offset_dict=offset_dict  # so we do not have to pass this around. It
                                  # does not change ever.

    if not suppress_print:
      print >>self.log,"Chains in this PDB file: ",chain_ids
    if self.verbose:
     for chain,chain_id,start in zip(chains,chain_ids,starting_residue_numbers):
       print >>self.log,"CHAIN: ",chain," ID: ",chain_id," START" ,start
     print >>self.log,"OFFSET LIST: "
     for id in self.offset_dict.keys():
       print >>self.log,id,self.offset_dict[id]


    # set up temp_dir if needed
    if params.simple_ncs_from_pdb.temp_dir:
      if os.path.isfile(params.simple_ncs_from_pdb.temp_dir):
        raise Sorry(
         "The directory "+str(params.simple_ncs_from_pdb.temp_dir)+" cannot be created...")
      if not os.path.isdir(params.simple_ncs_from_pdb.temp_dir):
        os.mkdir(params.simple_ncs_from_pdb.temp_dir)
      if not self.quiet:
        print >>self.log,"Working in ",params.simple_ncs_from_pdb.temp_dir

    groups=[]  # [ ['A','B'],['C','D','E]]
    #              group
    list_of_residue_range_list=[] # [ [      [       [1,120],[130-250] ]]]
    #                                 group  member  ranges
    [all_chains,all_chain_ids,all_starting_residue_numbers]=[chains,chain_ids,
      starting_residue_numbers]


    #========Get suggested NCS groups from phil object ======================

    #  If called with suggested_ncs_groups phil object, we pull out all those
    #  chains here...
    # initialize suggested NCS groups if any
    self.suggested_ncs_groups,self.suggested_group_list=\
          self.get_suggested_groups(suggested_ncs_groups,chain_ids)

    # NOTE: residue_range_list is in original residue numbers, with offsets
    used_ids=[]
    if self.suggested_group_list:  # see if we want to pull out these chains:
      print >> self.log, "Getting NCS from suggested chains: "
      self.suggested_ncs_groups=[] # not both
      for [group,residue_range_list] in self.suggested_group_list:
        rmsd_list,r_list,trans_list,center_list,residues_in_common_list=\
          self.get_rmsd(group,hierarchy,chains,chain_ids,
             starting_residue_numbers,residue_range_list)
        if not rmsd_list or len(rmsd_list)<2:
          pass
        elif rmsd_list[1]>params.simple_ncs_from_pdb.max_rmsd_user:
          print >>self.log,"Warning: requested alignment of ",group,\
               self.add_offsets(residue_range_list,group),\
               " \nrejected due \nto rmsd =",rmsd_list[1]," > ",\
                 params.simple_ncs_from_pdb.max_rmsd_user,\
           ". To keep it, set ncs.max_rmsd_user="+str(int(rmsd_list[1]+0.999))
        else:
          groups.append(group)
          list_of_residue_range_list.append(residue_range_list)
          print >>self.log,"RMSD for suggested group ",group,\
            self.add_offsets(residue_range_list,group)," is ",rmsd_list[1]
      # Now remove all chains in kept groups from list of chains so we only
      # look elsewhere
      used_ids=self.add_ids(groups,used_ids)
      [chains,chain_ids,starting_residue_numbers]=self.remove_used_chains(
                          chains,chain_ids,starting_residue_numbers,used_ids)
    #========End of suggested NCS groups from phil object ==============

    #======= Try to get possible NCS groups from direct comparison of ====
    #         all pairs of chains. Use a high njump to go quickly...
    #         This will work if all chains in a group are identical

    if allow_recursion:
      njump_use=params.simple_ncs_from_pdb.njump_recursion
      min_length_use=params.simple_ncs_from_pdb.min_length_recursion
      if args is not None:
       args_use=args
      else:
       args_use=[]
      args_use+=['njump='+str(njump_use),'min_length='+str(min_length_use)]
      if self.suggested_ncs_groups:   # present if NCS groups defined as "ACDE"
        args_use.append("suggested_ncs_groups="+str(suggested_ncs_groups))
      if self.verbose:
        logfile=log
      else:
        logfile=null_out()
      if self.exact_match_only:
         args_use.append('exact_match')
      args_use.append('no_recursion')

      quick_find_groups=simple_ncs_from_pdb (
                     args   = args_use,
                       # at end of arg list to overwrite..
                     pdb_inp = pdb_inp,  # 091608
                     params = params,
                     ignore_chains   = used_ids,
                     all_chain_proxies = all_chain_proxies,
                     hierarchy = hierarchy ,
                     suppress_print = True,
                     source_info    = source_info,
                     log=logfile,
                     quiet=True,
                     groups_only = True)

      self.suggested_ncs_groups=quick_find_groups.sequence_groups
      if not suppress_print:
        print >>self.log,"GROUPS BASED ON QUICK COMPARISON:",\
          quick_find_groups.sequence_groups

    # ======= End of getting NCS groups from direct comparison =========

    # ======= Getting groups quickly if this is called by simple_ncs ======
    #         This does the work for the recursive call above
    if groups_only:
      self.sequence_groups,self.sequence_list_of_residue_range_list=\
         self.find_groups(hierarchy,chains,chain_ids,
           starting_residue_numbers,
           min_length=params.simple_ncs_from_pdb.min_length,
           min_percent=params.simple_ncs_from_pdb.min_percent,
           max_rmsd=params.simple_ncs_from_pdb.max_rmsd,
           max_rmsd_user=params.simple_ncs_from_pdb.max_rmsd_user,
           called_by_self=True,
           exact_match_only=exact_match_only)
      return
    # ======= End of getting groups quickly if this is called by simple_ncs ======


    # ====== Get NCS groups using only sequence information and the
    #        suggested_ncs_groups found above, if any
    sequence_groups,sequence_list_of_residue_range_list=\
         self.find_groups(hierarchy,chains,chain_ids,
           starting_residue_numbers,
           min_length=params.simple_ncs_from_pdb.min_length,
           min_percent=params.simple_ncs_from_pdb.min_percent,
           max_rmsd=999999.,max_rmsd_user=15,called_by_self=True,
           exact_match_only=exact_match_only)
    if self.verbose:
      print >>self.log,"SEQUENCE-BASED GROUPS: ",sequence_groups
      print >>self.log,sequence_list_of_residue_range_list
    # ====== End of NCS groups using sequence information and suggested groups=======

    # JUST HERE problem sequence_list_of_residue_range_list has offsets in
    # some cases where it should not  021107...

    # =======Get new groups with domains on all sequence groups=======
    invariant_groups,invariant_list_of_residue_range_list=\
            self.find_invariant_groups(hierarchy,
             sequence_groups,sequence_list_of_residue_range_list,
             chains,chain_ids,
             starting_residue_numbers)
    if invariant_groups:
        if self.verbose:
          print >>self.log,"Invariant groups found:",\
              invariant_groups,invariant_list_of_residue_range_list
        groups+=invariant_groups
        list_of_residue_range_list+=invariant_list_of_residue_range_list

    # =======End of new groups with domains =============

    [chains,chain_ids,starting_residue_numbers]=[all_chains,  # restore these
        all_chain_ids,all_starting_residue_numbers]

    #  ======== Write out results and make an ncs object with them in it===
    ncs_object=ncs(exclude_h=self.exclude_h,exclude_d=self.exclude_d)
    count=0
    for group,residue_range_list in zip(groups,list_of_residue_range_list):
      count+=1
      # if necessary, add offsets from self.offset_dict to the values in
      #  residue_range_list
      residue_range_list_with_offsets=self.add_offsets(residue_range_list,group)
      if self.verbose:
        print >>self.log,"\nNCS GROUP",count,":",group  ,\
            residue_range_list_with_offsets
      # group is a list of chain_ids, with the reference one first
      # so get rmsd for members of the group from reference
      rmsd_list,r_list,trans_list,center_list,residues_in_common_list=\
        self.get_rmsd(group,hierarchy,chains,chain_ids,starting_residue_numbers,
        residue_range_list)  # NO OFFSET (I know, it's confusing)!

      if self.write_ncs_domain_pdb:
        ncs_domain_pdb=self.make_ncs_domain_pdb(
         stem=self.ncs_domain_pdb_stem,
         hierarchy=hierarchy,group_number=count,
         group=group,residue_range_list=residue_range_list_with_offsets,
         params=params)
      else:
        ncs_domain_pdb=None
      if not rmsd_list:
         print >>self.log,"\nNCS GROUP",count,":",group  ,\
               residue_range_list_with_offsets
         print >>self.log,"No rmsd found...giving up on this group"
      else:
        chain_residue_id=[group,residue_range_list_with_offsets]
        ncs_object.import_ncs_group(ncs_rota_matr=r_list,
         rmsd_list=rmsd_list,
         residues_in_common_list=residues_in_common_list,
         center_orth=center_list,
         trans_orth=trans_list,
         chain_residue_id=chain_residue_id,
         ncs_domain_pdb=ncs_domain_pdb,
         source_of_ncs_info=self.source_info)

    if len(ncs_object.ncs_groups()) >=1 and \
        self.too_few_residues_represented(ncs_object=ncs_object,
      total_residues=self.total_residues): # skip entirely
      print >>self.log,"Skipping NCS. Too few residues represented (< %6.1f percent of total)" %(100.*self.min_fraction_represented)
      ncs_object=ncs(exclude_h=self.exclude_h,exclude_d=self.exclude_d)

    if len(ncs_object.ncs_groups())<1:
      if not suppress_print:
        if self.source_info:
          print >>self.log,"\nNo NCS found from the chains in ",self.source_info
        else:
          print >>self.log,"\nNo NCS found"

    if not suppress_print:
      ncs_object.display_all(log=self.log)
      f=open("simple_ncs_from_pdb.resolve",'w')
      ncs_object.format_all_for_resolve(log=self.log,out=f)
      f.close()
      f=open("simple_ncs_from_pdb.ncs",'w')
      ncs_object.format_all_for_phenix_refine(log=self.log,out=f)
      f.close()
      f=open("simple_ncs_from_pdb.ncs_spec",'w')
      ncs_object.format_all_for_group_specification(log=self.log,out=f)
      f.close()

    self.ncs_object=ncs_object

    from phenix.utilities import citations
    citations.show(out=self.log, source='simple_ncs_from_pdb')

    #  ======== Done with writing out results and making an ncs object ====


  def too_few_residues_represented(self,ncs_object=None,
      total_residues=None):
    # count up how many residues are represented in the NCS and return True
    # if not a fraction self.min_fraction_represented
    total_represented=0
    for ncs_group in ncs_object.ncs_groups():
      for common in ncs_group.residues_in_common_list():
        total_represented+=common

    if total_residues is not None and total_residues>0:
      fraction=float(total_represented)/float(total_residues)
      if fraction < self.min_fraction_represented:
         return True # too few
      else:
         return False # not too few
    return False # do not know


  def make_ncs_domain_pdb(self,stem='',hierarchy=None,group_number=None,
         group=[],residue_range_list=[],params=None):
    # write out a file with the name of the group containing all atoms
    # associated with this group
    #if self.exclude_chains is set then EXCLUDE members of exclude_chains

    if hierarchy is None or group_number is None: return None
    if stem is None: stem=""
    file_name=stem+'group_'+str(group_number)+'.pdb'
    full_file_name=os.path.join(params.simple_ncs_from_pdb.temp_dir,file_name)
    f=open(full_file_name,'w')
    start_with_residue=1
    crystal_symmetry=None
    if crystal_symmetry:
      print >> out, pdb.format_cryst1_record(
         crystal_symmetry=crystal_symmetry)
    for model in hierarchy.models():
      for chain in model.chains():
        for conformer in chain.conformers():
          for residue in conformer.residues():
            resseq_int=residue.resseq_as_int()
            ok=False
            for id,residue_ranges in zip(group,residue_range_list):
                if id in self.exclude_chains: continue
                if chain.id==id:
                  for start,end in residue_ranges:
                    if resseq_int>=start and resseq_int<=end:
                      ok=True
                      break
                if ok: break
            if ok:
              for atom in residue.atoms():
                print >>f, atom.format_atom_record()

    f.close()
    return file_name

  def special_cases(self,args):
    # special cases for input files so user doesn't need to specify:
    new_args=[]
    if not args: return []
    for arg in args:
      if (os.path.isfile(arg)):
        if arg[-3:] in ['pdb','PDB']:
          arg_use='simple_ncs_from_pdb.pdb_in='+arg
        else:
          arg_use=arg
      else:
        arg_use=arg
      new_args.append(arg_use)
    return new_args

  def get_chain_length(self,chains):
    shortest_chain=9999999
    longest_chain=0
    for chain in chains:
      length=len(chain)
      if length<shortest_chain: shortest_chain=length
      if length>longest_chain: longest_chain=length
    return shortest_chain,longest_chain

  def remove_used_chains(self,chains,chain_ids,
            starting_residue_numbers,used_ids):
      [new_chains,new_chain_ids,new_starting_residue_numbers]=[[],[],[]]
      for chain,id,start in zip (chains,chain_ids,starting_residue_numbers):
        if not id in used_ids:
          new_chains.append(chain)
          new_chain_ids.append(id)
          new_starting_residue_numbers.append(start)
      return [new_chains,new_chain_ids,new_starting_residue_numbers]

  def add_ids(self,groups,used_ids):
      for group in groups:
          for id in group:
            if not id in used_ids: used_ids.append(id)
      return used_ids

  def max(self,x,y):
    if x>=y: return x
    return y

  def add_offsets(self,residue_range_list,group):
    # add offsets to all residue numbers in residue_range_list using
    # information in self.offset_dict and group:
    new_residue_range_list=[]
    id1=group[0]
    if self.verbose:print >>self.log,"OFFSETS:"
    for id2,residue_range in zip(group,residue_range_list):
      offset=self.offset_dict[id1+id2]
      if self.verbose: print >>self.log,residue_range,id1,id2,offset
      new_residue_range=[]
      for pair in residue_range:
        new_residue_range.append([pair[0]-offset,pair[1]-offset])
      new_residue_range_list.append(new_residue_range)
    if self.verbose:print >>self.log,"OLD, NEW: ",residue_range_list,new_residue_range_list
    return new_residue_range_list
  def subtract_offsets(self,residue_range_list,group):
    # remove offsets from all residue numbers in residue_range_list using
    # information in self.offset_dict and group:
    new_residue_range_list=[]
    id1=group[0]
    if self.verbose:print >>self.log,"OFFSETS:"
    for id2,residue_range in zip(group,residue_range_list):
      offset=self.offset_dict[id1+id2]
      if self.verbose: print >>self.log,residue_range,id1,id2,offset
      new_residue_range=[]
      for pair in residue_range:
        new_residue_range.append([pair[0]+offset,pair[1]+offset])
      new_residue_range_list.append(new_residue_range)
    if self.verbose:print >>self.log,"OLD, NEW: ",residue_range_list,new_residue_range_list
    return new_residue_range_list


  def ncs_object(self):
    # This is an instance of "ncs" as in "$PHENIX/phenix/phenix/autosol/ncs.py"
    if hasattr(self,'ncs_object'):
      return self.ncs_object

  def get_suggested_groups(self,suggested_ncs_groups,chain_ids):
    # expect a phil object or text...
    groups=[]
    group=[]
    suggested_group_list=[]
    if not suggested_ncs_groups: return groups,suggested_group_list
    if not type(suggested_ncs_groups)==type("ABCD"):
      suggested_group_list=self.get_suggested_group_list(suggested_ncs_groups,
           chain_ids)
      return groups,suggested_group_list
    else:
     input_text=suggested_ncs_groups
     # expect text: "[ABC][DE]"
     # return a list of lists: [["A","B","C"],["D","E"]]
     for char in input_text:
      if char=="[":  # start a group
        if group:
          raise Sorry(
            'Format for suggested_ncs_groups is "[ABCD][EFG]" or "ABCD"')
      elif char=="]":  # end a group
        if not group:
          raise Sorry(
            'Format for suggested_ncs_groups is "[ABCD][EFG]" or "ABCD"')
        group,groups=self.finish_suggested_group(group,groups)
      elif char in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
        if not char in chain_ids:
          raise Sorry ("The chain ID '"+str(char)+
         "' is not among the chain ids in this PDB file. \nAvailable chains: "+
         str(chain_ids))
        group.append(char)
     if group:
        group,groups=self.finish_suggested_group(group,groups)
     # check that each chain is used no more than once
     list_of_used_chains=[]
     for group in groups:
       for id in group:
         if id in list_of_used_chains:
            raise Sorry(
            "Each chain can be used only once in suggested_ncs_groups"+\
            '\nChain '+str(id)+' is listed more than once in '+str(input_text))
         list_of_used_chains.append(id)
     return groups,suggested_group_list

  def finish_suggested_group(self,group,groups):
        if len(group)<2:
          raise Sorry('The group '+str(group)+
            ' seems to have fewer than 2 members'+
            '\nFormat for suggested_ncs_groups is "[ABCD][EFG]" or "ABCD"')
        groups.append(group)
        group=[]
        return group,groups

  def get_suggested_group_list(self,suggested_ncs_groups,chain_ids):
      suggested_group_list=[]
      # 121706: make sure that all suggested ids are actually in chain_ids!
      # interpret suggested ncs groups and make a list of them in our format
      for i_seq, group in enumerate(suggested_ncs_groups):
        if(group.reference is not None):
           selections=[group.reference]+group.selection
           group=[]
           range_list=[]
           first=True
           for sel in selections:
             selection_object=self.atom_selection(sel)
             if selection_object is None:
                print >>self.log,"No matching residues in suggested NCS groups"
                return None
             id,ranges=self.get_id_and_ranges(selection_object,chain_ids)
             if id is None or ranges is None:
                print >>self.log,"No matching residues in suggested NCS groups"
                return None
             if not id in chain_ids:
                raise Sorry("\nThe suggested chain "+str(id)+" was not found "+
                  " in the input file...it could \nbe the chain is too "+
                  "short (try decreasing min_length) \nor that it is not "+
                  "recognized as a chain. Note that ligands\n"+
                  " and solvent are not recognized as chains for NCS\n")
             group.append(id)
             range_list.append(ranges)
           if group and range_list:
             range_list_without_offsets=self.subtract_offsets(range_list,group)
             suggested_group_list+=[[group,range_list_without_offsets]]
             # allow list of lists

      return suggested_group_list

  def get_id_and_ranges(self,selection_object,chain_ids):
    # we want to return the chain ID and list of residue ranges referred
    #  to by this selection object. If there are >1 ID: raise an exception as
    #   we cannot be sure that it is handled correctly.
    #  assume that all residue numbers go up sequentially!
    # In this routine there are no offsets.
    id_list=[]
    residue_range_list=[]
    hierarchy = self.hierarchy
    selection = selection_object
    assert selection.size() == self.pdb_inp.atoms().size()
    for atom,is_selected in zip(self.pdb_inp.atoms(), selection):
      atom.tmp = int(is_selected)
    for model in hierarchy.models():
      first_in_model = True
      for chain in model.chains():
        first_in_chain = True
        conformer_number=0
        for conformer in chain.conformers():
          conformer_number+=1
          if conformer_number>1: continue # take only first conformer
          first_in_conformer = True
          last_residue=None
          for residue in conformer.residues():
            if (self.is_selected(residue)):
              if (first_in_model):
                first_in_model = False
              if (first_in_chain):
                first_in_chain = False
              if (first_in_conformer):
                first_in_conformer = False

              if not chain.id in id_list: id_list.append(chain.id)
              resno=residue.resseq_as_int()
              if last_residue is None or resno != last_residue+1: # new segment
               if last_residue is not None:  # finished with a segment
                 residue_range_list.append([first_residue,last_residue])
               first_residue=resno
               last_residue=resno
              else:  # just increment last_residue
                last_residue=resno
          if last_residue is not None:
            residue_range_list.append([first_residue,last_residue])

    if len(id_list)==0 or not residue_range_list:
       return None,None

    if len(id_list)>1:
      raise Sorry("This version of simple_ncs_from_pdb requires any input NCS "+
       "\nspecifications to have only one chain in each reference or selection"+
       "\n('chain A or chain B' is not allowed)")

    return id_list[0],residue_range_list

  def is_selected(self,residue):
    for atom in residue.atoms():
      if (atom.tmp): return True
    return False

  def atom_selection(self, string):
    try:
      selection=self.all_chain_proxies.selection(string = string)
    except KeyboardInterrupt: raise
    except:
      selection=None
    return selection

  def find_chain_and_start(self,
       chain_id1,chains,chain_ids,starting_residue_numbers):
    # get chain and starting residue number for chain with id chain_id
    for chain,chain_id,start in zip(chains,chain_ids,starting_residue_numbers):
      if chain_id==chain_id1:
         return chain,start
    return None,None

  def get_rmsd(self,group,hierarchy,chains,chain_ids,starting_residue_numbers,
     residue_range_list):
    # get rmsd of all from the reference one, using only residues in
    # residue_range_list for each chain
    # residue_range_list is [residue_range1,residue_range2...]
    #  residue_range_1 is list of residue ranges: [[1,3],[7-9]...]

    # NOTE: residue ranges are absolute

    rmsd_list=[]
    r_list=[]
    t_list=[]
    center_list=[]
    residues_in_common_list=[]
    if self.verbose:
      print >>self.log,"GET_RMSD: residue_range_list: ",residue_range_list
    if len(group)<2:
      return rmsd_list,r_list,t_list,center_list,residues_in_common_list

    chain_id1=group[0]
    residue_range1=residue_range_list[0]
    chain1,start1=self.find_chain_and_start(
       chain_id1,chains,chain_ids,starting_residue_numbers)
    if not chain1:
      return rmsd_list,r_list,t_list,center_list,residues_in_common_list

    for chain_id2,residue_range2 in zip(group,residue_range_list):
     chain2,start2=self.find_chain_and_start(
       chain_id2,chains,chain_ids,starting_residue_numbers)
     if not chain2:
      return rmsd_list,r_list,t_list,center_list,residues_in_common_list
     if self.verbose: print >>self.log,"getting CA values for ",chain_id1," and ",chain_id2
     # pull out CA coords from chain1 and chain2 for all residues where both are
     #   present (as reflected in chains[])
     # residue range1 and 2 are absolute
     residue_list=self.get_residues_in_common(
          chain1,start1,chain2,start2,residue_range1,
          residue_range2,chain_id1,chain_id2)
     if self.verbose: print >>self.log,"Residues in common: ",residue_list
     offset1=0
     offset2=self.offset_dict[chain_id1+chain_id2]
     center1,xyz_list1=\
        self.get_xyz_list(chain_id1,hierarchy,residue_list,offset1)
     center2,xyz_list2=\
        self.get_xyz_list(chain_id2,hierarchy,residue_list,offset2)
     if self.verbose:print >>self.log,"Center: ",center1,center2
     if not xyz_list1 or len(xyz_list1)==0 or \
        not xyz_list2 or len(xyz_list2)==0 or \
        not len(xyz_list1)==len(xyz_list2):
       if self.verbose: print >>self.log,"Note: No overlap of segments in getting NCS rmsd"
       return [],[],[],[],[],

     else:  # all ok...
       lsq_fit = superpose.least_squares_fit(
        reference_sites=xyz_list1,
        other_sites=xyz_list2)
       rmsd = xyz_list1.rms_difference(lsq_fit.other_sites_best_fit())
       rmsd_list.append(rmsd)
       residues_in_common_list.append(len(xyz_list1))
       r_list.append(lsq_fit.r)
       t_list.append(lsq_fit.t)
       center_list.append(center2)

       if self.verbose:
         xyz=flex.vec3_double()
         xyz.append(lsq_fit.r * matrix.col(xyz_list1[0]) + lsq_fit.t)
         print >>self.log,"Getting atom 1 of chain ",chain_id2,\
              " from atom 1 of chain ",chain_id1
         print >>self.log,"XYZ of atom 1 chain ",chain_id1,\
                 " after application of RT: ",xyz[0]
         print >>self.log,"xyz of atom 1 chain ",chain_id2," : ",xyz_list2[0]

    return rmsd_list,r_list,t_list,center_list,residues_in_common_list

  def get_residues_in_common(self,chain1,start1,chain2,start2,residue_range1,
        residue_range2,id1,id2):
    # return list of residue numbers shared by chain1 and chain2
    # in range defined by intersection of residue_range1 and residue_range2
    #  residue_range_1 is list of residue ranges: [[1,3],[7-9]...]

    # NOTE: residue ranges are absolute

    id=id1+id2
    if id in self.offset_dict.keys():
      offset=self.offset_dict[id]
    else:
      print >>self.log,"NO ID : ",id
      offset=0
    start2_use=start2+offset

    if self.verbose:
       print >>self.log,"Getting residues in common"
       print >>self.log,"Chain1, start1: ",chain1,start1
       print >>self.log,"Chain2, start2: ",chain2,start2,start2_use
       print >>self.log,"Residue_range1",residue_range1
       print >>self.log,"Residue_range2",residue_range2
    residue_list=[]
    if not residue_range1 or not residue_range2:
      return residue_list

    all_residues1=self.list_all_integers_in_range(residue_range1)
    all_residues2=self.list_all_integers_in_range(residue_range2)
    all_to_use=[]
    for res in all_residues1:
       if res in all_residues2: all_to_use.append(res)
    for res in all_to_use:
      i1=(res-start1)//self.njump
      i2=(res-start2_use)//self.njump
      if i1 >=0 and i1 <len(chain1) and i2 >=0 and i2<len(chain2) and \
         chain1[i1] and chain2[i2] and chain1[i1]==chain2[i2]:  # Matching
          residue_list.append(res)
    return residue_list

  def list_all_integers_in_range(self,residue_range_list):
    # NOTE: residue ranges are absolute
    all_int=[]
    if not residue_range_list or len(residue_range_list)<1:return
    for range in residue_range_list:
       start=range[0]
       end=range[1]
       for i in xrange(start,end+1):
         if self.skip(i):continue
         all_int.append(i)
    return all_int

  def get_xyz_list(self,chain_id,hierarchy,residue_list,offset):
     if self.verbose: print >>self.log,"Getting xyz for ",chain_id,residue_list
     # NOTE: we are going to offset all the residue numbers in this chain
     # by offset...

     residue_found_list=[]
     xyz_list=flex.vec3_double()
     model_number=0
     center=flex.vec3_double()
     center.append( (0.,0.,0.))
     chain_ids_seen_already=[]
     for model in hierarchy.models():
      model_number+=1
      if model_number>1: continue  # take first model only
      for chain in model.chains():
        if chain.id in chain_ids_seen_already: continue
        chain_ids_seen_already.append(chain.id)
        if chain.id != chain_id: continue # we only want chain_id...
        conformer_number=0
        for conformer in chain.conformers():
          conformer_number+=1
          if conformer_number>1:continue  # take first conformer only
          residues_seen_already=[]
          for residue in conformer.residues():
            resseq_int = residue.resseq_as_int()
            if self.skip(resseq_int):
              continue  # based on self.njump
            if resseq_int in residues_seen_already:continue
            residues_seen_already.append(resseq_int)
            if (resseq_int+offset) in residue_list:  # take CA from this
              found_ca_or_p=False
              for atom in residue.atoms():
                 if (atom.name == " CA " or atom.name == " P  "):
                   xyz = atom.xyz
                   found_ca_or_p=True
                   break
              if found_ca_or_p:  # require CA: see get_chain_list
                xyz_list.append(xyz)
                residue_found_list.append(resseq_int+offset)
                center+=xyz
              else:
                if self.verbose:
                  print >>self.log,\
                  "NOTE: Missing CA/P for residue "+residue.resseq+ \
                  " of chain "+str(chain_id)
     if self.njump==1:   #  check only valid if njump==1
      for res in residue_list:
       if not res in residue_found_list:
          print >>self.log,"Residues expected: ",residue_list
          print >>self.log,"Residues found   : ",residue_found_list
          raise Sorry("Did not find expected residues in get_xyz_list"+
           "\nIt is possible that your PDB file has segments with one "+
           "chainID that are \ninterspersed with segments with another chainID")
     if len(xyz_list):
       center=center* (1./float(len(xyz_list)))
       return center[0],xyz_list
     else:
       return None,None


  def find_groups(self,hierarchy,chains,chain_ids,starting_residue_numbers,
     min_length=5,min_percent=100.,max_rmsd=5.,
     max_rmsd_user=15,called_by_self=False,exact_match_only=False):
    i1=-1
    match_list=[]
    for chain1,chain_id1,start1 in zip(
           chains,chain_ids,starting_residue_numbers):
      # if quick is set and all chains match to the first one we have enough
      if self.quick and self.have_enough(chain_ids,match_list):
        continue
      i1+=1
      i2=-1
      for chain2,chain_id2,start2 in zip(
           chains,chain_ids,starting_residue_numbers):
        i2+=1
        if i1>i2: continue
        if self.verbose:
           print >>self.log,"Comparing ",chain_id1,chain_id2
        matches,keep_range_list=self.matches_approx(
          chain1,start1,chain2,start2,chain_id1,chain_id2,
          min_length=min_length,min_percent=min_percent,
          exact_match_only=exact_match_only)
        if matches:
          # test for max_rmsd now and simply reject if too big

          group=[chain_id1,chain_id2]
          in_suggested_ncs_groups=self.are_in_suggested_ncs_groups(group)
          # reject if in >1 suggested NCS groups:
          if self.are_in_different_suggested_ncs_groups(group): continue

          residue_range=keep_range_list # list of ranges: [[1,3],[7-9]]
          residue_range_list=[residue_range,residue_range]
          # Look for Residue range used... it is 10-180
          rmsd_list,r_list,trans_list,center_list,residues_in_common_list=\
          self.get_rmsd(group,hierarchy,chains,chain_ids,
             starting_residue_numbers,residue_range_list)
          if not rmsd_list or len(rmsd_list)<2:
            pass
          elif not in_suggested_ncs_groups and rmsd_list[1]>max_rmsd:
            if self.verbose:
              print >>self.log,"Alignment ",group," rejected with rmsd =",rmsd_list[1]
          elif in_suggested_ncs_groups and rmsd_list[1]>max_rmsd_user and \
             rmsd_list[1]>max_rmsd:
            print >>self.log,"NOTE: User-suggested alignment ",group, \
                " rejected with rmsd =",rmsd_list[1]
            print >>self.log,"Set max_rmsd_user="+str(int(rmsd_list[1]+0.9999))+\
               " to keep it"
          else:
            match_list.append([matches,chain_id1,chain_id2,residue_range])
    if self.verbose: print >>self.log,"MATCH LIST: ",match_list

    # now figure out what are in clusters together, and which is the best
    #  reference for each group ...
    remaining_ids=chain_ids
    still_finding_groups=True
    groups=[]
    list_of_residue_range_list=[]

    while still_finding_groups:
       best_id=None
       best_score=-99999
       best_group=None
       best_residue_range_list=None
       for trial_id in remaining_ids: # use this as seed to pull out all others
         score,group,residue_range_list=self.get_score_and_group(trial_id,
             remaining_ids,match_list)
         if score>best_score:
            best_score=score
            best_id=trial_id
            best_residue_range_list=residue_range_list
            best_group=group

       if best_score>0:
         if self.verbose:
            print >>self.log,"Best new group:",best_group,best_score,\
              best_residue_range_list
         groups.append(best_group)
         list_of_residue_range_list.append(best_residue_range_list)
         new_remaining_ids=[]
         for id in remaining_ids:
           if not id in best_group:
             new_remaining_ids.append(id)
         remaining_ids=new_remaining_ids
         if self.verbose:
           print >>self.log,"Current remaining ids",remaining_ids
       else:
         if self.verbose:
            print >>self.log,"No new groups"
         still_finding_groups=False

    if not called_by_self and self.find_invariant_domains and remaining_ids:
       if self.verbose:
          print >>self.log,"Current groups: ",groups,list_of_residue_range_list
       # find groups from remaining_ids...this time taking ANY rmsd:
       remaining_chains=[]
       remaining_chain_ids=[]
       remaining_starting_residue_numbers=[]
       for chain,chain_id,starting_residue_number in zip (
          chains,chain_ids,starting_residue_numbers):
         if chain_id in remaining_ids:
          remaining_chains.append(chain)
          remaining_chain_ids.append(chain_id)
          remaining_starting_residue_numbers.append(starting_residue_number)
       additional_groups,additional_list_of_residue_range_list=\
         self.find_groups(hierarchy,remaining_chains,remaining_chain_ids,
           remaining_starting_residue_numbers,
           min_length=min_length,min_percent=min_percent,
           max_rmsd=999999.,max_rmsd_user=15,called_by_self=True)

       if additional_groups:
         invariant_groups,invariant_list_of_residue_range_list=\
            self.find_invariant_groups(hierarchy,
             additional_groups,additional_list_of_residue_range_list,
             remaining_chains,remaining_chain_ids,
             remaining_starting_residue_numbers)
         if invariant_groups:
           print >>self.log,"Invariant groups found:",\
              invariant_groups,invariant_list_of_residue_range_list
           return groups+invariant_groups, \
             list_of_residue_range_list+invariant_list_of_residue_range_list

    return groups,list_of_residue_range_list

  def have_enough(self,chain_ids,match_list):  # are all ids matched to a set
    # of id #1?
    if len(chain_ids)<2: return True
    for id in chain_ids[1:]:
      found=False
      for [matches,id1,id2,residue_range] in match_list:
        if id1==chain_ids[0] and id2==id: found=True
      if not found:
        if self.verbose: print >> self.log,"missing ",id,match_list
        return False
    if self.verbose: print >> self.log,"have ",chain_ids,match_list
    return True  # all the id's are matched to id #1


  def find_invariant_groups(self,hierarchy,
             additional_groups,additional_list_of_residue_range_list,
             remaining_chains,remaining_chain_ids,
             remaining_starting_residue_numbers):
    # for each group in additional groups, see if there is some way to break
    # it down into invariant sub-groups. Return list of these and new
    #  residue ranges...

    invariant_groups=[]
    invariant_list_of_residue_range_list=[]
    if not additional_groups or not additional_list_of_residue_range_list:
      return invariant_groups,invariant_list_of_residue_range_list

    for group,residue_range_list in zip(
      additional_groups,additional_list_of_residue_range_list):
      if self.verbose:
        print >>self.log,"Looking for invariant domains for ...:",group,\
          self.add_offsets(residue_range_list,group)
      # we are going to require that our domains be the same for all members
      # of this group....Find the domains in pairs with first chain in
      #  list...and choose the set of domains that best fits
      #  ALL comparisons to this chain
      id1 = group[0]
      chain1,start1=self.get_chain(id1,remaining_chains,remaining_chain_ids,
           remaining_starting_residue_numbers)
      residue_range1=residue_range_list[0]
      best_rmsd=None
      best_score=None
      best_matches=[]
      best_domain_list=[]
      best_domains=[]
      best_overall_rmsd_list=[]
      best_number_of_residues_used=[]
      for id2,residue_range2 in zip(group,residue_range_list):
          if id1==id2: continue
          if self.verbose: print >>self.log,"Using ",id1,id2," to find domains..."
          chain2,start2=self.get_chain(id2,remaining_chains,
             remaining_chain_ids,remaining_starting_residue_numbers)
          if not chain1 or not chain2:
             raise Sorry("Chain id "+str(id1)+" "+str(id2)+
               " are not both present in remaining chain_ids: "+
               str(remaining_chain_ids))

          residue_list=self.get_residues_in_common(
            chain1,start1,chain2,start2,residue_range1,residue_range2,
            id1,id2)
          if not residue_list: continue  # 021107
          if self.verbose:
            print >>self.log,"residues in common",residue_list
          # get sites from first and second...
          offset1=0
          offset2=self.offset_dict[id1+id2]
          center1,xyz_list1=\
            self.get_xyz_list(id1,hierarchy,residue_list,offset1)
          center2,xyz_list2=\
            self.get_xyz_list(id2,hierarchy,residue_list,offset2)
          if self.verbose: print >>self.log,"Centers: ",center1,center2
          min_size=int(
             self.min_fraction_domain*float(len(chain1))/float(self.njump)
             )  # 2011-03-30
          if xyz_list1 is None: continue  # 021107
          if min_size>len(xyz_list1):
             min_size=len(xyz_list1)

          if len(xyz_list1)<4: continue

          domains = find_domain(
            xyz_list1,
            xyz_list2,
            initial_rms=
               self.params.simple_ncs_from_pdb.domain_finding_parameters.initial_rms,
               match_radius =
                self.params.simple_ncs_from_pdb.domain_finding_parameters.match_radius,
            overlap_thres=
            self.params.simple_ncs_from_pdb.domain_finding_parameters.similarity_threshold,
            minimum_size=min_size)
          # make sure each residue is in only 1 domain...
          domain_list=[]
          used_list=[]
          for match in domains.matches:
            domain_residues=[]
            keep_list=[]
            iselect=match[0].as_1d()
            for i in iselect:
              if not i in used_list:
                keep_list.append(i)
                # save list of actual residue ranges
                domain_residues.append(residue_list[i])
            used_list+=keep_list
            domain_list.append(self.remove_small(self.connect_range(
             self.list_to_range(domain_residues)) ))
          # sum up and remove any domains that are too small or too large rmsd:
          if self.verbose:print >>self.log,"\nDomains in ",id1,id2
          final_domain_list=[]
          final_matches=[]
          for domain_residues,match in zip(domain_list,domains.matches):
            if not domain_residues or not match:
               if self.verbose:print >>self.log,"Domain completely removed..."
               continue
            if self.verbose:
              print >>self.log,"Domain: ",domain_residues
              print >>self.log,"RMSD: ",match[3]
              print >>self.log,"N: ",match[4]
            if match[3]>self.max_rmsd_domain:
               if self.verbose: print >>self.log,"RMSD too large...rejecting this domain"
               continue
            final_domain_list.append(domain_residues)
            final_matches.append(match)
          if not final_matches or not final_domain_list:
            if self.verbose:print >>self.log,"No domains from this pair"
            continue

          if self.verbose:
           print >>self.log,"Evaluating rmsd for all chain pairs using these final matches"
          overall_rmsd_list=[]
          # get rmsd of the whole group using residues in final_domain_list
          number_of_residues_used_list=[]
          for domain_residues,match in zip(final_domain_list,final_matches):
            if self.verbose: print >>self.log,"domain...",domain_residues
            list_of_domain_residues=[] # testing this definition of the domain..
            for id in group: list_of_domain_residues.append(domain_residues)

            rmsd_list,r_list,trans_list,center_list,residues_in_common_list=\
                self.get_rmsd(group,hierarchy,
                remaining_chains,remaining_chain_ids,
                remaining_starting_residue_numbers,list_of_domain_residues)
            if self.verbose:
              print >>self.log,"rmsd_list for ",domain_residues,':',rmsd_list

            if not rmsd_list or len(rmsd_list)<2:
               overall_rmsd_list=[]
               break
            overall_rmsd=0.
            mean_residues_used=0
            for rmsd,res in zip(rmsd_list[1:],residues_in_common_list[1:]):
              overall_rmsd+=rmsd
              mean_residues_used+=res
            overall_rmsd=overall_rmsd/float(len(rmsd_list)-1)  # skip self-rmsd
            mean_residues_used=mean_residues_used/float(len(rmsd_list)-1)
            overall_rmsd_list.append(overall_rmsd)
            number_of_residues_used_list.append(mean_residues_used)
          test_rmsd=0.
          score=None
          if overall_rmsd_list:
            test_mean_length=0.
            for rmsd,res in zip(overall_rmsd_list,number_of_residues_used_list):
              test_rmsd+=rmsd
              test_mean_length+=res
            test_rmsd=test_rmsd/float(len(overall_rmsd_list))
            test_mean_length=test_mean_length/float(len(overall_rmsd_list))
            score=self.max(0.,(self.max_rmsd_domain-test_rmsd))*test_mean_length
          if self.verbose:
           print >>self.log,"rmsd for this match: ",test_rmsd,score
          if score is not None and (best_score is None or score>best_score):
            best_score=score
            best_rmsd=test_rmsd
            best_matches=final_matches
            best_domain_list=final_domain_list
            best_domains=domains
            best_number_of_residues_used_list=number_of_residues_used_list
            best_overall_rmsd_list=overall_rmsd_list
      if best_domains is not None and best_score is not None:
        if self.verbose:
          print >>self.log,"Best rmsd for group ",group,": ", best_rmsd,score
          print >>self.log,"Rmsd and residues used by domain: "
          for rmsd,res in zip(best_overall_rmsd_list,
              best_number_of_residues_used_list): print >>self.log,rmsd,res
        if self.verbose: best_domains.show(out=self.log)
        for domain_list in best_domain_list:
          invariant_groups.append(group)
          list_of_residue_ranges_for_this_group=[]
          for i in group:
            list_of_residue_ranges_for_this_group.append(domain_list)
          invariant_list_of_residue_range_list.append(
            list_of_residue_ranges_for_this_group)

    return invariant_groups,invariant_list_of_residue_range_list

  def remove_small (self,range_list):
    # remove all ranges that are smaller than
    # self.params.simple_ncs_from_pdb.domain_finding_parameters.min_contig_length//self.njump
    contig_length=\
     self.params.simple_ncs_from_pdb.domain_finding_parameters.min_contig_length//self.njump
    new_range_list=[]
    for range in range_list:
       if range[1]-range[0]+1 < contig_length:
          if self.verbose: print >>self.log,"Rejecting short range ",range
       else:
          new_range_list.append(range)
    return new_range_list

  def connect_range(self,range_list):
    # If two ranges are within
    # self.params.simple_ncs_from_pdb.domain_finding_parameters.smooth_length, connect them
    smooth=\
     self.params.simple_ncs_from_pdb.domain_finding_parameters.smooth_length*self.njump
    new_range_list=[]
    work_range=[]
    for range in range_list:
      if not work_range:
        work_range=range
      elif range[0] <= work_range[1]+smooth+1:
        work_range[1]=range[1]
      else:
        new_range_list.append(work_range)
        work_range=range
    if work_range: new_range_list.append(work_range)
    return new_range_list


  def get_chain(self,id1,remaining_chains,remaining_chain_ids,
        remaining_starting_residue_numbers):
    for chain,id,start in zip(remaining_chains,remaining_chain_ids,
        remaining_starting_residue_numbers):
      if id1==id:
         return chain,start
    return [],0

  def are_in_suggested_ncs_groups(self,test_group):
    for group in self.suggested_ncs_groups:
      are_in=True
      for id in test_group:
        if not id in group:
          are_in=False
          break
      if are_in:
        return True  # all of test_group is contained in a suggested group
    return False

  def are_in_different_suggested_ncs_groups(self,test_group):
    # find a group that contains a member of test_group. If any other member
    #  is in another group, return True
    for group1 in self.suggested_ncs_groups:
     for id in test_group:
      if id in group1:  # group1 contains an id of test_group

       for group2 in self.suggested_ncs_groups:
        if group1 == group2: continue
        for id in test_group:
         if id in group2:
          return True  # members of test_group is contained in > 1 suggested ncs grp
    return False

  def a_suggested_ncs_group_is_in_test_group(self,test_group):
    for group in self.suggested_ncs_groups:
      are_in=True
      for id in group:
        if not id in test_group:
          are_in=False
          break
      if are_in:
        return True # all of a suggested_ncs_group are contained in test_group
    return False

  def get_score_and_group(self,trial_id,remaining_ids,match_list):
    # score as the MINIMUM number of matches of anything to the trial_id
    # if all the members of a user-specified group are present, add 10000
    #   to score
    # return residue ranges that are not necessarily the same for all members
    # of a group...
    # 090708: allow user to increase score for large groups with
    #  self.maximize_size_of_groups

    score=None
    group=[]
    residue_range_list=[]

    # 021107:
    # Solution: need to use only those with id1==trial_id to keep numbering
    #  consistent. Then we also need matches to include both directions.

    # 021207 this didn't do it completely...we want to return offsets
    # based not on trial_id but on the firstid in the list...

    # the real problem is that here we use a seed (C) to pull out the
    #   set of chains (A C E G) and the offset matches C here. Then later
    #   we assume that the offset matched A. One solution might be
    # to always use the quick option...using each set of chains
    # identified with sequence. HERE JUST this is the best plan?

    # 021207 better: make sure that trial_id is first in list!
    # 090808 allow calling program to require a particular set of
    # chains to always be present

    for [matches,id1,id2,residue_range] in match_list:
      if id1==trial_id and id2==trial_id:
             group.append(trial_id)
             residue_range_list.append(residue_range)
    for [matches,id1,id2,residue_range] in match_list:
          if id1==trial_id and id2 in remaining_ids and not id2 in group:
             group.append(id2)
             residue_range_list.append(residue_range)
             if score is None:
               score=matches
             elif score > matches:
               score=matches  # decrease to lowest number of matches
    if self.require_equal_start_match:
      residue_range_list= \
         self.require_matches_to_start_at_same_place(group,residue_range_list)

    if self.maximize_size_of_groups and score is not None:
       score=score*math.sqrt(float(len(group)))

    if self.required_chains:
      found=False
      for id in group:
        if id in self.required_chains:
          found=True
          break
      if not found:
        score=None
        group=[]
        residue_range_list=[]

    if len(group)<2:
       score=None
       group=[]
       residue_range_list=[]
    if score and self.a_suggested_ncs_group_is_in_test_group(group):
      score+=10000

    return score,group,residue_range_list

  def require_matches_to_start_at_same_place(self,group,residue_range_list):
    new_residue_range_list=[]
    # find last start for each range
    if len(residue_range_list)<1: return 0,[]
    number_of_ranges=len(residue_range_list[0])
    last_start_list=number_of_ranges*[None]
    for rr,g in zip(residue_range_list,group):
      if number_of_ranges is None:
         number_of_ranges = len(rr)
      #print "group, range: ",g,rr
      new_start_list=[]
      rr_mod=rr[:len(last_start_list)]+(len(last_start_list)-len(rr))*[None]
      for r,s in zip(rr_mod,last_start_list):
        if r is not None and (s is None or r[0] > s):
           s = r[0]
        new_start_list.append(s)
      last_start_list=new_start_list
    #print "start_list; ",last_start_list
    # now adjust all the starting to be last_start_list
    for rr,g in zip(residue_range_list,group):
      new_rr=[]
      rr_mod=rr[:len(last_start_list)]+(len(last_start_list)-len(rr))*[None]
      for r,s in zip(rr_mod,last_start_list):
        if r is not None and s is not None and r[0] != s:
           r[0] = s  # modifies r inside rr inside residue_range_list
    return residue_range_list


  def matches_approx(self,chain1,start1,chain2,start2,id1,id2,min_length=5,
        min_percent=95,exact_match_only=False):

    # return number of residues matching total, and a
    # list of residue ranges over which chain1 and chain2 match
    # within the tolerance min_percent and each with minimum length
    #  of min_length

    # NOTE: residue range is in terms of position in chain (relative to start,
    # and if self.njump>1 then the residue numbers are in  units of njump

    id=id1+id2
    if id in self.offset_dict.keys():
      offset=self.offset_dict[id]
    else:
      print >>self.log,"NO ID in matches_approx: ",id
      offset=0

    start2_use=start2+offset  # equivalent residue number for start2
    keep_range_list=[] # a list of ranges like [[1,3],[7-9]...]

    if not chain1 or not chain2: return 0,keep_range_list
    end1=start1+(len(chain1)-1)*self.njump
    end2=start2_use+(len(chain2)-1)*self.njump
    start=start1
    end=end1
    if start2_use>start1: start=start2_use
    if end2<end1: end=end2
    # we start,end at residue number start,end, numbered according to chain 1

    start_njump=start//self.njump
    end_njump=end//self.njump
    if self.verbose:
       print >>self.log,"Finding matching segments in range from ",\
       start," to ",end," for chains ",chain1,chain2
    # iteratively select the longest segment that has ends that match with
    # at least 2 consecutive residues matching at each end...
    # and >= min_percent and >= min_length until there is none...
    # then sort in order of first residue number

    # if exact_match_only then skip if first residues or last do not match...

    matched_residues=[]

    for i in xrange(start,end+self.njump,self.njump):
      res1=chain1[(i-start1)//self.njump]
      res2=chain2[(i-start2_use)//self.njump]
      if res1 and res2:
        if res1==res2:    #Matching
          matched_residues.append(True)
        else:
          matched_residues.append(False)
          if exact_match_only:
            return 0,keep_range_list
      else:
          matched_residues.append(False)

    segments=[]
    still_looking=True
    possible=len(matched_residues)
    matching=0
    while still_looking:
      best_segment,best_length=self.find_best_segment(
           matched_residues,segments,min_percent,min_length)
      if best_segment is None:
         still_looking=False
      else:
         segments.append(best_segment)
         matching+=best_length
    segments.sort()
    if self.verbose:
       print >>self.log,"List of segments: ",segments

    # sort and translate segments into residue numbers...
    for [pos_start,pos_end] in segments:
      keep_range_list.append(
         [pos_start*self.njump+start,pos_end*self.njump+start])

    if self.verbose:
     print >>self.log,"Possible: ",possible," Matching: ",matching
     print >>self.log,"Residue range used: ",segments
    return matching,keep_range_list

  def find_best_segment(self,matched_residues,segments,min_percent,min_length):
    # mask out all the residues that are used:
    available_residues=matched_residues
    for [start,end] in segments:
      for i in xrange(start,end+1):
        available_residues[i]=False
    if self.verbose:
        print >>self.log,"Available_residues: ",available_residues

    # and find longest segment in what remains, requiring matches on each
    # end and percent >= min_percent and length >=min_length
    # NOTE: min_length applies to absolute length; use min_length//self.njump
    best_start=0
    best_end=0
    best_length=0
    for start in xrange(0,len(available_residues)):
      if not available_residues[start]:continue
      for end in xrange(start,len(available_residues)):
        if not available_residues[end]:continue
        n_match=available_residues[start:end+1].count(True)
        n_not_match=available_residues[start:end+1].count(False)
        n_total=n_match+n_not_match
        if n_total <min_length//self.njump: continue
        if not n_total:
          percent=0.
        else:
          percent=100.*float(n_match)/float(n_total)
        if percent < min_percent: continue
        if n_match>best_length:
           best_length=n_match
           best_start=start
           best_end=end
    if self.verbose:
      print >>self.log,"Best segment: ",best_start,best_end,best_length
    if best_length:
      return [best_start,best_end],best_length
    else:
      return None,0





  def get_chain_list(self,hierarchy):
    # get list of chains. Each chain list is a list of (blanks or residue names)
    #  a chain list has a slot for each residue in a chain, it can be occupied
    #  or not.  Two chains are the same parent chain if all residues that are
    #  present match exactly.
    #  each chain starts at a specific residue number, stored as
    #    starting_residue_numbers for that chain

    # if self.njump is not 1, then we are only going to take residues with
    # numbers such that mod(resno,self.njump)=0
    # so if self.njump=10 then we take residue 0,10,20,30,40 and put them
    #  in position 1,2,3,4,5

    # ignore all in ignore_chains
    chains=[]
    chain_ids=[]
    starting_residue_numbers=[]
    model_number=0
    self.total_residues=0
    for model in hierarchy.models():
      model_number+=1
      if model_number>1: continue  # take first model only
      for chain in model.chains():
        if chain.id in self.ignore_chains: continue
        chain_id=chain.id
        if chain_id in chain_ids:
          if self.verbose:print >>self.log,\
                "NOTE: duplicate chain ID (ignoring): ",chain_id
          continue
        if self.verbose: print >>self.log,"Examining chain ",chain_id
        conformer_number=0
        for conformer in chain.conformers():
          conformer_number+=1
          if conformer_number>1:continue  # take first conformer only
          highest_residue_number=-999999
          lowest_residue_number=  9999999
          number_of_good_residues=0
          for residue in conformer.residues():  # get highest residue number
            resseq_int = residue.resseq_as_int()
            self.total_residues+=1
            if self.skip(resseq_int): continue  # based on self.njump
            if resseq_int>highest_residue_number:
               highest_residue_number=resseq_int
            if resseq_int<lowest_residue_number:
               lowest_residue_number=resseq_int
          if self.verbose: print >>self.log,"Residue range found: ",\
              lowest_residue_number,highest_residue_number
          seq=[]
          range=(highest_residue_number-lowest_residue_number)//self.njump+1
          for i in xrange(range):
             seq.append(None)
          residues_seen_already=[]
          for residue in conformer.residues():
            resseq_int = residue.resseq_as_int()
            if self.skip(resseq_int): continue  # based on self.njump
            if resseq_int in residues_seen_already:continue
            residues_seen_already.append(resseq_int)
            one_char=self.OneChar(residue.resname)
            if one_char:  # skip anything that is not protein and not DNA/RNA
              found_ca_or_p=False  # require ' CA ' or ' P  ': see get_xyz_list
              for atom in residue.atoms():
                 if (atom.name == " CA " or atom.name == " P  "):
                   found_ca_or_p=True
                   break
              if found_ca_or_p:
                number_of_good_residues+=1
                #    position of this residue
                i=(resseq_int-lowest_residue_number)//self.njump
                seq[i]=one_char
                if self.verbose: print >>self.log,"Residue",residue.resseq,one_char

        if number_of_good_residues >=self.min_length//self.njump:
          chains.append(seq)
          chain_ids.append(chain_id)
          starting_residue_numbers.append(lowest_residue_number)


    # Now if some chains are same as other chains but residue numbers are
    # offset, figure that out...
    offset_dict={}  # add offset_dict[id1+id2] to start2 to match seq of chain1
    for chain1,id1,start1 in zip(chains,chain_ids,starting_residue_numbers):
      for chain2,id2,start2 in zip(chains,chain_ids,starting_residue_numbers):
        offset=self.find_offset(chain1,start1,chain2,start2)
        if offset is None:
          offset=0
          if self.verbose:print >>self.log,"FAILED OFFSET ",offset,id1,id2,id1+id2
        else:
          if self.verbose:print >>self.log,"OFFSET ",offset,id1,id2,id1+id2
        offset_dict[id1+id2]=offset

    return chains,chain_ids,starting_residue_numbers,offset_dict

  def skip(self,resno):  # skip unless resno falls on n * self.njump
    if self.njump <=1:
      return False
    if int(0.5+math.fmod(100*self.njump+resno,self.njump))!=0:
      return True

  def find_offset(self,chain1,start1,chain2,start2):
    # return offset that makes chain1 same as chain2 (no mismatches allowed)
    # list chain1:
    # first residue of chain2 can match any one of residues of chain1.
    # offset can be:
    # start2+offset=start1-len(chain2)+1  to start2+offset=start1+len(chain1)
    best_overlap=0
    best_offset=None
    for offset_a in xrange(-len(chain2)+1,len(chain1)+1):
     offset=start1-start2+offset_a*self.njump
     overlap=self.get_overlap(chain1,start1,chain2,start2,offset)
     if overlap and overlap>best_overlap:
       best_overlap=overlap
       best_offset=offset
    return best_offset

  def get_overlap(self,chain1,start1,chain2,start2,offset):
    offset1=0
    offset2=offset
    n_match=0
    mismatches=0
    for i in xrange(len(chain1)):
      resno=i*self.njump+start1+offset1
      res1=chain1[i]
      #resno=j*self.njump+start2+offset2
      j=(resno-start2-offset2)//self.njump
      if j>=0 and j<len(chain2):
        res2=chain2[j]
      else:
        res2=None
      if res1!=None and res2!=None:
        if res1==res2:
          n_match+=1
        else:
          mismatches+=1
    if 100.*float(n_match)/float(self.max(1,n_match+mismatches)) < self.min_percent:
      return None
    if n_match >= self.min_length//self.njump:
      return n_match
    else:
      return None

  def OneChar(self,residue_name):  # return one char for DNA/RNA/Protein
    # 070508 catch +G or DG and call both of them G:
    # 2010-12-01 also Ar == A (ribo) call it A
    edited_residue_name=residue_name.lower().replace(" ","")
    if len(edited_residue_name)==2 and edited_residue_name[0] in ['d','+']:
      edited_residue_name=edited_residue_name[1]
    elif len(edited_residue_name)==2 and edited_residue_name[1] in ['r']:
      edited_residue_name=edited_residue_name[0] #2010-12-01
    elif edited_residue_name in ['ade','cyt','gua','thy','uri']:
       for one,three in zip(
           ['a','c','g','t','u'],
           ['ade','cyt','gua','thy','uri']):
         if edited_residue_name==three:
            edited_residue_name=one
            break
    aa_names= [
      'GLY','ALA','SER','VAL',
      'ILE','LEU','MET','CYS',
      'PHE','TYR','LYS','ARG',
      'TRP','HIS','GLU','ASP',
      'GLN','ASN','PRO','THR',
      'MSE', # treat MSE as MET for now
      'G','C','A','T','U']
    aa_letters=[
       'G','A','S','V',
       'I','L','M','C',
       'F','Y','K','R',
       'W','H','E','D',
       'Q','N','P','T',
       'M',                   # treat MSE as MET for now
       'G','C','A','T','U']

    for name,letter in zip(aa_names,aa_letters):
      if string.replace(string.lower(edited_residue_name)," ","") == \
         string.replace(string.lower(name)," ",""):
        return string.lower(letter)
    return None


  def is_pdb(self,file_name): # XXX REMOVE PDB 2013-08-20
    if not os.path.isfile(file_name): return False
    for line in open(file_name).readlines():
       if line[0:4]=='ATOM' or line[0:6]=='HETATM':
         return True
    return False




  def get_help(self,command_name,master_params,summary):
    print >>self.log,summary
    print >>self.log,"Default parameters:"
    master_params.format(python_object=
          master_params.fetch(sources=[]).extract()).show()

  def raise_missing(self,what):
      raise Sorry("""\
          Missing file name for %(what)s :
          Please add %(what)s=file_name
          to the command line to specify %(what)s .""" % vars())

  def get_summary_and_header(self,command_name):
    header="\n"
    header+="\n#                       "+str(command_name)
    header+="\n#"
    header+="\n# Find ncs among chains in a PDB file "
    header+="\n\n# type phenix.doc for help\n"

    summary= "usage: %s protein.pdb [parameter=value ...]" % command_name
    summary+="\n\nYou can set any parameter by specifying its path: "
    summary+="\nncs.max_rmsd =3 sets rmsd to 3"
    summary+="\nTo test use: %s exercise\n" % command_name
    return summary,header

  def exercise(self):  # run with a few ha sites
    text="""
ATOM      2  CA  GLY A   1      43.603 -11.488  24.325  1.00 35.57      2MLT 113
ATOM      6  CA  ILE A   2      44.200  -8.183  22.475  1.00 27.55      2MLT 117
ATOM     14  CA  GLY A   3      43.999 -10.264  19.329  1.00 21.05      2MLT 125
ATOM     18  CA  ALA A   4      40.378 -11.260  20.106  1.00 21.80      2MLT 129
ATOM     23  CA  VAL A   5      39.355  -7.658  21.083  1.00 19.34      2MLT 134
ATOM     30  CA  LEU A   6      41.062  -6.432  17.957  1.00 17.59      2MLT 141
ATOM     38  CA  LYS A   7      39.079  -8.646  15.636  1.00 22.55      2MLT 149
ATOM     47  CA  VAL A   8      35.792  -7.369  17.211  1.00 20.52      2MLT 158
ATOM     54  CA  LEU A   9      36.899  -3.759  16.725  1.00 16.83      2MLT 165
ATOM     62  CA  THR A  10      37.338  -4.508  13.084  1.00 19.41      2MLT 173
ATOM     69  CA  THR A  11      34.132  -6.405  12.343  1.00 24.14      2MLT 180
ATOM     76  CA  GLY A  12      31.584  -6.595  15.140  1.00 24.17      2MLT 187
ATOM     80  CA  LEU A  13      31.923  -2.919  16.364  1.00 23.24      2MLT 191
ATOM     88  CA  PRO A  14      31.026  -1.278  13.030  1.00 17.52      2MLT 199
ATOM     95  CA  ALA A  15      27.822  -3.418  12.724  1.00 17.10      2MLT 206
ATOM    100  CA  LEU A  16      26.958  -2.649  16.351  1.00 18.20      2MLT 211
ATOM    108  CA  ILE A  17      27.343   1.056  15.618  1.00 20.41      2MLT 219
ATOM    116  CA  SER A  18      24.825   0.827  12.744  1.00 19.98      2MLT 227
ATOM    122  CA  TRP A  19      22.492  -1.085  15.081  1.00 15.72      2MLT 233
ATOM    136  CA  ILE A  20      22.628   1.633  17.766  1.00 15.67      2MLT 247
ATOM    144  CA  LYS A  21      21.888   4.236  15.157  1.00 19.84      2MLT 255
ATOM    153  CA  ARG A  22      18.740   2.273  14.020  1.00 20.38      2MLT 264
ATOM    164  CA  LYS A  23      17.500   1.928  17.550  1.00 22.62      2MLT 275
ATOM    173  CA  ARG A  24      18.059   5.674  18.276  1.00 27.11      2MLT 284
ATOM    184  CA  GLN A  25      15.836   6.730  15.339  1.00 37.50      2MLT 295
ATOM    193  CA  GLN A  26      13.132   4.360  16.583  1.00 46.66      2MLT 304
ATOM    204  CA  GLY B 101      26.196  11.215  25.772  1.00 31.28      2MLT 315
ATOM    208  CA  ILE B 102      25.695   8.457  23.127  1.00 26.61      2MLT 319
ATOM    216  CA  GLY B 103      25.775  10.608  19.984  1.00 20.83      2MLT 327
ATOM    220  CA  ALA B 104      29.261  12.192  20.662  1.00 22.19      2MLT 331
ATOM    225  CA  VAL B 105      30.560   8.605  21.561  1.00 20.43      2MLT 336
ATOM    232  CA  LEU B 106      29.339   7.258  18.306  1.00 15.03      2MLT 343
ATOM    240  CA  LYS B 107      30.789  10.122  16.413  1.00 20.35      2MLT 351
ATOM    249  CA  VAL B 108      34.238   9.637  18.027  1.00 23.38      2MLT 360
ATOM    256  CA  LEU B 109      34.150   5.918  17.183  1.00 20.22      2MLT 367
ATOM    264  CA  THR B 110      33.081   6.344  13.562  1.00 23.85      2MLT 375
ATOM    271  CA  THR B 111      35.833   8.796  12.896  1.00 26.41      2MLT 382
ATOM    278  CA  GLY B 112      38.482   7.602  15.318  1.00 22.53      2MLT 389
ATOM    282  CA  LEU B 113      38.268   3.810  15.278  1.00 23.61      2MLT 393
ATOM    290  CA  PRO B 114      39.945   3.149  11.872  1.00 20.15      2MLT 401
ATOM    297  CA  ALA B 115      43.145   4.994  13.009  1.00 22.26      2MLT 408
ATOM    302  CA  LEU B 116      43.070   3.167  16.352  1.00 18.86      2MLT 413
ATOM    310  CA  ILE B 117      42.952  -0.212  14.601  1.00 14.88      2MLT 421
ATOM    318  CA  SER B 118      45.891   0.689  12.345  1.00 17.71      2MLT 429
ATOM    324  CA  TRP B 119      47.907   1.909  15.276  1.00 17.53      2MLT 435
ATOM    338  CA  ILE B 120      47.296  -1.149  17.359  1.00 14.07      2MLT 449
ATOM    346  CA  LYS B 121      48.314  -3.437  14.512  1.00 20.04      2MLT 457
ATOM    355  CA  ARG B 122      51.517  -1.513  14.068  1.00 25.56      2MLT 466
ATOM    366  CA  LYS B 123      52.367  -1.651  17.742  1.00 23.85      2MLT 477
ATOM    375  CA  ARG B 124      51.670  -5.431  17.781  1.00 28.09      2MLT 486
ATOM    386  CA  GLN B 125      54.153  -6.012  14.932  1.00 40.40      2MLT 497
ATOM    395  CA  GLN B 126      56.818  -4.124  16.883  1.00 45.23      2MLT 506
        """

    expected_result="""refinement.ncs.restraint_group {
reference = chain 'A' and (resseq 1:26 )
selection = chain 'B' and (resseq 101:126 )
}"""

    file_name='temp.pdb'
    f=open(file_name,'w')
    f.write(text)
    f.close()
    args=[file_name]
    args.append("min_length=1")
    args.append("min_percent=10")
    run_simple_ncs_from_pdb=simple_ncs_from_pdb(args=args,log=self.log)

    result=open('simple_ncs_from_pdb.ncs').read()
    if result and \
       string.replace(string.replace(result," ",""),"\n","")== \
       string.replace(string.replace(expected_result," ",""),"\n",""):
      print >>self.log,'OK'
    elif result:
      print >>self.log,"Output does not match. Result: ",result,\
        "\nExpected result: ",expected_result
    else:
      print >>self.log,"No result"

  def process_inputs(self,args_read):
   if not args_read: args_read=[]
   if self.log is None :
    self.log=sys.stdout
   self.args=[]
   sugg=None
   self.suggested_ncs_groups=None
   self.ignore_chains=[]
   self.allow_recursion=True
   self.exact_match_only=False

   for arg in args_read:
    if arg[:5]=='sugg=':
      sugg=arg[5:]
      print >>self.log,"Suggested groups:",sugg
    elif arg[:7]=='ignore=':
      self.ignore_chains=[]
      for char in arg[7:]:
        self.ignore_chains+=char
      print >>self.log,"Ignoring chains:",self.ignore_chains
    elif arg=='exact_match':
      self.exact_match_only=True
    elif arg=='no_recursion':
      self.allow_recursion=False
    else:
      self.args.append(arg)
   self.suggested_ncs_groups=sugg

# FIXME: can't just pass a phil file path for some reason
class launcher (object) :
  def __init__ (self, params, tmp_dir) :
    adopt_init_args(self, locals())

  def __call__ (self) :
    os.chdir(self.tmp_dir)
    ncs = simple_ncs_from_pdb(
      args=None,
      params=self.params,
      quiet=True,
      exclude_h=True,
      exclude_d=True,
      log=sys.stdout)
    return ncs.ncs_object

if (__name__ == "__main__"):
  argument_list=sys.argv[1:]
  if True or is_debug(argument_list).value:

    sys_stdout_sav=sys.stdout
    simple_ncs=simple_ncs_from_pdb(args=sys.argv[1:])
    sys.exit(0)
  try:
    sys_stdout_sav=sys.stdout
    simple_ncs=simple_ncs_from_pdb(args=sys.argv[1:])
  except KeyboardInterrupt:
    pass
  except Exception, e:
    print "\n************************************************"
    print e
    print "\n************************************************"
    if sys.stdout != sys_stdout_sav:
      sys.stdout=sys_stdout_sav # restore output stream so we can print error
      print "\n************************************************"
      print e
      print "\n************************************************"
    from phenix.utilities.is_raise_sorry import is_raise_sorry
    if is_raise_sorry(argument_list).value:
      from libtbx.utils import Sorry
      raise Sorry(e)

