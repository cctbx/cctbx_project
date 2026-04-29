
from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from libtbx.utils import null_out, Sorry
from libtbx.test_utils import approx_equal
from libtbx import easy_run
import libtbx.load_env
from cctbx import geometry_restraints
import iotbx.pdb
import iotbx.pdb.secondary_structure
import iotbx.phil
import mmtbx.secondary_structure
from mmtbx.secondary_structure import proteins, nucleic_acids
import sys, os
import time

from itertools import groupby
from operator import itemgetter
from six.moves import range

def contiguous_ss_selections(pdb_hierarchy):
  """
  Return list of atom selections that contiguous split molecule by SS and loops.
  """
  pdb_hierarchy.atoms().reset_i_seq()
  params = mmtbx.secondary_structure.sec_str_master_phil.extract()
  params.secondary_structure.protein.search_method="from_ca"
  ssm = mmtbx.secondary_structure.manager(
    pdb_hierarchy                = pdb_hierarchy,
    sec_str_from_pdb_file        = None,
    params                       = params.secondary_structure,
    log                          = null_out())
  n_atoms = pdb_hierarchy.atoms_size()
  ss_all = ssm.helix_selection().iselection()
  ss_all.extend(ssm.beta_selection().iselection())
  ss_all.extend(ssm.base_pair_selection().iselection())
  chain_selections = []
  for k, g in groupby(enumerate(ss_all), lambda i_x: i_x[0]-i_x[1]):
    chain_selections.append(flex.size_t([itemgetter(1) for _ in g]))
  #
  if(ss_all.size()>0):
    ss_all_as_one_array_bool = flex.bool(n_atoms, ss_all)
  else:
    ss_all_as_one_array_bool = flex.bool(n_atoms, False)
  #
  iseqs = []
  for atom in pdb_hierarchy.atoms():
    if(not ss_all_as_one_array_bool[atom.i_seq]):
      iseqs.append(atom.i_seq)
  for k, g in groupby(enumerate(iseqs), lambda i_x1: i_x1[0]-i_x1[1]):
    isels = [itemgetter(1) for _ in  g]
    chain_selections.append(flex.size_t(isels))
  # DEBUG
  C1 = flex.size_t()
  for s in chain_selections:
    C1.extend(s)
  s = flex.sort_permutation(C1)
  C1 = C1.select(s)
  C2 = flex.size_t(range(n_atoms))
  assert approx_equal(C1, C2)
  assert C1.size()==C2.size()
  #
  return chain_selections

#
# Notes. Now this class trying to keep updated SS information in form of
# phil parameters and as iotbx/secondary_structure/annotation. This is
# not the optimal way, it would be better to work only with annotation, but
# this requires massive rewriting of SS GUI which works only with phil...
#

sec_str_master_phil_str = """
secondary_structure
  .style = box auto_align noauto
{
  protein
    .style = box auto_align
  {
    enabled = True
      .type = bool
      .style = noauto
      .help = Turn on secondary structure restraints for protein
    search_method = *ksdssp from_ca cablam
      .type = choice
      .help = Particular method to search protein secondary structure.
    distance_ideal_n_o = 2.9
      .type = float
      .short_caption = Ideal N-O distance
      .help = Target length for N-O hydrogen bond
    distance_cut_n_o = 3.5
      .type = float
      .short_caption = N-O distance cutoff
      .help = Hydrogen bond with length exceeding this value will not be \
       established
    remove_outliers = True
      .type = bool
      .short_caption = Filter out h-bond outliers in SS
      .style = tribool
      .help = If true, h-bonds exceeding distance_cut_n_o length will not be \
       established
    restrain_hbond_angles = True
      .type = bool
      .short_caption = Restrain angles around hbonds in alpha helices
    %s
    %s
  }
  nucleic_acid
    .caption = If sigma and slack are not defined for nucleic acids, the \
      overall default settings for protein hydrogen bonds will be used \
      instead.
    .style = box auto_align
  {
    enabled = True
      .type = bool
      .style = noauto
      .help = Turn on secondary structure restraints for nucleic acids
    %s
  }
  ss_by_chain = True
    .type = bool
    .help = Only applies if search_method = from_ca. \
            Find secondary structure only within individual chains. \
            Alternative is to allow H-bonds between chains. Can be \
            much slower with ss_by_chain=False. If your model is complete \
            use ss_by_chain=True. If your model is many fragments, use \
            ss_by_chain=False.
    .short_caption = Secondary structure by chain
    .expert_level = 1
  from_ca_conservative = False
    .type = bool
    .help = various parameters changed to make from_ca method more \
      conservative, hopefully to closer resemble ksdssp.
    .short_caption = Conservative mode of from_ca
  max_rmsd = 1
    .type = float
    .help = Only applies if search_method = from_ca. \
            Maximum rmsd to consider two chains with identical sequences \
            as the same for ss identification
    .short_caption = Maximum rmsd
    .expert_level = 3
  use_representative_chains = True
    .type = bool
    .help = Only applies if search_method = from_ca. \
            Use a representative of all chains with the same sequence. \
            Alternative is to examine each chain individually. Can be \
            much slower with use_representative_of_chain=False if there \
            are many symmetry copies. Ignored unless ss_by_chain is True.
    .short_caption = Use representative chains
    .expert_level = 3
  max_representative_chains = 100
    .type = float
    .help = Only applies if search_method = from_ca. \
            Maximum number of representative chains
    .short_caption = Maximum representative chains
    .expert_level = 3

  enabled = False
    .short_caption = Use secondary structure restraints
    .type = bool
    .style = noauto bold
    .help = Turn on secondary structure restraints (main switch)
}
""" % (
       proteins.helix_group_params_str,
       proteins.sheet_group_params_str,
       nucleic_acids.dna_rna_params_str)

master_phil_str = sec_str_master_phil_str # for docs

sec_str_master_phil = iotbx.phil.parse(sec_str_master_phil_str)

class manager(object):
  def __init__(self,
                pdb_hierarchy,
                atom_selection_cache=None,
                geometry_restraints_manager=None,
                sec_str_from_pdb_file=None,
                params=None,
                was_initialized=False,
                mon_lib_srv=None,
                verbose=-1,
                show_summary_on=True,
                log=sys.stdout):
    self.pdb_hierarchy = pdb_hierarchy
    self.mon_lib_srv = mon_lib_srv
    self.actual_sec_str = None # Resulting SS either from pdb file header or
    # from automatic search. Instance of iotbx.pdb.secondary_structure.annotation
    self.grm = geometry_restraints_manager
    self.sec_str_from_pdb_file = sec_str_from_pdb_file
    self.params = sec_str_master_phil.extract()
    if params is not None:
      self.params.secondary_structure = params

    # checking params
    default_params = sec_str_master_phil.extract()
    ss_by_chain_set = self.params.secondary_structure.ss_by_chain != \
        default_params.secondary_structure.ss_by_chain
    max_rmsd_set = self.params.secondary_structure.max_rmsd != \
        default_params.secondary_structure.max_rmsd
    use_representative_chains_set = self.params.secondary_structure.use_representative_chains != \
        default_params.secondary_structure.use_representative_chains
    max_representative_chains_set = self.params.secondary_structure.max_representative_chains != \
        default_params.secondary_structure.max_representative_chains
    from_ca_conservative_set = self.params.secondary_structure.from_ca_conservative != \
        default_params.secondary_structure.from_ca_conservative
    from_ca_method_selected = self.params.secondary_structure.protein.search_method != "from_ca"
    for par, s in [(ss_by_chain_set, 'ss_by_chain'),
        (max_rmsd_set, 'max_rmsd'),
        (from_ca_conservative_set, 'from_ca_conservative'),
        (use_representative_chains_set, 'use_representative_chains'),
        (max_representative_chains_set, 'max_representative_chains')]:
      if par and from_ca_method_selected:
        raise Sorry("%s parameter is only used when search_method=from_ca.\n" % s +\
            "Please do not set it for other methods to avoid confusion." )

    self.verbose = verbose
    self.show_summary_on = show_summary_on
    self.log = log
    self.def_params = sec_str_master_phil.extract()
    self.stats = {'n_protein_hbonds':0, 'n_na_hbonds':0, 'n_na_hbond_angles':0,
        'n_na_basepairs':0, 'n_na_stacking_pairs':0}

    if not pdb_hierarchy:
      pdb_hierarchy = iotbx.pdb.hierarchy.root()
    atoms = pdb_hierarchy.atoms()
    i_seqs = atoms.extract_i_seq()
    if (i_seqs.all_eq(0)):
      atoms.reset_i_seq()
      i_seqs = atoms.extract_i_seq()
    self.n_atoms = atoms.size()
    self._was_initialized = was_initialized
    self.selection_cache = atom_selection_cache
    if self.selection_cache is None:
      self.selection_cache = pdb_hierarchy.atom_selection_cache()
    self.pdb_atoms = atoms
    self.initialize()

  @staticmethod
  def get_default_ss_params():
    return sec_str_master_phil.fetch().extract()

  def as_phil_str(self, master_phil=sec_str_master_phil):
    # used in secondary_structure_restraints...
    return master_phil.format(python_object=self.params)

  def initialize(self):
    if not self._was_initialized :
      self.find_automatically()
      if self.show_summary_on:
        self.show_summary()
      self._was_initialized = True

  def find_automatically(self):
    # find_automatically = self.params.secondary_structure.find_automatically
    protein_found = False
    use_segid = self.selection_cache.segid.size() > 1
    segids = {}
    if use_segid:
      for si in self.selection_cache.segid.keys():
        segids[si] = ""

    # XXX: check for presence of protein first?
    protein_ss_definition_present = False
    if (len(self.params.secondary_structure.protein.helix) > 1 or
        len(self.params.secondary_structure.protein.sheet) > 1):
      protein_ss_definition_present = True
    elif (len(self.params.secondary_structure.protein.helix) > 0 and
        self.params.secondary_structure.protein.helix[0].selection is not None):
      protein_ss_definition_present = True
    elif (len(self.params.secondary_structure.protein.sheet) > 0 and
        self.params.secondary_structure.protein.sheet[0].first_strand is not None):
      protein_ss_definition_present = True
    t0 = time.time()
    if protein_ss_definition_present:
      self.apply_phil_str(phil_string=None, phil_params=self.params, log=self.log)
    else:
      if (not protein_ss_definition_present and
          (self.sec_str_from_pdb_file is None or self.sec_str_from_pdb_file.is_empty())):
        if(self.verbose>0):
          print("No existing protein secondary structure definitions " + \
          "found in .pdb file or phil parameters.", file=self.log)
        # find_automatically = True
      else:
        # find_automatically = False
        protein_found = True
      if not protein_found:
        ss_params = []
        annot = None
        if use_segid:
          # Could get rid of this 'if' clause, but I want to avoid construction of
          # atom_selection_cache and selections when there is no segids in pdb
          # which is majority of cases.
          whole_annotation = iotbx.pdb.secondary_structure.annotation(helices=[], sheets=[])
          for segid in segids:
            isel = self.selection_cache.selection("segid '%s'" % segid).iselection()
            selected_pdb_h = self.pdb_hierarchy.select(isel)
            if selected_pdb_h.contains_protein():
              annot = self.find_sec_str(pdb_hierarchy=selected_pdb_h)
              if annot is not None:
                ss_phil = annot.as_restraint_groups(log=self.log,
                  prefix_scope="secondary_structure",
                  add_segid=segid)
                ss_params.append(ss_phil)
                whole_annotation.add_helices_and_sheets_simple(other_annot=annot)
          annot = whole_annotation
        else:
          if self.pdb_hierarchy.contains_protein():
            annot = self.find_sec_str(pdb_hierarchy=self.pdb_hierarchy)
            if annot is not None:
              ss_phil = annot.as_restraint_groups(log=self.log,
                prefix_scope="secondary_structure")
              ss_params.append(ss_phil)
              # self.actual_sec_str = annot
        ss_params_str = "\n".join(ss_params)
        self.apply_phil_str(ss_params_str, annot=annot, log=self.log)
      else:
        if (self.sec_str_from_pdb_file is not None):
          # self.actual_sec_str = self.sec_str_from_pdb_file
          removed_annot = self.sec_str_from_pdb_file.remove_empty_annotations(self.pdb_hierarchy)
          if not removed_annot.is_empty():
            print("\n  WARNING! Some SS annotations were removed because the model\n" +\
                "  does not have atoms for them:", file=self.log)
            rem_str = '  ' + removed_annot.as_pdb_str()
            rem_str = rem_str.replace('\n', '\n  ')
            print(rem_str+'\n', file=self.log)
          ss_params_str = self.sec_str_from_pdb_file.as_restraint_groups(
              log=self.log,
              prefix_scope="secondary_structure")
          self.apply_phil_str(ss_params_str, annot=self.sec_str_from_pdb_file, log=self.log)
        # else:
        #   # got phil SS, need to refactor later, when the class is fully
        #   # converted for operation with annotation object
        #   # phil_string = sec_str_master_phil.format(python_object=self.params)
        #   self.apply_phil_str(phil_string=None, phil_params=self.params, log=self.log)

    t1 = time.time()
    # print >> log, "    Time for finding protein SS: %f" % (t1-t0)
    # Step 2: nucleic acids
    if self.params.secondary_structure.nucleic_acid.enabled:
      na_found = False
      na_ss_definition_present = False
      if (len(self.params.secondary_structure.nucleic_acid.base_pair) > 1 or
          len(self.params.secondary_structure.nucleic_acid.stacking_pair) > 1):
        na_ss_definition_present = True
      elif (len(self.params.secondary_structure.nucleic_acid.base_pair) > 0 and
          self.params.secondary_structure.nucleic_acid.base_pair[0].base1 is not None):
        na_ss_definition_present = True
      elif (len(self.params.secondary_structure.nucleic_acid.stacking_pair) > 0 and
          self.params.secondary_structure.nucleic_acid.stacking_pair[0].base1 is not None):
        na_ss_definition_present = True

      if na_ss_definition_present:
        na_found = True

      if not na_found:
        pairs = self.find_base_pairs(log=self.log, use_segid=use_segid, segids=segids)
        if len(pairs) > 1:
          bp_phil = iotbx.phil.parse(pairs)
          bp_params = sec_str_master_phil.fetch(source=bp_phil).extract()
          self.params.secondary_structure.nucleic_acid.base_pair = \
            bp_params.secondary_structure.nucleic_acid.base_pair
          self.params.secondary_structure.nucleic_acid.stacking_pair = \
            bp_params.secondary_structure.nucleic_acid.stacking_pair
      t2 = time.time()
    # print >> log, "    Time for finding NA SS: %f" % (t2-t1)

  def records_for_pdb_file(self):
    if self.actual_sec_str is not None:
      return self.actual_sec_str.as_pdb_str()
    elif self.sec_str_from_pdb_file is not None:
      return self.sec_str_from_pdb_file.as_pdb_str()
    return None

  def find_sec_str(self, pdb_hierarchy):
    if ((pdb_hierarchy.atoms_size() > 99999 or not pdb_hierarchy.fits_in_pdb_format()) and
        self.params.secondary_structure.protein.search_method == "ksdssp"):
      print("\n".join([
          "Warning!!! ksdssp method is not applicable for",
          "structures that cannot fit in PDB format. Switching to from_ca."]), file=self.log)
      self.params.secondary_structure.protein.search_method = "from_ca"
    if self.params.secondary_structure.protein.search_method == "ksdssp":
      pdb_str = pdb_hierarchy.as_pdb_string() # PDB OK
      print("  running ksdssp...", file=self.log)
      (records, stderr) = run_ksdssp_direct(pdb_str)
      return iotbx.pdb.secondary_structure.annotation.from_records(
          records=records,
          log=self.log)
    elif self.params.secondary_structure.protein.search_method == "from_ca":
      from mmtbx.secondary_structure import find_ss_from_ca
      from_ca_args = []
      ca_label = "liberal"
      if self.params.secondary_structure.from_ca_conservative:
        from_ca_args = ["alpha.rise_tolerance=0.13",
                        "beta.rise_tolerance=0.3",
                        "alpha.dot_min=0.94",]
        ca_label = "conservative"
      print("  running find_ss_from_ca %s..." % ca_label, file=self.log)
      fss = find_ss_from_ca.find_secondary_structure(
          hierarchy=pdb_hierarchy,
          args = from_ca_args,
          ss_by_chain=self.params.secondary_structure.ss_by_chain,
          max_rmsd=self.params.secondary_structure.max_rmsd,
          use_representative_chains=\
           self.params.secondary_structure.use_representative_chains,
          max_representative_chains=\
           self.params.secondary_structure.max_representative_chains,
          out=null_out())
      return fss.get_annotation()
    elif self.params.secondary_structure.protein.search_method == "cablam":
      from mmtbx.validation import cablam
      print("  running cablam...", file=self.log)
      cablam_results = cablam.cablamalyze(
          pdb_hierarchy = pdb_hierarchy,
          outliers_only=False,
          out=null_out(),
          quiet=False)
      return cablam_results.as_secondary_structure()
    else:
      print("  WARNING: Unknown search method for SS. No SS found.", file=self.log)
      return iotbx.pdb.secondary_structure.annotation.from_records()


  def find_approximate_helices(self, log=sys.stdout):
    print("  Looking for approximately helical regions. . .", file=log)
    print("    warning: experimental, results not guaranteed to work!", file=log)
    find_helices = proteins.find_helices_simple(self.pdb_hierarchy)
    find_helices.show(out=log)
    restraint_groups = find_helices.as_restraint_groups()
    return restraint_groups

  def find_base_pairs(self, log=sys.stdout, use_segid=False, segids=None):
    if not use_segid:
      base_pairs = ""
      stacking_pairs = ""
      if self.pdb_hierarchy.contains_nucleic_acid(min_content=0.01):
        stacking_pairs = nucleic_acids.get_phil_stacking_pairs(
          pdb_hierarchy=self.pdb_hierarchy,
          prefix="secondary_structure.nucleic_acid",
          params=self.params.secondary_structure.nucleic_acid,
          log=log)
        if self.grm is not None:
          base_pairs = nucleic_acids.get_phil_base_pairs(
            pdb_hierarchy=self.pdb_hierarchy,
            nonbonded_proxies=self.grm.pair_proxies(
                sites_cart=self.pdb_hierarchy.atoms().extract_xyz()).\
                    nonbonded_proxies,
            prefix="secondary_structure.nucleic_acid",
            params=self.params.secondary_structure.nucleic_acid,
            log=log,
            verbose = self.verbose)
      return "\n".join([base_pairs, stacking_pairs])
    else:
      assert len(segids) > 1
      annotations = []
      for segid in segids:
        base_pairs = ""
        stacking_pairs = ""
        isel = self.selection_cache.selection("segid '%s'" % segid).iselection()
        selected_pdb_h = self.pdb_hierarchy.select(isel)
        if selected_pdb_h.contains_nucleic_acid(min_content=0.01):
          stacking_pairs = nucleic_acids.get_phil_stacking_pairs(
            pdb_hierarchy=selected_pdb_h,
            prefix="secondary_structure.nucleic_acid",
            log=log,
            add_segid=segid)
          if self.grm is not None:
            selected_grm = self.grm.select(iselection=isel)
            base_pairs = nucleic_acids.get_phil_base_pairs(
              pdb_hierarchy=selected_pdb_h,
              nonbonded_proxies=selected_grm.pair_proxies(
                  sites_cart=selected_pdb_h.atoms().extract_xyz()).\
                      nonbonded_proxies,
              prefix="secondary_structure.nucleic_acid",
              log=log,
              add_segid=segid,
              verbose = self.verbose)
          if (base_pairs is not None):
            annotations.append(base_pairs)
          if (stacking_pairs is not None):
            annotations.append(stacking_pairs)
      return "\n".join(annotations)

  def apply_phil_str(self,
      phil_string,
      phil_params=None,
      annot=None,
      log=None,
      verbose=False):
    assert [phil_string, phil_params].count(None) == 1
    if log is None:
      log = self.log
    new_ss_params = None
    if phil_string is not None:
      ss_phil = sec_str_master_phil.fetch(source=iotbx.phil.parse(phil_string))
      if verbose :
        ss_phil.show(out=log, prefix="    ")
      new_ss_params = ss_phil.extract()
    else:
      new_ss_params = phil_params
    self.params.secondary_structure.protein.helix = \
        new_ss_params.secondary_structure.protein.helix
    self.params.secondary_structure.protein.sheet = \
        new_ss_params.secondary_structure.protein.sheet
    self.actual_sec_str = annot
    if annot is None:
      self.actual_sec_str = iotbx.pdb.secondary_structure.annotation.from_phil(
          phil_helices=new_ss_params.secondary_structure.protein.helix,
          phil_sheets=new_ss_params.secondary_structure.protein.sheet,
          pdb_hierarchy=self.pdb_hierarchy,
          atom_selection_cache=self.selection_cache)

  def create_protein_hbond_proxies(self,
                            annotation=None,
                            log=sys.stdout):
    # assert as_regular_bond_proxies=True
    if annotation is None:
      annotation = self.actual_sec_str
    remove_outliers = self.params.secondary_structure.protein.remove_outliers

    from scitbx.array_family import flex
    atoms = self.pdb_hierarchy.atoms()
    hbond_counts = flex.int(atoms.size(), 0)

    distance_ideal = self.params.secondary_structure.protein.distance_ideal_n_o
    distance_cut = self.params.secondary_structure.protein.distance_cut_n_o
    if (distance_cut is None):
      distance_cut = -1
    generated_proxies = geometry_restraints.shared_bond_simple_proxy()
    hb_angle_proxies = []
    if self.params.secondary_structure.protein.enabled:
      for helix in self.params.secondary_structure.protein.helix :
        if helix.selection is not None:
          print("    Processing helix ", helix.selection, file=log)
          proxies, angle_proxies = proteins.create_helix_hydrogen_bond_proxies(
              params=helix,
              pdb_hierarchy=self.pdb_hierarchy,
              selection_cache=self.selection_cache,
              weight=1.0,
              hbond_counts=hbond_counts,
              distance_ideal=distance_ideal,
              distance_cut=distance_cut,
              remove_outliers=remove_outliers,
              restrain_hbond_angles=self.params.secondary_structure.protein.restrain_hbond_angles,
              log=log)
          if (proxies.size() == 0):
            print("      No H-bonds generated for '%s'" % helix.selection, file=log)
            continue
          else:
            generated_proxies.extend(proxies)
            hb_angle_proxies += angle_proxies
      for k, sheet in enumerate(self.params.secondary_structure.protein.sheet):
        print("    Processing sheet with id=%s, first strand: %s" % (
            sheet.sheet_id, sheet.first_strand), file=log)
        if sheet.first_strand is not None:
          proxies, angle_proxies = proteins.create_sheet_hydrogen_bond_proxies(
            sheet_params=sheet,
            pdb_hierarchy=self.pdb_hierarchy,
            selection_cache=self.selection_cache,
            weight=1.0,
            hbond_counts=hbond_counts,
            distance_ideal=distance_ideal,
            distance_cut=distance_cut,
            remove_outliers=remove_outliers,
            restrain_hbond_angles=self.params.secondary_structure.protein.restrain_hbond_angles,
            log=log)
          if (proxies.size() == 0):
            print("  No H-bonds generated for sheet with id=%s" % sheet.sheet_id, file=log)
            continue
          else:
            generated_proxies.extend(proxies)
            hb_angle_proxies += angle_proxies

    n_proxies = generated_proxies.size()
    print("", file=log)
    if (n_proxies == 0):
      print("    No hydrogen bonds defined for protein.", file=log)
    else :
      print("    %d hydrogen bonds defined for protein." % n_proxies, file=log)
      print("    %d hydrogen bond angles defined for protein." % len(hb_angle_proxies), file=log)
    return generated_proxies, geometry_restraints.shared_angle_proxy(hb_angle_proxies)

  def create_all_new_restraints(self,
      pdb_hierarchy,
      grm,
      log=sys.stdout):

    # initialize cache and monomer library server for underlying procedures
    if self.mon_lib_srv is None:
      from mmtbx.monomer_library import server
      self.mon_lib_srv = server.server()
    plane_cache = {}

    t0 = time.time()
    proteins_hbonds, prot_angle_proxies = self.create_protein_hbond_proxies(
        log=log,
        annotation=self.actual_sec_str)
    t1 = time.time()
    # print >> log, "    Time for creating protein proxies:%f" % (t1-t0)
    stacking_proxies = nucleic_acids.get_stacking_proxies(
        pdb_hierarchy=pdb_hierarchy,
        stacking_phil_params=self.params.secondary_structure.\
            nucleic_acid.stacking_pair,
        grm=grm,
        mon_lib_srv=self.mon_lib_srv,
        plane_cache=plane_cache)
    t2 = time.time()
    # print >> log, "    Time for creating stacking proxies:%f" % (t2-t1)
    (hb_bond_proxies, hb_angle_proxies, planarity_bp_proxies,
      parallelity_bp_proxies) = nucleic_acids.get_basepair_proxies(
          pdb_hierarchy=pdb_hierarchy,
          bp_phil_params=self.params.secondary_structure.nucleic_acid.base_pair,
          grm=grm,
          mon_lib_srv=self.mon_lib_srv,
          plane_cache=plane_cache,
          hbond_distance_cutoff=self.params.secondary_structure.\
            nucleic_acid.hbond_distance_cutoff,
          scale_bonds_sigma=self.params.secondary_structure.\
            nucleic_acid.scale_bonds_sigma)
    t4 = time.time()
    # print >> log, "    Time for creating basepair proxies (hbond, angle, planarity):%f" % (t4-t2)
    self.stats = {'n_protein_hbonds':0, 'n_na_hbonds':0, 'n_na_hbond_angles':0,
        'n_na_basepairs':0, 'n_na_stacking_pairs':0}
    print("    Restraints generated for nucleic acids:", file=log)
    print("      %d hydrogen bonds" % len(hb_bond_proxies), file=log)
    print("      %d hydrogen bond angles" % len(hb_angle_proxies), file=log)
    print("      %d basepair planarities" % len(planarity_bp_proxies), file=log)
    print("      %d basepair parallelities" % len(parallelity_bp_proxies), file=log)
    print("      %d stacking parallelities" % len(stacking_proxies), file=log)
    all_hbonds = proteins_hbonds.deep_copy()
    all_hbonds.extend(hb_bond_proxies)
    all_angle = prot_angle_proxies.deep_copy()
    all_angle.extend(hb_angle_proxies)
    return (all_hbonds,
        all_angle,
        planarity_bp_proxies, parallelity_bp_proxies+stacking_proxies)

  def get_simple_bonds(self, selection_phil=None):
    # assert 0 # used in GUI to draw bonds
    # this function wants
    # shared.stl_set_unsigned([(i_seq, j_seq),(i_seq, j_seq)])
    desired_annotation = None
    if (selection_phil is not None):
      if isinstance(selection_phil, str):
        selection_phil = iotbx.phil.parse(selection_phil)
      params = sec_str_master_phil.fetch(source=selection_phil).extract()
      desired_annotation = iotbx.pdb.secondary_structure.annotation.from_phil(
          phil_helices=params.secondary_structure.protein.helix,
          phil_sheets=params.secondary_structure.protein.sheet,
          pdb_hierarchy=self.pdb_hierarchy)
    else :
      desired_annotation = self.actual_sec_str
    hb_proxies, hb_angles = self.create_protein_hbond_proxies(
        annotation=desired_annotation,
        log=sys.stdout)
    from scitbx.array_family import shared
    simple_bonds = []
    for p in hb_proxies:
      simple_bonds.append(p.i_seqs)
    return shared.stl_set_unsigned(simple_bonds)

  def calculate_structure_content(self):
    if self.actual_sec_str is not None:
      oc = self.pdb_hierarchy.overall_counts()
      n_alpha = self.actual_sec_str.get_n_helix_residues()
      n_beta = self.actual_sec_str.get_n_sheet_residues()
      n_residues = oc.get_n_residues_of_classes(
          classes=['common_amino_acid', 'modified_amino_acid'])
    else:
      raise NotImplementedError("Should have defined secondary structure before calculating its content.")
    if n_residues == 0 :
      return (0.0, 0.0)
    return (n_alpha / n_residues, n_beta / n_residues)

  def show_summary(self, log=None):
    if log is None:
      log = self.log
    (frac_alpha, frac_beta) = self.calculate_structure_content()
    n_helices = len(self.params.secondary_structure.protein.helix)
    if n_helices == 1:
      if self.params.secondary_structure.protein.helix[0].selection is None:
        n_helices = 0
    n_sheets  = len(self.params.secondary_structure.protein.sheet)
    if n_sheets == 1:
      if self.params.secondary_structure.protein.sheet[0].first_strand is None:
        n_sheets = 0
    n_base_pairs = len(self.params.secondary_structure.nucleic_acid.base_pair)
    if n_base_pairs == 1:
      if self.params.secondary_structure.nucleic_acid.base_pair[0].base1 is None:
        n_base_pairs = 0
    n_stacking_pairs = len(self.params.secondary_structure.nucleic_acid.stacking_pair)
    if n_stacking_pairs == 1:
      if self.params.secondary_structure.nucleic_acid.stacking_pair[0].base1 is None:
        n_stacking_pairs = 0
    print("    Secondary structure from input PDB file:", file=log)
    print("      %d helices and %d sheets defined" % (n_helices,n_sheets), file=log)
    print("      %.1f%% alpha, %.1f%% beta" %(frac_alpha*100,frac_beta*100), file=log)
    print("      %d base pairs and %d stacking pairs defined." % (n_base_pairs,n_stacking_pairs), file=log)

  def helix_selections(self, limit=None, main_conf_only=False,
      alpha_only=False):
    sele = self.selection_cache.selection
    all_selections = []
    for helix in self.params.secondary_structure.protein.helix :
      if (helix.selection is not None):
        if (alpha_only) and (helix.helix_type != "alpha"):
          continue
        clauses = [ "(%s)" % helix.selection ]
        if (limit is not None):
          assert isinstance(limit, str)
          clauses.append("(%s)" % limit)
        if main_conf_only :
          clauses.append("(altloc ' ' or altloc 'A')")
        helix_sel = sele(" and ".join(clauses))
        all_selections.append(helix_sel)
    return all_selections

  def get_helix_types(self):
    return [ helix.helix_type for helix in self.params.secondary_structure.protein.helix ]

  def helix_selection(self, **kwds):
    whole_selection = flex.bool(self.n_atoms)
    for helix in self.helix_selections(**kwds):
      whole_selection |= helix
    return whole_selection

  def beta_selections(self, limit=None, main_conf_only=False):
    sele = self.selection_cache.selection
    all_selections = []
    for sheet in self.params.secondary_structure.protein.sheet :
      sheet_selection = flex.bool(self.n_atoms)
      clauses = []
      if (limit is not None):
        assert isinstance(limit, str)
        clauses.append("(%s)" % limit)
      if main_conf_only :
        clauses.append("(altloc ' ' or altloc 'A')")
      main_clause = [ "(%s)" % sheet.first_strand ]
      strand_sel = sele(" and ".join(main_clause+clauses))
      sheet_selection |= strand_sel
      for strand in sheet.strand :
        main_clause = [ "(%s)" % strand.selection ]
        strand_sel = sele(" and ".join(main_clause+clauses))
        sheet_selection |= strand_sel
      all_selections.append(sheet_selection)
    return all_selections

  def beta_selection(self, **kwds):
    whole_selection = flex.bool(self.n_atoms)
    for sheet in self.beta_selections(**kwds):
      whole_selection |= sheet
    return whole_selection

  def base_pair_selections(self, limit=None, main_conf_only=False):
    # assert 0, "Probably shouldn't be used"
    sele = self.selection_cache.selection
    all_selections = []
    for bp in self.params.secondary_structure.nucleic_acid.base_pair :
      if (bp.base1 is not None) and (bp.base2 is not None):
        clauses = [ "((%s) or (%s))" % (bp.base1, bp.base2) ]
        if (limit is not None):
          clauses.append("(%s)" % limit)
        if main_conf_only :
          clauses.append("(altloc ' ' or altloc 'A')")
        bp_sele = sele(" and ".join(clauses))
        all_selections.append(bp_sele)
    return all_selections

  def base_pair_selection(self, **kwds):
    # assert 0, "Probably shouldn't be used",
    # used in mmtbx/command_line/find_tls_groups.py
    whole_selection = flex.bool(self.n_atoms)
    for sheet in self.base_pair_selections(**kwds):
      whole_selection |= sheet
    return whole_selection

# =============================================================================
# General functions
# =============================================================================

def sec_str_from_phil(phil_str):
  # Used in wxGUI2/Programs/SecondaryStructure.py line 564, in update_from_strand_list
  ss_phil = iotbx.phil.parse(phil_str)
  ss_phil_fetched = sec_str_master_phil.fetch(source=ss_phil)
  # ss_phil_fetched.show()
  return ss_phil_fetched.extract()

def get_ksdssp_exe_path():
  if (not libtbx.env.has_module(name="ksdssp")):
    raise RuntimeError("ksdssp module is not configured")
  exe_path = libtbx.env.under_build("ksdssp/exe/ksdssp")
  if (os.name == "nt"):
    exe_path += ".exe"
  if (not os.path.isfile(exe_path)):
    raise RuntimeError("ksdssp executable is not available")
  return exe_path

def run_ksdssp(file_name, log=sys.stdout):
  if not os.path.isfile(file_name):
    raise RuntimeError("File %s not found.")
  exe_path = get_ksdssp_exe_path()
  print("  Running KSDSSP to generate HELIX and SHEET records", file=log)
  ksdssp_out = easy_run.fully_buffered(command="%s %s" % (exe_path, file_name))
#  if len(ksdssp_out.stderr_lines) > 0 :
#    print >> log, "\n".join(ksdssp_out.stderr_lines)
  return ksdssp_out.stdout_lines

def run_ksdssp_direct(pdb_str):
  exe_path = get_ksdssp_exe_path()
  ksdssp_out = easy_run.fully_buffered(command=exe_path, stdin_lines=pdb_str)
  return ( ksdssp_out.stdout_lines, ksdssp_out.stderr_lines )
