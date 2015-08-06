
from __future__ import division
from mmtbx.secondary_structure import proteins, nucleic_acids
from cctbx import geometry_restraints
import iotbx.pdb
import iotbx.pdb.secondary_structure
import iotbx.phil
from scitbx.array_family import flex
from libtbx.utils import null_out
from libtbx import easy_run
import libtbx.load_env
from math import sqrt
import sys, os
import time

#
# Notes. Now this class trying to keep updated SS information in form of
# phil parameters and as iotbx/secondary_structure/annotation. This is
# not the optimal way, it would be better to work only with annotation, but
# this requires massive rewriting of SS GUI which works only with phil...
#

sec_str_master_phil_str = """
secondary_structure
  .style = box auto_align hidden
{
  enabled = False
    .short_caption = Use secondary structure restraints (main switch)
    .type = bool
    .style = noauto bold
    .help = Turn on secondary structure restraints
  protein
    .style = box auto_align noauto
  {
    enabled = True
      .type = bool
      .help = Turn on secondary structure restraints for protein
    search_method = *ksdssp mmtbx_dssp from_ca
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
      .short_caption = Filter bond outliers
      .style = tribool
      .help = If true, h-bonds exceeding distance_cut_n_o length will not be \
       established
    %s
    %s
  }
  nucleic_acid
    .caption = If sigma and slack are not defined for nucleic acids, the \
      overall default settings for protein hydrogen bonds will be used \
      instead.
    .style = box auto_align noauto
  {
    enabled = True
      .type = bool
      .help = Turn on secondary structure restraints for nucleic acids
    %s
  }
}
""" % (
       proteins.helix_group_params_str,
       proteins.sheet_group_params_str,
       nucleic_acids.dna_rna_params_str)

master_phil_str = sec_str_master_phil_str # for docs

sec_str_master_phil = iotbx.phil.parse(sec_str_master_phil_str)
default_params = sec_str_master_phil.fetch().extract()

class manager (object) :
  def __init__ (self,
                pdb_hierarchy,
                geometry_restraints_manager=None,
                sec_str_from_pdb_file=None,
                params=None,
                mon_lib_srv=None,
                verbose=-1,
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

    self.verbose = verbose
    self.log = log
    self.def_params = sec_str_master_phil.extract()
    self.stats = {'n_protein_hbonds':0, 'n_na_hbonds':0, 'n_na_hbond_angles':0,
        'n_na_basepairs':0, 'n_na_stacking_pairs':0}

    atoms = pdb_hierarchy.atoms()
    i_seqs = atoms.extract_i_seq()
    if (i_seqs.all_eq(0)) :
      atoms.reset_i_seq()
      i_seqs = atoms.extract_i_seq()
    self.n_atoms = atoms.size()
    self._was_initialized = False
    self.selection_cache = pdb_hierarchy.atom_selection_cache()
    self.pdb_atoms = atoms
    self.initialize()

  def as_phil_str (self, master_phil=sec_str_master_phil) :
    # used in secondary_structure_restraints...
    return master_phil.format(python_object=self.params)

  def initialize(self):
    if not self._was_initialized :
      self.find_automatically()
      self.show_summary()
      self._was_initialized = True

  def find_automatically(self):
    # find_automatically = self.params.secondary_structure.find_automatically
    protein_found = False
    atom_labels = list(self.pdb_hierarchy.atoms_with_labels())
    segids = {}
    for a in atom_labels:
      segids[a.segid] = ""
    use_segid = len(segids) > 1
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
    if (not protein_ss_definition_present and
        self.sec_str_from_pdb_file is None):
      if(self.verbose>0):
        print >> log, "No existing protein secondary structure definitions " + \
        "found in .pdb file or phil parameters."
      # find_automatically = True
    else:
      # find_automatically = False
      protein_found = True
    if not protein_found:
      ss_params = []
      if use_segid:
        # Could get rid of this 'if' clause, but I want to avoid construction of
        # atom_selection_cache and selections when there is no segids in pdb
        # which is majority of cases.
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
      else:
        if self.pdb_hierarchy.contains_protein():
          annot = self.find_sec_str(pdb_hierarchy=self.pdb_hierarchy)
          if annot is not None:
            ss_phil = annot.as_restraint_groups(log=self.log,
              prefix_scope="secondary_structure")
            ss_params.append(ss_phil)
            # self.actual_sec_str = annot
      ss_params_str = "\n".join(ss_params)
      self.apply_phil_str(ss_params_str, log=self.log)
    else:
      if (self.sec_str_from_pdb_file is not None):
        # self.actual_sec_str = self.sec_str_from_pdb_file
        ss_params_str = self.sec_str_from_pdb_file.as_restraint_groups(
            log=self.log,
            prefix_scope="secondary_structure")
        self.apply_phil_str(ss_params_str, log=self.log)
      else:
        # got phil SS, need to refactor later, when the class is fully
        # converted for operation with annotation object
        # phil_string = sec_str_master_phil.format(python_object=self.params)
        self.apply_phil_str(phil_string=None, phil_params=self.params, log=self.log)

    t1 = time.time()
    # print >> log, "    Time for finding protein SS: %f" % (t1-t0)
    # Step 2: nucleic acids
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
    if self.params.secondary_structure.protein.search_method == "ksdssp":
      pdb_str = pdb_hierarchy.as_pdb_string()
      print >> self.log, "  running ksdssp..."
      (records, stderr) = run_ksdssp_direct(pdb_str)
      return iotbx.pdb.secondary_structure.annotation.from_records(
          records=records,
          log=self.log)
    elif self.params.secondary_structure.protein.search_method == "mmtbx_dssp":
      from mmtbx.secondary_structure import dssp
      print >> self.log, "  running mmtbx.dssp..."
      return dssp.dssp(
        pdb_hierarchy=pdb_hierarchy,
        pdb_atoms=self.pdb_atoms,
        out=null_out()).get_annotation()
    elif self.params.secondary_structure.protein.search_method == "from_ca":
      from mmtbx.secondary_structure import find_ss_from_ca
      print >> self.log, "  running find_ss_from_ca..."
      fss = find_ss_from_ca.find_secondary_structure(
          hierarchy=pdb_hierarchy,
          out=null_out())
      return fss.get_pdb_annotation()
    else:
      print >> self.log, "  WARNING: Unknown search method for SS. No SS found."
      return iotbx.pdb.secondary_structure.annotation.from_records()


  def find_approximate_helices (self, log=sys.stdout) :
    print >> log, "  Looking for approximately helical regions. . ."
    print >> log, "    warning: experimental, results not guaranteed to work!"
    find_helices = proteins.find_helices_simple(self.pdb_hierarchy)
    find_helices.show(out=log)
    restraint_groups = find_helices.as_restraint_groups()
    return restraint_groups

  def find_base_pairs (self, log=sys.stdout, use_segid=False, segids=None):
    if not use_segid:
      base_pairs = ""
      stacking_pairs = ""
      if self.pdb_hierarchy.contains_nucleic_acid():
        stacking_pairs = nucleic_acids.get_phil_stacking_pairs(
          pdb_hierarchy=self.pdb_hierarchy,
          prefix="secondary_structure.nucleic_acid",
          params=self.params.secondary_structure.nucleic_acid,
          log=log)
        if self.grm is not None:
          base_pairs = nucleic_acids.get_phil_base_pairs(
            pdb_hierarchy=self.pdb_hierarchy,
            nonbonded_proxies=self.grm.pair_proxies().nonbonded_proxies,
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
        if selected_pdb_h.contains_nucleic_acid():
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
          if (base_pairs is not None) :
            annotations.append(base_pairs)
          if (stacking_pairs is not None) :
            annotations.append(stacking_pairs)
      return "\n".join(annotations)

  def apply_phil_str(self,
      phil_string,
      phil_params=None,
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
    self.actual_sec_str = iotbx.pdb.secondary_structure.annotation.from_phil(
        phil_helices=new_ss_params.secondary_structure.protein.helix,
        phil_sheets=new_ss_params.secondary_structure.protein.sheet,
        pdb_hierarchy=self.pdb_hierarchy)

  def create_protein_hbond_proxies (self,
                            annotation=None,
                            log=sys.stdout):
    # assert as_regular_bond_proxies=True
    if annotation is None:
      annotation = self.actual_sec_str
    remove_outliers = self.params.secondary_structure.protein.remove_outliers

    from scitbx.array_family import flex
    atoms = self.pdb_hierarchy.atoms()
    hbond_counts = flex.int(atoms.size(), 0)
    selection_cache = self.pdb_hierarchy.atom_selection_cache()

    distance_ideal = self.params.secondary_structure.protein.distance_ideal_n_o
    distance_cut = self.params.secondary_structure.protein.distance_cut_n_o
    if (distance_cut is None) :
      distance_cut = -1
    generated_proxies = geometry_restraints.shared_bond_simple_proxy()
    if self.params.secondary_structure.protein.enabled:
      for helix in self.params.secondary_structure.protein.helix :
        if helix.selection is not None:
          print >> log, "    Processing helix ", helix.selection
          proxies = proteins.create_helix_hydrogen_bond_proxies(
              params=helix,
              pdb_hierarchy=self.pdb_hierarchy,
              selection_cache=selection_cache,
              weight=1.0,
              hbond_counts=hbond_counts,
              distance_ideal=distance_ideal,
              distance_cut=distance_cut,
              remove_outliers=remove_outliers,
              log=log)
          if (proxies.size() == 0) :
            print >> log, "      No H-bonds generated for '%s'" % helix.selection
            continue
          else:
            generated_proxies.extend(proxies)
      for k, sheet in enumerate(self.params.secondary_structure.protein.sheet) :
        print >> log, "    Processing sheet with id=%s, first strand: %s" % (
            sheet.sheet_id, sheet.first_strand)
        if sheet.first_strand is not None:
          proxies = proteins.create_sheet_hydrogen_bond_proxies(
            sheet_params=sheet,
            pdb_hierarchy=self.pdb_hierarchy,
            weight=1.0,
            hbond_counts=hbond_counts,
            distance_ideal=distance_ideal,
            distance_cut=distance_cut,
            remove_outliers=remove_outliers,
            log=log)
          if (proxies.size() == 0) :
            print >> log, \
                "  No H-bonds generated for sheet with id=%s" % sheet.sheet_id
            continue
          else:
            generated_proxies.extend(proxies)

    n_proxies = generated_proxies.size()
    print >> log, ""
    if (n_proxies == 0) :
      print >> log, "  No hydrogen bonds defined for protein."
    else :
      print >> log, "  %d hydrogen bonds defined for protein." % n_proxies
    # reg_proxies = []
    # for hb_p in build_proxies.proxies:
    #   reg_proxy = geometry_restraints.bond_simple_proxy(
    #     i_seqs=hb_p.i_seqs,
    #     distance_ideal=hb_p.distance_ideal,
    #     weight=hb_p.weight,
    #     slack=hb_p.slack,
    #     top_out=hb_p.top_out,
    #     limit=hb_p.limit,
    #     origin_id=1)
    #   reg_proxies.append(reg_proxy)
    # return geometry_restraints.shared_bond_simple_proxy(reg_proxies)
    return generated_proxies

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
    proteins_hbonds = self.create_protein_hbond_proxies(
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
    planarity_bp_proxies, parallelity_bp_proxies = nucleic_acids.\
        get_basepair_plane_proxies(
        pdb_hierarchy=pdb_hierarchy,
        bp_phil_params=self.params.secondary_structure.nucleic_acid.base_pair,
        grm=grm,
        mon_lib_srv=self.mon_lib_srv,
        plane_cache=plane_cache)
    t3 = time.time()
    # print >> log, "    Time for creating planar/parall proxies:%f" % (t3-t2)
    hb_bond_proxies, hb_angle_proxies = nucleic_acids.\
        get_basepair_hbond_proxies(
        pdb_hierarchy=pdb_hierarchy,
        bp_phil_params=self.params.secondary_structure.nucleic_acid.base_pair,
        hbond_distance_cutoff=self.params.secondary_structure.\
            nucleic_acid.hbond_distance_cutoff)
    t4 = time.time()
    # print >> log, "    Time for creating hbond-angle proxies:%f" % (t4-t3)
    self.stats = {'n_protein_hbonds':0, 'n_na_hbonds':0, 'n_na_hbond_angles':0,
        'n_na_basepairs':0, 'n_na_stacking_pairs':0}
    print >> log, "  Restraints generated for nucleic acids:"
    print >> log, "    %d hydrogen bonds" % len(hb_bond_proxies)
    print >> log, "    %d hydrogen bond angles" % len(hb_angle_proxies)
    print >> log, "    %d basepair planarities" % len(planarity_bp_proxies)
    print >> log, "    %d basepair parallelities" % len(parallelity_bp_proxies)
    print >> log, "    %d stacking parallelities" % len(stacking_proxies)
    all_hbonds = proteins_hbonds.deep_copy()
    all_hbonds.extend(hb_bond_proxies)
    return (all_hbonds, hb_angle_proxies,
        planarity_bp_proxies, parallelity_bp_proxies+stacking_proxies)

  def get_simple_bonds (self, selection_phil=None) :
    # assert 0 # used in GUI to draw bonds
    # this function wants
    # shared.stl_set_unsigned([(i_seq, j_seq),(i_seq, j_seq)])
    desired_annotation = None
    if (selection_phil is not None) :
      if isinstance(selection_phil, str) :
        selection_phil = iotbx.phil.parse(selection_phil)
      params = sec_str_master_phil.fetch(source=selection_phil).extract()
      desired_annotation = iotbx.pdb.secondary_structure.annotation.from_phil(
          phil_helices=params.secondary_structure.protein.helix,
          phil_sheets=params.secondary_structure.protein.sheet,
          pdb_hierarchy=self.pdb_hierarchy)
    else :
      desired_annotation = self.actual_sec_str
    hb_proxies = self.create_protein_hbond_proxies(
        annotation=desired_annotation,
        log=sys.stdout)
    from scitbx.array_family import shared
    simple_bonds = []
    for p in hb_proxies:
      simple_bonds.append(p.i_seqs)
    return shared.stl_set_unsigned(simple_bonds)

  def calculate_structure_content (self) :
    isel = self.selection_cache.iselection
    calpha = isel("name N and (altloc ' ' or altloc 'A')")
    alpha_sele = self.alpha_selection(limit="name N", main_conf_only=True)
    n_alpha = alpha_sele.count(True)
    beta_sele = self.beta_selection(limit="name N", main_conf_only=True)
    n_beta = beta_sele.count(True)
    if calpha.size() == 0 :
      return (0.0, 0.0)
    return (n_alpha / calpha.size(), n_beta / calpha.size())

  def show_summary (self, log=None) :
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
    print >> log, "Secondary structure from input PDB file:"
    print >> log, "  %d helices and %d sheets defined" % (n_helices,n_sheets)
    print >> log, "  %.1f%% alpha, %.1f%% beta" %(frac_alpha*100,frac_beta*100)
    print >> log, "  %d base pairs and %d stacking pairs defined." % (n_base_pairs,n_stacking_pairs)

  def helix_selections (self, limit=None, main_conf_only=False,
      alpha_only=False) :
    sele = self.selection_cache.selection
    all_selections = []
    for helix in self.params.secondary_structure.protein.helix :
      if (helix.selection is not None) :
        if (alpha_only) and (helix.helix_type != "alpha") :
          continue
        clauses = [ "(%s)" % helix.selection ]
        if (limit is not None) :
          assert isinstance(limit, str)
          clauses.append("(%s)" % limit)
        if main_conf_only :
          clauses.append("(altloc ' ' or altloc 'A')")
        helix_sel = sele(" and ".join(clauses))
        all_selections.append(helix_sel)
    return all_selections

  def get_helix_types (self) :
    return [ helix.helix_type for helix in self.params.secondary_structure.protein.helix ]

  def helix_selection (self, **kwds) :
    whole_selection = flex.bool(self.n_atoms)
    for helix in self.helix_selections(**kwds) :
      whole_selection |= helix
    return whole_selection

  # FIXME backwards compatibility
  def alpha_selection (self, **kwds) :
    return self.helix_selection(**kwds)

  def alpha_selections (self, **kwds) :
    return self.helix_selections(**kwds)

  def beta_selections (self, limit=None, main_conf_only=False) :
    sele = self.selection_cache.selection
    all_selections = []
    for sheet in self.params.secondary_structure.protein.sheet :
      sheet_selection = flex.bool(self.n_atoms)
      clauses = []
      if (limit is not None) :
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

  def beta_selection (self, **kwds) :
    whole_selection = flex.bool(self.n_atoms)
    for sheet in self.beta_selections(**kwds) :
      whole_selection |= sheet
    return whole_selection

  def base_pair_selections (self, limit=None, main_conf_only=False) :
    # assert 0, "Probably shouldn't be used"
    sele = self.selection_cache.selection
    all_selections = []
    for bp in self.params.secondary_structure.nucleic_acid.base_pair :
      if (bp.base1 is not None) and (bp.base2 is not None) :
        clauses = [ "((%s) or (%s))" % (bp.base1, bp.base2) ]
        if (limit is not None) :
          clauses.append("(%s)" % limit)
        if main_conf_only :
          clauses.append("(altloc ' ' or altloc 'A')")
        bp_sele = sele(" and ".join(clauses))
        all_selections.append(bp_sele)
    return all_selections

  def base_pair_selection (self, **kwds) :
    # assert 0, "Probably shouldn't be used",
    # used in mmtbx/command_line/find_tls_groups.py
    whole_selection = flex.bool(self.n_atoms)
    for sheet in self.base_pair_selections(**kwds) :
      whole_selection |= sheet
    return whole_selection

  def selections_as_ints (self) :
    assert 0, "Anybody using this?"
    sec_str = flex.int(self.n_atoms, 0)
    all_alpha = flex.int(self.n_atoms, 1)
    all_beta = flex.int(self.n_atoms, 2)
    helices = self.alpha_selection()
    sheets = self.beta_selection()
    sec_str.set_selected(helices, all_alpha.select(helices))
    sec_str.set_selected(sheets, all_beta.select(sheets))
    return sec_str

  def apply_params (self, params) :
    assert 0, "Not used anywhere?"
    self.params.helix = params.helix
    self.params.sheet = params.sheet

# =============================================================================
# General functions
# =============================================================================

def sec_str_from_phil (phil_str) :
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

def run_ksdssp (file_name, log=sys.stdout) :
  if not os.path.isfile(file_name) :
    raise RuntimeError("File %s not found.")
  exe_path = get_ksdssp_exe_path()
  print >> log, "  Running KSDSSP to generate HELIX and SHEET records"
  ksdssp_out = easy_run.fully_buffered(command="%s %s" % (exe_path, file_name))
#  if len(ksdssp_out.stderr_lines) > 0 :
#    print >> log, "\n".join(ksdssp_out.stderr_lines)
  return ksdssp_out.stdout_lines

def run_ksdssp_direct(pdb_str) :
  exe_path = get_ksdssp_exe_path()
  ksdssp_out = easy_run.fully_buffered(command=exe_path, stdin_lines=pdb_str)
  return ( ksdssp_out.stdout_lines, ksdssp_out.stderr_lines )




# =============================================================================
# Unused functions...
# =============================================================================

def _get_distances (bonds, sites_cart) :
  assert 0, "Hopefully is not used"
  distances = flex.double(bonds.size(), -1)
  for k, (i_seq, j_seq) in enumerate(bonds) :
    (x1, y1, z1) = sites_cart[i_seq]
    (x2, y2, z2) = sites_cart[j_seq]
    dist = sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
    distances[k] = dist
  return distances

def get_pdb_hierarchy (file_names) :
  assert 0, "Hopefully is not used"
  pdb_combined = iotbx.pdb.combine_unique_pdb_files(file_names=file_names)
  pdb_structure = iotbx.pdb.input(source_info=None,
    lines=flex.std_string(pdb_combined.raw_records))
  return pdb_structure.construct_hierarchy()

def analyze_distances (self, params, pdb_hierarchy=None, log=sys.stdout) :
  assert 0 # Not used anywhere?
  atoms = None
  if params.verbose :
    assert pdb_hierarchy is not None
    atoms = pdb_hierarchy.atoms()
  remove_outliers = params.secondary_structure.protein.remove_outliers
  atoms = pdb_hierarchy.atoms()
  hist =  flex.histogram(self.bond_lengths, 10)
  print >> log, "  Distribution of hydrogen bond lengths without filtering:"
  hist.show(f=log, prefix="    ", format_cutoffs="%.4f")
  print >> log, ""
  if not remove_outliers :
    return False
  for i, distance in enumerate(self.bond_lengths) :
    if distance > distance_max :
      self.flag_use_bond[i] = False
      if params.verbose :
        print >> log, "Excluding H-bond with length %.3fA" % distance
        i_seq, j_seq = self.bonds[i]
        print >> log, "  %s" % atoms[i_seq].fetch_labels().id_str()
        print >> log, "  %s" % atoms[j_seq].fetch_labels().id_str()
  print >> log, "  After filtering: %d bonds remaining." % \
    self.flag_use_bond.count(True)
  print >> log, "  Distribution of hydrogen bond lengths after applying cutoff:"
  hist = flex.histogram(self.bond_lengths.select(self.flag_use_bond), 10)
  hist.show(f=log, prefix="    ", format_cutoffs="%.4f")
  print >> log, ""
  return True

def find_nucleic_acids (pdb_hierarchy) :
  assert 0, "Anybody using this?"
  for model in pdb_hierarchy.models() :
    for chain in model.chains() :
      for conformer in chain.conformers() :
        if conformer.is_na() :
          return True
  return False

def manager_from_pdb_file (pdb_file) :
  assert 0, "will not work"
  from iotbx import file_reader
  assert os.path.isfile(pdb_file)
  pdb_in = file_reader.any_file(pdb_file, force_type="pdb")
  pdb_hierarchy = pdb_in.file_object.hierarchy
  pdb_hierarchy.atoms().reset_i_seq()
  ss_manager  = manager(pdb_hierarchy=pdb_hierarchy)
  return ss_manager

def calculate_structure_content (pdb_file) :
  assert 0, "will not work"
  ss_manager = manager_from_pdb_file(pdb_file)
  ss_manager.find_automatically()
  return ss_manager.calculate_structure_content()
