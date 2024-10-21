from __future__ import absolute_import, division, print_function

from mmtbx.secondary_structure import manager
import iotbx.pdb
import iotbx.phil
from libtbx.utils import Sorry, safe_div
from six.moves import cStringIO as StringIO
import sys
from libtbx import easy_mp
import mmtbx
from mmtbx.building.loop_closure.utils import get_phi_psi_atoms, get_pair_angles
from libtbx import group_args

import boost_adaptbx.boost.python as bp
from six.moves import zip
ext = bp.import_ext("mmtbx_validation_ramachandran_ext")
from mmtbx_validation_ramachandran_ext import rama_eval
from mmtbx.validation import ramalyze
from mmtbx.validation.phi_psi_2_data import phi_psi_2_mask_class

master_phil_str = '''
ss_validation {
  nproc = 1
    .type = int
  bad_hbond_cutoff = 3.5
    .type = float
  mediocre_hbond_cutoff = 3.0
    .type = float
  filter_annotation = False
    .type = bool
    .help = Output filtered annotations
}
'''

class gather_ss_stats(object):
  def __init__(
        self,
        pdb_h,
        mediocre_hbond_cutoff=3.0,
        bad_hbond_cutoff=3.5,
        rama_eval_manager=None):
    self.pdb_h = pdb_h
    self.bad_hbond_cutoff = bad_hbond_cutoff
    self.mediocre_hbond_cutoff = mediocre_hbond_cutoff
    self.asc = self.pdb_h.atom_selection_cache()
    self.atoms = pdb_h.atoms()
    self.r = rama_eval_manager
    if self.r is None:
      self.r = rama_eval()
    self.pp2m = phi_psi_2_mask_class()

  def __call__(self, hsh_tuple):
    temp_annot = iotbx.pdb.secondary_structure.annotation(
        helices = hsh_tuple[0],
        sheets = hsh_tuple[1])
    helix = len(hsh_tuple[0]) > 0
    # print temp_annot.as_pdb_str().replace('\n',' '),
    ss_params = mmtbx.secondary_structure.manager.get_default_ss_params()
    ss_params.secondary_structure.enabled=True
    ss_params.secondary_structure.protein.remove_outliers=False
    ss_params.secondary_structure.protein.helix=[]
    ss_params.secondary_structure.protein.sheet=[]
    ss_params.secondary_structure.nucleic_acid.enabled=False
    ss_m_log = StringIO()

    ss_manager = mmtbx.secondary_structure.manager(
        pdb_hierarchy=self.pdb_h,
        atom_selection_cache=self.asc,
        sec_str_from_pdb_file=temp_annot,
        params=ss_params.secondary_structure,
        show_summary_on=False,
        log = ss_m_log)
    h_bond_proxies, hb_angles = ss_manager.create_protein_hbond_proxies(log=ss_m_log)

    cutoff_bad = self.bad_hbond_cutoff
    cutoff_mediocre = self.mediocre_hbond_cutoff
    n_hbonds = 0
    n_bad_hbonds = 0
    n_mediocre_hbonds = 0
    hb_lens = []
    for hb_p in h_bond_proxies:
      # print dir(hb_p)
      n_hbonds += 1
      hb_len = self.atoms[hb_p.i_seqs[0]].distance(self.atoms[hb_p.i_seqs[1]])
      hb_lens.append(hb_len)
      if hb_len > cutoff_bad:
        n_bad_hbonds += 1
      elif hb_len > cutoff_mediocre:
        n_mediocre_hbonds += 1

    # Ramachandran outliers and wrong areas
    sele = self.asc.selection(temp_annot.as_atom_selections()[0])
    ss_h = self.pdb_h.select(sele)
    phi_psi_atoms = get_phi_psi_atoms(ss_h)

    n_outliers = 0
    n_wrong_region = 0
    for phi_psi_pair, rama_key in phi_psi_atoms:
      phi, psi = get_pair_angles(phi_psi_pair)
      r_eval = self.r.evaluate_angles(rama_key, phi, psi)
      if r_eval == ramalyze.RAMALYZE_OUTLIER:
        n_outliers += 1
      else:
        reg = pp2_helix_sheet_rama_region(phi, psi, self.pp2m)
        if (helix and reg != 'A') or (not helix and reg !='B'):
          n_wrong_region += 1
          # print "  Wrong region:", phi_psi_pair[0][2].id_str(), reg, helix
    del ss_manager
    del ss_params
    return n_hbonds, n_bad_hbonds, n_mediocre_hbonds, hb_lens, n_outliers, n_wrong_region

def pp2_helix_sheet_rama_region(phi, psi, pp2m=None):
  # result: A - helix, B - sheet.
  if pp2m is None:
    pp2m = phi_psi_2_mask_class()
  return pp2m.get_closest((phi,psi))

def is_ca_and_something(pdb_h):
  asc = pdb_h.atom_selection_cache()
  n_N = asc.iselection("name N").size()
  n_O = asc.iselection("name O").size()
  n_CA = asc.iselection("name CA").size()
  water = asc.iselection("water").size()
  n_O -= water
  assert n_CA > 0
  if abs(n_CA-n_N)/n_CA > 0.5:
    return True
  if abs(n_CA-n_O)/n_CA > 0.5:
    return True
  return False

def some_chains_are_ca(pdb_h):
  for chain in pdb_h.only_model().chains():
    if chain.is_ca_only():
      return True

class validate(object):
  def __init__(self, model, params=None, log=sys.stdout):
    self.model = model
    self.params = params
    if self.params is None:
      self.params = validate.get_default_params().ss_validation
    self.log = log
    self.results = None
    ss_log = StringIO()
    try:
      ss_annot = self.model.get_ss_annotation(log=ss_log)
    except Sorry as e:
      print(" Syntax error in SS: %s" % str(e), file=self.log)
      return
    ss_log_cont = ss_log.getvalue()
    n_bad_helices = ss_log_cont.count("Bad HELIX")
    n_bad_sheets = ss_log_cont.count("Bad SHEET")
    pdb_h = self.model.get_hierarchy()
    if ss_annot is None or ss_annot.is_empty():
      print("No SS annotation, nothing to analyze", file=self.log)
      return
    if n_bad_helices > 0:
      print("Number of helices with syntax error: %d" % n_bad_helices, file=self.log)
    if n_bad_helices > 0:
      print("Number of sheets with syntax error: %d" % n_bad_sheets, file=self.log)
    if model.get_number_of_models() != 1 :
      raise Sorry("Multiple models not supported.")
    if not pdb_h.contains_protein():
      print("Protein is not found in the model", file=self.log)
      return
    if pdb_h.is_ca_only():
      print("Error: CA-only model", file=self.log)
      return
    if is_ca_and_something(pdb_h):
      print("CA-only and something model", file=self.log)
      return
    if some_chains_are_ca(pdb_h):
      print("some chains are CA-only", file=self.log)
      return

    n_total_helix_sheet_records = ss_annot.get_n_helices()+ss_annot.get_n_sheets()
    n_bad_helix_sheet_records = 0
    # Empty stuff:
    empty_annots = ss_annot.remove_empty_annotations(pdb_h)
    number_of_empty_helices = empty_annots.get_n_helices()
    number_of_empty_sheets = empty_annots.get_n_sheets()
    n_bad_helix_sheet_records += (number_of_empty_helices+number_of_empty_sheets)
    if number_of_empty_helices > 0:
      print("Helices without corresponding atoms in the model (%d):" % number_of_empty_helices, file=self.log)
      for h in empty_annots.helices:
        print("  ", h.as_pdb_str(), file=self.log)
    if number_of_empty_sheets > 0:
      print("Sheets without corresponding atoms in the model (%d):" % number_of_empty_sheets, file=self.log)
      for sh in empty_annots.sheets:
        print("  ", sh.as_pdb_str(), file=self.log)

    print("Checking annotations thoroughly, use nproc=<number> if it is too slow...", file=self.log)

    hsh_tuples = []
    for h in ss_annot.helices:
      hsh_tuples.append(([h],[]))
    for sh in ss_annot.sheets:
      hsh_tuples.append(([],[sh]))
    calc_ss_stats = gather_ss_stats(
        pdb_h,
        mediocre_hbond_cutoff=self.params.mediocre_hbond_cutoff,
        bad_hbond_cutoff=self.params.bad_hbond_cutoff)
    results = []
    if len(hsh_tuples) > 0:
      results = easy_mp.pool_map(
          processes=self.params.nproc,
          fixed_func=calc_ss_stats,
          args=hsh_tuples)

    cumm_n_hbonds = 0
    cumm_n_bad_hbonds = 0
    cumm_n_mediocre_hbonds = 0
    cumm_n_rama_out = 0
    cumm_n_wrong_reg = 0

    n_elem_with_wrong_rama = 0
    n_elem_with_rama_out = 0
    n_elem_with_bad_hbond = 0
    #
    # Hydrogen Bonds in Proteins: Role and Strength
    # Roderick E Hubbard, Muhammad Kamran Haider
    # ENCYCLOPEDIA OF LIFE SCIENCES & 2010, John Wiley & Sons, Ltd. www.els.net
    #
    # See also: http://proteopedia.org/wiki/index.php/Hydrogen_bonds
    #
    for ss_elem, r in zip(ss_annot.helices+ss_annot.sheets, results):
      if r is not None:
        n_hbonds, n_bad_hbonds, n_mediocre_hbonds, hb_lens, n_outliers, n_wrong_region = r
        cumm_n_hbonds += n_hbonds
        cumm_n_bad_hbonds += n_bad_hbonds
        cumm_n_mediocre_hbonds += n_mediocre_hbonds
        cumm_n_rama_out += n_outliers
        cumm_n_wrong_reg += n_wrong_region
        if n_wrong_region > 0:
          n_elem_with_wrong_rama += 1
        if n_outliers > 0:
          n_elem_with_rama_out += 1
        if n_bad_hbonds > 0:
          n_elem_with_bad_hbond += 1
        if n_bad_hbonds + n_outliers + n_wrong_region > 0:
          n_bad_helix_sheet_records += 1
        if n_bad_hbonds + n_mediocre_hbonds + n_outliers + n_wrong_region > 0:
          # this is bad annotation, printing it to log with separate stats:
          print("Bad annotation found:", file=self.log)
          print("%s" % ss_elem.as_pdb_str(), file=self.log)
          print("  Total hb: %d, mediocre: %d, bad: %d, Rama outliers: %d, Rama wrong %d" % (
              n_hbonds, n_mediocre_hbonds, n_bad_hbonds, n_outliers, n_wrong_region), file=self.log)
          print("-"*80, file=self.log)

    # n1 = percentage of bad SS elements (per given model);
    # bad here means: n_bad_hbonds + n_outliers + n_wrong_region > 0
    n1 = safe_div(n_bad_helix_sheet_records,n_total_helix_sheet_records)*100.
    # n2 = percentage of SS elements that have at least one residue belonging to a wrong region of Ramachandran plot (per given model);
    n2 = safe_div(n_elem_with_wrong_rama,n_total_helix_sheet_records)*100.
    # n3 = percentage of SS elements that have at least one residue being a Ramachandran plot outlier (per given model);
    n3 = safe_div(n_elem_with_rama_out,n_total_helix_sheet_records)*100.
    # n4 = percentage of bad H bonds (per given model).
    n4 = safe_div(cumm_n_bad_hbonds,cumm_n_hbonds)*100. # No per SS element separation
    # percentage of SS elements that have at least one bad H bond (per given model)
    n5 = safe_div(n_elem_with_bad_hbond,n_total_helix_sheet_records)*100.
    print("Overall info:", file=self.log)
    print("  Total HELIX+SHEET recods       :", n_total_helix_sheet_records, file=self.log)
    print("  Total bad HELIX+SHEET recods   :", n_bad_helix_sheet_records, file=self.log)
    print("  Total declared H-bonds         :", cumm_n_hbonds, file=self.log)
    print("  Total mediocre H-bonds (%.1f-%.1fA):" % (
        self.params.mediocre_hbond_cutoff, self.params.bad_hbond_cutoff), \
        cumm_n_mediocre_hbonds, file=self.log)
    print("  Total bad H-bonds (>%.1fA)      :" % self.params.bad_hbond_cutoff, \
        cumm_n_bad_hbonds, file=self.log)
    print("  Total Ramachandran outliers    :", cumm_n_rama_out, file=self.log)
    print("  Total wrong Ramachandrans      :", cumm_n_wrong_reg, file=self.log)
    print("All done.", file=self.log)
    help_string = """\
  Total bad HELIX+SHEET recods does not include records with syntax mistakes
    (they are outputted separately in the beginning of the log),
    but includes empty records (without corresponding atoms in the model)
    and records with any deviations in geometry (bad/mediocre bonds,
    Ramachandran angles are outliers or wrong).
  Ramachandran outliers - residues in disallowed region of Ramachandran plot.
  Wrong Ramachandrans - residues in favored and allowed regions of Ramachandran
    plot, but don't belong to region of annotated secondary structure element.
    For example, residue annotated as HELIX has phi-psi angles in beta-strand
    region and vice versa.
"""
    print(help_string, file=self.log)

    if self.params.filter_annotation:
      filtered_ann = ss_annot.filter_annotation(hierarchy=pdb_h)
      print("Filtered annotation:", file=self.log)
      print(filtered_ann.as_pdb_str(), file=self.log)
    self.results = group_args(
      n_total_helix_sheet_records = n_total_helix_sheet_records,
      n_bad_helix_sheet_records   = n_bad_helix_sheet_records,
      n_hbonds                    = cumm_n_hbonds,
      n_mediocre_hbonds           = cumm_n_mediocre_hbonds,
      n_bad_hbonds                = cumm_n_bad_hbonds,
      n_rama_out                  = cumm_n_rama_out,
      n_wrong_reg                 = cumm_n_wrong_reg,
      n1                          = n1,
      n2                          = n2,
      n3                          = n3,
      n4                          = n4,
      n5                          = n5,
      # Number of helices with syntax error. Specifically, those producing
      # ValueError on converting the field to a number.
      n_bad_helices               = n_bad_helices,
      n_bad_sheets                = n_bad_sheets)

  @staticmethod
  def get_default_params():
    """
    Get extracted params
    """
    return iotbx.phil.parse(
          input_string=master_phil_str,
          process_includes=True).extract()

  def get_results(self):
    return self.results
