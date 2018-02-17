from __future__ import division

from mmtbx.secondary_structure import manager
import iotbx.pdb
import iotbx.phil
from libtbx.utils import Sorry, safe_div
import cStringIO
import sys
from libtbx import easy_mp
import mmtbx
from mmtbx.building.loop_closure.utils import get_phi_psi_atoms, get_rama_score, \
    rama_evaluate, get_pair_angles
from libtbx import group_args

import boost.python
ext = boost.python.import_ext("mmtbx_validation_ramachandran_ext")
from mmtbx_validation_ramachandran_ext import rama_eval
from mmtbx.validation import ramalyze

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

  def __call__(self, hsh_tuple):
    temp_annot = iotbx.pdb.secondary_structure.annotation(
        helices = hsh_tuple[0],
        sheets = hsh_tuple[1])
    helix = len(hsh_tuple[0]) > 0
    # print temp_annot.as_pdb_str().replace('\n',' '),
    ss_params = mmtbx.secondary_structure.default_params
    ss_params.secondary_structure.enabled=True
    ss_params.secondary_structure.protein.remove_outliers=False
    ss_params.secondary_structure.protein.helix=[]
    ss_params.secondary_structure.protein.sheet=[]
    ss_params.secondary_structure.nucleic_acid.enabled=False
    ss_m_log = cStringIO.StringIO()

    ss_manager = mmtbx.secondary_structure.manager(
        pdb_hierarchy=self.pdb_h,
        atom_selection_cache=self.asc,
        sec_str_from_pdb_file=temp_annot,
        params=ss_params.secondary_structure,
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
    sele = ss_manager.selection_cache.selection(temp_annot.as_atom_selections()[0])
    ss_h = self.pdb_h.select(sele)
    phi_psi_atoms = get_phi_psi_atoms(ss_h)

    n_outliers = 0
    n_wrong_region = 0
    for phi_psi_pair, rama_key in phi_psi_atoms:
      rama_score = get_rama_score(phi_psi_pair, self.r, rama_key)
      if rama_evaluate(phi_psi_pair, self.r, rama_key) == ramalyze.RAMALYZE_OUTLIER:
        n_outliers += 1
      else:
        reg = gather_ss_stats.helix_sheet_rama_region(phi_psi_pair)
        if (reg == 1 and not helix) or (reg == 2 and helix):
          n_wrong_region += 1
          # print "  Wrong region:", phi_psi_pair[0][2].id_str(), reg, helix

    del ss_manager
    del ss_params
    return n_hbonds, n_bad_hbonds, n_mediocre_hbonds, hb_lens, n_outliers, n_wrong_region

  @classmethod
  def helix_sheet_rama_region(cls,phi_psi_pair):
    # result 1 - helix, 2 - sheet, 0 - other
    # cutoff: phi < 70 - helix, phi>=70 - sheet
    phi_psi_angles = get_pair_angles(phi_psi_pair, round_coords=False)
    if phi_psi_angles[1] < 70:
      return 1
    else:
      return 2

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
    ss_log = cStringIO.StringIO()
    try:
      ss_annot = self.model.get_ss_annotation(log=ss_log)
    except Sorry as e:
      print >> out, " Syntax error in SS: %s" % e.message
      return
    ss_log_cont = ss_log.getvalue()
    n_bad_helices = ss_log_cont.count("Bad HELIX")
    n_bad_sheets = ss_log_cont.count("Bad SHEET")
    pdb_h = self.model.get_hierarchy()
    if ss_annot is None or ss_annot.is_empty():
      print >> self.log, "No SS annotation, nothing to analyze"
      return
    if n_bad_helices > 0:
      print >> self.log, "Number of bad helices: %d" % n_bad_helices
    if n_bad_helices > 0:
      print >> self.log, "Number of bad sheets: %d" % n_bad_sheets
    if model.get_number_of_models() != 1 :
      raise Sorry("Multiple models not supported.")
    if not pdb_h.contains_protein():
      print >> self.log, "Protein is not found in the model"
      return
    if pdb_h.is_ca_only():
      print >> self.log, "Error: CA-only model"
      return
    if is_ca_and_something(pdb_h):
      print >> self.log, "CA-only and something model"
      return
    if some_chains_are_ca(pdb_h):
      print >> self.log, "some chains are CA-only"
      return

    n_total_helix_sheet_records = ss_annot.get_n_helices()+ss_annot.get_n_sheets()
    n_bad_helix_sheet_records = 0
    # Empty stuff:
    empty_annots = ss_annot.remove_empty_annotations(pdb_h)
    number_of_empty_helices = empty_annots.get_n_helices()
    number_of_empty_sheets = empty_annots.get_n_sheets()
    n_bad_helix_sheet_records += (number_of_empty_helices+number_of_empty_sheets)
    if number_of_empty_helices > 0:
      print >> self.log, "Helices without corresponding atoms in the model (%d):" % number_of_empty_helices
      for h in empty_annots.helices:
        print >> self.log, "  ", h.as_pdb_str()
    if number_of_empty_sheets > 0:
      print >> self.log, "Sheets without corresponding atoms in the model (%d):" % number_of_empty_sheets
      for sh in empty_annots.sheets:
        print >> self.log, "  ", sh.as_pdb_str()

    print >> self.log, "Checking annotations thoroughly, use nproc=<number> if it is too slow..."

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
          print >> self.log, "Bad annotation found:"
          print >> self.log, "%s" % ss_elem.as_pdb_str()
          print >> self.log, "  Total hb: %d, mediocre: %d, bad: %d, Rama outliers: %d, Rama wrong %d" % (
              n_hbonds, n_mediocre_hbonds, n_bad_hbonds, n_outliers, n_wrong_region)
          print >> self.log, "-"*80

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
    print >> self.log, "Overall info:"
    print >> self.log, "  Total HELIX+SHEET recods       :", n_total_helix_sheet_records
    print >> self.log, "  Total bad HELIX+SHEET recods   :", n_bad_helix_sheet_records
    print >> self.log, "  Total declared H-bonds         :", cumm_n_hbonds
    print >> self.log, "  Total mediocre H-bonds (%.1f-%.1fA):" % (
        self.params.mediocre_hbond_cutoff, self.params.bad_hbond_cutoff), \
        cumm_n_mediocre_hbonds
    print >> self.log, "  Total bad H-bonds (>%.1fA)      :" % self.params.bad_hbond_cutoff, \
        cumm_n_bad_hbonds
    print >> self.log, "  Total Ramachandran outliers    :", cumm_n_rama_out
    print >> self.log, "  Total wrong Ramachandrans      :", cumm_n_wrong_reg
    print >> self.log, "All done."

    if self.params.filter_annotation:
      filtered_ann = ss_annot.filter_annotation(hierarchy=pdb_h)
      print >> self.log, "Filtered annotation:"
      print >> self.log, filtered_ann.as_pdb_str()
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
