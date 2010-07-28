
from __future__ import division
import mmtbx.base_pairing
import iotbx.pdb
import iotbx.pdb.secondary_structure
from scitbx.array_family import shared, flex
import libtbx.phil
import libtbx.object_oriented_patterns as oop
from libtbx import smart_open, easy_run
from libtbx.utils import Sorry, Usage
from libtbx import adopt_init_args, group_args
import libtbx.load_env
from math import sqrt
import cStringIO
import sys, os

default_sigma = 0.05
default_slack = 0.00
h_o_distance_ideal = 1.975
h_o_distance_max = 2.5
n_o_distance_ideal = 2.9
n_o_distance_max = 3.5

ss_restraint_params_str = """
  verbose = False
    .type = bool
  substitute_n_for_h = None
    .type = bool
    .short_caption = Substitute N for H
    .style = tribool
  restrain_helices = True
    .type = bool
  alpha_only = False
    .type = bool
    .short_caption = Use alpha helices only
  restrain_sheets = True
    .type = bool
  restrain_base_pairs = True
    .type = bool
  remove_outliers = True
    .type = bool
    .short_caption = Filter bond outliers
  restrain_initial_values = False
    .type = bool
    .style = hidden
  sigma = %.3f
    .type = float
    .short_caption = Restraint sigma
    .style = bold
  slack = %.3f
    .type = float
    .short_caption = Restraint slack
    .style = bold
  h_o_distance_ideal = %.3f
    .type = float
    .short_caption = Ideal H-O distance
  h_o_distance_max = %.3f
    .type = float
    .short_caption = Max. allowable H-O distance
  n_o_distance_ideal = %.3f
    .type = float
    .short_caption = Ideal N-O distance
  n_o_distance_max = %.3f
    .type = float
    .short_caption = Max. allowable N-O distance
""" % (default_sigma, default_slack, h_o_distance_ideal, h_o_distance_max,
       n_o_distance_ideal, n_o_distance_max)

ss_tardy_params_str = """\
  group_helix_backbone = False
    .style = bool
"""

helix_classes = ["unknown"] * 10
helix_classes[0] = "alpha"
helix_classes[2] = "pi"
helix_classes[4] = "3_10"

def get_helix_class (helix_class) :
  class_id = int(helix_class) - 1
  return helix_classes[class_id]

helix_group_params_str = """
helix
  .multiple = True
  .optional = True
  .style = noauto
{
  selection = None
    .type = str
    .style = bold selection
  helix_type = *alpha pi 3_10 unknown
    .type = choice
    .help = Type of helix, defaults to alpha.  Only alpha, pi, and 3_10 \
      helices are used for hydrogen-bond restraints.
    .style = bold
  restraint_sigma = None
    .type = float
  restraint_slack = None
    .type = float
  backbone_only = False
    .type = bool
    .help = Only applies to rigid-body groupings, and not H-bond restraints \
      which are already backbone-only.
}"""

sheet_group_params_str = """
sheet
  .multiple = True
  .optional = True
  .style = noauto
{
  first_strand = None
    .type = str
    .style = bold selection
  strand
    .multiple = True
    .optional = True
  {
    selection = None
      .type = str
      .style = bold selection
    sense = parallel antiparallel *unknown
      .type = choice
      .style = bold
    bond_start_current = None
      .type = str
      .style = bold selection
    bond_start_previous = None
      .type = str
      .style = bold selection
  }
  restraint_sigma = None
    .type = float
  restraint_slack = None
    .type = float
  backbone_only = False
    .type = bool
    .help = Only applies to rigid-body groupings, and not H-bond restraints \
      which are already backbone-only.
}
"""

ss_group_params_str = """%s\n%s""" % (helix_group_params_str,
  sheet_group_params_str)
ss_tardy_params_str = "" # XXX: remove this later

sec_str_master_phil_str = """
input
  .style = box auto_align
{
%s
  find_automatically = None
    .type = bool
    .style = bold tribool
  preserve_protein_segid = False
    .type = bool
    .style = bold
  preserve_nucleic_acid_segid = False
    .type = bool
    .style = bold
  force_nucleic_acids = False
    .type = bool
    .short_caption = Force base pair detection
    .help = This will ignore the automatic chain type detection and run \
      the base pair detection using PROBE even if no nucleic acids are found. \
      Useful for tRNAs which have a large number of modified bases.
}
h_bond_restraints
  .short_caption = Hydrogen bonding restraints
  .style = box auto_align
{
%s
}
tardy
  .help = UNUSED
  .style = box auto_align noauto
{
%s
}
%s
nucleic_acids
  .caption = If sigma and slack are not defined for nucleic acids, the \
    overall default settings for protein hydrogen bonds will be used \
    instead.
  .style = box auto_align noauto
{
  sigma = None
    .type = float
    .short_caption = Sigma for nucleic acid H-bonds
    .help = Defaults to global setting
  slack = None
    .type = float
    .short_caption = Slack for nucleic acid H-bonds
    .help = Defaults to global setting
#  planar = False
#    .type = bool
#    .short_caption = Use planar restraints by default
#    .style = noauto
  %s
}
""" % (iotbx.pdb.secondary_structure.ss_input_params_str,
       ss_restraint_params_str, ss_tardy_params_str, ss_group_params_str,
       mmtbx.base_pairing.dna_rna_params_str)

sec_str_master_phil = libtbx.phil.parse(sec_str_master_phil_str)
default_params = sec_str_master_phil.fetch().extract()

def sec_str_from_phil (phil_str) :
  ss_phil = libtbx.phil.parse(phil_str)
  return sec_str_master_phil.fetch(source=ss_phil).extract()

use_resids = False # XXX: for debugging purposes only

class _annotation (oop.injector, iotbx.pdb.secondary_structure.annotation) :
  def as_restraint_groups (self, log=sys.stderr, prefix_scope="",
      add_segid=None) :
    phil_strs = []
    for helix in self.helices :
      helix_phil = helix.as_restraint_group(log, prefix_scope, add_segid)
      if helix_phil is not None :
        phil_strs.append(helix_phil)
    for sheet in self.sheets :
      sheet_phil = sheet.as_restraint_group(log, prefix_scope, add_segid)
      if sheet_phil is not None :
        phil_strs.append(sheet_phil)
    return "\n".join(phil_strs)

class _pdb_helix (oop.injector, iotbx.pdb.secondary_structure.pdb_helix) :
  def as_restraint_group (self, log=sys.stderr, prefix_scope="",
      add_segid=None) :
    if self.start_chain_id != self.end_chain_id :
      print >> log, "Helix chain ID mismatch: starts in %s, ends in %s" % (
        self.start_chain_id, self.end_chain_id)
      return None
    segid_extra = ""
    if add_segid is not None :
      segid_extra = "and segid '%s' " % add_segid
    if use_resids :
      resid_start = "%d%s" % (self.start_resseq, self.start_icode)
      resid_end = "%d%s" % (self.end_resseq, self.end_icode)
      sele = "chain '%s' %sand resid %s through %s" % (self.start_chain_id,
        segid_extra, resid_start, resid_end)
    else :
      sele = "chain '%s' %sand resseq %d:%d" % (self.start_chain_id,
        segid_extra, self.start_resseq, self.end_resseq)
    if prefix_scope != "" and not prefix_scope.endswith(".") :
      prefix_scope += "."
    rg = """\
%shelix {
  selection = "%s"
  helix_type = %s
}""" % (prefix_scope, sele, helix_classes[self.helix_class - 1])
    return rg

class _pdb_sheet (oop.injector, iotbx.pdb.secondary_structure.pdb_sheet) :
  def as_restraint_group (self, log=sys.stderr, prefix_scope="",
      add_segid=None) :
    if len(self.strands) == 0 :
      return None
    selections = []
    senses = []
    reg_curr = []
    reg_prev = []
    segid_extra = ""
    if add_segid is not None :
      segid_extra = "and segid '%s' " % add_segid
    for (strand,registration) in zip(self.strands, self.registrations) :
      if use_resids :
        resid_start = "%d%s" % (strand.start_resseq, strand.start_icode)
        resid_end = "%d%s" % (strand.end_resseq, strand.end_icode)
        sele = "chain '%s' %sand resid %s through %s" % (strand.start_chain_id,
          segid_extra, resid_start, resid_end)
      else :
        sele = "chain '%s' %sand resseq %d:%d" % (strand.start_chain_id,
          segid_extra, strand.start_resseq, strand.end_resseq)
      selections.append(sele)
      if strand.sense == 0 :
        senses.append("unknown")
      elif strand.sense == -1 :
        senses.append("antiparallel")
      elif strand.sense == 1 :
        senses.append("parallel")
      else :
        raise Sorry("Sense must be 0, 1, or -1.")
      if registration is not None :
        if use_resids :
          sele_base = "chain '%s' %sand resid %s"
          resid_curr = "%d%s" % (registration.cur_resseq,registration.cur_icode)
          resid_prev = "%d%s" % (registration.prev_resseq,
            registration.prev_icode)
          reg_curr.append(sele_base % (registration.cur_chain_id,segid_extra,
            resid_curr))
          reg_prev.append(sele_base % (registration.prev_chain_id,segid_extra,
            resid_prev))
        else :
          reg_curr.append("chain '%s' %sand resseq %d" % (
            registration.cur_chain_id, segid_extra, registration.cur_resseq))
          reg_prev.append("chain '%s' %sand resseq %d" % (
            registration.prev_chain_id, segid_extra, registration.prev_resseq))
      else :
        reg_curr.append(None)
        reg_prev.append(None)
    n = 0
    first_strand = None
    strands = []
    for (sele, sense, curr, prev) in zip(selections,senses,reg_curr,reg_prev) :
      if n == 0 :
        first_strand = sele
      else :
        strands.append("""\
  strand {
    selection = "%s"
    sense = %s
    bond_start_current = "%s"
    bond_start_previous = "%s"
  }""" % (sele, sense, curr, prev))
      n += 1
    assert first_strand is not None
    if prefix_scope != "" and not prefix_scope.endswith(".") :
      prefix_scope += "."
    phil_str = """
%ssheet {
  first_strand = "%s"
%s
}""" % (prefix_scope, first_strand, "\n".join(strands))
    return phil_str

class hydrogen_bond_table (object) :
  def __init__ (self, bonds, distance, sigma, slack, bond_lengths) :
    assert (bonds.size() == distance.size() == sigma.size() == slack.size() ==
      bond_lengths.size())
    self.bonds = bonds
    self.distance = distance
    self.sigma = sigma
    self.slack = slack
    self.bond_lengths = bond_lengths
    self.flag_use_bond = flex.bool(bonds.size(), True)

  def analyze_distances (self, params, pdb_hierarchy=None, log=sys.stderr) :
    atoms = None
    if params.verbose :
      assert pdb_hierarchy is not None
      atoms = pdb_hierarchy.atoms()
    remove_outliers = params.remove_outliers
    distance_max = params.h_o_distance_max
    distance_ideal = params.h_o_distance_ideal
    if params.substitute_n_for_h :
      distance_max = params.n_o_distance_max
      distance_ideal = params.n_o_distance_ideal
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

  def get_bond_restraint_data (self, filter=True) :
    for i, (donor_i_seq, acceptor_i_seq) in enumerate(self.bonds) :
      if not filter or self.flag_use_bond[i] :
        yield group_args(donor_i_seq=donor_i_seq,
          acceptor_i_seq=acceptor_i_seq,
          sigma=self.sigma[i],
          slack=self.slack[i],
          distance_ideal=self.distance[i])

  def get_simple_bonds (self, filter=True) :
    if filter :
      for i, bond in enumerate(self.bonds) :
        if self.flag_use_bond[i] :
          yield bond
    else :
      for i, bond in enumerate(self.bonds) :
        yield bond

  def as_pymol_dashes (self, pdb_hierarchy, filter=True, out=sys.stdout) :
    atoms = pdb_hierarchy.atoms()
    for (i_seq, j_seq) in self.get_simple_bonds(filter=filter) :
      atom1 = atoms[i_seq].fetch_labels()
      atom2 = atoms[j_seq].fetch_labels()
      base_sele = """chain "%s" and resi %s and name %s"""
      sele1 = base_sele % (atom1.chain_id, atom1.resseq, atom1.name)
      sele2 = base_sele % (atom2.chain_id, atom2.resseq, atom2.name)
      print >>out, "dist %s, %s" % (sele1, sele2)

  def as_refmac_restraints (self, pdb_hierarchy, filter=True, out=sys.stdout) :
    atoms = pdb_hierarchy.atoms()
    for bond in self.get_bond_restraint_data(filter=filter) :
      donor = atoms[bond.donor_i_seq].fetch_labels()
      acceptor = atoms[bond.acceptor_i_seq].fetch_labels()
      cmd = (("exte dist first chain %s residue %s atom %s " +
              "second chain %s residue %s atom %s value %.3f sigma %.2f") %
        (donor.chain_id, donor.resseq, donor.name, acceptor.chain_id,
         acceptor.resseq, acceptor.name, bond.distance_ideal, bond.sigma))
      print >> out, cmd

def hydrogen_bonds_from_selections (
    pdb_hierarchy,
    params,
    log=sys.stderr) :
  bond_i_seqs = shared.stl_set_unsigned()
  sigmas = flex.double()
  slacks = flex.double()
  atoms = pdb_hierarchy.atoms()
  n_atoms = atoms.size()
  sites = atoms.extract_xyz()
  has_bond = flex.bool(sites.size(), False)
  selection_cache = pdb_hierarchy.atom_selection_cache()
  isel = selection_cache.iselection
  donor_name = "H"
  distance_ideal = params.h_bond_restraints.h_o_distance_ideal
  if params.h_bond_restraints.substitute_n_for_h :
    distance_ideal = params.h_bond_restraints.n_o_distance_ideal
    donor_name = "N"
  if params.h_bond_restraints.restrain_helices :
    for helix in params.helix :
      helix_class = helix.helix_type
      if helix_class != "alpha" and params.h_bond_restraints.alpha_only :
        print >> log, "Skipping non-alpha helix (class %s):" % helix_class
        print >> log, "  %s" % helix.selection
        continue
      try :
        donor_isel, acceptor_isel = _donors_and_acceptors(
          base_sele=helix.selection,
          selection_cache=selection_cache,
          atoms=atoms,
          donor_name=donor_name,
          ss_type="helix")
      except RuntimeError, e :
        print >> log, str(e)
        continue
      else :
        if helix_class == "alpha" :
          j = 4
        elif helix_class == "pi" :
          j = 5
        elif helix_class == "3_10" :
          j = 3
        else :
          print >> log, "Don't know bonding for helix class %s." % helix_class
          continue
        sigma = params.h_bond_restraints.sigma
        if helix.restraint_sigma is not None :
          sigma = helix.restraint_sigma
        elif sigma is None :
          raise Sorry(("Please either set the global sigma for hydrogen bond "+
            "restraints, or set the sigma for helix '%s'.") % helix.selection)
        slack = params.h_bond_restraints.slack
        if helix.restraint_slack is not None :
          slack = helix.restraint_slack
        elif slack is None :
          raise Sorry(("Please either set the global slack for hydrogen bond "+
            "restraints, or set the slack for helix '%s'.") % helix.selection)
        n_donors = acceptor_isel.size()
        for n, i_seq in enumerate(acceptor_isel) :
          if (n + j) < n_donors :
            j_seq =  donor_isel[n+j]
            if j_seq == n_atoms :
              continue # dummy index - missing H from PRO
            elif atoms[j_seq].fetch_labels().resname == "PRO" :
              print >> log, "Skipping proline residue in middle of helix:"
              print >> log, "  %s" % atoms[j_seq].fetch_labels().id_str()
              continue
            if has_bond[i_seq] or has_bond[j_seq] :
              print >> log, "One or more atoms already bonded:"
              print >> log, "  %s" % atoms[i_seq].fetch_labels().id_str()
              print >> log, "  %s" % atoms[j_seq].fetch_labels().id_str()
              continue
            has_bond[i_seq] = True
            has_bond[j_seq] = True
            bond_i_seqs.append((j_seq, i_seq))
            sigmas.append(sigma)
            slacks.append(slack)
          else :
            break
  if params.h_bond_restraints.restrain_sheets :
    for i_sheet, sheet in enumerate(params.sheet) :
      if sheet.first_strand is None :
        raise Sorry("First strand must be a valid atom selection.")
      sheet_bonds = []
      prev_strand_sele = sheet.first_strand
      for curr_strand in sheet.strand :
        if curr_strand.selection is None :
          raise Sorry("All strands must have a valid atom selection.")
        elif curr_strand.bond_start_current is None :
          raise Sorry("Missing start of bonding for strand '%s'." %
            curr_strand.selection)
        elif curr_strand.bond_start_previous is None :
          raise Sorry("Missing start of bonding for strand previous to '%s'." %
            curr_strand.selection)
        try :
          try :
            if curr_strand.sense == "unknown" :
              raise RuntimeError(("Skipping strand of unknown sense:\n" +
                "%s") % curr_strand.selection)
            curr_donors, curr_acceptors = _donors_and_acceptors(
              base_sele=curr_strand.selection,
              selection_cache=selection_cache,
              atoms=atoms,
              donor_name=donor_name,
              ss_type="sheet")
            prev_donors, prev_acceptors = _donors_and_acceptors(
              base_sele=prev_strand_sele,
              selection_cache=selection_cache,
              atoms=atoms,
              donor_name=donor_name,
              ss_type="sheet")
            curr_donor_start, curr_acceptor_start = _donors_and_acceptors(
              base_sele=curr_strand.bond_start_current,
              selection_cache=selection_cache,
              atoms=atoms,
              donor_name=donor_name,
              ss_type="sheet")
            prev_donor_start, prev_acceptor_start = _donors_and_acceptors(
              base_sele=curr_strand.bond_start_previous,
              selection_cache=selection_cache,
              atoms=atoms,
              donor_name=donor_name,
              ss_type="sheet")
            new_bonds = _hydrogen_bonds_from_strand_pair(atoms=atoms,
              prev_strand_donors=prev_donors,
              prev_strand_acceptors=prev_acceptors,
              prev_strand_start=prev_donor_start[0],
              curr_strand_donors=curr_donors,
              curr_strand_acceptors=curr_acceptors,
              curr_strand_start=curr_acceptor_start[0],
              sense=curr_strand.sense)
            if new_bonds is None :
              raise RuntimeError(
                ("Can't determine start of bonding for strand pair:\n" +
                "  %s\n  %s\n") % (prev_strand_sele, curr_strand.selection))
            elif len(new_bonds) == 0 :
              raise RuntimeError(("No bonds found for strand pair:\n"+
                "  %s\n  %s\n") % (prev_strand_sele, curr_strand.selection))
            sheet_bonds.extend(new_bonds)
          except RuntimeError, e :
            print >> log, str(e)
        finally :
          prev_strand_sele = curr_strand.selection
      sigma = params.h_bond_restraints.sigma
      if sheet.restraint_sigma is not None :
        sigma = sheet.restraint_sigma
      elif sigma is None :
        raise Sorry(("Please either set the global sigma for hydrogen bond "+
          "restraints, or set the sigma for sheet #%d.") % i_sheet)
      slack = params.h_bond_restraints.slack
      if sheet.restraint_slack is not None :
        slack = sheet.restraint_slack
      elif slack is None :
        raise Sorry(("Please either set the global slack for hydrogen bond "+
          "restraints, or set the slack for sheet #%d.") % i_sheet)
      for (i_seq, j_seq) in sheet_bonds :
        if has_bond[i_seq] or has_bond[j_seq] :
          print >> log, "One or more atoms already bonded:"
          print >> log, "  %s" % atoms[i_seq].fetch_labels().id_str()
          print >> log, "  %s" % atoms[j_seq].fetch_labels().id_str()
          continue
        has_bond[i_seq] = True
        has_bond[j_seq] = True
        bond_i_seqs.append((i_seq, j_seq))
        sigmas.append(sigma)
        slacks.append(slack)
  if params.h_bond_restraints.restrain_base_pairs :
    sigma = params.nucleic_acids.sigma
    if sigma is None :
      sigma = params.h_bond_restraints.sigma
    slack = params.nucleic_acids.slack
    if slack is None :
      slack = params.h_bond_restraints.slack
    for i, base_pair in enumerate(params.nucleic_acids.base_pair) :
      try :
        resname1 = _get_residue_name_from_selection(
          resi_sele=base_pair.base1,
          selection_cache=selection_cache,
          atoms=atoms)
        resname2 = _get_residue_name_from_selection(
          resi_sele=base_pair.base2,
          selection_cache=selection_cache,
          atoms=atoms)
        atom_pairs = mmtbx.base_pairing.get_h_bond_atoms(
          residues=(resname1,resname2),
          pair_type=base_pair.pair_type.upper(),
          use_hydrogens=(not params.h_bond_restraints.substitute_n_for_h))
        for (name1, name2) in atom_pairs :
          sele1 = """name %s and %s""" % (name1, base_pair.base1)
          sele2 = """name %s and %s""" % (name2, base_pair.base2)
          # XXX these aren't necessarily in donor/acceptor order, but it
          # doesn't really matter here.
          (i_seq,j_seq) = _hydrogen_bond_from_selection_pair(sele1, sele2,
            selection_cache)
          if has_bond[i_seq] or has_bond[j_seq] :
            print >> log, "One or more atoms already bonded:"
            print >> log, "  %s" % atoms[i_seq].fetch_labels().id_str()
            print >> log, "  %s" % atoms[j_seq].fetch_labels().id_str()
            continue
          has_bond[i_seq] = True
          has_bond[j_seq] = True
          bond_i_seqs.append((i_seq, j_seq))
          sigmas.append(sigma)
          slacks.append(slack)
      except RuntimeError, e :
        print >> log, str(e)
  return hydrogen_bond_table(bonds=bond_i_seqs,
    distance=flex.double(bond_i_seqs.size(), distance_ideal),
    sigma=sigmas,
    slack=slacks,
    bond_lengths=_get_distances(bond_i_seqs, sites))

def _get_distances (bonds, sites_cart) :
  distances = flex.double(bonds.size(), -1)
  for k, (i_seq, j_seq) in enumerate(bonds) :
    (x1, y1, z1) = sites_cart[i_seq]
    (x2, y2, z2) = sites_cart[j_seq]
    dist = sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
    distances[k] = dist
  return distances

def _donors_and_acceptors (base_sele, selection_cache, atoms, donor_name,
    ss_type) :
  isel = selection_cache.iselection
  donor_sele = "(%s) and (altloc 'A' or altloc ' ') and name %s" % (
    base_sele, donor_name)
  acceptor_sele = "(%s) and (altloc 'A' or altloc ' ') and name O"% base_sele
  donor_isel = isel(donor_sele)
  acceptor_isel = isel(acceptor_sele)
  n_donors = donor_isel.size()
  n_acceptors = acceptor_isel.size()
  n_atoms = atoms.size()
  if n_acceptors == 0 :
    raise RuntimeError("No atoms for selection %s." % acceptor_sele)
  elif n_donors != n_acceptors :
    n_pro = 0
    for k, i_seq in enumerate(acceptor_isel) :
      acceptor_atom = atoms[i_seq].fetch_labels()
      if acceptor_atom.resname.strip() == "PRO" :
        donor_isel.insert(k, n_atoms)
        n_pro += 1
    if (n_donors + n_pro) != n_acceptors :
      raise RuntimeError("""\
hydrogen_bonds_from_selections: incomplete non-PRO residues in %s.
  \"%s\" => %d donors
  \"%s\" => %d acceptors""" % (ss_type, donor_sele, donor_isel.size(),
      acceptor_sele, acceptor_isel.size()))
  return donor_isel, acceptor_isel

def _get_residue_name_from_selection (resi_sele, selection_cache, atoms) :
  i_seqs = selection_cache.iselection(resi_sele)
  if len(i_seqs) == 0 :
    raise RuntimeError("Empty selection '%s'" % resi_sele)
  resnames = []
  for i_seq in i_seqs :
    resnames.append(atoms[i_seq].fetch_labels().resname)
  unique_resnames = set(resnames)
  assert len(unique_resnames) != 0
  if len(unique_resnames) > 1 :
    raise RuntimeError("More than one residue name in selection '%s' (%s)" %
      (resi_sele, ", ".join(unique_resnames)))
  else :
    return resnames[0]

def _hydrogen_bonds_from_strand_pair (atoms,
    prev_strand_donors,
    prev_strand_acceptors,
    prev_strand_start,
    curr_strand_donors,
    curr_strand_acceptors,
    curr_strand_start,
    sense) :
  n_atoms = atoms.size()
  assert sense != "unknown"
  assert prev_strand_donors.size() == prev_strand_acceptors.size()
  assert curr_strand_donors.size() == curr_strand_acceptors.size()
  start_bonding = False
  bonds = []
  n_prev_strand = prev_strand_donors.size()
  n_curr_strand = curr_strand_donors.size()
  i = j = None
  for k, donor_i_seq in enumerate(prev_strand_donors) :
    if donor_i_seq == prev_strand_start :
      i = k
      break
  #print curr_strand_start, curr_strand_acceptors
  for k, acceptor_i_seq in enumerate(curr_strand_acceptors) :
    if acceptor_i_seq == curr_strand_start :
      j = k
      break
  if None in [i, j] :
    return None
  if sense == "antiparallel" :
    while i < n_prev_strand and j > 0 :
      donor1_i_seq = prev_strand_donors[i]
      acceptor1_i_seq = curr_strand_acceptors[j]
      if (donor1_i_seq != n_atoms and
          atoms[donor1_i_seq].fetch_labels().resname.strip() != "PRO") :
        bonds.append((donor1_i_seq, acceptor1_i_seq))
      donor2_i_seq = curr_strand_donors[j]
      acceptor2_i_seq = prev_strand_acceptors[i]
      if (donor2_i_seq != n_atoms and
          atoms[donor2_i_seq].fetch_labels().resname.strip() != "PRO") :
        bonds.append((donor2_i_seq, acceptor2_i_seq))
      i += 2
      j -= 2
  else :
    while i < n_prev_strand and j < n_curr_strand :
      donor1_i_seq = prev_strand_donors[i]
      acceptor1_i_seq = curr_strand_acceptors[j]
      if (donor1_i_seq != n_atoms and
          atoms[donor1_i_seq].fetch_labels().resname.strip() != "PRO") :
        bonds.append((donor1_i_seq, acceptor1_i_seq))
      if (j + 2) >= n_curr_strand :
        break
      donor2_i_seq = curr_strand_donors[j+2]
      acceptor2_i_seq = prev_strand_acceptors[i]
      if (donor2_i_seq != n_atoms and
          atoms[donor2_i_seq].fetch_labels().resname.strip() != "PRO") :
        bonds.append((donor2_i_seq, acceptor2_i_seq))
      i += 2
      j += 2
  return bonds

def _hydrogen_bond_from_selection_pair (donor_sele, acceptor_sele,
    selection_cache) :
  isel = selection_cache.iselection
  donor_i_seqs = isel(donor_sele)
  acceptor_i_seqs = isel(acceptor_sele)
  n_donor_sel = donor_i_seqs.size()
  n_acceptor_sel = acceptor_i_seqs.size()
  if n_donor_sel == 0 or n_acceptor_sel == 0 :
    raise RuntimeError("""\
analyze_h_bonds: one or more atoms missing
  %s (%d atoms)
  %s (%d atoms)""" % (donor_sele, n_donor_sel, acceptor_sele, n_acceptor_sel))
  elif n_donor_sel > 1 or n_acceptor_sel > 1 :
    raise RuntimeError("""\
analyze_h_bonds: multiple atoms matching a selection
  %s (%d atoms)
  %s (%d atoms)""" % (donor_sele, n_donor_sel, acceptor_sele, n_acceptor_sel))
  return (donor_i_seqs[0], acceptor_i_seqs[0])

def get_pdb_hierarchy (file_names) :
  pdb_combined = iotbx.pdb.combine_unique_pdb_files(file_names=file_names)
  pdb_structure = iotbx.pdb.input(source_info=None,
    lines=flex.std_string(pdb_combined.raw_records))
  return pdb_structure.construct_hierarchy()

def restraint_groups_as_pdb_helices (pdb_hierarchy, helices, log=sys.stderr) :
  isel = pdb_hierarchy.atom_selection_cache().iselection
  atoms = [ a for a in pdb_hierarchy.atoms_with_labels() ]
  pdb_helices = []
  for i, helix_params in enumerate(helices) :
    if helix_params.selection is None :
      print >> log, "Empty helix at serial %d." % (i+1)
      continue
    sele_str = ("(%s) and name N and (altloc 'A' or altloc ' ')" %
                helix_params.selection)
    amide_isel = isel(sele_str)
    start_atom = atoms[amide_isel[0]]
    end_atom = atoms[amide_isel[-1]]
    if helix_params.helix_type == "unknown" :
      helix_class = 2
    else :
      helix_class = helix_classes.index(helix_params.helix_type)
    current_helix = iotbx.pdb.secondary_structure.pdb_helix(
      serial=i+1,
      helix_id=i+1,
      start_resname=start_atom.resname,
      start_chain_id=start_atom.chain_id,
      start_resseq=start_atom.resseq,
      start_icode=start_atom.icode,
      end_resname=end_atom.resname,
      end_chain_id=end_atom.chain_id,
      end_resseq=end_atom.resseq,
      end_icode=end_atom.icode,
      helix_class=helix_class,
      comment="",
      length=amide_isel.size())
    pdb_helices.append(current_helix)
  return pdb_helices

def restraint_groups_as_pdb_sheets (pdb_hierarchy, sheets, log=sys.stderr) :
  isel = pdb_hierarchy.atom_selection_cache().iselection
  atoms = [ a for a in pdb_hierarchy.atoms_with_labels() ]
  pdb_sheets = []
  for i, sheet in enumerate(sheets) :
    sheet_id = string.uppercase[i]
    if sheet.first_strand is None :
      print >> log, "Missing first strand in sheet %s" % sheet_id
    current_sheet = iotbx.pdb.secondary_structure.pdb_sheet(
      sheet_id=sheet_id,
      n_strands=1+len(sheet.strand),
      strands=[],
      registrations=[])
    first_strand = __strand_group_as_pdb_strand(isel=isel,
      selection=sheet.first_strand,
      atoms=atoms,
      log=log,
      sense=None)
    current_sheet.add_strand(first_strand)
    current_sheet.add_registration(None)
    base_sele = "(%s) and name N and (altloc 'A' or altloc ' ')"
    for strand in sheet.strand :
      pdb_strand = __strand_group_as_pdb_strand(isel=isel,
        selection=strand.selection,
        atoms=atoms,
        log=log,
        sense=strand.sense)
      current_sheet.add_strand(pdb_strand)
      s1 = base_sele % strand.bond_start_current
      s2 = base_sele % strand.bond_start_previous
      reg_curr_isel = isel(s1)
      reg_prev_isel = isel(s2)
      if reg_curr_isel.size() == 0 or reg_prev_isel.size() == 0 :
        current_sheet.add_registration(None)
        continue
      reg_curr_atom = atoms[reg_curr_isel[0]]
      reg_prev_atom = atoms[reg_prev_isel[0]]
      registration = group_args(
        cur_atom=donor_name, #reg_curr_atom.name,
        cur_resname=reg_curr_atom.resname,
        cur_chain_id=reg_curr_atom.chain_id,
        cur_resseq=reg_curr_atom.resseq,
        cur_icode=reg_curr_atom.icode,
        prev_atom=reg_prev_atom.name,
        prev_resname="O", #reg_prev_atom.resname,
        prev_chain_id=reg_prev_atom.chain_id,
        prev_resseq=reg_prev_atom.resseq,
        prev_icode=reg_prev_atom.icode)
      current_sheet.add_registration(registration)
    pdb_sheets.append(current_sheet)
  return pdb_sheets

def __strand_group_as_pdb_strand (isel, selection, atoms, log, sense) :
  if sense is None or sense == "unknown" :
    int_sense = 0
  elif sense == "parallel" :
    int_sense = 1
  elif sense == "antiparallel" :
    int_sense = -1
  strand_isel = isel("(%s) and name N and (altloc 'A' or altloc ' ')" % (
    selection))
  start_atom = atoms[strand_isel[0]]
  end_atom = atoms[strand_isel[-1]]
  pdb_strand = group_args(
    sheet_id=sheet_id,
    strand_id=i+1,
    start_resname=start_atom.resname,
    start_chain_id=start_atom.chain_id,
    start_resseq=start_atom.resseq,
    start_icode=start_atom.icode,
    end_resname=end_atom.resname,
    end_chain_id=end_atom.chain_id,
    end_resseq=end_atom.resseq,
    end_icode=end_atom.icode,
    sense=int_sense)
  return pdb_strand

########################################################################
class manager (object) :
  def __init__ (self,
                pdb_hierarchy,
                xray_structure,
                sec_str_from_pdb_file=None,
                params=None,
                assume_hydrogens_all_missing=None,
                tmp_dir=None) :
    adopt_init_args(self, locals())
    self._was_initialized = False
    if self.params is None :
      self.params = sec_str_master_phil.fetch().extract()
    if self.tmp_dir is None :
      self.tmp_dir = os.getcwd()
    if self.xray_structure is None :
      self.xray_structure = pdb_hierarchy.extract_xray_structure()
    if self.assume_hydrogens_all_missing is None :
      xrs = self.xray_structure
      sctr_keys = xrs.scattering_type_registry().type_count_dict().keys()
      self.assume_hydrogens_all_missing = not ("H" in sctr_keys or
        "D" in sctr_keys)
    self.selection_cache = pdb_hierarchy.atom_selection_cache()
    if self.params.h_bond_restraints.substitute_n_for_h is None :
      self.params.h_bond_restraints.substitute_n_for_h = \
        self.assume_hydrogens_all_missing

  def as_phil_str (self, master_phil=sec_str_master_phil) :
    return master_phil.format(python_object=self.params)

  def initialize (self, log=sys.stderr) :
    if not self._was_initialized :
      self.find_automatically(log=log)
      self.show_summary(out=log)
      self._was_initialized = True

  def find_automatically (self, log=sys.stderr) :
    params = self.params
    find_automatically = params.input.find_automatically
    # XXX: check for presence of protein first?
    if len(params.helix) == 0 and len(params.sheet) == 0 :
      print >> log, "No existing secondary structure definitions found."
      if self.sec_str_from_pdb_file is None and find_automatically != False :
        print >> log, "No HELIX or SHEET records found in PDB file."
        find_automatically = True
    if find_automatically :
      if params.input.preserve_protein_segid :
        self.sec_str_from_pdb_file = self.find_sec_str_with_segids(log=log)
      else :
        self.sec_str_from_pdb_file = self.find_sec_str(log=log)
    if self.sec_str_from_pdb_file is not None :
      if isinstance(self.sec_str_from_pdb_file, list) :
        print >> log, "  Interpreting HELIX and SHEET records for individual chains"
        ss_params = []
        for annotation, segid in self.sec_str_from_pdb_file :
          ss_phil = annotation.as_restraint_groups(log=log,
            prefix_scope="",
            add_segid=segid)
          ss_params.append(ss_phil)
        ss_params_str = "\n".join(ss_params)
      else :
        print >> log, "  Interpreting HELIX and SHEET records from PDB file"
        ss_params_str = self.sec_str_from_pdb_file.as_restraint_groups(log=log,
          prefix_scope="")
      self.apply_phil_str(ss_params_str, log=log)
    # Step 2: nucleic acids
    if (find_nucleic_acids(self.pdb_hierarchy) or
        params.input.force_nucleic_acids) :
      find_automatically = params.input.find_automatically
      if len(params.nucleic_acids.base_pair) == 0 :
        if find_automatically != False :
          find_automatically = True
      if find_automatically :
        if params.input.preserve_nucleic_acid_segid :
          base_pairs = self.find_base_pairs_with_segids(log=log,
            force=params.input.force_nucleic_acids)
        else :
          base_pairs = self.find_base_pairs(log=log)
        if base_pairs is not None :
          bp_phil = libtbx.phil.parse(base_pairs)
          bp_params = sec_str_master_phil.fetch(source=bp_phil).extract()
          self.params.nucleic_acids.base_pair = \
            bp_params.nucleic_acids.base_pair

  def find_sec_str (self, log=sys.stderr) :
    pdb_str = self.pdb_hierarchy.as_pdb_string()
    (records, stderr) = run_ksdssp_direct(pdb_str)
    sec_str_from_pdb_file = iotbx.pdb.secondary_structure.process_records(
      records=records,
      allow_none=True)
    return sec_str_from_pdb_file

  def find_sec_str_with_segids (self, log=sys.stderr) :
    annotations = []
    for chain in self.pdb_hierarchy.models()[0].chains() :
      if not chain.conformers()[0].is_protein() :
        continue
      segid = chain.atoms()[0].segid
      detached_hierarchy = iotbx.pdb.hierarchy.new_hierarchy_from_chain(chain)
      pdb_str = detached_hierarchy.as_pdb_string()
      (records, stderr) = run_ksdssp_direct(pdb_str)
      sec_str_from_pdb_file = iotbx.pdb.secondary_structure.process_records(
        records=records,
        allow_none=True)
      if sec_str_from_pdb_file is not None :
        annotations.append((sec_str_from_pdb_file, segid))
    return annotations

  def find_base_pairs (self, log=sys.stderr) :
    base_pairs = mmtbx.base_pairing.get_phil_base_pairs(
      pdb_hierarchy=self.pdb_hierarchy,
      prefix=None,
      log=log)
    return base_pairs

  def find_base_pairs_with_segids (self, log=sys.stderr, force=False) :
    annotations = []
    for chain in self.pdb_hierarchy.models()[0].chains() :
      if not force and not chain.conformers()[0].is_na() :
        continue
      segid = chain.atoms()[0].segid
      detached_hierarchy = iotbx.pdb.hierarchy.new_hierarchy_from_chain(chain)
      pdb_str = detached_hierarchy.as_pdb_string()
      base_pairs = mmtbx.base_pairing.get_phil_base_pairs(
        pdb_hierarchy=detached_hierarchy,
        prefix=None,
        log=log,
        add_segid=segid)
      annotations.append(base_pairs)
    return "\n".join(annotations)

  def apply_phil_str (self, phil_string, log=sys.stderr, verbose=False) :
    ss_phil = sec_str_master_phil.fetch(source=libtbx.phil.parse(phil_string))
    if verbose :
      ss_phil.show(out=log, prefix="    ")
    new_ss_params = ss_phil.extract()
    self.params.helix = new_ss_params.helix
    self.params.sheet = new_ss_params.sheet

  def apply_params (self, params) :
    self.params.helix = params.helix
    self.params.sheet = params.sheet
    if self.params.h_bond_restraints.substitute_n_for_h is None :
      self.params.h_bond_restraints.substitute_n_for_h = \
        self.assume_hydrogens_all_missing

  def get_bonds_table (self, log=sys.stderr, verbose=True) :
    params = self.params
    bonds_table = hydrogen_bonds_from_selections(
      pdb_hierarchy=self.pdb_hierarchy,
      params=params,
      log=log)
    if verbose :
      print >> log, ""
      print >> log, "  %d hydrogen bonds defined." % bonds_table.bonds.size()
    bonds_table.analyze_distances(params=params.h_bond_restraints,
      pdb_hierarchy=self.pdb_hierarchy,
      log=log)
    return bonds_table

  def calculate_structure_content (self) :
    isel = self.selection_cache.iselection
    calpha = isel("name N and (altloc ' ' or altloc 'A')")
    n_alpha = self.alpha_selection(limit="and name N").count(True)
    n_beta = self.beta_selection(limit="and name N").count(True)
    if calpha.size() == 0 :
      return (0.0, 0.0)
    return (n_alpha / calpha.size(), n_beta / calpha.size())

  def show_summary (self, out=sys.stdout) :
    (frac_alpha, frac_beta) = self.calculate_structure_content()
    n_helices = len(self.params.helix)
    n_sheets  = len(self.params.sheet)
    print >> out, ""
    print >> out, "  %d helices and %d sheets defined" % (n_helices,n_sheets)
    print >> out, "  %.1f%% alpha, %.1f%% beta" %(frac_alpha*100,frac_beta*100)
    print >> out, ""

  def alpha_selections (self, limit="") :
    sele = self.selection_cache.selection
    all_selections = []
    for helix in self.params.helix :
      if helix.selection is not None :
        helix_sel = sele("(%s) %s and (altloc ' ' or altloc 'A')" %
          (helix.selection, limit))
        all_selections.append(helix_sel)
    return all_selections

  def get_helix_types (self) :
    return [ helix.helix_type for helix in self.params.helix ]

  def alpha_selection (self, limit="") :
    whole_selection = flex.bool(self.xray_structure.sites_cart().size())
    for helix in self.alpha_selections(limit=limit) :
      whole_selection |= helix
    return whole_selection

  def beta_selections (self, limit="") :
    sele = self.selection_cache.selection
    all_selections = []
    for sheet in self.params.sheet :
      sheet_selection = flex.bool(self.xray_structure.sites_cart().size())
      strand_sel = sele("(%s) %s and (altloc ' ' or altloc 'A')" %
        (sheet.first_strand, limit))
      sheet_selection |= strand_sel
      for strand in sheet.strand :
        strand_sel = sele("(%s) %s and (altloc ' ' or altloc 'A')" %
          (strand.selection, limit))
        sheet_selection |= strand_sel
      all_selections.append(sheet_selection)
    return all_selections

  def beta_selection (self, limit="") :
    whole_selection = flex.bool(self.xray_structure.sites_cart().size())
    for sheet in self.beta_selections(limit=limit) :
      whole_selection |= sheet
    return whole_selection

def process_structure (params, processed_pdb_file, tmp_dir, log,
    assume_hydrogens_all_missing=None) :
  acp = processed_pdb_file.all_chain_proxies
  try :
    sec_str_from_pdb_file = acp.extract_secondary_structure()
  except Exception :
    sec_str_from_pdb_file = None
  pdb_hierarchy = acp.pdb_hierarchy
  xray_structure = acp.extract_xray_structure()
  structure_manager = manager(
    pdb_hierarchy=pdb_hierarchy,
    xray_structure=xray_structure,
    sec_str_from_pdb_file=sec_str_from_pdb_file,
    params=params,
    assume_hydrogens_all_missing=assume_hydrogens_all_missing,
    tmp_dir=tmp_dir)
  return structure_manager

def find_nucleic_acids (pdb_hierarchy) :
  for model in pdb_hierarchy.models() :
    for chain in model.chains() :
      for conformer in chain.conformers() :
        if conformer.is_na() :
          return True
  return False

def get_ksdssp_exe_path():
  if (not libtbx.env.has_module(name="ksdssp")):
    raise RuntimeError("ksdssp module is not configured")
  exe_path = libtbx.env.under_build("ksdssp/exe/ksdssp")
  if (os.name == "nt"):
    exe_path += ".exe"
  if (not os.path.isfile(exe_path)):
    raise RuntimeError("ksdssp executable is not available")
  return exe_path

def run_ksdssp (file_name, log=sys.stderr) :
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

def manager_from_pdb_file (pdb_file) :
  assert os.path.isfile(pdb_file)
  pdb_in = file_reader.any_file(pdb_file, force_type="pdb")
  pdb_hierarchy = pdb_in.file_object.construct_hierarchy()
  xray_structure = pdb_in.file_object.xray_structure_simple()
  ss_manager  = manager(pdb_hierarchy=pdb_hierarchy,
    xray_structure=xray_structure)
  return ss_manager

def calculate_structure_content (pdb_file) :
  ss_manager = manager_from_pdb_file(pdb_file)
  ss_manager.find_automatically()
  return ss_manager.calculate_structure_content()

def run (args, out=sys.stdout, log=sys.stderr) :
  pdb_files = []
  sources = []
  force_new_annotation = False
  master_phil = libtbx.phil.parse("""
    show_all_params = False
      .type = bool
    show_histograms = False
      .type = bool
    format = *phenix phenix_bonds pymol refmac
      .type = choice
%s""" % sec_str_master_phil_str)
  parameter_interpreter = libtbx.phil.command_line.argument_interpreter(
    master_phil=master_phil,
    home_scope="")
  for arg in args :
    if os.path.isfile(arg) :
      if iotbx.pdb.is_pdb_file(arg) :
        pdb_files.append(os.path.abspath(arg))
      else :
        try :
          user_phil = libtbx.phil.parse(file_name=arg)
        except RuntimeError :
          print "Unrecognizable file format for %s" % arg
        else :
          sources.append(user_phil)
    else :
      if arg.startswith("--") :
        arg = arg[2:] + "=True"
      try :
        user_phil = parameter_interpreter.process(arg=arg)
        sources.append(user_phil)
      except RuntimeError :
        print "Unrecognizable parameter %s" % arg
  params = master_phil.fetch(sources=sources).extract()
  if len(pdb_files) > 0 :
    params.input.file_name.extend(pdb_files)
  pdb_files = params.input.file_name
  if len(pdb_files) == 0 :
    raise Sorry("No PDB files specified.")
  pdb_combined = iotbx.pdb.combine_unique_pdb_files(file_names=pdb_files)
  pdb_structure = iotbx.pdb.input(source_info=None,
    lines=flex.std_string(pdb_combined.raw_records))
  pdb_hierarchy = pdb_structure.construct_hierarchy()
  xray_structure = pdb_structure.xray_structure_simple()
  if len(pdb_hierarchy.models()) != 1 :
    raise Sorry("Multiple models not supported.")
  m = manager(pdb_hierarchy=pdb_hierarchy,
    xray_structure=xray_structure,
    sec_str_from_pdb_file=pdb_structure.extract_secondary_structure(),
    params=params)
  m.find_automatically(log=log)
  prefix_scope="refinement.secondary_structure"
  if params.show_histograms or params.format != "phenix" :
    prefix_scope = ""
  ss_phil = None
  working_phil = m.as_phil_str(master_phil=master_phil)
  phil_diff = sec_str_master_phil.fetch_diff(source=working_phil)
  #params = working_phil.extract()
  if params.show_histograms :
    #working_phil.show()
    phil_diff.show()
    print >> out, ""
    print >> out, "========== Analyzing hydrogen bonding distances =========="
    print >> out, ""
    bonds_table = m.get_bonds_table(log=log)
  elif params.format == "phenix_bonds" :
    raise Sorry("Not yet implemented.")
  elif params.format in ["pymol", "refmac"] :
    bonds_table = m.get_bonds_table(log=log)
    if params.format == "pymol" :
      bonds_table.as_pymol_dashes(pdb_hierarchy, filter=True, out=out)
    else :
      bonds_table.as_refmac_restraints(pdb_hierarchy, filter=True, out=out)
  else :
    #working_phil.show(out=out)
    print "# These parameters are suitable for use in phenix.refine."
    print "refinement.secondary_structure {"
    if params.show_all_params :
      working_phil.show(prefix="  ")
    else :
      phil_diff.show(prefix="  ")
    print "}"
    #print >> out, ss_params_str
    return working_phil.as_str()

########################################################################
def exercise () :
  from iotbx import file_reader
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1ywf.pdb",
    test=os.path.isfile)
  pdb_file_h = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1ywf_h.pdb",
    test=os.path.isfile)
  if pdb_file is None :
    print "Skipping exercise(): input file not available."
    return False
  log = cStringIO.StringIO()
  pdb_in = file_reader.any_file(pdb_file_h, force_type="pdb").file_object
  pdb_hierarchy = pdb_in.construct_hierarchy()
  xray_structure = pdb_in.xray_structure_simple()
  sec_str_from_pdb_file = pdb_in.extract_secondary_structure()
  m = manager(pdb_hierarchy=pdb_hierarchy,
    xray_structure=xray_structure,
    sec_str_from_pdb_file=sec_str_from_pdb_file)
  m.find_automatically(log=log)
  bonds_table = m.get_bonds_table(log=log)
  assert bonds_table.bonds.size() == 109
  m.params.h_bond_restraints.substitute_n_for_h = True
  bonds_table = m.get_bonds_table(log=log)
  assert bonds_table.flag_use_bond.count(True) == 106
  (frac_alpha, frac_beta) = m.calculate_structure_content()
  assert ("%.3f" % frac_alpha) == "0.643"
  assert ("%.3f" % frac_beta) == "0.075"
  del m
  # using KSDSSP
  if (not libtbx.env.has_module(name="ksdssp")):
    print "Skipping KSDSSP tests: ksdssp module not available."
  else:
    m = manager(pdb_hierarchy=pdb_hierarchy,
      xray_structure=xray_structure,
      sec_str_from_pdb_file=None)
    m.find_automatically(log=log)
    bonds_table = m.get_bonds_table(log=log)
    assert bonds_table.bonds.size() == 93
    m.params.h_bond_restraints.substitute_n_for_h = True
    bonds_table = m.get_bonds_table(log=log)
    assert bonds_table.flag_use_bond.count(True) == 86
    (frac_alpha, frac_beta) = m.calculate_structure_content()
    assert ("%.3f" % frac_alpha) == "0.552"
    assert ("%.3f" % frac_beta) == "0.066"
    del m
    del pdb_hierarchy
    del xray_structure
  # without hydrogens
  pdb_in = file_reader.any_file(pdb_file, force_type="pdb").file_object
  pdb_hierarchy = pdb_in.construct_hierarchy()
  xray_structure = pdb_in.xray_structure_simple()
  sec_str_from_pdb_file = pdb_in.extract_secondary_structure()
  m = manager(pdb_hierarchy=pdb_hierarchy,
    xray_structure=xray_structure,
    sec_str_from_pdb_file=sec_str_from_pdb_file)
  m.find_automatically(log=log)
  bonds_table = m.get_bonds_table(log=log)
  assert bonds_table.bonds.size() == 109
  del m
  # using KSDSSP
  if (libtbx.env.has_module(name="ksdssp")):
    m = manager(pdb_hierarchy=pdb_hierarchy,
      xray_structure=xray_structure,
      sec_str_from_pdb_file=None)
    m.find_automatically(log=log)
    bonds_table = m.get_bonds_table(log=log)
    assert bonds_table.bonds.size() == 93
  # Nucleic acids (requires REDUCE and PROBE)
  if (libtbx.env.has_module(name="reduce") and
      libtbx.env.has_module(name="probe")):
    pass # TODO
  print "OK"

if __name__ == "__main__" :
  exercise()
