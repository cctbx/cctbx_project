
import iotbx.pdb.secondary_structure
from scitbx.array_family import shared, flex
import libtbx.phil
import libtbx.object_oriented_patterns as oop
from libtbx import smart_open
from libtbx.utils import Sorry
from libtbx import adopt_init_args, group_args
from math import sqrt
import sys, os

ss_restraint_params_str = """
  substitute_n_for_h = False
    .type = bool
  restrain_helices = True
    .type = bool
  alpha_only = False
    .type = bool
  restrain_sheets = True
    .type = bool
  remove_outliers = True
    .type = bool
  restrain_initial_values = False
    .type = bool
  slack = 0.1
    .type = float
  sigma = 0.05
    .type = float
  n_o_distance_ideal = 3.0
    .type = float
  n_o_outlier_cutoff = 3.5
    .type = float
  h_o_distance_ideal = 2.0
    .type = float
  h_o_outlier_cutoff = 2.5
    .type = float
"""

ss_group_params_str = """
helix
  .multiple = True
  .optional = True
{
  selection = None
    .type = str
    .style = bold selection
  helix_class = *1 2 3 4 5 6 7 8 9 10
    .type = choice
    .caption = alpha other pi other 3_10 other other other other other
    .help = Type of helix, defaults to alpha (1).  Only alpha, pi, and 3_10 \
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
}
sheet
  .multiple = True
  .optional = True
{
  first_strand = None
    .type = str
    .style = bold selection
  strand
    .multiple = True
    .optional = False
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

sec_str_master_phil_str = """
input
  .style = box auto_align
{
%s
  find_automatically = None
    .type = bool
    .style = bold tribool
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
  group_helix_backbone = False
    .style = bool
}
%s
""" % (iotbx.pdb.secondary_structure.ss_input_params_str,
       ss_restraint_params_str, ss_group_params_str)

sec_str_master_phil = libtbx.phil.parse(sec_str_master_phil_str)

use_resid_range = False # XXX: for debugging purposes only

class _annotation (oop.injector, iotbx.pdb.secondary_structure.annotation) :
  def as_restraint_groups (self, log=sys.stderr, prefix_scope="") :
    phil_strs = []
    for helix in self.helices :
      helix_phil = helix.as_restraint_group(log, prefix_scope)
      if helix_phil is not None :
        phil_strs.append(helix_phil)
    for sheet in self.sheets :
      sheet_phil = sheet.as_restraint_group(log, prefix_scope)
      if sheet_phil is not None :
        phil_strs.append(sheet_phil)
    return "\n".join(phil_strs)

class _pdb_helix (oop.injector, iotbx.pdb.secondary_structure.pdb_helix) :
  def as_restraint_group (self, log=sys.stderr, prefix_scope="") :
    if self.start_chain_id != self.end_chain_id :
      print >> log, "Helix chain ID mismatch: starts in %s, ends in %s" % (
        self.start_chain_id, self.end_chain_id)
      return None
    if use_resid_range :
      resid_start = "%d%s" % (self.start_resseq, self.start_icode)
      resid_end = "%d%s" % (self.end_resseq, self.end_icode)
      sele = "chain '%s' and resid %s through %s" % (self.start_chain_id,
        resid_start, resid_end)
    else :
      sele = "chain '%s' and resseq %d:%d" % (self.start_chain_id,
        self.start_resseq, self.end_resseq)
    if prefix_scope != "" and not prefix_scope.endswith(".") :
      prefix_scope += "."
    rg = """\
%shelix {
  selection = %s
  helix_class = %d
}""" % (prefix_scope, sele, self.helix_class)
    return rg

class _pdb_sheet (oop.injector, iotbx.pdb.secondary_structure.pdb_sheet) :
  def as_restraint_group (self, log=sys.stderr, prefix_scope="") :
    if len(self.strands) == 0 :
      return None
    selections = []
    senses = []
    reg_curr = []
    reg_prev = []
    for (strand,registration) in zip(self.strands, self.registrations) :
      if use_resid_range :
        resid_start = "%d%s" % (strand.start_resseq, strand.start_icode)
        resid_end = "%d%s" % (strand.end_resseq, strand.end_icode)
        sele = "chain '%s' and resid %s through %s" % (strand.start_chain_id,
          resid_start, resid_end)
      else :
        sele = "chain '%s' and resseq %d:%d" % (strand.start_chain_id,
          strand.start_resseq, strand.end_resseq)
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
        sele_base = "chain '%s' and resid %s"
        resid_curr = "%d%s" % (registration.cur_resseq,registration.cur_icode)
        resid_prev = "%d%s" % (registration.prev_resseq,registration.prev_icode)
        reg_curr.append(sele_base % (registration.cur_chain_id,resid_curr))
        reg_prev.append(sele_base % (registration.prev_chain_id,resid_prev))
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
    selection = %s
    sense = %s
    bond_start_current = %s
    bond_start_previous = %s
  }""" % (sele, sense, curr, prev))
      n += 1
    assert first_strand is not None
    if prefix_scope != "" and not prefix_scope.endswith(".") :
      prefix_scope += "."
    phil_str = """
%ssheet {
  first_strand = %s
%s
}""" % (prefix_scope, first_strand, "\n".join(strands))
    return phil_str

class hydrogen_bond_table (object) :
  def __init__ (self, bonds, distance, sigma, slack) :
    adopt_init_args(self, locals())
    self.flag_use_bond = flex.bool(bonds.size(), True)

  def filter_h_bonds (self, params, pdb_hierarchy, log=sys.stderr) :
    atoms = pdb_hierarchy.atoms()
    if params.remove_outliers and self.distance[i] > distance_ideal :
      pass

  def get_final_bonds (self) :
    if self.flag_use_bond.count(True) == 0 :
      raise Sorry("No bonds meeting usability criteria.")
    for i, (donor_i_seq, acceptor_i_seq) in enumerate(self.bonds) :
      yield group_args(donor_i_seq=donor_i_seq,
        acceptor_i_seq=acceptor_i_seq,
        sigma=self.sigma[i],
        slack=self.slack[i],
        distance=self.distance[i])

def hydrogen_bonds_from_selections (
    pdb_hierarchy,
    params,
    log=sys.stderr) :
  bond_i_seqs = shared.stl_set_unsigned()
  sigmas = flex.double()
  slacks = flex.double()
  atoms = pdb_hierarchy.atoms()
  sites = atoms.extract_xyz()
  has_bond = flex.bool(sites.size(), False)
  selection_cache = pdb_hierarchy.atom_selection_cache()
  isel = selection_cache.iselection
  donor_name = "H"
  if params.h_bond_restraints.substitute_n_for_h :
    donor_name = "N"
  if params.h_bond_restraints.restrain_helices :
    for helix in params.helix :
      helix_class = int(helix.helix_class)
      if helix_class != 1 and params.h_bond_restraints.alpha_only :
        print >> log, "Skipping non-alpha helix (class %d):" % helix_class
        print >> log, "  %s" % helix.selection
        continue
      try :
        donor_isel, acceptor_isel = _donors_and_acceptors(
          base_sele=helix.selection,
          selection_cache=selection_cache,
          donor_name=donor_name,
          ss_type="helix")
      except RuntimeError, e :
        print >> log, str(e)
        continue
      else :
        if helix_class == 1 :
          j = 4
        elif helix_class == 3 :
          j = 5
        elif helix_class == 5 :
          j = 3
        else :
          print >> log, "Don't know bonding for helix class %d." % helix_class
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
      sheet_bonds = []
      prev_strand_sele = sheet.first_strand
      for curr_strand in sheet.strand :
        try :
          if curr_strand.sense == "unknown" :
            raise RuntimeError(("Skipping sheet of unknown sense:\n" +
              "%s") % curr_strand.selection)
          curr_donors, curr_acceptors = _donors_and_acceptors(
            base_sele=curr_strand.selection,
            selection_cache=selection_cache,
            donor_name=donor_name,
            ss_type="sheet")
          prev_donors, prev_acceptors = _donors_and_acceptors(
            base_sele=prev_strand_sele,
            selection_cache=selection_cache,
            donor_name=donor_name,
            ss_type="sheet")
          curr_donor_start, curr_acceptor_start = _donors_and_acceptors(
            base_sele=curr_strand.bond_start_current,
            selection_cache=selection_cache,
            donor_name=donor_name,
            ss_type="sheet")
          prev_donor_start, prev_acceptor_start = _donors_and_acceptors(
            base_sele=curr_strand.bond_start_previous,
            selection_cache=selection_cache,
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
          if len(new_bonds) == 0 :
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
        sigmas.append(sheet.sigma)
        slacks.append(sheet.slack)
  distances = _get_distances(bond_i_seqs, sites)
  return hydrogen_bond_table(bonds=bond_i_seqs,
    distance=distances,
    sigma=sigmas,
    slack=slacks)

def _get_distances (bonds, sites_cart) :
  distances = flex.double(bonds.size(), -1)
  for k, (i_seq, j_seq) in enumerate(bonds) :
    (x1, y1, z1) = sites_cart[i_seq]
    (x2, y2, z2) = sites_cart[j_seq]
    dist = sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
    distances[k] = dist
  return distances

def _donors_and_acceptors (base_sele, selection_cache, donor_name, ss_type) :
  isel = selection_cache.iselection
  donor_sele = "(%s) and (altloc 'A' or altloc ' ') and name %s" % (
    base_sele, donor_name)
  acceptor_sele = "(%s) and (altloc 'A' or altloc ' ') and name O"% base_sele
  donor_isel = isel(donor_sele)
  acceptor_isel = isel(acceptor_sele)
  if donor_isel.size() != acceptor_isel.size() :
    raise RuntimeError("""\
hydrogen_bonds_from_selections: incomplete residues in %s.
  \"%s\" => %d donors
  \"%s\" => %d acceptors""" % (ss_type, donor_sele, donor_isel.size(),
      acceptor_sele, acceptor_isel.size()))
  return donor_isel, acceptor_isel

def _hydrogen_bonds_from_strand_pair (atoms,
    prev_strand_donors,
    prev_strand_acceptors,
    prev_strand_start,
    curr_strand_donors,
    curr_strand_acceptors,
    curr_strand_start,
    sense) :
  assert sense != "unknown"
  assert prev_strand_donors.size() == prev_strand_acceptors.size()
  assert curr_strand_donors.size() == curr_strand_acceptors.size()
  start_bonding = False
  bonds = []
  n_prev_strand = prev_strand_donors.size()
  n_curr_strand = curr_strand_donors.size()
  i = j = None
  for k in prev_strand_donors :
    if k == prev_strand_start :
      i = k
  for k in curr_strand_acceptors :
    if k == curr_strand_start :
      j = k
  if None in [i, j] :
    raise RuntimeError("Can't find start of bonding.")
  if sense == "parallel" :
    while i < n_prev_strand and j > 0 :
      donor1_i_seq = prev_strand_donors[i]
      acceptor1_i_seq = curr_strand_acceptors[j]
      bonds.append(donor1_i_seq, acceptor1_i_seq)
      donor2_i_seq = curr_strand_donors[i]
      acceptor2_i_seq = prev_strand_acceptors[i]
      bonds.append(donor2_i_seq, acceptor2_i_seq)
      i += 2
      j -= 2
  else :
    while i < n_prev_strand and j < n_curr_strand :
      donor1_i_seq = curr_strand_donors[i]
      acceptor1_i_seq = curr_strand_acceptors[j]
      bonds.append(donor1_i_seq, acceptor1_i_seq)
      if (j + 2) < n_curr_strand :
        break
      donor2_i_seq = curr_strand_donors[j+2]
      acceptor2_i_seq = prev_strand_acceptors[i]
      bonds.append(donor2_i_seq, acceptor2_i_seq)
      i += 2
      j += 2
  return bonds

def _hydrogen_bond_from_selection_pair (donor_sele, acceptor_sele,
    selection_cache) :
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
  import iotbx.pdb
  from scitbx.array_family import flex
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
      helix_class=int(helix_params.helix_class),
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
        cur_atom="N", #reg_curr_atom.name,
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

# TODO
def run_ksdssp (file_name, log=sys.stderr) :
  return []

def run (args, out=sys.stdout, log=sys.stderr) :
  pdb_files = []
  sources = []
  for arg in args :
    if os.path.isfile(arg) :
      if iotbx.pdb.is_pdb_file(arg) :
        pdb_files.append(os.path.abspath(arg))
    else :
      sources.append(libtbx.phil.parse(arg))
  secondary_structure = iotbx.pdb.secondary_structure.process_files(pdb_files)
  if secondary_structure is None :
    records = run_ksdssp(pdb_files[0], log=log)
    secondary_structure = iotbx.pdb.secondary_structure.process_files(
      pdb_files=[], records=records, allow_none=False)
  ss_params_str = secondary_structure.as_restraint_groups(log=log,
    prefix_scope="refinement.secondary_structure")
  print >> out, ss_params_str

def get_bonds (file_name, out=sys.stdout, log=sys.stderr) :
  secondary_structure = iotbx.pdb.secondary_structure.process_files([file_name])
  assert secondary_structure is not None
  ss_params_str = secondary_structure.as_restraint_groups(log=sys.stderr)
  pdb_hierarchy = get_pdb_hierarchy([file_name])
  atoms = pdb_hierarchy.atoms()
  sources = [libtbx.phil.parse(ss_params_str)]
  working_phil = sec_str_master_phil.fetch(sources=sources)
  params = working_phil.extract()
  params.h_bond_restraints.substitute_n_for_h = True
  bonds_table = hydrogen_bonds_from_selections (
    pdb_hierarchy,
    params=params,
    log=sys.stderr)
  return bonds_table

def exercise () :
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1ywf.pdb",
    test=os.path.isfile)
  if pdb_file is None :
    print "Skipping"
    return False
  bonds_table = get_bonds(pdb_file)
  for bond in bonds_table.get_final_bonds() :
    donor = atoms[bond.donor_i_seq].fetch_labels()
    acceptor = atoms[bond.acceptor_i_seq].fetch_labels()
    sele1 = "(chain '%s' and resi %s and name %s)" % (donor.chain_id,
      donor.resseq, donor.name)
    sele2 = "(chain '%s' and resi %s and name %s)" % (acceptor.chain_id,
      acceptor.resseq, acceptor.name)
    print "dist %s, %s" % (sele1, sele2)

if __name__ == "__main__" :
  if "--test" in sys.argv :
    exercise()
  else :
    (sys.argv[1:])
