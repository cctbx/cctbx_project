from __future__ import division
import iotbx.phil
from libtbx.utils import Sorry
from math import sqrt
from cctbx import geometry_restraints
from cctbx.array_family import flex
import sys

helix_group_params_str = """
helix
  .multiple = True
  .optional = True
  .style = noauto
{
  serial_number = None
    .type = int
    .optional = True
  helix_identifier = None
    .type = str
    .optional = True
  enabled = True
    .type = bool
    .help = Restrain this particular helix
  selection = None
    .type = atom_selection
    .style = bold
  helix_type = *alpha pi 3_10 unknown
    .type = choice
    .help = Type of helix, defaults to alpha.  Only alpha, pi, and 3_10 \
      helices are used for hydrogen-bond restraints.
    .style = bold
  sigma = 0.05
    .type = float
  slack = 0
    .type = float
  top_out = False
    .type = bool
  hbond
    .multiple = True
    .optional = True
    .style = noauto
  {
    donor = None
      .type = atom_selection
    acceptor = None
      .type = atom_selection
  }
}"""

sheet_group_params_str = """
sheet
  .multiple = True
  .optional = True
  .style = noauto
{
  enabled = True
    .type = bool
    .help = Restrain this particular sheet
  first_strand = None
    .type = atom_selection
    .style = bold
  sheet_id = None
    .type = str
  strand
    .multiple = True
    .optional = True
  {
    selection = None
      .type = atom_selection
      .style = bold
    sense = parallel antiparallel *unknown
      .type = choice
      .style = bold
    bond_start_current = None
      .type = atom_selection
      .style = bold
    bond_start_previous = None
      .type = atom_selection
      .style = bold
  }
  sigma = 0.05
    .type = float
  slack = 0
    .type = float
  top_out = False
    .type = bool
  backbone_only = False
    .type = bool
    .help = Only applies to rigid-body groupings, and not H-bond restraints \
      which are already backbone-only.
  hbond
    .multiple = True
    .optional = True
    .style = noauto
  {
    donor = None
      .type = atom_selection
    acceptor = None
      .type = atom_selection
  }
}
"""

master_helix_phil = iotbx.phil.parse(helix_group_params_str)

def _create_hbond_proxy (
    acceptor_atoms,
    donor_atoms,
    hbond_counts,
    distance_ideal,
    distance_cut,
    remove_outliers,
    weight=1.0,
    sigma=None,
    slack=None,
    top_out=False,
    log=sys.stdout) :
  assert sigma is not None
  assert slack is not None
  donors = []
  acceptors = []
  for atom in acceptor_atoms :
    if (atom.name == ' O  ') :
      acceptors.append(atom)
  for atom in donor_atoms :
    if (atom.name == ' N  ') :
      donors.append(atom)
  if len(donors) > 0 and len(acceptors) > 0:
    # make pairs of connecting atoms
    donor_acceptor_pairs = []
    if len(donors) == 1:
      for acc in acceptors:
        donor_acceptor_pairs.append((donors[0], acc))
    elif len(donors) == 2 and len(acceptors) == 2:
      donor_acceptor_pairs.append((donors[0], acceptors[0]))
      donor_acceptor_pairs.append((donors[1], acceptors[1]))
    elif len(donors) == 2 and len(acceptors) == 1:
      for donor in donors:
        donor_acceptor_pairs.append((donor, acceptors[0]))
    result = []
    for donor, acceptor in donor_acceptor_pairs:
      # print "  linking:", donor.id_str(), acceptor.id_str()
      if (hbond_counts[donor.i_seq] > 0) :
        print >> log, "      WARNING: donor atom is already bonded, skipping"
        print >> log, "    %s" % donor_labels.id_str()
        return None
      elif (hbond_counts[acceptor.i_seq] > 0) :
        print >> log, "      WARNING: acceptor atom is already bonded, skipping"
        print >> log, "    %s" % acceptor_labels.id_str()
        return None
      if (remove_outliers) and (distance_cut > 0) :
        dist = donor.distance(acceptor)
        if (dist > distance_cut) :
          print >> log, "      removed outlier: %.3fA  %s --> %s (cutoff:%.3fA)"%(
              dist, donor.id_str(), acceptor.id_str(), distance_cut)
          return None
      limit = -1
      if (top_out) :
        limit = (distance_cut - distance_ideal)**2 * weight/(sigma**2)
        print "limit: %.2f" % limit
      proxy = geometry_restraints.bond_simple_proxy(
        i_seqs=(donor.i_seq, acceptor.i_seq),
        distance_ideal=distance_ideal,
        weight=weight/(sigma ** 2),
        slack=slack,
        top_out=top_out,
        limit=limit,
        origin_id=1)
      result.append(proxy)
    return result
  else :
    print >> log, "WARNING: missing atoms!"
    return None

def create_helix_hydrogen_bond_proxies (
    params,
    pdb_hierarchy,
    selection_cache,
    weight,
    hbond_counts,
    distance_ideal,
    distance_cut,
    remove_outliers,
    log=sys.stdout) :
  assert (not None in [distance_ideal, distance_cut])
  generated_proxies = geometry_restraints.shared_bond_simple_proxy()
  helix_class = params.helix_type
  if helix_class == "alpha" :
    helix_step = 4
  elif helix_class == "pi" :
    helix_step = 5
  elif helix_class == "3_10" :
    helix_step = 3
  else :
    print >> log, "  Don't know bonding for helix class %s." % helix_class
    return generated_proxies
  try :
    helix_selection = selection_cache.selection(params.selection)
  except Exception, e :
    print >> log, str(e)
    return generated_proxies
  assert (helix_step in [3, 4, 5])
  helix_rgs = _get_residue_groups_from_selection(pdb_hierarchy, helix_selection)
  i = 0
  while i < len(helix_rgs)-helix_step:
    if helix_rgs[i+helix_step].atom_groups()[0].resname.strip() == "PRO":
      print >> log, "      Proline residue: %s - end of helix" % \
        (helix_rgs[i+helix_step].id_str())
      i += 3
      continue # XXX is this safe?
    proxies = _create_hbond_proxy(
      acceptor_atoms=helix_rgs[i].atoms(),
      donor_atoms=helix_rgs[i+helix_step].atoms(),
      hbond_counts=hbond_counts,
      distance_ideal=distance_ideal,
      distance_cut=distance_cut,
      remove_outliers=remove_outliers,
      weight=weight,
      sigma=params.sigma,
      slack=params.slack,
      log=log)
    if proxies is not None:
      for proxy in proxies:
        generated_proxies.append(proxy)
    i += 1
  return generated_proxies

def create_sheet_hydrogen_bond_proxies (
    sheet_params,
    pdb_hierarchy,
    weight,
    hbond_counts,
    distance_ideal,
    distance_cut,
    remove_outliers,
    log=sys.stdout) :
  assert (not None in [distance_ideal, distance_cut])
  cache = pdb_hierarchy.atom_selection_cache()
  prev_strand = sheet_params.first_strand
  prev_selection = cache.selection(prev_strand)
  prev_rgs = _get_residue_groups_from_selection(
      pdb_hierarchy=pdb_hierarchy,
      bool_selection=prev_selection)
  n_proxies = 0
  k = 0
  generated_proxies = geometry_restraints.shared_bond_simple_proxy()
  while k < len(sheet_params.strand) :
    curr_strand = sheet_params.strand[k]
    curr_selection = cache.selection(curr_strand.selection)
    curr_start = None
    prev_start = None
    if curr_strand.bond_start_current is not None:
      curr_start = cache.selection(curr_strand.bond_start_current)
    if curr_strand.bond_start_previous is not None:
      prev_start = cache.selection(curr_strand.bond_start_previous)
    curr_rgs = _get_residue_groups_from_selection(
        pdb_hierarchy=pdb_hierarchy,
        bool_selection=curr_selection)
    i = j = 0
    len_prev_residues = len(prev_rgs)
    len_curr_residues = len(curr_rgs)
    if curr_start is not None and prev_start is not None:
      if curr_start.count(True) < 1 or prev_start.count(True) < 1:
        error_msg = """\
Wrong registration in SHEET record. One of these selections
"%s" or "%s"
yielded zero or several atoms. Possible reason for it is the presence of
insertion codes or alternative conformations for one of these residues or
the .pdb file was edited without updating SHEET records.""" \
% (curr_strand.bond_start_current, curr_strand.bond_start_previous)
        raise Sorry(error_msg)

      current_start_res_is_donor = pdb_hierarchy.atoms().select(curr_start)[0].name.strip() == 'N'
      if (len_curr_residues > 0) and (len_prev_residues > 0) :
        i = _find_start_residue(
          residues=prev_rgs,
          start_selection=prev_start)
        j = _find_start_residue(
          residues=curr_rgs,
          start_selection=curr_start)
        if (i >= 0) and (j >= 0) :
          # move i,j pointers from registration residues to the beginning of
          # beta-strands
          while (1 < i and
              ((1 < j and curr_strand.sense == "parallel") or
              (j < len_curr_residues-2 and curr_strand.sense == "antiparallel"))):
            if curr_strand.sense == "parallel":
              i -= 2
              j -= 2
            elif curr_strand.sense == "antiparallel":
              i -= 2
              j += 2
          if (curr_strand.sense == "parallel") :
            # some tweaking for ensure correct donor assignment
            if i >= 2 and not current_start_res_is_donor:
              i -= 2
              current_start_res_is_donor = not current_start_res_is_donor
            if j >= 2 and current_start_res_is_donor:
              j -= 2
              current_start_res_is_donor = not current_start_res_is_donor
            while (i < len_prev_residues) and (j < len_curr_residues):
              if current_start_res_is_donor:
                donor_residue = curr_rgs[j]
                acceptor_residue = prev_rgs[i]
                i += 2
              else:
                donor_residue = prev_rgs[i]
                acceptor_residue = curr_rgs[j]
                j += 2
              current_start_res_is_donor = not current_start_res_is_donor
              if donor_residue.atom_groups()[0].resname.strip() != "PRO":
                proxies = _create_hbond_proxy(
                    acceptor_atoms=acceptor_residue.atoms(),
                    donor_atoms=donor_residue.atoms(),
                    hbond_counts=hbond_counts,
                    distance_ideal=distance_ideal,
                    distance_cut=distance_cut,
                    remove_outliers=remove_outliers,
                    weight=weight,
                    sigma=sheet_params.sigma,
                    slack=sheet_params.slack,
                    top_out=sheet_params.top_out,
                    log=log)
                if proxies is not None:
                  for proxy in proxies:
                    generated_proxies.append(proxy)
          elif (curr_strand.sense == "antiparallel") :
            while(i < len_prev_residues and j >= 0):
              if (prev_rgs[i].atom_groups()[0].resname.strip() != "PRO") :
                proxies = _create_hbond_proxy(
                  acceptor_atoms=curr_rgs[j].atoms(),
                  donor_atoms=prev_rgs[i].atoms(),
                  hbond_counts=hbond_counts,
                  distance_ideal=distance_ideal,
                  distance_cut=distance_cut,
                  remove_outliers=remove_outliers,
                  weight=weight,
                  sigma=sheet_params.sigma,
                  slack=sheet_params.slack,
                  top_out=sheet_params.top_out,
                  log=log)
                if proxies is not None:
                  for proxy in proxies:
                    generated_proxies.append(proxy)

              if (curr_rgs[j].atom_groups()[0].resname.strip() != "PRO") :
                proxies = _create_hbond_proxy(
                  acceptor_atoms=prev_rgs[i].atoms(),
                  donor_atoms=curr_rgs[j].atoms(),
                  hbond_counts=hbond_counts,
                  distance_ideal=distance_ideal,
                  distance_cut=distance_cut,
                  remove_outliers=remove_outliers,
                  weight=weight,
                  sigma=sheet_params.sigma,
                  slack=sheet_params.slack,
                  top_out=sheet_params.top_out,
                  log=log)
                if proxies is not None:
                  for proxy in proxies:
                    generated_proxies.append(proxy)
              i += 2;
              j -= 2;
          else :
            print >> log, "  WARNING: strand direction not defined!"
            print >> log, "    previous: %s" % prev_strand
            print >> log, "    current: %s" % curr_strand.selection
        else :
          print >> log, "  WARNING: can't find start of bonding for strands!"
          print >> log, "    previous: %s" % prev_strand
          print >> log, "    current: %s" % curr_strand.selection
      else :
        print >> log, "  WARNING: can't find one or more strands!"
        print >> log, "    previous: %s" % prev_strand
        print >> log, "    current: %s" % curr_strand.selection
    k += 1
    prev_strand = curr_strand.selection
    prev_selection = curr_selection
    prev_rgs = curr_rgs
  return generated_proxies

def _find_start_residue (
    residues,
    start_selection) :
  start_i_seqs = start_selection.iselection()
  for i, residue in enumerate(residues) :
    atom_i_seqs = residue.atoms().extract_i_seq()
    if (atom_i_seqs.intersection(start_i_seqs).size() > 0) :
      return i
  return -1

def _get_residue_groups_from_selection(pdb_hierarchy, bool_selection):
  assert isinstance(bool_selection, flex.bool)
  i_seqs = bool_selection.iselection()
  rgs = []
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      chain_i_seqs = flex.size_t(sorted(chain.atoms().extract_i_seq()))
      if (chain_i_seqs.intersection(i_seqs).size() == 0) :
        continue
      for rg in chain.residue_groups():
        if bool_selection[rg.atoms()[0].i_seq]:
          rgs.append(rg)
  return rgs

########################################################################
# ANNOTATION
#
# Won't be used soon
#
class find_helices_simple (object) :
  """
  Identify helical regions, defined as any three or more contiguous residues
  with phi and psi within specified limits:
    -120 < phi < -20
    -80 < psi < -10
  This is much more tolerant of distorted models than KSDSSP, but will still
  miss helices in poor structures.
  """
  def __init__ (self, pdb_hierarchy) :
    self.pdb_hierarchy = pdb_hierarchy
    self._helices = []
    self._current_helix = []
    self.run()

  def process_break (self) :
    if (len(self._current_helix) > 3) :
      self._helices.append(self._current_helix)
    self._current_helix = []

  def push_back (self, chain_id, resseq) :
    if (len(self._current_helix) > 0) :
      last_chain, last_resseq = self._current_helix[-1]
      if (last_chain == chain_id) :
        self._current_helix.append((chain_id,resseq))
    else :
      self._current_helix = [(chain_id, resseq)]

  def run (self) :
    import mmtbx.rotamer
    import iotbx.pdb
    get_class = iotbx.pdb.common_residue_names_get_class
    atoms = self.pdb_hierarchy.atoms()
    sites_cart = atoms.extract_xyz()
    current_helix = []
    for model in self.pdb_hierarchy.models() :
      for chain in model.chains() :
        main_conf = chain.conformers()[0]
        residues = main_conf.residues()
        for i_res in range(len(residues) - 1) :
          residue1 = residues[i_res]
          if (get_class(residue1.resname) == "common_amino_acid") :
            residue0 = None
            if (i_res > 0) :
              residue0 = residues[i_res - 1]
            residue2 = residues[i_res+1]
            resseq1 = residue1.resseq_as_int()
            resseq2 = residue2.resseq_as_int()
            if (residue0 is not None) :
              if ((resseq2 == (resseq1 + 1)) or
                  ((resseq2 == resseq1) and
                   (residue1.icode != residue2.icode))) :
                resseq0 = residue0.resseq_as_int()
                if ((resseq0 == (resseq1 - 1)) or ((resseq0 == resseq1) and
                    (residue0.icode != residue1.icode))) :
                  phi_psi_i_seqs = mmtbx.rotamer.get_phi_psi_indices(
                    prev_res=residue0,
                    residue=residue1,
                    next_res=residue2)
                  if (phi_psi_i_seqs.count(None) > 0) :
                    continue
                  (phi, psi) = mmtbx.rotamer.phi_psi_from_sites(
                    i_seqs=phi_psi_i_seqs,
                    sites_cart=sites_cart)
                  if is_approximately_helical(phi, psi) :
                    self.push_back(chain.id, resseq1)
                    continue
                  else :
                    pass
            self.process_break()

  def build_selections (self) :
    atom_selections = []
    for helix in self._helices :
      chain_id = helix[0][0]
      selection = """chain '%s' and resseq %d:%d""" % (helix[0][0],
        helix[0][1], helix[-1][1])
      atom_selections.append(selection)
    return atom_selections

  def as_restraint_group_phil (self) :
    phil_strs = []
    for selection in self.build_selections() :
      helix_str = """helix {\n  selection = "%s"\n}""" % selection
      phil_strs.append(helix_str)
    if (len(phil_strs) > 0) :
      master_phil = iotbx.phil.parse(helix_group_params_str)
      helix_phil = iotbx.phil.parse("\n".join(phil_strs))
      return master_phil.fetch(source=helix_phil)
    return None

  def as_restraint_groups (self) :
    helix_phil = self.as_restraint_group_phil()
    if (helix_phil is not None) :
      return helix_phil.extract()
    return None

  def as_pdb_records (self) :
    pass

  def show (self, out=sys.stdout) :
    if (len(self._helices) == 0) :
      print >> out, "No recognizable helices."
    else :
      print >> out, "%d helix-like regions found:" % len(self._helices)
    for selection in self.build_selections() :
      print >> out, "  %s" % selection

def is_approximately_helical (phi, psi) :
  if (-120 < phi < -20) and (-80 < psi < -10) :
    return True
  return False

# FIXME
def _find_strand_bonding_start (atoms,
    prev_strand_donors,
    prev_strand_acceptors,
    curr_strand_donors,
    curr_strand_acceptors,
    sense,
    max_distance_cutoff=4.5) :
  assert sense != "unknown"
  assert prev_strand_donors.size() == prev_strand_acceptors.size()
  assert curr_strand_donors.size() == curr_strand_acceptors.size()
  sites_cart = atoms.extract_xyz()
  min_dist = max_distance_cutoff
  best_pair = (None, None)
  for donor_i_seq in prev_strand_donors :
    for acceptor_j_seq in curr_strand_acceptors :
      (x1, y1, z1) = sites_cart[donor_i_seq]
      (x2, y2, z2) = sites_cart[acceptor_j_seq]
      dist = sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
      if (dist < min_dist) :
        best_pair = (donor_i_seq, acceptor_i_seq)
  return best_pair
