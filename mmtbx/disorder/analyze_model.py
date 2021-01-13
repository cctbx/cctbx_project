
from __future__ import absolute_import, division, print_function
from mmtbx.disorder import backbone
from scitbx.array_family import flex
from scitbx.matrix import col
from libtbx.str_utils import format_value as fv
from libtbx import Auto, slots_getstate_setstate
import math
import sys
from six.moves import range

# XXX in order to make this run in parallel over many PDB IDs, I need to cheat
# slightly and substitute pickle-able objects for the original classes in
# iotbx.pdb.hierarchy.  Note that parent relationships will be lost in the
# process.
class residue_group_proxy(slots_getstate_setstate):
  """Pickle-able stand-in for iotbx.pdb.hierarchy.residue_group."""
  __slots__ = ["resseq", "icode", "_atom_groups", "_id_str", ]
  def __init__(self, residue_group):
    self.resseq = residue_group.resseq
    self.icode = residue_group.icode
    self._id_str = residue_group.id_str()
    self._atom_groups = [ ]
    for ag in residue_group.atom_groups():
      self._atom_groups.append(atom_group_proxy(ag))

  def id_str(self):
    return self._id_str

  def atom_groups(self):
    return self._atom_groups

class atom_group_proxy(slots_getstate_setstate):
  """Pickle-able stand-in for iotbx.pdb.hierarchy.atom_group."""
  __slots__ = [ "resname", "altloc", "_atoms", ]
  def __init__(self, atom_group):
    self.resname = atom_group.resname
    self.altloc = atom_group.altloc
    self._atoms = atoms_proxy(atom_group.atoms())

  def atoms(self):
    return self._atoms

class atoms_proxy(slots_getstate_setstate):
  """
  Pickle-able stand-in for af::shared<atom> array, using the atom_with_labels
  objects as elements.
  """
  __slots__ = [ "_atoms" ]
  def __init__(self, atoms):
    self._atoms = [ a.fetch_labels() for a in atoms ]

  def __getitem__(self, idx):
    return self._atoms[idx]

  def extract_occ(self):
    return flex.double([ a.occ for a in self._atoms ])

  def extract_b(self):
    return flex.double([ a.b for a in self._atoms ])

class disordered_segment(object):
  """
  A group of one or more adjacent residues presumed to form continuous
  alternate conformations.
  """
  def __init__(self, residue_group):
    self.residue_groups = [ ]
    self.outliers = {}
    self.rotamers = {}
    self.ramachandran = {}
    self.backrubs = []
    self.append_residue_group(residue_group)

  def __str__(self):
    if (self.n_residues() == 1):
      return self.residue_groups[0].id_str()
    else :
      return "%s --> %s" % (self.residue_groups[0].id_str(),
         self.residue_groups[-1].id_str())

  def show(self, prefix="", out=sys.stdout):
    if (self.n_residues() == 1):
      print(prefix + "Segment: 1 residue (%s), %d conformers" % \
        (self.residue_groups[0].id_str(), self.n_confs()), file=out)
    else :
      print(prefix+"Segment: %d residues (%s --> %s), %d conformers" %\
        (self.n_residues(), self.residue_groups[0].id_str(),
         self.residue_groups[-1].id_str(), self.n_confs()), file=out)
    for i_res, rg in enumerate(self.residue_groups):
      print(prefix+"  residue_group=%s" % rg.id_str(), file=out)
      for ag in rg.atom_groups():
        rama = rota = None
        for o in self.ramachandran.get(rg.id_str(), []):
          if (o.altloc == ag.altloc):
            rama = o
            break
        for o in self.rotamers.get(rg.id_str(), []):
          if (o.altloc == ag.altloc):
            rota = o
            break
        print(prefix + "    " + \
          "atom_group=%1s %3s  occ=%.2f phi=%-6s psi=%-6s rot=%-7s" %\
          (ag.altloc,
           ag.resname,
           flex.mean(ag.atoms().extract_occ()),
           fv("%.1f", getattr(rama, "phi", None)),
           fv("%.1f", getattr(rama, "psi", None)),
           getattr(rota, "rotamer_name", None)), file=out)
      if (len(self.backrubs[i_res]) > 0):
        for backrub in self.backrubs[i_res] :
          backrub.show(out=out, prefix=prefix+"    ")
      outliers = self.outliers[rg.id_str()]
      if (len(outliers) > 0):
        print(prefix+"     MolProbity outliers:", file=out)
      for outlier in outliers :
        print(prefix+"       %s: %s" % (type(outlier).__name__,
          str(outlier)), file=out)

  def get_previous_conformer(self, index=0):
    rg = self.residue_groups[-1]
    i_group = 0
    for atom_group in rg.atom_groups():
      if (atom_group.altloc.strip() != ''):
        if (i_group == index):
          return atom_group
        else :
          i_group += 1
    return None

  def is_part_of_segment(self, other,
      ignore_inconsistent_occupancy=False,
      ignore_inconsistent_n_conformers=False,
      max_peptide_bond_distance_within_conformer=2.0):
    """
    Determine whether a residue_group object is part of the same continuous
    disordered segment.  The precise meaning of this can be adjusted depending
    on user preferences; by default a continuous segment must have the same
    number of conformers for each residue, and occupancies must be constrained
    for each conformation.  The latter assumption will probably be violated
    most often.
    """
    other_groups = other.atom_groups()
    assert len(other_groups) >= 2
    if (len(other_groups) != len(self.residue_groups[-1].atom_groups())):
      if (not ignore_inconsistent_n_conformers):
        return False
    i_group = 0
    for atom_group in other_groups :
      if (atom_group.altloc != ''):
        other_atoms = atom_group.atoms()
        prev_group = self.get_previous_conformer(index=i_group)
        if (prev_group is None):
          assert ignore_inconsistent_n_conformers
          break
        i_group += 1
        if (prev_group.altloc != atom_group.altloc):
          return False
        prev_atoms = prev_group.atoms()
        if (prev_atoms[0].occ != other_atoms[0].occ):
          if (not ignore_inconsistent_occupancy):
            return False
        curr_n, prev_c = None, None
        for atom in prev_atoms :
          if (atom.name == " C  "):
            prev_c = atom.xyz
            break
        for atom in other_atoms :
          if (atom.name == " N  "):
            curr_n = atom.xyz
            break
        if (curr_n is None) or (prev_c is None):
          return False
        dist = abs(col(curr_n) - col(prev_c))
        if (dist > max_peptide_bond_distance_within_conformer):
          return False
    return True

  def append_residue_group(self, rg):
    self.residue_groups.append(residue_group_proxy(rg))
    rg_backrubs = backbone.find_backrubs(residue_group=rg)
    self.backrubs.append(rg_backrubs)

  def detect_sequence_disorder(self):
    """
    Find any residue groups with heterogeneous chemical identity.
    """
    disordered = []
    for rg in self.residue_groups :
      resnames = set([ ag.resname.upper() for ag in rg.atom_groups() ])
      if (len(resnames) > 1):
        disordered.append((rg.id_str(), sorted(list(resnames))))
    return disordered

  def n_residues(self):
    return len(self.residue_groups)

  def n_partial_splits(self, join_at_calpha=False):
    """
    Count the number of residues where not all atoms have alternates.
    """
    n_partial = 0
    for residue_group in self.residue_groups :
      for atom_group in residue_group.atom_groups():
        if (atom_group.altloc.strip() == ''):
          if (join_at_calpha):
            for atom in atom_group.atoms():
              if (atom.name == " CA "):
                n_partial += 1
                break
          else :
            n_partial += 1
          break
    return n_partial

  def n_confs(self):
    """
    Count the number of alternate conformations.  Sometimes this may not be
    the same for all residue groups, in which case a list is returned.
    """
    all_n_confs = []
    for residue_group in self.residue_groups :
      all_n_confs.append(0)
      for atom_group in residue_group.atom_groups():
        if (atom_group.altloc.strip() != ''):
          all_n_confs[-1] += 1
    all_n_confs_uniq = set(all_n_confs)
    if (len(all_n_confs_uniq) != 1):
      return sorted(list(all_n_confs_uniq))
    return all_n_confs_uniq.pop()

  def n_confs_max(self):
    n_confs = self.n_confs()
    if isinstance(n_confs, int):
      return n_confs
    return max(n_confs)

  def minimum_atom_group_occupancy(self):
    occ_min = 1.
    for rg in self.residue_groups :
      for ag in rg.atom_groups():
        ag_atoms = ag.atoms()
        total = 0
        n_non_hd = 0
        for atom in ag.atoms():
          if (atom.element.strip() not in ["H", "D"]) and (atom.occ != 0):
            total += atom.occ
            n_non_hd += 1
        if (total != 0):
          occ_mean_ag = total / n_non_hd
          occ_min = min(occ_min, occ_mean_ag)
    return occ_min

  def get_all_conformer_distances(self, backbone=None):
    n_confs = self.n_confs()
    assert isinstance(n_confs, int)
    pairwise_distances = []
    for i_conf in range(n_confs - 1):
      indices = [i_conf, i_conf + 1]
      pairwise_distances.append(self.get_conformer_distances(
        conformer_indices=indices,
        backbone=backbone))
    return pairwise_distances

  def get_conformer_distances(self,
      conformer_indices=Auto,
      backbone=None):
    """
    Calculate the distances between atoms in the specified pair of conformers
    (must be present for all residue groups).
    """
    # XXX the way this is handled is somewhat clumsy, but necessary because
    # there is no requirement that atom groups have the same number of atoms or
    # even the same chemical identity (although they are assumed to be amino
    # acids)
    distances = []
    for rg in self.residue_groups :
      i_ag = 0
      atom_groups = rg.atom_groups()
      if (conformer_indices is Auto):
        if (atom_groups[0].altloc.strip() == ''):
          if (len(atom_groups) <= 2):
            continue
          else :
            conformer_indices = (1,2)
        else :
          conformer_indices = (0,1)
      else :
        assert (len(conformer_indices) == 2)
      ag1 = rg.atom_groups()[conformer_indices[0]]
      ag2 = rg.atom_groups()[conformer_indices[1]]
      if ((ag1.altloc.strip() == '') and
          (conformer_indices[0] == 0) and
          (len(atom_groups) >= 3)):
        ag1 = rg.atom_groups()[conformer_indices[0]+1]
        ag2 = rg.atom_groups()[conformer_indices[1]+1]
      for atom1 in ag1.atoms():
        name = atom1.name.strip()
        element = atom1.element.upper().strip()
        if (element in ["H","D"]):
          continue
        if (backbone is not None):
          if (((backbone) and (not name in ["C","CA","N","O"])) or
              ((not backbone) and (name in ["C","CA","N","O"]))):
            continue
        for atom2 in ag2.atoms():
          if (atom1.name == atom2.name):
            distances.append(abs(col(atom1.xyz) - col(atom2.xyz)))
    return distances

  def max_distance_between_conformers(self, backbone=None):
    paired_distances = self.get_all_conformer_distances(backbone=backbone)
    paired_max = []
    for distances in paired_distances :
      if (len(distances) > 0):
        paired_max.append(max(distances))
    if (len(paired_max) > 0):
      return max(paired_max)
    return None

  def max_rmsd_between_conformers(self, backbone=None):
    paired_distances = self.get_all_conformer_distances(backbone=backbone)
    rmsd_max = None
    for distances in paired_distances :
      if (len(distances) == 0):
        continue
      rmsd = math.sqrt(sum([ dxyz**2 for dxyz in distances]) / len(distances))
      if (rmsd_max is None) or (rmsd > rmsd_max):
        rmsd_max = rmsd
    return rmsd_max

  def extract_validation_results(self, multi_criterion):
    """
    Find the matching validation result objects from the multi-criterion
    object (see mmtbx/validation/molprobity/__init__.py).
    """
    for rg in self.residue_groups :
      self.outliers[rg.id_str()] = []
      self.rotamers[rg.id_str()] = []
      self.ramachandran[rg.id_str()] = []
      results = multi_criterion.get_residue_group_data(rg)
      for result in results.outliers :
        if result.is_outlier():
          self.outliers[rg.id_str()].append(result)
        if type(result).__name__ == "rotamer" :
          self.rotamers[rg.id_str()].append(result)
        elif type(result).__name__ == "ramachandran" :
          self.ramachandran[rg.id_str()].append(result)

  def n_rotamer_changes(self, resname=None):
    n_changes = 0
    for rg in self.residue_groups :
      resnames = set([ ag.resname.upper() for ag in rg.atom_groups() ])
      if (len(resnames) > 1):
        continue
      elif (resname is not None) and (resnames.pop() != resname.upper()):
        continue
      rotamers = set(self.rotamers.get(rg.id_str(), []))
      if (len(rotamers) > 1):
        n_changes += 1
    return n_changes

  def find_peptide_flips(self, angle_cutoff=150):
    residues_and_angles = []
    for rg in self.residue_groups :
      peptide_angle = carbonyl_oxygen_angle(rg)
      if peptide_angle is not None and peptide_angle >= angle_cutoff:
        residues_and_angles.append((rg.id_str(), peptide_angle))
    return residues_and_angles

  def n_cbeta_outliers(self):
    return self.n_outliers_of_type(analysis_type='cbeta')

  def n_outliers_of_type(self, analysis_type):
    n_outliers = 0
    for rg in self.residue_groups :
      results = self.outliers.get(rg.id_str(), [])
      for result in results :
        if (type(result).__name__ == analysis_type) and result.is_outlier():
          n_outliers += 1
    return n_outliers

#-----------------------------------------------------------------------
# utility methods
def is_joined_at_calpha(residue_group):
  for atom_group in residue_group.atom_groups():
    if (atom_group.altloc.strip() == ''):
      for atom in atom_group.atoms():
        if (atom.name == " CA "):
          return True
  return False

def carbonyl_oxygen_angle(residue_group):
  """
  Calculate angles between carbonyl oxygen (C=O) bonds in each pair of atom
  groups, and return the maximum value (or None if fewer than two such bonds
  are found).
  """
  c_o_vectors = []
  for atom_group in residue_group.atom_groups():
    c_xyz = o_xyz = None
    for atom in atom_group.atoms():
      if (atom.name.strip() == "O"):
        o_xyz = col(atom.xyz)
      elif (atom.name.strip() == "C"):
        c_xyz = col(atom.xyz)
    if (not None in [c_xyz, o_xyz]):
      c_o_vectors.append(c_xyz - o_xyz)
  if (len(c_o_vectors) >= 2):
    angles = []
    i_ag = 0
    while (i_ag < len(c_o_vectors) - 1):
      angles.append(c_o_vectors[i_ag].angle(c_o_vectors[i_ag+1], deg=True))
      i_ag += 1
    return max(angles)
  return None

def only_amide_hydrogen_split(residue_group):
  """
  Detect cases where the only alternate conformation is for the amide hydrogen,
  presumably because the previous residue was split and Reduce was used to
  add hydrogens.  These residues are ignored in our analyses.
  """
  for atom in residue_group.atoms():
    labels = atom.fetch_labels()
    if (labels.altloc.strip() != '') and (atom.name != " H  "):
      return False
  return True

# XXX unused?
def get_nconfs(pdb_hierarchy):
  """
  Count the number of conformers in a structure.
  """
  if (len(pdb_hierarchy.models()) > 1):
    n_confs = -1 # multiple MODELs aren't handled
  else :
    for chain in pdb_hierarchy.only_model().chains():
      if (chain.is_protein()):
        confs = chain.conformers()
        if (len(confs) > n_confs):
          n_confs = len(confs)
  return n_confs

#-----------------------------------------------------------------------
class process_residue_groups(object):
  def __init__(self, chain,
      multi_criterion_validation=None,
      ignore_inconsistent_occupancy=False,
      log=sys.stdout):
    self.segments = []
    self.chain_id = chain.id
    self.n_residue_groups = 0
    self.n_disordered = 0
    self.residue_counts = {}
    self.disordered_residue_counts = {}
    assert chain.is_protein()
    segment = None
    for residue_group in chain.residue_groups():
      self.n_residue_groups += 1
      atom_groups = residue_group.atom_groups()
      resname_1 = atom_groups[0].resname
      if (not resname_1 in self.residue_counts):
        self.residue_counts[resname_1] = 0
      self.residue_counts[resname_1] += 1
      if (len(atom_groups) > 1):
        self.n_disordered += 1
        if only_amide_hydrogen_split(residue_group):
          print("    residue %s only has alt. confs. for H" % \
            residue_group.id_str(), file=log)
          segment = None
          continue
        else :
          if (not resname_1 in self.disordered_residue_counts):
            self.disordered_residue_counts[resname_1] = 0
          self.disordered_residue_counts[resname_1] += 1
          if (segment is None):
            segment = disordered_segment(residue_group)
            self.segments.append(segment)
          else :
            if segment.is_part_of_segment(other=residue_group,
                ignore_inconsistent_occupancy=ignore_inconsistent_occupancy):
              segment.append_residue_group(residue_group)
            else :
              segment = disordered_segment(residue_group)
              self.segments.append(segment)
      else :
        segment = None
    if (multi_criterion_validation is not None):
      for segment in self.segments :
        segment.extract_validation_results(multi_criterion_validation)

  def show(self, prefix="", out=sys.stdout):
    print(prefix+"Chain '%s': %d residues, %d disordered" % (
      self.chain_id, self.n_residue_groups, self.n_disordered), file=out)
    for segment in self.segments :
      segment.show(out=out, prefix=prefix+"  ")

class process_pdb_hierarchy(object):
  def __init__(self, pdb_hierarchy,
      validation,
      ignore_inconsistent_occupancy=False,
      log=sys.stdout):
    self.chains = []
    self.n_residue_groups = 0
    self.n_disordered = 0
    self.sequence_disorder = []
    self.n_rama_outliers = validation.ramalyze.n_outliers
    self.n_rota_outliers = validation.rotalyze.n_outliers
    self.n_cbeta_outliers = validation.cbetadev.n_outliers
    multi_criterion_validation = None
    if (validation is not None):
      multi_criterion_validation = validation.as_multi_criterion_view()
    for chain in pdb_hierarchy.only_model().chains():
      if (chain.is_protein()):
        print("  processing chain '%s'" % chain.id, file=log)
        chain_info = process_residue_groups(chain=chain,
          multi_criterion_validation=multi_criterion_validation,
          ignore_inconsistent_occupancy=ignore_inconsistent_occupancy,
          log=log)
        self.chains.append(chain_info)
        self.n_residue_groups += chain_info.n_residue_groups
        self.n_disordered += chain_info.n_disordered
        for segment in chain_info.segments :
          self.sequence_disorder.extend(segment.detect_sequence_disorder())
      else :
        print("  skipping non-protein chain '%s'" % chain.id, file=log)
    # TODO post-analysis

  @property
  def segments(self):
    for chain in self.chains :
      for segment in chain.segments :
        yield segment

  def max_rmsd_between_conformers(self, backbone=None):
    rmsd_max = segment_max = None
    for segment in self.segments :
      rmsd = segment.max_rmsd_between_conformers(backbone=backbone)
      if (rmsd_max is None) or (rmsd > rmsd_max):
        rmsd_max = rmsd
        segment_max = segment
    return rmsd_max, segment_max

  def max_distance_between_conformers(self, backbone=None):
    dist_max = segment_max = None
    for segment in self.segments :
      dist = segment.max_distance_between_conformers(backbone=backbone)
      if (dist_max is None) or (dist > dist_max):
        dist_max = dist
        segment_max = segment
    return dist_max, segment_max

  def show(self, out=sys.stdout, verbose=True):
    print("", file=out)
    print("Overall: %d protein chain(s)" % len(self.chains), file=out)
    print("         %d residues" % self.n_residue_groups, file=out)
    print("         %d disorered in %d segments" % (self.n_disordered,
      sum([ len(c.segments) for c in self.chains ])), file=out)
    if (len(self.sequence_disorder) > 0):
      print("%d heterogeneous residues:" % len(self.sequence_disorder), file=out)
      for rg_id, resnames in self.sequence_disorder :
        print("  %s (%s)" % (rg_id, ",".join(resnames)))
    n_rotamer_changes = n_cbeta_dev = n_partial_splits = 0
    peptide_flips = []
    for segment in self.segments :
      n_rotamer_changes += segment.n_rotamer_changes()
      n_cbeta_dev += segment.n_cbeta_outliers()
      n_partial_splits += segment.n_partial_splits(join_at_calpha=True)
      peptide_flips.extend(segment.find_peptide_flips())
    print("%d disordered residues have multiple rotamers" % \
      n_rotamer_changes, file=out)
    if (n_partial_splits > 0):
      print("%d disordered residues have a single C-alpha atom" % \
        n_partial_splits, file=out)
    if (n_cbeta_dev > 0):
      print("%d disordered residues have C-beta deviations" % \
        n_cbeta_dev, file=out)
    if (len(peptide_flips) > 0):
      print("%d apparent peptide flips:", file=out)
      for residue_id_str, angle in peptide_flips :
        print("  %s (angle=%.1f)" % (residue_id_str, angle), file=out)
    # distances and RMSDs
    rmsd_max, segment_max = self.max_rmsd_between_conformers()
    rmsd_mc_max, segment_mc_max = self.max_rmsd_between_conformers(
      backbone=True)
    assert (rmsd_max is not None)
    print("Max. RMSD between conformers:", file=out)
    print("  %6.3f (%s) [all non-H atoms]" % (rmsd_max, segment_max), file=out)
    if (rmsd_mc_max is not None):
      print("  %6.3f (%s) [backbone only]" %(rmsd_mc_max,
        segment_mc_max), file=out)
    dist_max, segment_max = self.max_distance_between_conformers()
    dist_mc_max, segment_mc_max = self.max_distance_between_conformers(
      backbone=True)
    assert (dist_max is not None)
    print("Max. distance between conformers:", file=out)
    print("  %6.3f (%s) [all non-H atoms]" % (dist_max, segment_max), file=out)
    if (dist_mc_max is not None):
      print("  %6.3f (%s) [backbone only]" %(dist_mc_max,
        segment_mc_max), file=out)
    # verbose output
    if (verbose):
      for chain in self.chains :
        chain.show(out=out)
    else :
      print("Run with --verbose to show per-residue results.", file=out)
    print("", file=out)
