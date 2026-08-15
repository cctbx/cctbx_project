
from __future__ import absolute_import, division, print_function
from libtbx import easy_pickle
from libtbx import group_args
from libtbx.utils import Sorry
from collections import defaultdict
from iotbx.pdb import common_residue_names_get_class
import os.path
import math
import sys
from six.moves import range

def export_ramachandran_distribution(n_dim_table, scale_factor=0.25):
  """
  Convert a MolProbity Ramachandran distribution to a format suitable for
  display using matplotlib (see wxtbx.plots).
  """
  import numpy
  z = n_dim_table.lookupTable
  z_grid = [ [ z[i + (j * 180)] for j in range(180) ]
                          for i in range(180) ]
  npz = numpy.array(z_grid)
  return npz ** scale_factor

def export_rotamer_distribution(n_dim_table, scale_factor=0.5):
  """
  Convert a MolProbity rotamer distribution to a format suitable for
  display using matplotlib (see wxtbx.plots).  Will reduce dimensionality to
  2 if necessary.
  """
  import numpy
  z = n_dim_table.lookupTable
  n_dim = n_dim_table.nDim
  assert n_dim >= 2
  x_offset = 1
  for nbins in n_dim_table.nBins[1:] :
    x_offset *= nbins
  y_width = 1
  if n_dim > 2 :
    for nbins in n_dim_table.nBins[2:] :
      y_width *= nbins
  z_grid = [ [] for i in range(n_dim_table.nBins[1]) ]
  for i in range(n_dim_table.nBins[0]):
    for j in range(n_dim_table.nBins[1]):
      z_total = 0
      for k in range(y_width):
        z_total += z[(i * x_offset) + (j * y_width) + k]
      z_grid[j].append(z_total)
  npz = numpy.array(z_grid)
  return npz ** scale_factor

def get_rotarama_data(residue_type=None, pos_type=None, db="rama",
    convert_to_numpy_array=False):
  from mmtbx.rotamer import ramachandran_eval
  from mmtbx.rotamer.rotamer_eval import find_rotarama_data_dir
  # backwards compatibility
  if (pos_type == "proline") : pos_type = "trans-proline"
  if (pos_type == "prepro") : pos_type = "pre-proline"
  assert (pos_type in ["general", "cis-proline", "trans-proline", "glycine",
    "isoleucine or valine", "pre-proline",None])
  assert (db in ["rama", "rota"])
  assert (residue_type is not None) or (pos_type is not None)
  if pos_type is not None :
    residue_type = ramachandran_eval.aminoAcids_8000[pos_type]
  if residue_type.lower() in ["phe", "tyr"] :
    residue_type = "phetyr"
  assert (residue_type is not None)
  rama_data_dir = find_rotarama_data_dir()
  if (db == "rama"):
    pkl_file = "%s8000-%s.pickle" % (db, residue_type)
  else :
    pkl_file = "%s8000-%s.pickle" % (db, residue_type.lower())
  ndt = easy_pickle.load(os.path.join(rama_data_dir, pkl_file))
  if convert_to_numpy_array :
    if (db == "rama"):
      return export_ramachandran_distribution(ndt)
    else :
      return export_rotamer_distribution(ndt)
  else :
    return ndt

def decode_atom_str(atom_id):
  chain_id = atom_id[8:10].strip()
  if (chain_id == ""):
    chain_id = " "
  return group_args(
    name = atom_id[0:4],
    altloc = atom_id[4],
    resname = atom_id[5:8],
    chain_id = chain_id,
    resid = atom_id[10:],
    resseq = atom_id[10:-1].strip())

def find_sequence_mismatches(pdb_hierarchy,
                              sequences,
                              assume_same_order=True,
                              expected_sequence_identity=0.8,
                              log=sys.stdout):
  chains = pdb_hierarchy.models()[0].chains()
  chain_ids = []
  actual_seqs = []
  expected_seqs = []
  if (len(chains) != len(sequences)) or (not assume_same_order):
    print("Can't determine sequence->chain mapping autoamtically", file=log)
    print("Running sequence alignments. . .", file=log)
    from mmtbx.alignment import pairwise_global_wrapper
    for chain in chains :
      chain_seq = chain.as_padded_sequence()
      actual_seqs.append(chain_seq)
      chain_ids.append(chain.id)
      best_identity = 0
      best_sequence = None
      for sequence in sequences :
        pg = pairwise_global_wrapper(chain_seq, sequence)
        identity = pg.calculate_sequence_identity()
        if (identity >= expected_sequence_identity):
          if (identity >= best_identity):
            best_identity = identity
            best_sequence = sequence
      expected_seqs.append(best_sequence)
  mismatches = []

def molprobity_score(clashscore, rota_out, rama_fav):
  """
  Calculate the overall Molprobity score, as described here:
    http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2877634/?tool=pubmed
    http://kinemage.biochem.duke.edu/suppinfo/CASP8/methods.html
  """
  if clashscore is not None and rota_out is not None and rama_fav is not None \
       and (clashscore >= 0) and (rota_out >= 0) and (rama_fav >= 0):
    rama_iffy = 100. - rama_fav
    mpscore = (( 0.426 * math.log(1 + clashscore) ) +
             ( 0.33 * math.log(1 + max(0, rota_out - 1)) ) +
             ( 0.25 * math.log(1 + max(0, rama_iffy - 2)) )) + 0.5
  else :
    return -1 # FIXME prevents crashing on RNA and None in inputs
  return mpscore

from typing import Dict, Any, Optional, List

def _clash_severity(abs_overlap, num_clashes):
    """Map clash overlap magnitude and count to a continuous severity.

    Uses a non-linear power curve on the worst overlap so that small
    clashes near the 0.4 A threshold are mild while large overlaps
    grow sublinearly, staying below the twisted-peptide tier.  The
    count bonus uses log2 with a cap at 4.0, because high clash counts
    typically reflect many atom-atom contacts from a single bad
    pairwise interaction rather than independent problems.

    Severity scale (single clash, approximate):
      0.4 A overlap -> 1.4    (minor)
      0.5 A         -> 2.0    (moderate)
      0.7 A         -> 3.4    (moderate)
      0.9 A         -> 5.0    (significant)
      1.5 A         -> 10.3   (severe, still below twisted peptide)

    Count bonus (log2, capped):
      2 clashes  -> +1.0
      5 clashes  -> +2.3
      16 clashes -> +4.0  (cap reached)
      50 clashes -> +4.0
    """
    severity = max(0.0, (abs_overlap - 0.1) * 10.0) ** 1.3 / 3.0
    if num_clashes > 1:
        severity += min(math.log2(num_clashes), 4.0) * 1.0
    return severity

CBETA_THRESHOLD = 0.25          # Angstroms, MolProbity's C-beta outlier cut
# Anchored so 0.25 A -> 1.0 and 0.50 A -> 4.0, the values the linear form documented.
CBETA_K = 3.0 / math.log(0.50 / 0.25)


def _cbeta_severity(deviation):
    """Map C-beta deviation to a continuous severity.

    Outlier threshold is 0.25 A.  Logarithmic in the deviation, which is the form that
    makes severity proportional to -log(population frequency):

      0.25 A -> 1.0
      0.50 A -> 4.0
      1.00 A -> 7.0
      2.00 A -> 10.0

    The previous form was linear, (deviation - 0.13) * 12.0, and uncapped.  Measured over
    139,418 C-beta evaluations in 204 PDB-REDO structures it reached 33.6, more than twice
    the ceiling of every discrete tier in the metric, which let one C-beta deviation
    outrank a residue carrying five independent problems.  Fitting the survival function
    P(dev >= x) over 0.25-1.99 A, a power law beat both exponential and Gaussian
    (R^2 0.770 against 0.529 and 0.381), so the logarithm is the surprisal-calibrated
    shape rather than a clamp picked to make the numbers behave.  It self-limits as a
    result: the largest deviation in that sample, 2.93 A, scores 11.6.

    That R^2 is the weakest of the three graded populations, and the reason is known: the
    fitted range spans a bimodal distribution (see the suppression note in
    calculate_overall_residue_quality_score), so a single power law is being asked to
    cover two populations.  The shape is right for the first mode; the second is handled
    by suppression rather than by the curve.

    Note the linear form never hit its own documented anchors; it returned 1.44 at 0.25 A,
    not 1.0.  This one does.

    Deviations at or below the threshold return 1.0, because the caller only invokes this
    for residues already flagged as outliers.
    """
    if deviation is None or deviation <= CBETA_THRESHOLD:
        return 1.0
    return 1.0 + CBETA_K * math.log(deviation / CBETA_THRESHOLD)

def _omega_twist_severity(omega):
    """Map a twisted peptide's omega to a continuous severity.

    Every twisted peptide previously scored a flat 15.0, the highest severity in the metric,
    above a chirality inversion. Measured over 1.40 million peptides in 2168 PDB-REDO
    structures, that flat tier was spending the top of the scale on a quantity that barely
    varies: the median twist is 36.8 degrees in deposited models and 36.4 in re-refined ones,
    unchanged across every resolution from 0.8 to 4.0 A, with 43% of all twists inside the
    first five degrees past the threshold and only 8% reaching 60. A flag that fires on
    something that constant carries no information beyond "a twist exists", and discards both
    the 30-90 degree range of the angle itself and the 40-fold variation in how often it
    happens.

    The shape is not fitted, it is the physics. A peptide bond resists rotation because the
    amide is conjugated, and that conjugation is lost as sin^2 of the rotation: zero when
    planar, maximal when perpendicular. So severity follows sin^2 of the twist, which also
    means it will not need re-tuning when the reference population changes.

    Anchored at the two ends that matter:

      twist 30 deg (MolProbity's Twisted border) ->  3.0, the rotamer tier
      twist 45 deg                               ->  7.0
      twist 60 deg                               -> 11.0, above a chirality inversion
      twist 90 deg (perpendicular)               -> 15.0, the old flat value

    The floor is 3.0 rather than 0.0 deliberately. The threshold itself is arbitrary, a bare
    constant in omegalyze.find_omega_type with no derivation, and the population shows no
    feature at 30 degrees: the cumulative share runs 0.145 / 0.102 / 0.075 / 0.060 percent at
    25 / 28 / 30 / 32 degrees, a smooth decay straight through. That argues for continuity
    across the border, which a zero floor would give. But a 31-degree twist is still in the
    top 0.1% of all peptides; it only looks borderline against the 90-degree maximum. Scoring
    the 99.9th percentile of non-planarity as 0.00 would drop half of all twisted peptides
    below the high-triage cut of 2.0 and out of hotspot counts entirely.

    omega is the dihedral in degrees. If it is unavailable the old flat 15.0 is returned, so
    a caller that cannot supply the angle is no worse off than before.
    """
    if omega is None:
        return 15.0
    twist = min(abs(omega), 180.0 - abs(omega))
    frac = (math.sin(math.radians(twist)) ** 2 - 0.25) / 0.75   # 0.25 = sin^2(30 deg)
    return max(3.0, min(15.0, 3.0 + 12.0 * frac))


BOND_ANGLE_THRESHOLD = 4.0      # sigma, MolProbity's bond/angle outlier cut
# Anchored so 4 sigma -> 1.0 and 10 sigma -> 4.0, the values the linear form documented.
BOND_ANGLE_K = 3.0 / math.log(10.0 / 4.0)


def _bond_angle_severity(num_outliers, worst_sigma):
    """Map bond/angle outlier count and worst sigma to a continuous severity.

    Uses the worst deviation in sigma units as the primary signal, logarithmically, with
    the 4-sigma outlier threshold as the floor.

      4 sigma   -> 1.0
      10 sigma  -> 4.0
      30 sigma  -> 7.6
      100 sigma -> 11.5
      164 sigma -> 13.1   (the worst seen in 1.2M bond restraints; still under 15)

    The previous form was a straight line, (|sigma| - 2) * 0.5, whose docstring only
    anchored values to 10 sigma and which had no ceiling past there.  Across the scored
    benchmark it reached 86.6, nearly six times the ceiling of every discrete tier, and
    8% of structures ended up with a residue whose *only* problem was covalent geometry
    sorted above residues carrying five or six independent problems.  That inverts the
    one guarantee the max-plus-quarter-sum combination rule exists to provide.

    The shape is not a clamp, it is what the population says.  Over every restraint in 204
    PDB-REDO structures (1.2M bonds, 1.7M angles), fitting the survival function
    P(sigma >= x) from the 4-sigma threshold out to where the data runs out, a power law
    beat both exponential and Gaussian:

      bond    fitted 4.0 - 52.5 sigma   power law R^2 0.867, exp 0.708, Gaussian 0.605
      angle   fitted 4.0 - 44.8 sigma   power law R^2 0.931, exp 0.719, Gaussian 0.537

    so -ln(frequency) is linear in ln(sigma) and a surprisal-proportional severity goes as
    ln(sigma).  Note this is emphatically NOT the physics: a harmonic restraint says
    quadratic and so does Gaussian surprisal, and both lose to the data, because past
    roughly 10 sigma the number stops measuring strain and starts measuring how badly
    something is wrong.  A 31-sigma and a 170-sigma deviation are the same finding with
    the same fix.

    Bond and angle keep one shared function because both are anchored to the same two
    documented points, not because their populations are identical: the fitted survival
    exponents are 1.40 and 2.19, implying density exponents near 2.4 and 3.2.  Angle
    outliers really do thin out faster.  A per-metric k is the obvious refinement if the
    tiers are ever recalibrated as a set.

    The fit deliberately stops where fewer than 10 observations remain above a point.
    Beyond that the sample is a handful of restraints, and some of them are not geometry
    at all: PDB-REDO 0cyc files occasionally swap a backbone N with a side-chain O or N
    (seen in 4f35 residue 151 across three chains, and 6baq H176), which manufactures a
    150+ sigma "bond" out of an atom-naming error.  Those must not be allowed to set the
    shape of a severity curve.

    The count bonus now matches _clash_severity: log2, capped at 4.0.  It was +0.5 per
    extra outlier with no cap, the same unbounded-growth problem in a second term.

    A caller that cannot supply the sigma gets the floor of 1.0 plus the count bonus.
    """
    if worst_sigma is None or worst_sigma == 0:
        severity = 1.0                       # sigma unavailable; count is all we have
    else:
        s = abs(worst_sigma)
        severity = 1.0 if s <= BOND_ANGLE_THRESHOLD else \
            1.0 + BOND_ANGLE_K * math.log(s / BOND_ANGLE_THRESHOLD)
    if num_outliers > 1:
        severity += min(math.log2(num_outliers), 4.0)
    return severity

def _rna_suite_severity(is_outlier, suiteness):
    """Map an RNA backbone suite to a graded severity.

    A suite is the backbone unit between two riboses, classified by suitename
    into named rotamer-like conformers.  A full outlier ("!!") matches no known
    cluster and is the serious case.  A suite that is assigned a cluster but fits
    it poorly (low suiteness "wannabe") gets a mild disfavored-tier severity, so
    the metric has dynamic range instead of a single binary flag.  suiteness is
    0..1 (higher = better fit).

    Suites are the RNA analogue of protein rotamers (a few % outliers in good
    structures), so the full-outlier tier sits just below the Ramachandran-
    equivalent sugar-pucker tier (see _rna_pucker_severity).
      full outlier             -> 4.0
      assigned, suiteness < 0.3 -> 1.5   (disfavored-tier "wannabe")
      otherwise                -> 0.0
    """
    if is_outlier:
        return 4.0
    if suiteness is not None and suiteness < 0.3:
        return 1.5
    return 0.0

def _rna_pucker_severity(is_delta_outlier, is_epsilon_outlier,
                         is_pucker_outlier=False):
    """Map an RNA sugar-pucker problem to a graded severity.

    A *delta* outlier means the base-phosphate perpendicular (P-perp) distance
    and the delta torsion disagree about the sugar pucker (C2'- vs C3'-endo) -- a
    genuinely wrong, rare, chemically well-defined error, so it sits at the
    Ramachandran tier.  An *epsilon-only* outlier (epsilon torsion out of the
    valid range while the sugar pucker itself is consistent) is a milder backbone
    problem.  Sugar-pucker outliers are rarer than suite outliers, hence weighted
    above them.
      delta outlier (wrong sugar pucker) -> 5.0
      epsilon-only outlier               -> 3.0
      otherwise                          -> 0.0
    is_pucker_outlier is a fallback for the (delta OR epsilon) flag when the typed
    sub-flags are unavailable; treated as the delta (wrong-pucker) tier.
    """
    if is_delta_outlier:
        return 5.0
    if is_epsilon_outlier:
        return 3.0
    if is_pucker_outlier:
        return 5.0
    return 0.0

def calculate_overall_residue_quality_score(
    residue_data: Dict[str, Any],
) -> Optional[float]:
    """
    Calculates a triage priority score for a single residue based on
    validation outliers.

    The score is designed to rank residues by how urgently they need
    attention: higher scores indicate problems that are more likely to
    be genuine modeling errors.  The combination rule ensures that a
    single severe problem (e.g. a twisted peptide) always outranks an
    accumulation of minor issues.

    Combination rule:
      score = worst_severity + 0.25 * sum(remaining_severities)

    Where possible, continuous values (clash overlap in A, C-beta
    deviation in A, bond/angle deviation in sigma) are used rather
    than discrete outlier/not-outlier bins.

    Args:
        residue_data: Dictionary of validation metrics for one residue.
            Recognized keys:
              num_bad_clashes_res (int), worst_clash_overlap (float, negative),
              ramalyze_type (str), ramalyze_category (str),
              rotalyze_category (str),
              is_glycine (bool), is_cbeta_outlier (bool), cbeta_deviation (float),
              cablam_outlier_type (str),
              omega_type (str), omega_dihedral (float, degrees), is_proline (bool),
              num_bond_outliers_res (int), worst_bond_deviation (float, sigma),
              num_angle_outliers_res (int), worst_angle_deviation (float, sigma),
              num_chiral_handedness_res (int), num_chiral_tetrahedral_res (int),
              num_chiral_pseudochiral_res (int),
              num_chiral_outliers_res (int, fallback if typed counts unavailable),
              is_rna_residue (bool),
              is_rna_suite_outlier (bool), rna_suite_suiteness (float, 0..1),
              is_rna_pucker_outlier (bool), is_delta_outlier (bool),
              is_epsilon_outlier (bool).

    Returns:
        The triage score rounded to 1 decimal place (0.0 = no issues,
        higher = worse), or None if no applicable metrics were found.
    """

    severities: List[float] = []
    has_any_metric = False

    def get(key, default=None):
        return residue_data.get(key, default)

    # --- 1. Steric Clashes ---
    num_bad_clashes = get('num_bad_clashes_res', 0)
    if num_bad_clashes > 0:
        has_any_metric = True
        worst_overlap = get('worst_clash_overlap', 0.0)
        abs_overlap = abs(worst_overlap) if worst_overlap else 0.4
        severities.append(_clash_severity(abs_overlap, num_bad_clashes))

    # --- 2. Ramachandran ---
    if get('ramalyze_type') not in ['not_applicable', 'not_evaluated']:
        has_any_metric = True
        if get('ramalyze_category') == 'outlier':
            severities.append(5.0)

    # --- 3. Rotamer ---
    if get('rotalyze_category') not in ['not_evaluated', None]:
        has_any_metric = True
        if get('rotalyze_category') == 'outlier':
            severities.append(3.0)

    # --- 4. C-beta Deviation ---
    # Suppressed when a chirality handedness swap is already flagged on this residue,
    # because then the two metrics are reporting one physical fact, not two findings.
    # C-beta deviation is bimodal: a first mode below 1.0 A (4419 outliers measured, 4.2%
    # chirality-flagged) and, past a gap at 1.0-1.8 A holding 1.0% of the population, a
    # sharp second mode at 1.9-2.3 A of which 98.9% already carries a chirality flag.
    # That second mode is the D-amino-acid population; cctbx names the same phenomenon at
    # the same place, in cbetadev's own exclude_d_peptides guard at dev >= 2.0.  Scoring
    # it twice gave those residues 10.0 for the handedness swap plus 25-33 on top.
    #
    # Only the typed handedness count suppresses.  A tetrahedral-geometry outlier is a
    # genuinely different problem, and the untyped num_chiral_outliers_res fallback cannot
    # distinguish the two, so neither is allowed to silence C-beta.
    if not get('is_glycine', False):
        has_any_metric = True
        if get('is_cbeta_outlier') and not get('num_chiral_handedness_res', 0):
            deviation = get('cbeta_deviation', 0.0) or 0.0
            severities.append(_cbeta_severity(deviation))

    # --- 5. CaBLAM ---
    cablam_type = get('cablam_outlier_type')
    if cablam_type not in ['not_evaluated', None]:
        has_any_metric = True
        if cablam_type in ('outlier', 'ca_geom_outlier'):
            severities.append(3.0)
        elif cablam_type == 'disfavored':
            severities.append(1.0)

    # --- 6. Omega Angle ---
    omega_type = get('omega_type')
    if omega_type not in ['not_applicable', 'not_evaluated', None]:
        has_any_metric = True
        if omega_type == 'twisted':
            severities.append(_omega_twist_severity(get('omega_dihedral')))
        elif omega_type == 'cis' and not get('is_proline'):
            severities.append(8.0)

    # --- 7. Bond Lengths ---
    num_bond_outliers = get('num_bond_outliers_res', 0)
    if num_bond_outliers > 0:
        has_any_metric = True
        worst_bond = get('worst_bond_deviation')
        severities.append(_bond_angle_severity(num_bond_outliers, worst_bond))

    # --- 8. Bond Angles ---
    num_angle_outliers = get('num_angle_outliers_res', 0)
    if num_angle_outliers > 0:
        has_any_metric = True
        worst_angle = get('worst_angle_deviation')
        severities.append(_bond_angle_severity(num_angle_outliers, worst_angle))

    # --- 9. Chirality ---
    # Three distinct types with very different implications:
    #   - Handedness swap: true chirality inversion (e.g. L->D), serious error
    #   - Tetrahedral geometry: distorted geometry, needs attention
    #   - Pseudochiral naming: swapped names on chemically identical atoms
    #     (e.g. VAL CG1/CG2), cosmetic fix only
    n_handedness = get('num_chiral_handedness_res', 0)
    n_tetrahedral = get('num_chiral_tetrahedral_res', 0)
    n_pseudochiral = get('num_chiral_pseudochiral_res', 0)
    n_chiral_total = n_handedness + n_tetrahedral + n_pseudochiral
    if n_chiral_total > 0:
        has_any_metric = True
        if n_handedness > 0:
            severities.append(10.0)
        if n_tetrahedral > 0:
            severities.append(5.0)
        if n_pseudochiral > 0:
            severities.append(1.5)
    elif get('num_chiral_outliers_res', 0) > 0:
        # Fallback if only total count is available
        has_any_metric = True
        severities.append(5.0)

    # --- 10. RNA Suite (graded: full outlier + low-suiteness wannabe) ---
    if get('is_rna_residue', False):
        has_any_metric = True
        suite_sev = _rna_suite_severity(
            get('is_rna_suite_outlier', False),
            get('rna_suite_suiteness'))
        if suite_sev > 0:
            severities.append(suite_sev)

    # --- 11. RNA Pucker (graded: wrong-pucker delta outranks epsilon-only) ---
    if get('is_rna_residue', False):
        has_any_metric = True
        pucker_sev = _rna_pucker_severity(
            get('is_delta_outlier', False),
            get('is_epsilon_outlier', False),
            get('is_rna_pucker_outlier', False))
        if pucker_sev > 0:
            severities.append(pucker_sev)

    if not has_any_metric:
        return None
    if not severities:
        return 0.0

    # Combination rule: worst finding dominates, additional findings
    # contribute at 25% so accumulation of minor issues cannot outrank
    # a single severe problem.
    severities.sort(reverse=True)
    score = severities[0] + 0.25 * sum(severities[1:])
    return round(score, 1)

def use_segids_in_place_of_chainids(hierarchy, strict=False):
  use_segids = False
  for model in hierarchy.models():
    for chain in model.chains():
      if chain.id in [' ', '  ']:
        cur_segid = None
        for atom in chain.atoms():
          # new as of 20150203
          if atom.segid not in ['    ', '']:
            return True
          # It makes no sense to require indentical segID for
          # Chains with blank chainID. This was commented out by BJH on 20150203
          #if cur_segid is None:
          #  cur_segid = atom.segid
          #if atom.segid not in ['    ', '']:
          #  if atom.segid != cur_segid:
          #    if strict:
          #      raise Sorry("Chains with blank chainID may not have multiple"+
          #                  " segID values")
          #    else:
          #      return False
        #if len(cur_segid.strip()) > 0:
        #  use_segids = True
        #else:
        #  return False
  return use_segids

#this function assumes that use_segids_in_place_of_chainids() is True
def get_segid_as_chainid(chain):
  for atom in chain.atoms():
    return atom.segid

def get_rna_backbone_dihedrals(processed_pdb_file,
      geometry=None, pdb_hierarchy=None):
  # at present, this will only return measurements for angles arising from
  # atoms with altloc ' ' or altloc 'A'
  # TO-DO: extend to more alternates JJH 140108
  from cctbx import geometry_restraints
  bb_dihedrals = defaultdict(dict)
  formatted_out = []
  alt_tracker = {}
  if (processed_pdb_file is not None):
    sites_cart = processed_pdb_file.all_chain_proxies.sites_cart
    geometry = processed_pdb_file.geometry_restraints_manager()
    pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  else :
    assert (not None in [geometry, pdb_hierarchy])
    sites_cart = pdb_hierarchy.atoms().extract_xyz()
  dihedral_proxies = geometry.dihedral_proxies
  i_seq_name_hash = build_name_hash(pdb_hierarchy=pdb_hierarchy)

  def is_blank_or_alt_a(proxy):
    for i in proxy.i_seqs:
       alt = i_seq_name_hash[i][4:5]
       if alt not in [' ', 'A']:
         return False
    return True

  for dp in dihedral_proxies:
    atoms = []
    debug_key = ""
    invert_sign = False
    dp.sort_i_seqs()
    for i in dp.i_seqs:
      atoms.append(i_seq_name_hash[i][0:4].strip())
      debug_key = debug_key+i_seq_name_hash[i]
    if len(atoms) != 4:
      continue
    name = match_dihedral_to_name(atoms=atoms)
    #handle dihedral equivalences
    if name == None:
      inverted_atoms = get_inverted_atoms(atoms=atoms, improper=False)
      name = match_dihedral_to_name(atoms=inverted_atoms)
      if name == None:
        inverted_atoms = get_inverted_atoms(atoms=atoms, improper=True)
        name = match_dihedral_to_name(atoms=inverted_atoms)
        if name is not None:
          invert_sign = True
    if (name is not None) and (is_blank_or_alt_a(dp)):
      restraint = geometry_restraints.dihedral(
                                               sites_cart=sites_cart,
                                               proxy=dp)
      key = i_seq_name_hash[dp.i_seqs[1]][4:]
      if alt_tracker.get(key[1:]) is None:
        alt_tracker[key[1:]] = []
      if key[0:1] not in alt_tracker[key[1:]]:
        alt_tracker[key[1:]].append(key[0:1])
      bb_dihedrals[key][name] = restraint.angle_model
      if invert_sign:
        bb_dihedrals[key][name] = bb_dihedrals[key][name] * -1.0
  for key in list(bb_dihedrals.keys()):
    altloc = key[0:1]
    resname = key[1:4]
    chainID = key[4:6]
    resnum = key[6:10]
    i_code = key[10:]
    if 'A' in alt_tracker[key[1:]]:
      if altloc != 'A':
        continue
    if bb_dihedrals[key].get('alpha') is not None:
      alpha = "%.3f" % bb_dihedrals[key]['alpha']
    # FIXME will the lookup below ever actually work?
    elif altloc == 'A' and \
         bb_dihedrals[' '+key[1:]].get('alpha') is not None:
      alpha = "%.3f" % bb_dihedrals[' '+key[1:]]['alpha']
    else:
      alpha = '__?__'
    if bb_dihedrals[key].get('beta') is not None:
      beta = "%.3f" % bb_dihedrals[key]['beta']
    elif altloc == 'A' and \
         bb_dihedrals[' '+key[1:]].get('beta') is not None:
      beta = "%.3f" % bb_dihedrals[' '+key[1:]]['beta']
    else:
      beta = '__?__'
    if bb_dihedrals[key].get('gamma') is not None:
      gamma = "%.3f" % bb_dihedrals[key]['gamma']
    elif altloc == 'A' and \
         bb_dihedrals[' '+key[1:]].get('gamma') is not None:
      gamma = "%.3f" % bb_dihedrals[' '+key[1:]]['gamma']
    else:
      gamma = '__?__'
    if bb_dihedrals[key].get('delta'):
      delta = "%.3f" % bb_dihedrals[key]['delta']
    elif altloc == 'A' and \
         bb_dihedrals[' '+key[1:]].get('delta') is not None:
      delta = "%.3f" % bb_dihedrals[' '+key[1:]]['delta']
    else:
      delta = '__?__'
    if bb_dihedrals[key].get('epsilon'):
      epsilon = "%.3f" % bb_dihedrals[key]['epsilon']
    elif altloc == 'A' and \
         bb_dihedrals[' '+key[1:]].get('epsilon') is not None:
      epsilon = "%.3f" % bb_dihedrals[' '+key[1:]]['epsilon']
    else:
      epsilon = '__?__'
    if bb_dihedrals[key].get('zeta'):
      zeta = "%.3f" % bb_dihedrals[key]['zeta']
    elif altloc == 'A' and \
         bb_dihedrals[' '+key[1:]].get('zeta') is not None:
      zeta = "%.3f" % bb_dihedrals[' '+key[1:]]['zeta']
    else:
      zeta = '__?__'
    eval = "%s:%s:%s:%s:%s:%s:%s:%s:%s:%s:%s:%s:%s" \
           % (" ",
              "1",
              chainID,
              resnum,
              i_code,
              altloc,
              resname,
              alpha,
              beta,
              gamma,
              delta,
              epsilon,
              zeta)
    formatted_out.append(eval)
  formatted_out.sort()
  backbone_dihedrals = ""
  for line in formatted_out:
    backbone_dihedrals += line+'\n'
  return backbone_dihedrals

def get_inverted_atoms(atoms, improper=False):
  temp = []
  if not improper:
    temp.append(atoms[3])
    temp.append(atoms[2])
    temp.append(atoms[1])
    temp.append(atoms[0])
  else:
    temp.append(atoms[3])
    temp.append(atoms[1])
    temp.append(atoms[2])
    temp.append(atoms[0])
  return temp

def match_dihedral_to_name(atoms):
  name = None
  alpha = ["O3'","P","O5'","C5'"]
  beta = ["P","O5'","C5'","C4'"]
  gamma = ["O5'","C5'","C4'","C3'"]
  delta = ["C5'","C4'","C3'","O3'"]
  epsilon = ["C4'","C3'","O3'","P"]
  zeta = ["C3'","O3'","P","O5'"]
  if atoms == alpha:
    name = "alpha"
  elif atoms == beta:
    name = "beta"
  elif atoms == gamma:
    name = "gamma"
  elif atoms == delta:
    name = "delta"
  elif atoms == epsilon:
    name = "epsilon"
  elif atoms == zeta:
    name = "zeta"
  return name

def build_name_hash(pdb_hierarchy):
  i_seq_name_hash = dict()
  for atom in pdb_hierarchy.atoms():
    i_seq_name_hash[atom.i_seq]=atom.pdb_label_columns()
  return i_seq_name_hash

def get_mmtype_from_resname(resname):
  class_string = common_residue_names_get_class(resname)
  if "amino_acid" in class_string:
    return "PROTEIN"
  elif "rna_dna" in class_string:
    return "NA"
  else:
    return "OTHER"

def exercise():
  from libtbx.test_utils import approx_equal
  try :
    import numpy
  except ImportError :
    test_numpy = False
    print("Numpy not installed, will skip array conversion.")
  else :
    test_numpy = True
  # ramachandran
  z_data = get_rotarama_data(pos_type="general",
    convert_to_numpy_array=test_numpy)
  z_data = get_rotarama_data(pos_type="pre-proline",
    convert_to_numpy_array=test_numpy)
  # rotamer
  z_data = get_rotarama_data(residue_type="arg",
    db="rota",
    convert_to_numpy_array=test_numpy)
  z_data = get_rotarama_data(residue_type="phe",
    db="rota",
    convert_to_numpy_array=test_numpy)
  atom_info = decode_atom_str(" OD2 ASP A  14L")
  assert (atom_info.name == " OD2") and (atom_info.resname == "ASP")
  assert (atom_info.altloc == " ") and (atom_info.chain_id == "A")
  assert (atom_info.resid == "  14L") and (atom_info.resseq == "14")
  mpscore = molprobity_score(48.0, 9.95, 86.44) # 2hr0
  assert approx_equal(mpscore, 3.55, eps=0.01)
  mpscore = molprobity_score(215.8, 17.99, 52.18) # 3mku
  assert approx_equal(mpscore, 4.71, eps=0.01)

if __name__ == "__main__" :
  exercise()
  print("OK")
