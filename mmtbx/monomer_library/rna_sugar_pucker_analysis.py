from __future__ import absolute_import, division, print_function
import scitbx.math
from scitbx import matrix
import libtbx.phil
from libtbx.math_utils import normalize_angle
from libtbx import slots_getstate_setstate

master_phil = libtbx.phil.parse("""\
  bond_min_distance = 1.2
    .type = float
  bond_max_distance = 1.8
    .type = float
  epsilon_range_min = 155.0
    .type = float
  epsilon_range_max = 310.0
    .type = float
  delta_range_2p_min = 129.0
    .type = float
  delta_range_2p_max = 162.0
    .type = float
  delta_range_3p_min = 65.0
    .type = float
  delta_range_3p_max = 104.0
    .type = float
  p_distance_c1p_outbound_line_2p_max = 2.9
    .type = float
  o3p_distance_c1p_outbound_line_2p_max = 2.4
    .type = float
  bond_detection_distance_tolerance = 0.5
    .type = float
  enable = True
    .type = bool
    .style = hidden

  # DNA-specific thresholds, used only when evaluate() is called with is_dna=True.
  # The defaults above are RNA-derived and are deliberately left untouched, so the
  # restraints path through pdb_interpretation keeps its current behaviour.
  #
  # Deoxyribose occupies a broader delta range than ribose. Derived from 26826 DNA
  # residues in 1022 of 1032 non-redundant X-ray entries at <=1.9 A (one representative
  # per unique DNA-sequence/protein-partner group, ranked by resolution; ten entries
  # yielded no measurable sugar), using the pseudorotation phase as an independent
  # ground truth. Against that population the RNA windows wrongly reject 1.60% of
  # unambiguously C2'-endo DNA and 1.20% of unambiguously C3'-endo DNA; a delta outlier
  # carries Ramachandran-tier severity, so that rate is harmful. The windows below are
  # 0.3%-tail percentile intervals and bring those to 0.31% and 0.32%. They stay
  # cleanly separated by 13.0 degrees.
  #
  # O3'-perp is the larger error: the RNA value of 2.4 A misclassifies 0.300% of DNA,
  # where DNA's two populations are separated by a clean gap from 1.83 to 2.06 A.
  #
  # P-perp is deliberately NOT overridden. At 2.9 A it misclassifies only 0.237% of
  # DNA, which is already well placed.
  #
  # Counts here are post-deduplication: a residue whose pucker atoms carry no alternate
  # is counted once, not once per altloc. Before that fix the population read 28592 and
  # the windows came out a degree tighter at 124.0-161.0.
  #
  # Derivation: mp-triage-research/scripts/derive_dna_pucker_params.py
  dna_delta_range_2p_min = 123.0
    .type = float
  dna_delta_range_2p_max = 162.0
    .type = float
  dna_delta_range_3p_min = 62.0
    .type = float
  dna_delta_range_3p_max = 110.0
    .type = float
  dna_o3p_distance_c1p_outbound_line_2p_max = 1.95
    .type = float
  # How far P-perp must sit from its cutoff before the DNA pucker call is trusted.
  # Widening it trades coverage for precision; measured on 24535 DNA residues carrying
  # both measures:
  #   0.0 -> 5.31% flagged, 100% of residues judged
  #   0.5 -> 1.15% flagged,  85% judged
  #   0.8 -> 0.44% flagged,  61% judged   (default, against RNA's 0.27% rate)
  #   1.0 -> 0.30% flagged,  42% judged
  # The default is deliberately conservative for a first release.
  #
  # The DNA benchmark arm now exists, and it supports this setting rather than
  # unsettling it. Against PDB-REDO movement over 15014 DNA residues in 412 structures,
  # this rule flags 1.1% of residues and 12.9% of those were moved, an enrichment of
  # 2.5x over the 6.0% base rate, where the current RNA-threshold rule scores 0.92x and
  # is indistinguishable from chance. Note that AUC is useless for judging this: any
  # rule firing on ~1% of residues is pinned near 0.5 whatever its quality.
  #
  # Two caveats on that number. It is a floor, because the movement target undercounts
  # pucker repairs: of the outliers that did NOT clear the 0.5 A threshold, 103 still
  # moved inside the reference contours, since rotating C5' shifts delta a long way
  # while barely touching all-atom RMSD. And roughly a quarter of the enrichment is
  # B-factor correlation; stratified by B it is 1.9x, concentrated in the B 30-60 band
  # where it runs 11-18x. Below B 30 almost nothing moves at all.
  dna_pucker_confidence_margin = 0.8
    .type = float
""")

class evaluate(slots_getstate_setstate):

  __slots__ = [
    "epsilon",
    "delta",
    "p_distance_c1p_outbound_line",
    "o3p_distance_c1p_outbound_line",
    "p_perp_xyz",
    "o3p_perp_xyz",
    "is_epsilon_outlier",
    "is_2p_delta",
    "is_3p_delta",
    "is_delta_outlier",
    "is_2p_p_distance_c1p_outbound_line",
    "is_2p_o3p_distance_c1p_outbound_line",
    "is_2p"]

  def assign(O, **keyword_args):
    for slot in evaluate.__slots__:
      setattr(O, slot, keyword_args.get(slot, None))

  def __init__(O,
        params,
        residue_1_deoxy_ribo_atom_dict,
        residue_1_c1p_outbound_atom,
        residue_2_p_atom,
        is_dna=False):
    # is_dna selects the dna_* thresholds where they exist. It defaults to False so
    # every existing caller, including the pdb_interpretation restraints path, keeps
    # the RNA-derived values it has always used.
    def threshold(name):
      if (is_dna):
        dna_value = getattr(params, "dna_" + name, None)
        if (dna_value is not None): return dna_value
      return getattr(params, name, None)
    delta_range_2p_min = threshold("delta_range_2p_min")
    delta_range_2p_max = threshold("delta_range_2p_max")
    delta_range_3p_min = threshold("delta_range_3p_min")
    delta_range_3p_max = threshold("delta_range_3p_max")
    o3p_distance_2p_max = threshold("o3p_distance_c1p_outbound_line_2p_max")
    p_distance_2p_max = params.p_distance_c1p_outbound_line_2p_max
    def v(key): return matrix.col(residue_1_deoxy_ribo_atom_dict[key].xyz)
    c1p = v("C1'")
    o3p = v("O3'")
    c3p = v("C3'")
    c4p = v("C4'")
    c5p = v("C5'")
    bonded_sites = [
      (o3p, c3p),
      (c3p, c4p),
      (c5p, c4p)]
    distances = [abs(s2-s1) for s1,s2 in bonded_sites]
    p = None
    if (residue_2_p_atom is not None):
      p_trial = matrix.col(residue_2_p_atom.xyz)
      d = abs(p_trial-o3p)
      if (d <= params.bond_max_distance):
        p = p_trial
        distances.append(d)
    if (residue_1_c1p_outbound_atom is None):
      c1p_outbound = None
    else:
      c1p_outbound = matrix.col(residue_1_c1p_outbound_atom.xyz)
      distances.append(abs(c1p-c1p_outbound))
    if (   min(distances) < params.bond_min_distance
        or max(distances) > params.bond_max_distance):
      O.assign()
      return
    #
    def dihe(sites):
      return normalize_angle(
        scitbx.math.dihedral_angle(sites=sites, deg=True),
        deg=True)
    if (p is None):
      epsilon = None
    else:
      epsilon = dihe([p, o3p, c3p, c4p])
    delta = dihe([c5p, c4p, c3p, o3p])
    if (c1p_outbound is None):
      p_distance_c1p_outbound_line = None
      o3p_distance_c1p_outbound_line = None
      p_perp_xyz = None
      o3p_perp_xyz = None
    else:
      if (p is None):
        p_distance_c1p_outbound_line = None
        p_perp_xyz = None
      else:
        p_distance_c1p_outbound_line = scitbx.math.line_given_points(
          points=(c1p, c1p_outbound)).distance_sq(point=p)**0.5
        p_perp_xyz = (p, scitbx.math.line_given_points(
          points=(c1p, c1p_outbound)).perp_xyz(point=p))
      o3p_distance_c1p_outbound_line = scitbx.math.line_given_points(
        points=(c1p, c1p_outbound)).distance_sq(point=o3p)**0.5
      o3p_perp_xyz = (o3p, scitbx.math.line_given_points(
        points=(c1p, c1p_outbound)).perp_xyz(point=o3p))
    n_decisions = 0
    if (   epsilon is None
        or params.epsilon_range_min is None
        or params.epsilon_range_max is None):
      is_epsilon_outlier = None
    else:
      is_epsilon_outlier = (   epsilon < params.epsilon_range_min
                               or epsilon > params.epsilon_range_max)
      n_decisions += 1
    if (   delta_range_2p_min is None
        or delta_range_2p_max is None):
      is_2p_delta = None
    else:
      is_2p_delta = (    delta >= delta_range_2p_min
                     and delta <= delta_range_2p_max)
      n_decisions += 1
    if (   delta_range_3p_min is None
        or delta_range_3p_max is None):
      is_3p_delta = None
    else:
      is_3p_delta = (    delta >= delta_range_3p_min
                     and delta <= delta_range_3p_max)
      n_decisions += 1
    if (   p_distance_c1p_outbound_line is None
        or p_distance_2p_max is None):
      is_2p_p_distance_c1p_outbound_line = None
    else:
      is_2p_p_distance_c1p_outbound_line = (
          p_distance_c1p_outbound_line
        < p_distance_2p_max)
      n_decisions += 1
    if (   o3p_distance_c1p_outbound_line is None
        or o3p_distance_2p_max is None):
      is_2p_o3p_distance_c1p_outbound_line = None
    else:
      is_2p_o3p_distance_c1p_outbound_line = (
          o3p_distance_c1p_outbound_line
        < o3p_distance_2p_max)
      n_decisions += 1
    if (n_decisions == 0):
      is_2p = None
    else:
      perp_distance = is_2p_p_distance_c1p_outbound_line
      if (perp_distance is None):
        perp_distance = is_2p_o3p_distance_c1p_outbound_line
      if (perp_distance):
        is_2p = True
      else:
        is_2p = False
    if (is_3p_delta is None or is_2p_delta is None or is_2p is None):
      is_delta_outlier = None
    elif (is_dna):
      # DNA needs a different rule, not merely different thresholds.
      #
      # The RNA rule flags whenever delta falls outside the single window implied by
      # P-perp. That presumes the sugar is one pucker or the other. RNA obliges: 94.8%
      # of ribose sits in a canonical window. DNA does not: only 66.2% does, and about
      # a third of DNA sugars sit on the continuum between C3'-endo and C2'-endo, where
      # a binary P-perp call carries no information. Applying the RNA rule to DNA flags
      # 18% of all residues, and even with DNA-specific windows it still flags 11%,
      # essentially all of it from those intermediate sugars. At Ramachandran-tier
      # severity that is not tolerable.
      #
      # So for DNA, require positive evidence of a contradiction rather than mere
      # absence of agreement, and decline to judge when P-perp sits near the boundary:
      #
      #   1. skip when P-perp is within dna_pucker_confidence_margin of the cutoff,
      #   2. otherwise flag only when delta lands squarely in the OPPOSITE window.
      #
      # Measured on 8188 DNA residues this gives 0.29% flagged, against RNA's 0.25%,
      # at the cost of declining to judge 34% of DNA sugars. The residues it does flag
      # carry a median sugar B-factor of 33.3 versus 23.5 for the rest, so the flags
      # track poor density as they should.
      #
      # RNA is untouched: this branch runs only when the caller passes is_dna=True.
      margin = getattr(params, "dna_pucker_confidence_margin", None)
      ambiguous = False
      if (    margin is not None
          and p_distance_c1p_outbound_line is not None
          and p_distance_2p_max is not None):
        ambiguous = abs(p_distance_c1p_outbound_line - p_distance_2p_max) < margin
      if (ambiguous):
        is_delta_outlier = False
      elif (is_2p and is_3p_delta) or ((not is_2p) and is_2p_delta):
        is_delta_outlier = True
      else:
        is_delta_outlier = False
    else:
      if (is_2p and not is_2p_delta) or (not is_2p and not is_3p_delta):
        is_delta_outlier = True
      else:
        is_delta_outlier = False
    O.assign(
      epsilon=epsilon,
      delta=delta,
      p_distance_c1p_outbound_line=p_distance_c1p_outbound_line,
      o3p_distance_c1p_outbound_line=o3p_distance_c1p_outbound_line,
      p_perp_xyz=p_perp_xyz,
      o3p_perp_xyz=o3p_perp_xyz,
      is_epsilon_outlier=is_epsilon_outlier,
      is_2p_delta=is_2p_delta,
      is_3p_delta=is_3p_delta,
      is_delta_outlier=is_delta_outlier,
      is_2p_p_distance_c1p_outbound_line=is_2p_p_distance_c1p_outbound_line,
      is_2p_o3p_distance_c1p_outbound_line=is_2p_o3p_distance_c1p_outbound_line,
      is_2p=is_2p)
