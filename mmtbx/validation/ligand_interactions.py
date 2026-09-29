"""
Ligand interaction profile from geometric tools (H-bonds, clashes, vdW contacts,
salt bridges).

Input: an mmtbx.model.manager with electron-cloud H, as the ligand validation tool
prepares it (this module does not place H), processed with restraints, the
ligand's iselection and its selection string. The selection must be one molecule:
connected by the restraints' bonds (it may span residues, e.g. a glycan; one atom,
e.g. an ion, is connected; alternate conformers of an atom count as connected).
Sources:
  pnp     cctbx process_nonbonded_proxies on the ligand and the residues within
          3 A (ligand_overlaps; what validate_ligands reports): clashes and H-bonds
          with at least one ligand atom; H-bond criteria from the hbond scope
          (defaults: pnp.h_bond())
  probe2  library call (mmtbx.programs.probe2, approach=both, ligand selection as
          source, the rest as target, raw output, no files written); per atom
          pair the dots per class (both directions; the pair class is the first
          class in pair_class_order with dots), the minimum atom-pair gap, and per
          side the dots and contact area (dots / density) per class and in total:
          on the ligand atom's surface (ligand -> environment dots) and on the
          environment atom's surface (environment -> ligand dots). Classes: wc cc
          so hb bo; wh with allow_weak_hydrogen_bonds; wo with
          separate_worse_clashes; any other class is an error
  charged groups  standard amino acids by template, nucleotides by name, other
          residues from mmtbx.ligands.rdkit_utils.residue_molecule (below)
Entries (type, subtype, atoms, operators, labels, residue, symop, geometry,
sources, cross_check, model_support), ligand-environment only:
  hbond        one per D-H...A: pnp's H-bonds and probe2's hb and wh pairs
               (subtype "weak (probe2)" for wh)
  clash        pnp's clashes and probe2's bo and wo pairs; an atom clashing with
               two bonded atoms in line with it is one clash (pnp's rule in
               _process_clashes: the partners bonded, |cos| > 0.707, pnp.cos_vec;
               not applied to symmetry pairs, as in pnp): one entry listing every
               pair (pairs), the pair pnp kept (else the shortest) as its atoms
  vdw          probe2's wc, cc or so pairs (subtype), dots and area per class,
               minimum gap
  salt_bridge  oppositely charged groups (below); subtype from Kumar & Nussinov
               (2002); lists the H-bond entries between the same groups
cross_check (hbond, clash): "pnp and probe2", "pnp only", "probe2 only" (the last
two listed as disagreements), or "symmetry, probe2 not applicable" (probe2 does not
see symmetry-related copies). Ligand-internal pnp records and salt bridges are
listed apart (internal), not as entries. model_support is left unset.
The donor of an H in a probe2 H-bond is its bonded heavy atom from the restraints'
1-2 connectivity (same conformer); an H without exactly one bonded heavy atom is
listed in unresolved_donors and the donor left unset. probe2's hb geometry adds
d_HA and a_DHA from the model (probe2's hb class is overlap-based, without angles).
Symmetry convention: the selected ligand stays in place; an entry's operators (one
per atom, "x,y,z" for untransformed atoms) move its partner atoms (the
non-ligand side; for the ligand with its own symmetry copy, the acceptor of an
H-bond, the second atom of a clash, the partner group of a salt bridge). symop is
the partner's operator. pnp's operator belongs to its own atom pair (rt_mx_ji);
the operator or its inverse is taken, whichever applied to the partner reproduces
pnp's distance. The partner residue is labelled with the operator
("A 856 301 (x+1,y,z)"), also when it is the ligand itself.
Charged groups (find_charged_groups, charged_group_pairs; independent of the
ligand): examined on "(ligand) or residues_within(search, ligand)" (symmetry
included; search = max(atom_pair_cutoff, charge_centre_cutoff + 2.5 A)), full-model
i_seqs and bonds, per conformer (blank-altloc atoms plus one altloc; resname from
that conformer's atom_group). A group found identically in every conformer is
reported once with blank altloc, else with its atoms' or conformer's altloc; groups
with different non-blank altlocs are not paired (as pnp). Metal ions are not
salt-bridge partners (bonds to metals ignored).
Standard amino acids (iotbx.pdb class common_amino_acid), by template: Asp (CG;
OD1, OD2), Glu (CD; OE1, OE2), Lys (NZ), Arg (CZ; NE, NH1, NH2), His (CE1; ND1,
NE2), the N-terminal N (no bond to another residue), the C-terminal carboxylate
(C; O, OXT) only with OXT (a C without OXT and without a following residue is a
chain break: no group, nothing reported). Charge as modelled, from H/D bonded to
the group atoms in that conformer: Asp, Glu, C-terminus neutral if an O carries
H; His charged only with H on ND1 and NE2; Lys and N-terminus charged with four
bonded atoms on N; Arg charged unless a guanidinium H is missing. usual_charge: at
pH 7 (Asp, Glu, C-terminus -1; Lys, Arg, N-terminus +1; His 0). A residue
conformer without H gets usual_charge, state "assumed (no H)" (the N-terminus
only for the first residue of its chain). A missing template atom is reported
(charged_group_missing_atoms) and the group skipped.
Nucleotides (common_rna_dna), by name: the phosphate (P; OP1, OP2, and OP3 when
present), one negative charge per terminal O without H beyond the first;
usual_charge -1, -2 with OP3; without H in the residue conformer usual_charge,
"assumed (no H)"; no P (5' end): no group.
Other residues (the ligand, modified residues, cofactors, other het groups):
rdkit_utils.residue_molecule per conformer (builder_groups). Groups: a charged
heteroatom of the builder's molecule and its resonance partners (same element,
sharing a heavy neighbour; O and S partners terminal), overlapping sets merged.
Charged atoms bonded to an opposite charge form clusters (connected through such
bonds): net 0 (nitro, N-oxide, organic azide) is not a group and is listed with
the dropped groups ("charge-separated, net 0"); a nonzero net (nitrate, azide ion)
makes a group from the cluster's atoms of the net charge's sign, with the net
charge. Kind by SMARTS at the group's centre: carboxylate, guanidinium,
imidazolium, amidinium, ammonium (quaternary included), phosphate/phosphonate,
sulfate/sulfonate; any other charged heteroatom (pyridinium, tetrazolate,
phenolate, thiolate, nitrate, ...) "other". Charge: the sum of the builder's
formal charges over the group (a cluster: its net charge). A group with a capped atom (linked, or next to a missing
atom) is dropped and listed (charged_groups_dropped); a metal-bound atom flags the
group (metal_bound). Uncertain (certain False): the builder's total by search (no
formal charges), or, with H in the model, H completed from the restraint file on
a charged atom or an atom bonded to the centre; uncertain groups are listed but
make no salt bridges. Without H in the model: state "assumed (no H)". Per group:
source "builder", charge_source (the builder's total_charge_source), hydrogens,
notes. A residue whose residue_molecule fails has no groups and is listed
(charged_group_failures, with the reason; its formal charges are not compared).
A group's charge centre is the centroid of its charged atoms.
Formal charges in the restraint file (parsed, as the model's monomer library
server resolves it: files supplied with the model, else GeoStd, else the monomer
library; chem_comp_atom.formal_charge(), 0 if absent) and in the CCD are compared
for the ligand and the residues near it: template residues per group, and per
dictionary-charged atom outside the groups (against 0); builder residues atom by
atom, the builder's formal charges against the dictionary's, resonance-equivalent
atoms (same element sharing a heavy neighbour) summed. Where the dictionary's H on the heavy
atoms within two bonds match the model's, a different charge is a conflict
(reported, not resolved); where they differ, "protonation differs".
Salt bridges: oppositely charged groups, criterion atom_pair (default: at least one
pair of charged atoms within 4.0 A; Barlow & Thornton, J. Mol. Biol. 168, 867
(1983): <= 4 A between charged groups) or charge_centre (charge centres within
5.5 A; PLIP, Salentin et al., Nucleic Acids Res. 43, W443 (2015): Barlow &
Thornton's 4 A + 1.5 A). Subtype (Kumar & Nussinov, Biophys. J. 83, 1595 (2002)):
"salt bridge (K&N)" if the charge centres and at least one N-O pair are within
4 A, "N-O bridge (K&N)" if only an N-O pair is, "longer-range (K&N)" otherwise.
Contact patches: probe2's ligand -> environment dots, at their location on the
ligand's surface, grouped by spatial connectivity (dots closer than
patch_link_distance are linked); per patch the dots and area per class, the mean
dot location (center), atom pairs, ligand atoms and residues. No patch-level type
and no rule for counting patches.
"""
from __future__ import absolute_import, division, print_function
import math
import os
import sys
from six.moves import cStringIO as StringIO
import iotbx.phil
import cctbx.geometry_restraints.process_nonbonded_proxies as pnp
from libtbx import group_args
from libtbx.utils import Sorry, null_out
from scitbx.array_family import flex
from cctbx import crystal, sgtbx

# probe2's classes with the default options (in this order in reports)
probe_classes = ("wc", "cc", "so", "hb", "bo")
# every class probe2 writes: wh with allow_weak_hydrogen_bonds, wo with
# separate_worse_clashes
probe_all_classes = ("wc", "cc", "wh", "so", "hb", "bo", "wo")
# pair_class_order to use when wh or wo can occur
extended_pair_class_order = ("wo", "bo", "hb", "wh", "so", "cc", "wc")
clash_classes = ("bo", "wo")
hbond_classes = ("hb", "wh")
vdw_classes = ("so", "cc", "wc")

# pnp.manager._process_clashes: two clashes of one atom are one if the other two
# atoms are bonded and |pnp.cos_vec| exceeds this (hard-coded there)
inline_clash_cos_min = 0.707

metal_elements = frozenset("""LI BE NA MG AL K CA SC TI V CR MN FE CO NI CU ZN GA RB
  SR Y ZR NB MO TC RU RH PD AG CD IN SN CS BA LA CE PR ND PM SM EU GD TB DY HO ER
  TM YB LU HF TA W RE OS IR PT AU HG TL PB BI PO FR RA AC TH PA U NP PU""".split())

def _hbond_phil_str():
  h = pnp.h_bond()
  return """
  hbond
    .help = "process_nonbonded_proxies (pnp) H-bond criteria, passed to pnp.manager \\
as h_bond_params. Defaults are pnp.h_bond()'s values, which validate_ligands uses."
  {
    d_HA_cutoff = %s %s
      .type = floats(size=2)
      .help = "H...A distance range (A)."
    d_DA_cutoff = %s %s
      .type = floats(size=2)
      .help = "D...A distance range (A)."
    a_DHA_cutoff = %s
      .type = float
      .help = "Minimum D-H...A angle (deg)."
    a_YAH_cutoff = %s %s
      .type = floats(size=2)
      .help = "Y-A...H angle range (deg); stored by pnp but not applied."
    min_bonds_H_A = %s
      .type = int
      .help = "H and A at least this many bonds apart (same copy)."
    hydrogen_elements = %s
      .type = strings
    donor_elements = %s
      .type = strings
    acceptor_elements = %s
      .type = strings
  }
""" % (h.d_HA_cutoff[0], h.d_HA_cutoff[1], h.d_DA_cutoff[0], h.d_DA_cutoff[1],
    h.a_DHA_cutoff, h.a_YAH_cutoff[0], h.a_YAH_cutoff[1], h.min_bonds_H_A,
    " ".join(h.Hs), " ".join(h.Ds), " ".join(h.As))

master_phil_str = """
ligand_interactions
{
  pair_class_order = bo hb so cc wc
    .type = strings
    .help = "Pair class for an atom pair with dots of several probe2 classes: the \\
first class in this list that has dots. Must list every class that can occur: \\
with allow_weak_hydrogen_bonds (wh) or separate_worse_clashes (wo) use e.g. \\
wo bo hb wh so cc wc."
  patch_link_distance = 0.5
    .type = float
    .help = "Contact patches: ligand -> environment dots closer than this (A) \\
are linked. At probe density 16 dots/A^2 the nearest-neighbour spacing is \\
0.256 A (median; 90th percentile 0.261); 0.5 A bridges the seams between \\
neighbouring atoms' dot sets. 5XH3 ligand 856: 11 patches at 0.4 A, 9 at \\
0.45-0.5 A, 6 at 0.55-0.7 A."
  use_neutron_distances = False
    .type = bool
    .help = "probe2's use_neutron_distances (electron-cloud H is the convention)."
  separate_worse_clashes = False
    .type = bool
    .help = "probe2's output.separate_worse_clashes (class wo)."
  include scope mmtbx.probe.Helpers.probe_phil_parameters
""" + _hbond_phil_str() + """
  salt_bridge
    .help = "Salt bridges between oppositely charged groups."
  {
    criterion = *atom_pair charge_centre
      .type = choice
      .help = "atom_pair: at least one pair of charged atoms of the two groups \\
within atom_pair_cutoff (Barlow & Thornton, J. Mol. Biol. 168, 867 (1983)). \\
charge_centre: charge centres within charge_centre_cutoff (PLIP, Salentin et al., \\
Nucleic Acids Res. 43, W443 (2015)), a looser variant."
    atom_pair_cutoff = 4.0
      .type = float
      .help = "A. Barlow & Thornton (1983): ion pair if <= 4 A between charged groups."
    charge_centre_cutoff = 5.5
      .type = float
      .help = "A. PLIP's SALTBRIDGE_DIST_MAX: Barlow & Thornton's 4 A + 1.5 A."
    kumar_nussinov_cutoff = 4.0
      .type = float
      .help = "A. Subtype after Kumar & Nussinov, Biophys. J. 83, 1595 (2002): salt \\
bridge if the charge centres and at least one N-O pair are within this, N-O bridge \\
if only an N-O pair is, longer-range otherwise."
  }
}
"""

def master_params():
  return iotbx.phil.parse(master_phil_str, process_includes=True)

def h_bond_params(p=None):
  """pnp.h_bond() with the values of a ligand_interactions.hbond scope (None: pnp's defaults)."""
  h = pnp.h_bond()
  if p is not None:
    h.d_HA_cutoff = list(p.d_HA_cutoff)
    h.d_DA_cutoff = list(p.d_DA_cutoff)
    h.a_DHA_cutoff = p.a_DHA_cutoff
    h.a_YAH_cutoff = list(p.a_YAH_cutoff)
    h.min_bonds_H_A = p.min_bonds_H_A
    h.Hs = [e.upper() for e in p.hydrogen_elements]
    h.Ds = [e.upper() for e in p.donor_elements]
    h.As = [e.upper() for e in p.acceptor_elements]
  return h

# ------------------------------------------------------------------------------
# atoms

def atom_key(atom):
  """(chain, resid+icode, resname, name, altloc) of a hierarchy atom, as probe2 names it."""
  ag = atom.parent()
  rg = ag.parent()
  return (rg.parent().id.strip(), str(rg.resseq_as_int()) + rg.icode.strip(),
    ag.resname.strip().upper(), atom.name.strip(), ag.altloc.strip())

def atom_label(atom):
  ag = atom.parent()
  rg = ag.parent()
  alt = ag.altloc.strip()
  return "%s %s %s %s%s" % (rg.parent().id.strip(), ag.resname.strip(),
    (rg.resseq + rg.icode).strip(), atom.name.strip(), (" alt " + alt) if alt else "")

def residue_label(atom, resname=None):
  """Chain, resname (of the atom's atom_group unless given), resseq+icode."""
  ag = atom.parent()
  rg = ag.parent()
  return "%s %s %s" % (rg.parent().id.strip(), resname or ag.resname.strip(),
    (rg.resseq + rg.icode).strip())

def _element(atom):
  return atom.element.strip().upper()

# ------------------------------------------------------------------------------
# probe2

def probe_atom_field(atom):
  """The atom field of probe2's raw output for a hierarchy atom (probe2's format)."""
  ag = atom.parent()
  rg = ag.parent()
  icode = rg.icode if rg.icode != "" else " "
  return "{:>2s}{:>4s}{}{:>3s} {:<3s}{:1s}".format(rg.parent().id,
    str(rg.resseq_as_int()), icode, ag.resname.strip().upper(), atom.name, ag.altloc)

def probe_atom_table(atoms):
  """{probe2 atom field: i_seq}; Sorry if two atoms share a field."""
  result = {}
  for a in atoms:
    f = probe_atom_field(a)
    if f in result:
      raise Sorry("probe2 atom field %r is not unique." % f)
    result[f] = a.i_seq
  return result

def probe_key(field):
  """
  atom_key from a probe2 raw atom field, by probe2's columns
  "{:>2s}{:>4s}{}{:>3s} {:<3s}{:1s}": chain 2, resid 4, icode 1, resname (3 or
  more), a space, atom name, altloc 1. A resid wider than 4 characters shifts the
  columns; mapping to i_seqs uses probe_atom_table instead.
  """
  if len(field) < 15:
    raise Sorry("probe2 atom field %r not understood." % field)
  t = field[7:-1].lstrip().split(" ", 1)
  if len(t) != 2:
    raise Sorry("probe2 atom field %r not understood." % field)
  return (field[0:2].strip(), field[2:6].strip() + field[6].strip(), t[0],
    t[1].strip(), field[-1].strip())

def pair_class(dots, order=("bo", "hb", "so", "cc", "wc")):
  """The first class in order with dots (dots: {class: count})."""
  for c in order:
    if dots.get(c, 0):
      return c
  return None

def parse_probe_raw(text, order=("bo", "hb", "so", "cc", "wc")):
  """
  probe2 raw output (':direction:class:source:target:gap:dotGap:spike x:y:z:
  spikeLen:score:srcClass:targetClass:loc x:y:z:srcB:targetB'). Returns
  group_args(pairs, dots): pairs {frozenset((source field, target field)):
  dict(dots={class: n}, min_gap, pair_class)} over both directions (min_gap: the
  atom-pair gap; dots keyed by the classes in order); dots [dict(direction, cls,
  source, target, loc, spike, gap)] (loc: the dot on the source atom's surface;
  gap: the dot's gap). Sorry for a class probe2 does not write or not in order.
  """
  pairs, dots = {}, []
  for line in text.splitlines():
    f = line.split(":")
    if len(f) < 17:
      continue
    if f[2] not in probe_all_classes:
      raise Sorry("Unknown probe2 class %r in: %s" % (f[2], line))
    if f[2] not in order:
      raise Sorry("probe2 class %r is not in pair_class_order (%s)." % (f[2],
        " ".join(order)))
    source, target = f[3], f[4]
    p = pairs.setdefault(frozenset([source, target]),
      dict(dots=dict([(c, 0) for c in order]), min_gap=None))
    p["dots"][f[2]] += 1
    gap = float(f[5])
    p["min_gap"] = gap if p["min_gap"] is None else min(p["min_gap"], gap)
    dots.append(dict(direction=f[1], cls=f[2], source=source, target=target,
      spike=(float(f[7]), float(f[8]), float(f[9])),
      loc=(float(f[14]), float(f[15]), float(f[16])), gap=float(f[6])))
  for p in pairs.values():
    p["pair_class"] = pair_class(p["dots"], order)
  return group_args(pairs=pairs, dots=dots)

def probe2_parameters(source_selection, target_selection, probe=None,
                      use_neutron_distances=False, separate_worse_clashes=False):
  """probe2's master PHIL and params: approach=both, raw, no files; probe scope copied."""
  import iotbx.cli_parser
  from mmtbx.programs import probe2
  parser = iotbx.cli_parser.CCTBXParser(program_class=probe2.Program, logger=null_out())
  parser.parse_args([
    "source_selection=%s" % source_selection,
    "target_selection=%s" % target_selection,
    "approach=both", "output.format=raw", "output.write_files=False",
    "output.filename=probe2_ligand_interactions.txt",
    "output.separate_worse_clashes=%s" % separate_worse_clashes,
    "use_neutron_distances=%s" % use_neutron_distances])
  params = parser.working_phil.extract()
  if probe is not None:
    for name in [n for n in dir(probe) if not n.startswith("_")]:
      if hasattr(params.probe, name):
        setattr(params.probe, name, getattr(probe, name))
  return parser.master_phil, params

def run_probe2(model, source_selection, target_selection, probe=None,
               use_neutron_distances=False, separate_worse_clashes=False):
  """
  probe2 as a library call on a deep copy of model (run() can add phantom H to
  waters): returns the raw output string. The model must carry its H. Raw, not
  JSON: probe2's approach=both writes no JSON.
  """
  from iotbx.data_manager import DataManager
  from mmtbx.programs import probe2
  master_phil, params = probe2_parameters(source_selection, target_selection, probe,
    use_neutron_distances, separate_worse_clashes)
  dm = DataManager(["model"])
  dm.add_model("ligand_interactions_model", model)
  p2 = probe2.Program(dm, params, master_phil=master_phil, logger=null_out())
  p2.overrideModel(model.deep_copy(), processed=False)
  results, output = p2.run()
  return output

# ------------------------------------------------------------------------------
# process_nonbonded_proxies

def ligand_overlaps(model, sel_str, within_radius=3.0, h_bond_params=None):
  """
  Clashes and H-bonds with at least one ligand atom (cctbx
  process_nonbonded_proxies on the ligand and the residues within within_radius;
  h_bond_params: a pnp.h_bond(), None for pnp's defaults).
  Returns group_args: n_clashes, clashscore, n_clashes_sym, clashes_str,
  n_hbonds (validate_ligands' report), clash_records [dict(i_seq, j_seq,
  distance, sum_vdw_radii, overlap, symop)], hbond_records [dict(d_seq, h_seq,
  a_seq, d_HA, d_DA, a_DHA, symop)] (full-model i_seqs), hbond_criteria and
  clash_criteria (the values applied). The model must have H.
  """
  sel_within_str = '%s or (residues_within (%s, %s))' \
    % (sel_str, within_radius, sel_str)
  sel_within = model.selection(sel_within_str)
  model_within = model.select(sel_within)
  isel_ligand_within = model_within.iselection(sel_str)

  processed_nbps = pnp.manager(model = model_within, h_bond_params = h_bond_params)
  clashes = processed_nbps.get_clashes()
  hbonds = processed_nbps.get_hbonds()

  clashes_dict = clashes._clashes_dict
  hbonds_dict = hbonds._hbonds_dict

  ligand_clashes_dict = {}
  for iseq_tuple, record in clashes_dict.items():
    if (iseq_tuple[0] in isel_ligand_within or
        iseq_tuple[1] in isel_ligand_within):
      ligand_clashes_dict[iseq_tuple] = record

  ligand_clashes = pnp.clashes(
                  clashes_dict = ligand_clashes_dict,
                  model        = model_within)

  ligand_hbonds_dict = {}
  # iseq_tuple is (donor, H, acceptor)
  for iseq_tuple, record in hbonds_dict.items():
    if any(i_seq in isel_ligand_within for i_seq in iseq_tuple):
      ligand_hbonds_dict[iseq_tuple] = record

  ligand_hbonds = pnp.hbonds(
                  hbonds_dict  = ligand_hbonds_dict,
                  model        = model_within)

  results_hbonds = ligand_hbonds.get_results()

  string_io = StringIO()
  ligand_clashes.show(log=string_io, show_clashscore=False,
    show_header=False)

  results = ligand_clashes.get_results()

  # i_seqs from model_within to the full model; symop from the rt_mx (r[4]; r[3]
  # is a display flag)
  isel_within = sel_within.iselection()
  def xyz(rt_mx):
    return "" if rt_mx is None else rt_mx.as_xyz()
  clash_records = [dict(
    i_seq = int(isel_within[i]), j_seq = int(isel_within[j]),
    distance = r[0], sum_vdw_radii = r[1], overlap = r[2], symop = xyz(r[4]))
    for (i, j), r in ligand_clashes_dict.items()]
  hbond_records = [dict(
    d_seq = int(isel_within[d]), h_seq = int(isel_within[h]),
    a_seq = int(isel_within[a]), d_HA = r[0], d_DA = r[1], a_DHA = r[2],
    symop = xyz(r[4]))
    for (d, h, a), r in ligand_hbonds_dict.items()]
  p = processed_nbps
  # a_YAH_cutoff is stored by pnp but not applied
  hbond_criteria = dict(Hs=list(p.Hs), As=list(p.As), Ds=list(p.Ds),
    d_HA_cutoff=list(p.d_HA_cutoff), d_DA_cutoff=list(p.d_DA_cutoff),
    a_DHA_cutoff=p.a_DHA_cutoff, min_bonds_H_A=p.min_bonds_H_A)
  # hard-coded in pnp.manager: model_distance - vdw_sum < -0.40
  clash_criteria = dict(min_overlap=0.4, within_radius=within_radius)
  hbond_criteria["within_radius"] = within_radius

  return group_args(
    n_clashes      = results.n_clashes,
    clashscore     = results.clashscore,
    n_clashes_sym  = results.n_clashes_sym,
    clashes_str    = string_io.getvalue(),
    n_hbonds       = results_hbonds.n_hbonds,
    clash_records  = clash_records,
    hbond_records  = hbond_records,
    hbond_criteria = hbond_criteria,
    clash_criteria = clash_criteria)

def _identity(symop):
  return symop in (None, "", "x,y,z")

# ------------------------------------------------------------------------------
# charged groups and formal charges

# Standard amino acids (iotbx.pdb class common_amino_acid): groups by atom name.
# kind, centre, charged atoms, charge at pH 7
amino_acid_templates = {
  "ASP": ("carboxylate", "CG", ("OD1", "OD2"), -1),
  "GLU": ("carboxylate", "CD", ("OE1", "OE2"), -1),
  "LYS": ("ammonium", "NZ", ("NZ",), 1),
  "ARG": ("guanidinium", "CZ", ("NE", "NH1", "NH2"), 1),
  "HIS": ("imidazolium", "CE1", ("ND1", "NE2"), 0),
}
n_terminus_template = ("ammonium", "N", ("N",), 1)
c_terminus_template = ("carboxylate", "C", ("O", "OXT"), -1)

# Groups on residue_molecule's molecule: (kind, SMARTS); the first atom is the
# group centre (explicit H count in X)
builder_group_smarts = (
  ("carboxylate", "[CX3](~[OX1])~[OX1]"),
  ("guanidinium", "[CX3](~[NX3])(~[NX3])~[NX3]"),
  ("imidazolium", "[#6X3]1~[#7X3]~*~*~[#7X3]~1"),
  ("amidinium", "[CX3](~[NX3])~[NX3]"),
  ("ammonium", "[NX4+]"),
  ("phosphate", "[PX4](~[OX1])~[OX1]"),
  ("sulfate", "[SX4](~[OX1])(~[OX1])~[OX1]"),
)

def _altloc(atom):
  return atom.parent().altloc.strip()

def builder_groups(r, atoms):
  """
  Charged groups on a residue_molecule result r (ok): charged heteroatoms of r.mol
  with their resonance partners (same element, sharing a heavy neighbour; O and S
  partners terminal), overlapping sets merged; a charge balanced by a directly
  bonded opposite charge is resolved by clusters: charged atoms joined through bonds
  to an opposite charge; net 0 (nitro, N-oxide, organic azide) is not a group
  (dropped, "charge-separated, net 0"); otherwise (nitrate, azide ion) the cluster's
  atoms carrying the net charge's sign seed the group, with the net charge. Kind by
  builder_group_smarts at the group's centre (phosphate/sulfate: "phosphonate",
  "sulfonate" with a non-O heavy neighbour), else "other"; charge: the sum of the
  seeds' charges (a cluster counts its net charge). Returns (groups, dropped):
  groups as in find_charged_groups (i_seqs; altloc and resname unset); dropped:
  net-0 clusters and groups with a capped atom (linked, or next to a missing atom)
  [dict(kind, charge, center, charged, reason)].
  """
  from rdkit import Chem
  mol = r.mol
  iseq = r.rdkit_to_iseq
  def heavy_nb(a):
    return [n for n in a.GetNeighbors() if n.GetAtomicNum() > 1]
  charged = [a for a in mol.GetAtoms() if a.GetFormalCharge() and
    a.GetAtomicNum() not in (1, 6) and a.GetIdx() in iseq]
  q_of = dict([(a.GetIdx(), a.GetFormalCharge()) for a in charged])
  # charge-separated clusters: charged atoms joined through bonds to an opposite charge
  parent = dict([(k, k) for k in q_of])
  def find(k):
    while parent[k] != k:
      k = parent[k]
    return k
  for a in charged:
    for n in a.GetNeighbors():
      if n.GetIdx() in q_of and q_of[n.GetIdx()] * q_of[a.GetIdx()] < 0:
        parent[find(a.GetIdx())] = find(n.GetIdx())
  clusters = {}
  for k in q_of:
    clusters.setdefault(find(k), []).append(k)
  separated, net_zero = {}, []
  for members in clusters.values():
    if len(members) == 1:
      continue
    net = sum([q_of[k] for k in members])
    if net == 0:
      net_zero.append(sorted(members))
    else:
      for k in members:
        separated[k] = (net, members)
  # seeds: (atoms carrying the charge, charge); partners: same element sharing a heavy
  # neighbour, O and S partners terminal; overlapping sets merged, charges summed
  seeds, done = [], set()
  for a in charged:
    k = a.GetIdx()
    if k in done or [m for m in net_zero if k in m]:
      continue
    if k in separated:
      net, members = separated[k]
      done.update(members)
      seeds.append(([m for m in members if q_of[m] * net > 0], net))
    else:
      seeds.append(([k], q_of[k]))
  opposite = set([m for m in separated]) | set([m for z in net_zero for m in z])
  sets = []
  for atoms_q, q in seeds:
    members = set(atoms_q)
    for k in atoms_q:
      a = mol.GetAtomWithIdx(k)
      for c in heavy_nb(a):
        for p in heavy_nb(c):
          if p.GetAtomicNum() == a.GetAtomicNum() and p.GetIdx() in iseq and \
              (p.GetIdx() not in opposite or p.GetIdx() in atoms_q) and \
              (p.GetAtomicNum() not in (8, 16) or len(heavy_nb(p)) == 1):
            members.add(p.GetIdx())
    for s in [s for s in sets if s[0] & members]:
      members |= s[0]
      q += s[1]
      sets.remove(s)
    sets.append((members, q))
  patterns = [(k, Chem.MolFromSmarts(p)) for k, p in builder_group_smarts]
  capped = {}
  for c in r.caps:
    capped.setdefault(c["on"], []).append(c)
  uncertain = set(r.uncertain_iseqs)
  no_h = r.hydrogens == "restraint file (no H in the model)"
  groups, dropped = [], []
  for members in net_zero:
    seqs = sorted([iseq[k] for k in members])
    positive = [k for k in members if q_of[k] > 0]
    dropped.append(dict(kind="charge-separated", charge=0,
      center=iseq[positive[0]] if positive else seqs[0], charged=seqs,
      reason="charge-separated, net 0"))
  for members, q in sets:
    if q == 0:
      continue
    if len(members) == 1:
      centre = list(members)[0]
    else:
      counts = {}
      for k in members:
        for n in heavy_nb(mol.GetAtomWithIdx(k)):
          counts[n.GetIdx()] = counts.get(n.GetIdx(), 0) + 1
      centre = sorted(counts, key=lambda k: (-counts[k], k))[0]
    kind = "other"
    for k, patt in patterns:
      hits = [m for m in mol.GetSubstructMatches(patt) if m[0] == centre and
        (len(m) == 1 and centre in members or set(m[1:]) & members)]
      if hits:
        kind = k
        break
    c_atom = mol.GetAtomWithIdx(centre)
    if kind in ("phosphate", "sulfate") and [n for n in heavy_nb(c_atom)
        if n.GetAtomicNum() != 8]:
      kind = {"phosphate": "phosphonate", "sulfate": "sulfonate"}[kind]
    seqs = sorted([iseq[k] for k in members])
    center = iseq.get(centre, seqs[0])
    caps = [c for i in seqs for c in capped.get(i, [])]
    if caps:
      why = []
      for c in caps:
        if c["kind"] == "linked":
          why.append("%s linked to %s" % (atoms[c["on"]].name.strip(),
            atom_label(atoms[c["partner"]])))
        else:
          why.append("%s next to missing %s" % (atoms[c["on"]].name.strip(), c["partner"]))
      dropped.append(dict(kind=kind, charge=q, center=center, charged=seqs,
        reason="capped: %s" % "; ".join(why)))
      continue
    near = set(seqs) | set([iseq[n.GetIdx()] for n in heavy_nb(c_atom)
      if n.GetIdx() in iseq])
    notes = []
    if not r.charge_certain:
      notes.append("total charge by search (no formal charges)")
    unsure = [] if no_h else sorted(near & uncertain)
    if unsure:
      notes.append("H completed from the restraint file on %s" % " ".join(
        [atoms[i].name.strip() for i in unsure]))
    metal = [i for i in seqs if i in r.metal_bound]
    if metal:
      notes.append("metal-bound: %s" % " ".join([atoms[i].name.strip() for i in metal]))
    groups.append(dict(kind=kind, charge=q, usual_charge=None,
      state="assumed (no H)" if no_h else "modelled", center=center, charged=seqs,
      source="builder", certain=bool(r.charge_certain) and not unsure,
      metal_bound=bool(metal), charge_source=r.total_charge_source,
      hydrogens=r.hydrogens, notes=notes))
  return groups, dropped

def find_charged_groups(model, selection, fsc0=None, use_templates=True):
  """
  Charged groups among the selected atoms (flex.bool or flex.size_t; full-model
  i_seqs and fsc0), per conformer (blank-altloc atoms plus one altloc; resname
  from that conformer's atom_group). Standard amino acids by template (Asp, Glu,
  C-terminus: neutral if an O carries H; Lys, N-terminus: charged with four bonded
  atoms on N; Arg: charged unless a guanidinium H is missing; His: charged with H
  on ND1 and NE2; a residue conformer without H gets the charge at pH 7, state
  "assumed (no H)"; a missing template atom is reported and the group skipped; the
  C-terminal carboxylate only with OXT (a C without OXT and without a following
  residue is a chain break: no group, no report); without H the N-terminus only
  for the first residue of its chain). Nucleotides (common_rna_dna) by name: the
  phosphate (P; OP1, OP2, and OP3 when present), one negative charge per terminal
  O without H beyond the first; at pH 7 -1, -2 with OP3; without H in the residue
  conformer the charge at pH 7, "assumed (no H)"; no P: no group. All other
  residues (with use_templates=False all residues): builder_groups on
  mmtbx.ligands.rdkit_utils.residue_molecule per conformer; a failure gives no
  groups and is listed. A group found identically in every conformer has blank
  altloc, else the altloc of its atoms or its conformer. Returns
  group_args(groups=[dict(kind, charge, usual_charge, state, center, charged,
  altloc, resname, source, certain, metal_bound, charge_source, hydrogens,
  notes)], missing=[dict(residue, altloc, kind, atoms)], dropped=[dict(residue,
  altloc, kind, charge, atoms, reason)], failures=[dict(residue, altloc,
  reason)], builder_charges={(chain, resseq, icode, altloc): {i_seq: formal
  charge}} for the residue conformers built).
  """
  import iotbx.pdb
  from mmtbx.ligands import rdkit_utils
  atoms = model.get_hierarchy().atoms()
  if fsc0 is None:
    fsc0 = model.get_restraints_manager().geometry.shell_sym_tables[0] \
      .full_simple_connectivity()
  if isinstance(selection, flex.bool):
    selection = selection.iselection()
  sel = set(selection)
  el = [_element(a) for a in atoms]
  alt = [_altloc(a) for a in atoms]
  alts = sorted(set([alt[i] for i in sel if alt[i]])) or [""]
  def residue_key(i):
    rg = atoms[i].parent().parent()
    return (rg.parent().id, rg.resseq, rg.icode)
  found = {}
  missing, dropped, failures = [], [], []
  built = {}
  def freeze(g):
    return tuple(sorted([(k, tuple(v) if isinstance(v, list) else v)
      for k, v in g.items()]))
  for conf in alts:
    def visible(k):
      return alt[k] in ("", conf)
    def nb(i):
      return [k for k in fsc0[i] if visible(k) and el[k] not in metal_elements]
    def h_on(i):
      return [k for k in nb(i) if el[k] in ("H", "D")]
    def add(g):
      found.setdefault(freeze(g), set()).add(conf)
    def template_group(kind, charge, usual, state, center, charged, resname):
      add(dict(kind=kind, charge=charge, usual_charge=usual, state=state, center=center,
        charged=charged, resname=resname, source="template", certain=True,
        metal_bound=False, charge_source=None, hydrogens=None, notes=[]))
    # residue conformers in the selection
    residues = {}
    for i in sel:
      if visible(i):
        residues.setdefault(residue_key(i), []).append(i)
    for key, seqs in sorted(residues.items()):
      ags = [atoms[i].parent() for i in seqs]
      ag = ([a for a in ags if a.altloc.strip()] or ags)[0]
      rg = ag.parent()
      resname = ag.resname.strip().upper()
      rclass = iotbx.pdb.common_residue_names_get_class(resname)
      # a residue can be split into residue groups (e.g. Asp A / Asn B)
      rgs = [x for x in rg.parent().residue_groups() if (x.resseq, x.icode) ==
        (rg.resseq, rg.icode)]
      rg_atoms = [a for x in rgs for a in x.atoms() if visible(a.i_seq)]
      names = dict([(a.name.strip(), a.i_seq) for a in rg_atoms])
      has_h = len([a for a in rg_atoms if el[a.i_seq] in ("H", "D")]) > 0
      label = residue_label(rg_atoms[0], resname)
      if use_templates and rclass == "common_amino_acid":
        todo = []
        if resname in amino_acid_templates:
          todo.append(amino_acid_templates[resname])
        n = names.get("N")
        first = rg.parent().residue_groups()[0]
        if n is not None and not [k for k in fsc0[n] if residue_key(k) != key] and (
            has_h or (first.resseq, first.icode) == (rg.resseq, rg.icode)):
          todo.append(n_terminus_template)
        c = names.get("C")
        if c is not None and "OXT" in names and \
            not [k for k in fsc0[c] if residue_key(k) != key]:
          todo.append(c_terminus_template)
        for kind, center, charged, usual in todo:
          absent = [x for x in (center,) + charged if x not in names]
          if absent:
            missing.append((label, kind, tuple(absent), conf))
            continue
          ci, qi = names[center], [names[x] for x in charged]
          if not has_h:
            q, state = usual, "assumed (no H)"
          else:
            state = "modelled"
            if kind == "carboxylate":
              q = 0 if [k for k in qi if h_on(k)] else -1
            elif kind == "ammonium":
              q = 1 if len(nb(qi[0])) == 4 else 0
            elif kind == "guanidinium":
              q = 1 if sum([len(h_on(k)) for k in qi]) == 5 else 0
            else:
              q = 1 if not [k for k in qi if not h_on(k)] else 0
          template_group(kind, q, usual, state, ci, qi, resname)
        continue
      if use_templates and rclass == "common_rna_dna":
        if "P" not in names:
          continue
        charged = [x for x in ("OP1", "OP2") if x in names]
        absent = [x for x in ("OP1", "OP2") if x not in names]
        if absent:
          missing.append((label, "phosphate", tuple(absent), conf))
          continue
        if "OP3" in names:
          charged.append("OP3")
        qi = [names[x] for x in charged]
        usual = -(len(qi) - 1)
        if not has_h:
          q, state = usual, "assumed (no H)"
        else:
          free = [k for k in qi if not h_on(k)]
          q, state = -max(0, len(free) - 1), "modelled"
        template_group("phosphate", q, usual, state, names["P"], qi, resname)
        continue
      # all other residues: the builder, per conformer of the residue
      own = sorted(set([ag_.altloc.strip() for ag_ in rg.atom_groups() if ag_.altloc.strip()]))
      eff = conf if conf in own else (own[0] if own else "")
      bkey = (key, eff)
      if bkey not in built:
        try:
          r = rdkit_utils.residue_molecule(model, rg, altloc=eff, fsc0=fsc0)
        except Exception as e:
          from libtbx import group_args as _ga
          r = _ga(ok=False, reason="%s: %s" % (type(e).__name__, e))
        built[bkey] = (r, builder_groups(r, atoms) if r.ok else ([], []))
      r, (bgroups, bdropped) = built[bkey]
      rg_label = residue_label(atoms[rg.atoms()[0].i_seq], resname)
      if not r.ok:
        failures.append((rg_label, r.reason, conf))
        continue
      for g in bgroups:
        if [i for i in g["charged"] if i in sel]:
          add(dict(g, resname=resname))
      for d in bdropped:
        dropped.append((rg_label, d["kind"], d["charge"],
          tuple([atoms[i].name.strip() for i in d["charged"]]), d["reason"], conf))
  groups = []
  for key, confs in sorted(found.items(), key=lambda x: (dict(x[0])["center"],
      dict(x[0])["kind"], sorted(x[1]))):
    g = dict([(k, list(v) if isinstance(v, tuple) else v) for k, v in key])
    own = sorted(set([alt[k] for k in [g["center"]] + g["charged"] if alt[k]]))
    if own:
      altlocs = own
    elif set(confs) == set(alts):
      altlocs = [""]
    else:
      altlocs = sorted(confs)
    for a in altlocs:
      groups.append(dict(g, altloc=a))
  def per_conformer(items):
    confs = {}
    for x in items:
      confs.setdefault(x[:-1], set()).add(x[-1])
    for k, cs in sorted(confs.items()):
      for a in ([""] if set(cs) == set(alts) else sorted(cs)):
        yield k, a
  return group_args(groups=groups,
    missing=[dict(residue=k[0], altloc=a, kind=k[1], atoms=list(k[2]))
      for k, a in per_conformer(missing)],
    dropped=[dict(residue=k[0], altloc=a, kind=k[1], charge=k[2], atoms=list(k[3]),
      reason=k[4]) for k, a in per_conformer(dropped)],
    failures=[dict(residue=k[0], altloc=a, reason=k[1]) for k, a in per_conformer(failures)],
    builder_charges=dict([(k + (eff,), dict([(i, r.mol.GetAtomWithIdx(x).GetFormalCharge())
      for x, i in r.rdkit_to_iseq.items()])) for (k, eff), (r, g) in built.items() if r.ok]))

def moved_site(unit_cell, xyz, rt_mx):
  return unit_cell.orthogonalize(rt_mx * unit_cell.fractionalize(xyz))

def charged_group_pairs(model, groups, cutoff):
  """
  [(k1, k2, op)]: oppositely charged groups (indices into groups) with a pair of
  charged atoms within cutoff, op moving group k2 next to group k1 (both
  directions listed; crystal symmetry included). Groups with different non-blank
  altlocs, and uncertain groups (certain False), are not paired.
  """
  atoms = model.get_hierarchy().atoms()
  sites = [(k, i) for k, g in enumerate(groups) if g["charge"] and g.get("certain", True)
    for i in g["charged"]]
  xyz = flex.vec3_double([atoms[i].xyz for k, i in sites])
  cs = model.crystal_symmetry()
  raw = []
  if cs is None or cs.unit_cell() is None or cs.space_group() is None:
    for a in range(len(sites)):
      for b in range(a + 1, len(sites)):
        if atoms[sites[a][1]].distance(atoms[sites[b][1]]) <= cutoff:
          raw.append((a, b, "x,y,z"))
  else:
    uc = cs.unit_cell()
    sps = cs.special_position_settings()
    asu = sps.asu_mappings(buffer_thickness=cutoff)
    asu.process_sites_cart(original_sites=xyz,
      site_symmetry_table=sps.site_symmetry_table(sites_cart=xyz))
    for p in crystal.neighbors_fast_pair_generator(asu, distance_cutoff=cutoff):
      rt = asu.get_rt_mx_i(p).inverse().multiply(asu.get_rt_mx_j(p))
      d = math.sqrt(p.dist_sq)
      best = None
      for c in (rt, rt.inverse()):
        x = moved_site(uc, atoms[sites[p.j_seq][1]].xyz, c)
        err = abs(atoms[sites[p.i_seq][1]].distance(x) - d)
        if best is None or err < best[0]:
          best = (err, c.as_xyz())
      raw.append((p.i_seq, p.j_seq, best[1]))
  result = set()
  for a, b, op in raw:
    ka, kb = sites[a][0], sites[b][0]
    ga, gb = groups[ka], groups[kb]
    if ga["charge"] * gb["charge"] >= 0:
      continue
    if ga["altloc"] and gb["altloc"] and ga["altloc"] != gb["altloc"]:
      continue
    op = "x,y,z" if _identity(op) else op
    inv = sgtbx.rt_mx(op).inverse().as_xyz()
    result.add((ka, kb, op))
    result.add((kb, ka, "x,y,z" if _identity(inv) else inv))
  return sorted(result)

def _with_h(atoms, bonds):
  h = dict([(n, set()) for n in atoms])
  for a, b in bonds:
    if a in atoms and b in atoms:
      if atoms[b][0] in ("H", "D"):
        h[a].add(b)
      if atoms[a][0] in ("H", "D"):
        h[b].add(a)
  return dict(atoms=atoms, h=h)

def restraints_formal_charges(model, resname):
  """
  (dictionary, file name) from the parsed restraint file as the model's monomer
  library server resolves it (files supplied with the model, else GeoStd, else
  the monomer library): {atoms: {name: (element, formal charge; 0 if absent)},
  h: {name: set(H names)}}; (None, file name) without a charge column, (None,
  None) without a file.
  """
  comp = model.get_mon_lib_srv().get_comp_comp_id_direct(resname)
  if comp is None:
    return None, None
  info = comp.source_info or ""
  file_name = info[len("file: "):] if info.startswith("file: ") else None
  # a file supplied with the model: its name as supplied
  for name, cif_object in (model.get_restraint_objects() or []):
    if file_name and os.path.basename(name) == os.path.basename(file_name):
      file_name = name
  if not [a for a in comp.atom_list if getattr(a, "charge", None) is not None]:
    return None, file_name
  atoms = dict([(a.atom_id.strip(), ((a.type_symbol or "").strip().upper(),
    a.formal_charge() or 0)) for a in comp.atom_list])
  bonds = [(b.atom_id_1.strip(), b.atom_id_2.strip()) for b in comp.bond_list]
  return _with_h(atoms, bonds), file_name

def ccd_formal_charges(resname):
  """The CCD entry's formal charges and H (chem_data), None if absent."""
  import mmtbx.chemical_components
  cif = mmtbx.chemical_components.get_cif_dictionary(resname)
  if not cif or "_chem_comp_atom" not in cif:
    return None
  atoms = {}
  for a in cif["_chem_comp_atom"]:
    try:
      q = int(a.charge)
    except (TypeError, ValueError):
      q = 0
    atoms[a.atom_id.strip()] = (a.type_symbol.strip().upper(), q)
  bonds = [(b.atom_id_1.strip(), b.atom_id_2.strip()) for b in cif.get("_chem_comp_bond", [])]
  return _with_h(atoms, bonds)

# ------------------------------------------------------------------------------
# patches

def dot_patches(xyz, link_distance):
  """Connected components of points (flex.vec3_double) linked below link_distance: list of index lists."""
  n = xyz.size()
  parent = list(range(n))
  def find(i):
    while parent[i] != i:
      parent[i] = parent[parent[i]]
      i = parent[i]
    return i
  cells = {}
  for i, p in enumerate(xyz):
    cells.setdefault(tuple([int(c // link_distance) for c in p]), []).append(i)
  d2 = link_distance * link_distance
  for (cx, cy, cz), members in cells.items():
    for dx in (-1, 0, 1):
      for dy in (-1, 0, 1):
        for dz in (-1, 0, 1):
          other = cells.get((cx + dx, cy + dy, cz + dz))
          if other is None:
            continue
          for i in members:
            for j in other:
              if j <= i:
                continue
              p, q = xyz[i], xyz[j]
              if (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2 < d2:
                ri, rj = find(i), find(j)
                if ri != rj:
                  parent[ri] = rj
  groups = {}
  for i in range(n):
    groups.setdefault(find(i), []).append(i)
  return sorted(groups.values(), key=lambda g: (-len(g), g[0]))

# ------------------------------------------------------------------------------

class manager(object):
  """
  The ligand interaction profile. model: mmtbx.model.manager with H, processed
  with restraints; ligand_isel: flex.size_t; sel_str: the ligand's selection
  string (selects exactly ligand_isel); params: extract of master_phil_str (the
  ligand_interactions scope), default master values.
  """
  def __init__(self, model, ligand_isel, sel_str, params=None, log=None):
    if params is None:
      params = master_params().extract().ligand_interactions
    self.model = model
    self.ligand_isel = flex.size_t(list(ligand_isel))
    self.sel_str = sel_str
    self.params = params
    self.log = log if log is not None else null_out()
    if not model.has_hd():
      raise Sorry("ligand_interactions needs a model with H (electron-cloud positions).")
    possible = list(probe_classes)
    if params.probe.allow_weak_hydrogen_bonds:
      possible.append("wh")
    if params.separate_worse_clashes:
      possible.append("wo")
    order = list(params.pair_class_order)
    unknown = [c for c in order if c not in probe_all_classes]
    missing = [c for c in possible if c not in order]
    if unknown or missing or len(set(order)) != len(order):
      raise Sorry("pair_class_order must list each probe2 class that can occur (%s) "
        "once, got %s. With allow_weak_hydrogen_bonds or separate_worse_clashes use "
        "e.g. %s." % (" ".join(possible), " ".join(order),
        " ".join(extended_pair_class_order)))
    self.order = order
    self.classes = [c for c in probe_all_classes if c in possible]
    self.entries = []
    self.internal = []
    self.unresolved_donors = []
    self.patches = []
    self.disagreements = []
    self.formal_charges = []

  def run(self):
    atoms = self.model.get_hierarchy().atoms()
    self._atoms = atoms
    lig = set(self.ligand_isel)
    self._lig = lig
    self._fsc0 = self.model.get_restraints_manager().geometry.shell_sym_tables[0] \
      .full_simple_connectivity()
    if set(self.model.selection(self.sel_str).iselection()) != lig:
      raise Sorry("sel_str %r does not select the ligand_isel atoms." % self.sel_str)
    self._check_one_molecule()
    table = probe_atom_table(atoms)
    # probe2; by selection, not resname: other copies of the ligand are environment
    self.probe_output = run_probe2(self.model, "(%s)" % self.sel_str, "not (%s)" % self.sel_str,
      probe=self.params.probe, use_neutron_distances=self.params.use_neutron_distances,
      separate_worse_clashes=self.params.separate_worse_clashes)
    parsed = parse_probe_raw(self.probe_output, self.order)
    # fields that are not model atoms (e.g. probe2's phantom water H) are skipped
    self.probe_pairs = {}
    self.probe_unmapped = set()
    for k, v in parsed.pairs.items():
      seqs = [table.get(x) for x in k]
      if None in seqs:
        self.probe_unmapped.update([x for x in k if x not in table])
        continue
      a, b = seqs if len(seqs) == 2 else seqs * 2
      i, j = (a, b) if a in lig else (b, a)
      if (i in lig) != (j in lig):
        self.probe_pairs[(i, j)] = v
    self.probe_dots = [dict(d, source=table[d["source"]], target=table[d["target"]])
      for d in parsed.dots if d["source"] in table and d["target"] in table]
    # dots per (source atom, target atom) and class: the source atom's surface
    self._side_dots = {}
    for d in self.probe_dots:
      c = self._side_dots.setdefault((d["source"], d["target"]), {})
      c[d["cls"]] = c.get(d["cls"], 0) + 1
    self.overlaps = ligand_overlaps(self.model, self.sel_str,
      h_bond_params=h_bond_params(self.params.hbond))
    self.clash_criteria = dict(self.overlaps.clash_criteria, inline_merge=dict(
      rule="an atom clashing with two atoms that are bonded to each other and in "
        "line with it is one clash (pnp.manager._process_clashes); not applied to "
        "symmetry pairs", cos_min=inline_clash_cos_min, function="pnp.cos_vec"))
    self._build_entries()
    self._build_salt_bridges()
    self._build_patches()
    return self

  def _check_one_molecule(self):
    """Sorry unless the ligand is one molecule by the restraints' bonds (altlocs of an atom joined)."""
    atoms = self._atoms
    parent = dict([(i, i) for i in self._lig])
    def find(i):
      while parent[i] != i:
        parent[i] = parent[parent[i]]
        i = parent[i]
      return i
    def union(i, j):
      parent[find(i)] = find(j)
    same_atom = {}
    for i in self._lig:
      for j in self._fsc0[i]:
        if j in parent:
          union(i, j)
      rg = atoms[i].parent().parent()
      k = (rg.parent().id, rg.resseq, rg.icode, atoms[i].name)
      if k in same_atom:
        union(i, same_atom[k])
      same_atom[k] = i
    parts = {}
    for i in sorted(self._lig):
      parts.setdefault(find(i), []).append(i)
    if len(parts) > 1:
      raise Sorry("The ligand selection %r is not one molecule: %d fragments not "
        "connected by the restraints' bonds (%s). Select one molecule." % (self.sel_str,
        len(parts), "; ".join([residue_label(atoms[g[0]]) + " " + atoms[g[0]].name.strip()
        for g in sorted(parts.values())])))

  # -- entries -----------------------------------------------------------------

  def _donor_of(self, h):
    """The heavy atom bonded to H (restraints' 1-2 connectivity); None, listed, if not exactly one."""
    heavy = [k for k in self._fsc0[h] if not self._is_h(k)]
    if len(heavy) == 1:
      return heavy[0]
    self.unresolved_donors.append(dict(h=atom_label(self._atoms[h]),
      bonded_heavy=[atom_label(self._atoms[k]) for k in heavy]))
    return None

  def _moved_xyz(self, xyz, rt_mx):
    return moved_site(self.model.crystal_symmetry().unit_cell(), xyz, rt_mx)

  def _moved(self, i, rt_mx):
    return self._moved_xyz(self._atoms[i].xyz, rt_mx)

  def _site(self, i, op):
    return self._atoms[i].xyz if _identity(op) else self._moved(i, sgtbx.rt_mx(op))

  def _partner_operator(self, fixed, partner, symop, distance):
    """pnp's operator or its inverse: the one that, applied to partner, gives distance from fixed."""
    rt = sgtbx.rt_mx(symop)
    best = None
    for c in (rt, rt.inverse()):
      x = self._moved(partner, c)
      err = abs(self._atoms[fixed].distance(x) - distance)
      if best is None or err < best[0]:
        best = (err, c.as_xyz())
    if best[0] > 1.e-3:
      raise Sorry("symmetry operator %s does not reproduce the distance %.3f A of %s ... %s."
        % (symop, distance, atom_label(self._atoms[fixed]), atom_label(self._atoms[partner])))
    return best[1]

  def _entry(self, type_, subtype, atom_seqs, geometry, sources, operators=None,
             cross_check=None):
    """operators: one per atom (None: all "x,y,z"); non-identity ones mark the partner."""
    atoms = self._atoms
    lig = self._lig
    if operators is None:
      operators = ["x,y,z"] * len(atom_seqs)
    moved = [k for k, i in enumerate(atom_seqs) if i is not None and
      not _identity(operators[k])]
    partner = moved or [k for k, i in enumerate(atom_seqs) if i is not None and
      i not in lig]
    ligand_atoms = [i for k, i in enumerate(atom_seqs) if i is not None and i in lig
      and k not in moved]
    def op(k):
      return "" if k not in moved else " (%s)" % operators[k]
    residue = None
    if partner:
      residue = residue_label(atoms[atom_seqs[partner[0]]]) + op(partner[0])
    return dict(type=type_, subtype=subtype,
      atoms=[None if i is None else int(i) for i in atom_seqs],
      operators=list(operators),
      labels=[None if i is None else atom_label(atoms[i]) + op(k)
        for k, i in enumerate(atom_seqs)],
      ligand_atoms=[atom_label(atoms[i]) for i in ligand_atoms],
      residue=residue, symop=(operators[moved[0]] if moved else None),
      geometry=geometry, sources=sorted(sources),
      cross_check=cross_check, model_support=None)

  def _is_h(self, i):
    return self._atoms[i].element.strip().upper() in ("H", "D")

  def _probe_geometry(self, i, j, p):
    """
    Pair (ligand atom i, environment atom j): dots of both directions per class
    (the pair class's basis), minimum gap, and per side the dots and the contact
    area (dots / density) per class and in total: ligand (dots on i's surface,
    ligand -> environment) and environment (on j's surface).
    """
    density = self.params.probe.density
    result = dict(dots=dict(p["dots"]), min_gap=p["min_gap"])
    for side, key in (("ligand", (i, j)), ("environment", (j, i))):
      n = self._side_dots.get(key, {})
      dots = dict([(c, n.get(c, 0)) for c in p["dots"]])
      result["dots_%s" % side] = dots
      result["area_%s" % side] = dict([(c, k / density) for c, k in dots.items()])
      result["area_%s_total" % side] = sum(dots.values()) / density
    return result

  def _add_checked(self, type_, atom_seqs, e, operators, i, j, subtype=None, extra=None):
    """hbond/clash entry with its pnp-probe2 cross-check; disagreements listed."""
    if operators is not None:
      check = "symmetry, probe2 not applicable"
    elif e["sources"] == set(["pnp", "probe2"]):
      check = "pnp and probe2"
    else:
      check = "%s only" % list(e["sources"])[0]
    entry = self._entry(type_, subtype, atom_seqs, e["geometry"], e["sources"],
      operators, check)
    if extra:
      entry.update(extra)
    self.entries.append(entry)
    if check.endswith(" only"):
      self.disagreements.append(dict(type=type_, labels=entry["labels"],
        sources=entry["sources"], missing=sorted(set(["pnp", "probe2"]) - e["sources"]),
        probe_class=self._probe_class_of(i, j)))

  def _build_entries(self):
    lig = self._lig
    atoms = self._atoms
    # H-bonds: key (H, A, operators of D, H, A)
    hb = {}
    for r in self.overlaps.hbond_records:
      d, h, a = r["d_seq"], r["h_seq"], r["a_seq"]
      g = dict(d_HA=r["d_HA"], d_DA=r["d_DA"], a_DHA=r["a_DHA"])
      if _identity(r["symop"]):
        ops = ""
        if d in lig and h in lig and a in lig:
          self.internal.append(dict(type="hbond", labels=[atom_label(atoms[k])
            for k in (d, h, a)], geometry=dict(pnp=g)))
          continue
      elif h in lig:
        op = self._partner_operator(h, a, r["symop"], r["d_HA"])
        ops = ("x,y,z", "x,y,z", op)
      else:
        op = self._partner_operator(a, h, r["symop"], r["d_HA"])
        ops = (op, op, "x,y,z")
      e = hb.setdefault((h, a, ops), dict(d=d, geometry={}, sources=set(), weak=False))
      e["sources"].add("pnp")
      e["geometry"]["pnp"] = g
    for (i, j), p in self.probe_pairs.items():
      if p["pair_class"] not in hbond_classes:
        continue
      g = self._probe_geometry(i, j, p)
      if self._is_h(i) != self._is_h(j):
        h, a = (i, j) if self._is_h(i) else (j, i)
        d = self._donor_of(h)
        # from the model: probe2's hb class does not use angles
        g["d_HA"] = atoms[h].distance(atoms[a])
        g["a_DHA"] = None if d is None else atoms[h].angle(atoms[a], atoms[d], deg=True)
      else:
        h, a, d = None, j, i
      e = hb.setdefault((h, a, ""), dict(d=d, geometry={}, sources=set(), weak=False))
      e["sources"].add("probe2")
      e["geometry"]["probe2"] = g
      e["weak"] = p["pair_class"] == "wh"
    for (h, a, ops), e in sorted(hb.items(), key=lambda x: (str(x[0][0]), x[0][1], str(x[0][2]))):
      self._add_checked("hbond", [e["d"], h, a], e, list(ops) if ops else None, h, a,
        subtype="weak (probe2)" if e["weak"] else None)
    self._build_clashes()
    # vdW contacts: probe2's wc, cc, so pairs
    for (i, j), p in sorted(self.probe_pairs.items()):
      if p["pair_class"] in vdw_classes:
        self.entries.append(self._entry("vdw", p["pair_class"], [i, j],
          dict(probe2=self._probe_geometry(i, j, p)), ["probe2"]))
    self.pair_class_order = list(self.order)

  def _build_clashes(self):
    """pnp's clashes and probe2's bo/wo pairs; inline pairs of one atom merged (pnp's rule)."""
    lig = self._lig
    atoms = self._atoms
    # key (ligand atom, partner atom, partner operator); for the ligand with its
    # own copy the lower i_seq stays
    cl = {}
    for r in self.overlaps.clash_records:
      i, j = sorted([r["i_seq"], r["j_seq"]])
      if j in lig and i not in lig:
        i, j = j, i
      g = dict(distance=r["distance"], sum_vdw_radii=r["sum_vdw_radii"],
        overlap=r["overlap"])
      if _identity(r["symop"]):
        op = ""
        if i in lig and j in lig:
          self.internal.append(dict(type="clash", labels=[atom_label(atoms[i]),
            atom_label(atoms[j])], geometry=dict(pnp=g)))
          continue
      else:
        op = self._partner_operator(i, j, r["symop"], r["distance"])
      e = cl.setdefault((i, j, op), dict(geometry={}, sources=set()))
      e["sources"].add("pnp")
      e["geometry"]["pnp"] = g
    for (i, j), p in self.probe_pairs.items():
      if p["pair_class"] not in clash_classes:
        continue
      e = cl.setdefault((i, j, ""), dict(geometry={}, sources=set()))
      e["sources"].add("probe2")
      e["geometry"]["probe2"] = self._probe_geometry(i, j, p)
    # merge: shared atom c, partners bonded and in line (pnp.cos_vec), no operator
    keys = sorted(cl)
    parent = dict([(k, k) for k in keys])
    def find(k):
      while parent[k] != k:
        k = parent[k]
      return k
    by_atom = {}
    for k in keys:
      if k[2]:
        continue
      by_atom.setdefault(k[0], []).append(k)
      by_atom.setdefault(k[1], []).append(k)
    for c, ks in sorted(by_atom.items()):
      for x in range(len(ks)):
        for y in range(x + 1, len(ks)):
          o1 = ks[x][1] if ks[x][0] == c else ks[x][0]
          o2 = ks[y][1] if ks[y][0] == c else ks[y][0]
          if o2 not in self._fsc0[o1] or atoms[o1].xyz == atoms[o2].xyz:
            continue
          cos = pnp.cos_vec(atoms[o1].xyz, atoms[o2].xyz, atoms[c].xyz)
          if abs(cos) > inline_clash_cos_min:
            parent[find(ks[y])] = find(ks[x])
    groups = {}
    for k in keys:
      groups.setdefault(find(k), []).append(k)
    for group in sorted(groups.values()):
      pairs = []
      for (i, j, op) in group:
        e = cl[(i, j, op)]
        pairs.append(dict(key=(i, j, op), atoms=[i, j],
          labels=[atom_label(atoms[i]), atom_label(atoms[j]) + (" (%s)" % op if op else "")],
          distance=atoms[i].distance(self._site(j, op)), sources=sorted(e["sources"]),
          pnp_kept="pnp" in e["sources"], geometry=e["geometry"]))
      kept = [q for q in pairs if q["pnp_kept"]]
      rep = min(kept or pairs, key=lambda q: q["distance"])
      for q in pairs:
        q["representative"] = q is rep
      i, j, op = rep["key"]
      sources = set()
      for q in pairs:
        sources.update(q["sources"])
      extra = None
      if len(pairs) > 1:
        extra = dict(pairs=[dict((k, v) for k, v in q.items() if k != "key")
          for q in sorted(pairs, key=lambda q: q["distance"])])
      self._add_checked("clash", [i, j], dict(sources=sources,
        geometry=dict(cl[rep["key"]]["geometry"])), ["x,y,z", op] if op else None,
        i, j, extra=extra)

  def _probe_class_of(self, i, j):
    if i is None or j is None:
      return None
    k = (i, j) if i in self._lig else (j, i)
    p = self.probe_pairs.get(k)
    return p["pair_class"] if p else None

  # -- salt bridges ------------------------------------------------------------

  def _group_info(self, g, op="x,y,z"):
    atoms = self._atoms
    suffix = "" if _identity(op) else " (%s)" % op
    return dict(kind=g["kind"], charge=g["charge"], usual_charge=g["usual_charge"],
      state=g["state"], source=g["source"], altloc=g["altloc"],
      atoms=[atom_label(atoms[i]) + suffix for i in g["charged"]],
      residue=residue_label(atoms[g["center"]], g["resname"]) + suffix,
      certain=g["certain"], metal_bound=g["metal_bound"],
      charge_source=g["charge_source"], hydrogens=g["hydrogens"], notes=list(g["notes"]))

  def _build_salt_bridges(self):
    """Salt bridges between the ligand's charged groups and the environment's (see the module docstring)."""
    sp = self.params.salt_bridge
    atoms = self._atoms
    lig = self._lig
    self.salt_bridge_criteria = dict(criterion=sp.criterion,
      atom_pair_cutoff=sp.atom_pair_cutoff, charge_centre_cutoff=sp.charge_centre_cutoff,
      kumar_nussinov_cutoff=sp.kumar_nussinov_cutoff)
    search = max(sp.atom_pair_cutoff, sp.charge_centre_cutoff + 2.5)
    # the ligand and the residues within search (symmetry included)
    region = self.model.selection("(%s) or (residues_within (%s, %s))" % (
      self.sel_str, search, self.sel_str))
    found = find_charged_groups(self.model, region, self._fsc0)
    groups = found.groups
    for k, g in enumerate(groups):
      g["id"] = k
      g["ligand"] = g["center"] in lig
    self.charged_groups = [self._group_info(g) for g in groups if g["ligand"]]
    self.examined_groups = [self._group_info(g) for g in groups]
    self.charged_group_missing_atoms = found.missing
    self.charged_groups_dropped = found.dropped
    self.charged_group_failures = found.failures
    self._builder_charges = found.builder_charges
    candidates = [c for c in charged_group_pairs(self.model, groups, search)
      if groups[c[0]]["ligand"]]
    partner_residues = set()
    for s, p, op in candidates:
      stay, partner = groups[s], groups[p]
      partner_residues.add(partner["center"])
      xs = [atoms[i].xyz for i in stay["charged"]]
      xp = [self._site(i, op) for i in partner["charged"]]
      pairs = [math.sqrt(sum([(x[c] - y[c]) ** 2 for c in range(3)])) for x in xs for y in xp]
      k = min(range(len(pairs)), key=lambda k: pairs[k])
      d_min = pairs[k]
      cs_ = [sum([x[c] for x in xs]) / len(xs) for c in range(3)]
      cp_ = [sum([x[c] for x in xp]) / len(xp) for c in range(3)]
      d_centre = math.sqrt(sum([(cs_[c] - cp_[c]) ** 2 for c in range(3)]))
      if sp.criterion == "atom_pair":
        ok = d_min <= sp.atom_pair_cutoff
      else:
        ok = d_centre <= sp.charge_centre_cutoff
      if not ok:
        continue
      kn = sp.kumar_nussinov_cutoff
      if d_min <= kn and d_centre <= kn:
        subtype = "salt bridge (K&N)"
      elif d_min <= kn:
        subtype = "N-O bridge (K&N)"
      else:
        subtype = "longer-range (K&N)"
      a_stay = stay["charged"][k // len(xp)]
      a_part = partner["charged"][k % len(xp)]
      geometry = dict(charged_groups=dict(min_atom_distance=d_min,
        charge_centre_distance=d_centre,
        closest_pair=[atom_label(atoms[a_stay]), atom_label(atoms[a_part]) +
          ("" if _identity(op) else " (%s)" % op)],
        ligand_group=self._group_info(stay), partner_group=self._group_info(partner, op)))
      if partner["ligand"] and _identity(op):
        self.internal.append(dict(type="salt_bridge", subtype=subtype,
          labels=geometry["charged_groups"]["ligand_group"]["atoms"] +
            geometry["charged_groups"]["partner_group"]["atoms"], geometry=geometry))
        continue
      seqs = stay["charged"] + partner["charged"]
      ops = ["x,y,z"] * len(stay["charged"]) + [op] * len(partner["charged"])
      entry = self._entry("salt_bridge", subtype, seqs, geometry, ["charged groups"],
        None if _identity(op) else ops)
      entry["hbonds"] = self._overlapping_hbonds(stay, partner, op)
      self.entries.append(entry)
    self._compare_formal_charges(groups, region, partner_residues)

  def _overlapping_hbonds(self, stay, partner, op):
    """H-bond entries between the two groups' charged atoms (charge-assisted H-bonds)."""
    result = []
    s, p = set(stay["charged"]), set(partner["charged"])
    for k, e in enumerate(self.entries):
      if e["type"] != "hbond":
        continue
      d, h, a = e["atoms"]
      od, oh, oa = e["operators"]
      for x, ox, y, oy in ((d, od, a, oa), (a, oa, d, od)):
        if (x in s and _identity(ox) and y in p and
            (_identity(oy) if _identity(op) else oy == op)):
          result.append(dict(index=k, labels=e["labels"]))
          break
    return result

  def _compare_formal_charges(self, groups, region, partner_centers):
    """
    Formal charges (restraint dictionary, CCD) against the perception, per residue
    and conformer: template residues per group (and dictionary-charged atoms outside
    the groups against 0), builder residues atom by atom (_atomic_charge_checks).
    """
    atoms = self._atoms
    near = set(region.iselection())
    near.update(partner_centers)
    fsc0 = self._fsc0
    done = set()
    failed = dict([((d["residue"], d["altloc"]), d["reason"])
      for d in self.charged_group_failures])
    for rg in self.model.get_hierarchy().residue_groups():
      rg_atoms = rg.atoms()
      if not [a for a in rg_atoms if a.i_seq in near]:
        continue
      alts = sorted(set([ag.altloc for ag in rg.atom_groups() if ag.altloc.strip()])) or [""]
      for alt in alts:
        ags = [ag for ag in rg.atom_groups() if ag.altloc in ("", alt)]
        conf = [a for ag in ags for a in ag.atoms()]
        resname = ([ag for ag in ags if ag.altloc.strip()] or ags)[0].resname.strip().upper()
        key = (rg.memory_id() if hasattr(rg, "memory_id") else id(rg), alt)
        if key in done:
          continue
        done.add(key)
        names = dict([(a.name.strip(), a.i_seq) for a in conf])
        conf_seqs = set(names.values())
        model_h = dict([(n, set([atoms[k].name.strip() for k in fsc0[i] if self._is_h(k)]))
          for n, i in names.items()])
        label = residue_label(conf[0]) + ((" alt " + alt) if alt.strip() else "")
        conf_groups = [g for g in groups if g["center"] in conf_seqs and
          g["altloc"] in ("", alt.strip())]
        restraints, file_name = restraints_formal_charges(self.model, resname)
        why = failed.get((residue_label(conf[0], resname), alt.strip()),
          failed.get((residue_label(conf[0], resname), "")))
        for source, d, where in (("restraints", restraints, file_name),
                                 ("CCD", ccd_formal_charges(resname), None)):
          if d is None:
            continue
          if why is not None:
            self.formal_charges.append(dict(residue=label, source=source, file=where,
              group=None, atoms=[], perceived_charge=None, dictionary_charge=None,
              status="not compared (residue_molecule failed: %s)" % why))
            continue
          bc = self._builder_charges.get((rg.parent().id, rg.resseq, rg.icode,
            alt.strip()), self._builder_charges.get((rg.parent().id, rg.resseq,
            rg.icode, "")))
          if bc is not None:
            self._atomic_charge_checks(label, source, where, d, bc, names, conf_seqs,
              model_h)
            continue
          covered = set()
          for g in conf_groups:
            seqs = [g["center"]] + [k for k in g["charged"] if k != g["center"]]
            gn = [atoms[i].name.strip() for i in seqs]
            covered.update(gn)
            self.formal_charges.append(self._charge_check(label, source, where, g["kind"],
              gn, self._local_names(seqs, conf_seqs), g["charge"], d, model_h))
          for n, (el, q) in sorted(d["atoms"].items()):
            if q != 0 and n not in covered and n in names:
              self.formal_charges.append(self._charge_check(label, source, where, None,
                [n], self._local_names([names[n]], conf_seqs), 0, d, model_h))

  def _atomic_charge_checks(self, label, source, file_name, d, charges, names, conf_seqs,
                            model_h):
    """
    A builder residue: the builder's atomic formal charges (charges {i_seq: q})
    against the dictionary's, per resonance group (same element sharing a heavy
    neighbour, as residue_molecule's differences["charges"]); groups where either
    side is charged are recorded, with the protonation rule of _charge_check.
    """
    atoms = self._atoms
    heavy = sorted([i for i in names.values() if not self._is_h(i) and
      _element(atoms[i]) not in metal_elements])
    group = dict([(i, i) for i in heavy])
    def find(i):
      while group[i] != i:
        i = group[i]
      return i
    for c in heavy:
      nbs = [k for k in self._fsc0[c] if k in group]
      for a in nbs:
        for b in nbs:
          if a < b and _element(atoms[a]) == _element(atoms[b]):
            group[find(a)] = find(b)
    sets = {}
    for i in heavy:
      sets.setdefault(find(i), []).append(i)
    for members in sorted(sets.values()):
      gn = [atoms[i].name.strip() for i in members]
      perceived = sum([charges.get(i, 0) for i in members])
      dq = sum([d["atoms"][n][1] for n in gn if n in d["atoms"]])
      if perceived == 0 and dq == 0:
        continue
      self.formal_charges.append(self._charge_check(label, source, file_name, None, gn,
        self._local_names(members, conf_seqs), perceived, d, model_h))

  def _local_names(self, seqs, conf_seqs):
    """Names of the heavy atoms within two bonds of seqs (same residue and conformer)."""
    atoms = self._atoms
    seen = set(seqs)
    shell = set(seqs)
    for step in range(2):
      shell = set([k for i in shell for k in self._fsc0[i]
        if k in conf_seqs and k not in seen and not self._is_h(k)])
      seen.update(shell)
    return sorted([atoms[i].name.strip() for i in seen])

  def _charge_check(self, residue, source, file_name, kind, names, local, perceived, d,
                    model_h):
    """
    Dictionary charge on names against the perceived charge; compared only where the
    H on the heavy atoms within two bonds (local) are those of the dictionary.
    """
    missing = [n for n in names if n not in d["atoms"]]
    if missing:
      status, q = "atoms not in dictionary: %s" % " ".join(missing), None
    else:
      q = sum([d["atoms"][n][1] for n in names])
      if [n for n in local if n in d["atoms"] and d["h"][n] != model_h[n]]:
        status = "protonation differs"
      else:
        status = "agrees" if q == perceived else "conflict"
    return dict(residue=residue, source=source, file=file_name, group=kind,
      atoms=names, perceived_charge=perceived, dictionary_charge=q, status=status)

  def formal_charge_conflicts(self):
    return [c for c in self.formal_charges if c["status"] == "conflict"]

  # -- patches -----------------------------------------------------------------

  def _build_patches(self):
    atoms = self._atoms
    lig = self._lig
    density = self.params.probe.density
    dots = [d for d in self.probe_dots if d["source"] in lig and d["target"] not in lig]
    xyz = flex.vec3_double([d["loc"] for d in dots])
    for group in dot_patches(xyz, self.params.patch_link_distance):
      counts = dict([(c, 0) for c in self.classes])
      pairs, lig_atoms, residues = {}, set(), set()
      for k in group:
        d = dots[k]
        counts[d["cls"]] += 1
        i, j = d["source"], d["target"]
        pairs[(i, j)] = pairs.get((i, j), 0) + 1
        lig_atoms.add(atom_label(atoms[i]))
        residues.add(residue_label(atoms[j]))
      center = xyz.select(flex.size_t(group)).mean()
      self.patches.append(dict(n_dots=len(group), dots=counts,
        area=dict([(c, n / density) for c, n in counts.items()]),
        area_total=len(group) / density, center=center,
        pairs=[dict(ligand=atom_label(atoms[i]), environment=atom_label(atoms[j]), dots=n)
          for (i, j), n in sorted(pairs.items(), key=lambda x: -x[1])],
        ligand_atoms=sorted(lig_atoms), residues=sorted(residues)))

  # -- summaries ---------------------------------------------------------------

  def counts(self):
    """Counts per type (subtype where set), per residue and per ligand atom."""
    per_type, per_residue, per_atom = {}, {}, {}
    for e in self.entries:
      t = e["type"] if e["subtype"] is None else "%s:%s" % (e["type"], e["subtype"])
      per_type[t] = per_type.get(t, 0) + 1
      r = e["residue"]
      per_residue.setdefault(r, {})
      per_residue[r][t] = per_residue[r].get(t, 0) + 1
      for a in e["ligand_atoms"]:
        per_atom.setdefault(a, {})
        per_atom[a][t] = per_atom[a].get(t, 0) + 1
    return group_args(per_type=per_type, per_residue=per_residue, per_ligand_atom=per_atom)

  def probe_parameters(self):
    p = self.params.probe
    return dict([(n, getattr(p, n)) for n in dir(p) if not n.startswith("_")])

  def as_dict(self):
    c = self.counts()
    return dict(ligand=self.sel_str, pair_class_order=list(self.order),
      patch_link_distance=self.params.patch_link_distance,
      use_neutron_distances=self.params.use_neutron_distances,
      separate_worse_clashes=self.params.separate_worse_clashes,
      probe=self.probe_parameters(), hbond_criteria=self.overlaps.hbond_criteria,
      clash_criteria=self.clash_criteria, salt_bridge_criteria=self.salt_bridge_criteria,
      entries=self.entries,
      counts=dict(per_type=c.per_type, per_residue=c.per_residue,
        per_ligand_atom=c.per_ligand_atom),
      disagreements=self.disagreements, internal=self.internal, patches=self.patches,
      unresolved_donors=self.unresolved_donors,
      probe_unmapped=sorted(self.probe_unmapped),
      charged_groups=self.charged_groups,
      charged_groups_dropped=self.charged_groups_dropped,
      charged_group_failures=self.charged_group_failures,
      formal_charges=self.formal_charges,
      formal_charge_conflicts=self.formal_charge_conflicts())

  def show(self, log=None):
    if log is None:
      log = sys.stdout
    c = self.counts()
    print("Ligand interactions: %s" % self.sel_str, file=log)
    print("  pair class order: %s; patch link distance %.2f A" % (
      " ".join(self.order), self.params.patch_link_distance), file=log)
    print("  pnp H-bond criteria: %s" % ", ".join(["%s=%s" % (k, v) for k, v in
      sorted(self.overlaps.hbond_criteria.items())]), file=log)
    print("  salt bridges: %s" % ", ".join(["%s=%s" % (k, v) for k, v in
      sorted(self.salt_bridge_criteria.items())]), file=log)
    print("  counts: %s" % ", ".join(["%s %d" % (k, v) for k, v in sorted(c.per_type.items())]),
      file=log)
    for e in self.entries:
      g = e["geometry"]
      detail = []
      for s in sorted(g):
        detail.append("%s(%s)" % (s, ", ".join(["%s=%s" % (k, ("%.2f" % v)
          if isinstance(v, float) else v) for k, v in sorted(g[s].items())
          if not isinstance(v, (dict, list))])))
      print("  %-6s %-3s %s  [%s] %s" % (e["type"], e["subtype"] or "",
        " ... ".join([l for l in e["labels"] if l]),
        e["cross_check"] or ", ".join(e["sources"]), " ".join(detail)), file=log)
      if e.get("pairs"):
        for q in e["pairs"]:
          print("         pair %s ... %s %.2f A [%s]%s" % (q["labels"][0], q["labels"][1],
            q["distance"], ", ".join(q["sources"]), " (pnp kept)" if q["pnp_kept"] else ""),
            file=log)
    if self.disagreements:
      print("  disagreements:", file=log)
      for d in self.disagreements:
        print("    %s %s: reported by %s, not by %s (probe class %s)" % (d["type"],
          " ... ".join([l for l in d["labels"] if l]), ", ".join(d["sources"]),
          ", ".join(d["missing"]), d["probe_class"]), file=log)
    if self.unresolved_donors:
      print("  H without exactly one bonded heavy atom (donor unset):", file=log)
      for d in self.unresolved_donors:
        print("    %s: %s" % (d["h"], ", ".join(d["bonded_heavy"]) or "none"), file=log)
    if self.internal:
      print("  ligand-internal (not entries):", file=log)
      for d in self.internal:
        print("    %s %s" % (d["type"], " ... ".join(d["labels"])), file=log)
    uncertain = [g for g in self.examined_groups if not g["certain"]]
    if uncertain:
      print("  charged groups not used for salt bridges (uncertain):", file=log)
      for g in uncertain:
        print("    %s %s %+d: %s" % (g["residue"], g["kind"], g["charge"],
          "; ".join(g["notes"])), file=log)
    if self.charged_groups_dropped:
      print("  charged groups dropped:", file=log)
      for d in self.charged_groups_dropped:
        print("    %s%s %s %+d (%s): %s" % (d["residue"], (" alt " + d["altloc"])
          if d["altloc"] else "", d["kind"], d["charge"], " ".join(d["atoms"]),
          d["reason"]), file=log)
    if self.charged_group_failures:
      print("  residues without charged groups (residue_molecule failed):", file=log)
      for d in self.charged_group_failures:
        print("    %s%s: %s" % (d["residue"], (" alt " + d["altloc"]) if d["altloc"]
          else "", d["reason"]), file=log)
    conflicts = self.formal_charge_conflicts()
    if conflicts:
      print("  formal-charge conflicts (dictionary vs perceived, same protonation):", file=log)
      for d in conflicts:
        print("    %s %s %s: %s %s, perceived %s" % (d["residue"], " ".join(d["atoms"]),
          d["group"] or "", d["source"], d["dictionary_charge"], d["perceived_charge"]),
          file=log)
    print("  contact patches (ligand -> environment dots):", file=log)
    for k, p in enumerate(self.patches):
      print("    %d: %d dots, %.1f A^2 (%s); ligand %s; residues %s" % (k + 1, p["n_dots"],
        p["area_total"], " ".join(["%s %d" % (c, p["dots"][c]) for c in self.classes
        if p["dots"][c]]), " ".join([a.split()[3] for a in p["ligand_atoms"]]),
        ", ".join(p["residues"])), file=log)
