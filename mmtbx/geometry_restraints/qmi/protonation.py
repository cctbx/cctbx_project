"""In-pipeline protonation for the qmi metal-restraint Program.

The QM cluster needs hydrogens, and metal-coordinating sidechain donors need
the *right* protonation (reduce2 protonates every sidechain to its neutral
biological-pH form, which is wrong at a metal site). This module ports the
protonation prep from the calibration ``prepare_qm_calcs.py``:

- ``add_hydrogens``: run ``mmtbx.reduce2`` on the model to add H.
- ``deprotonate_metal_donors``: fix the protonation of donors coordinating any
  metal. Cys/Tyr/Asp/Glu carboxylate donors are stripped to their anionic
  form; His is handled tautomer-aware (it stays NEUTRAL, not an imidazolate):
  the coordinating ring N is made H-free while the partner ring N keeps a
  proton. A coordinating Cu/Ni N-terminus is deprotonated to its ATCUN form: a
  backbone amide becomes an amidate and the free alpha-amino terminus becomes
  -NH2. The RULES for all of this are
  :mod:`mmtbx.geometry_restraints.qmi.metal_protonation`; what stays here is
  qmi's policy (which metals are in scope, which cutoff) and its reporting.
- ``prepare_model``: ``add_hydrogens``, then ``deprotonate_metal_donors``,
  returning a fresh protonated model.

The Program's validate() guard uses ``mmtbx.model.manager.has_hd()`` to decide
whether the model is already protonated.
"""

from libtbx.utils import null_out
from iotbx.pdb import common_residue_names_get_class

from mmtbx.geometry_restraints.qmi import metal_protonation
from mmtbx.geometry_restraints.qmi import metals as qmi_metals
from mmtbx.geometry_restraints.qmi import util


# The metal-donor protonation RULES themselves live in metal_protonation, not
# here: they are chemistry, not qmi pipeline, and the split keeps them callable
# by a protonation program that wants to place the right hydrogens in the first
# place instead of qmi removing hydrogens that were just added (that module
# depends on nothing else in qmi). What stays in this module is the
# orchestration: which metals qmi scopes to, the coordination cutoff it shares
# with the rest of the pipeline, the log format, and the atom re-sort its
# caller depends on.


def add_hydrogens(model, water=None, use_neutron_distances=True, log=None):
  """Add hydrogens to *model* with reduce2 and return the H-completed
  ``mmtbx.model.manager``; *model* itself is left as it was. Neutron X-H
  distances by default; *water* sets reduce2's water H placement (None: its
  defaults).
  """
  from iotbx.cli_parser import CCTBXParser
  from mmtbx.programs import reduce2
  if log is None:
    log = null_out()
  print(f"  adding hydrogens: reduce2 {reduce2.version}, "
        f"use_neutron_distances={use_neutron_distances}", file=log)
  parser = CCTBXParser(program_class=reduce2.Program, logger=null_out())
  args = [f"use_neutron_distances={use_neutron_distances}",
          "add_h_to_water=True", "output.write_files=False"]
  if water is not None:
    args += [f"water_hydrogens.refine.max_sweeps={water.refine.max_sweeps}",
             f"water_hydrogens.refine.tolerance={water.refine.tolerance}",
             f"water_hydrogens.basin.rounds={water.basin.rounds}"]
  parser.parse_args(args)
  parser.data_manager.add_model("qmi", model)
  program = reduce2.Program(
    parser.data_manager, parser.working_phil.extract(), logger=null_out())
  program.validate()
  program.run()
  return program.get_results().model


# ---------------------------------------------------------------------------
def _describe(change):
  """The qmi log sentence for one metal-protonation change.

  Presentation lives here, not in the rule module: a rule that returns English
  cannot be consumed by anything but a log."""
  donor = change.donor_atom.name.strip()
  h_name = change.h_atom.name.strip()
  if change.action == metal_protonation.MOVE_TO_PARTNER:
    partner = metal_protonation.HIS_PARTNER[donor]
    return (f"relocated {h_name}({donor}) to {change.new_name}({partner}) "
            f"wrong tautomer, d={change.distance:.2f} A")
  if change.rule == "his-hip":
    return (f"removed {h_name} from {donor} (HIP to neutral, "
            f"d={change.distance:.2f} A)")
  return f"removed {h_name}"


# ---------------------------------------------------------------------------
# metal-donor deprotonation
# ---------------------------------------------------------------------------
def deprotonate_metal_donors(hier, metals=None,
                             coord_cutoff=None,
                             log=None):
  """Fix, IN PLACE, the protonation of donors in *hier* coordinating any metal
  atom (within ``coord_cutoff(metal_element, donor_element)``).

  The rules are :mod:`mmtbx.geometry_restraints.qmi.metal_protonation`; this is
  the wrapper that supplies the pipeline's own policy and reporting. It decides
  three things the rule module deliberately refuses to assume:

  - which metals count (*metals*, default all METALS). The in-pipeline caller
    passes the metals qmi may model, so a donor near an entirely out-of-scope
    metal is left untouched. That is a scoping choice for a tool that models
    only some metals, not a statement about the chemistry.
  - what counts as coordinating (*coord_cutoff*, default the shared qmi
    cutoff), so protonation, first-shell detection and the QM region cannot
    silently disagree about the same shell.
  - that the hierarchy is re-sorted afterwards. qmi's caller immediately builds
    a fresh model from it; a caller holding i_seq-indexed structures must not
    have its atoms reordered, so the rule module never does this itself.

  Returns ``(n_donors_corrected, n_his_corrected)``. A no-op (0, 0) when the
  model has no in-scope metals.
  """
  if log is None:
    log = null_out()
  if coord_cutoff is None:
    coord_cutoff = util.make_coordination_cutoff()
  metal_atoms = qmi_metals.metal_atoms(hier.atoms(), metals=metals)
  if not metal_atoms:
    return 0, 0

  def _is_standard_aa(resname):
    return common_residue_names_get_class(resname) in util.AA_RESIDUE_CLASSES

  changes = metal_protonation.propose(
    hier, metal_atoms, coord_cutoff,
    is_standard_aa=_is_standard_aa,
    normalize_element=qmi_metals.normalize_element)
  metal_protonation.apply_changes(hier, changes)

  # Canonicalise atom order (a relocated His proton is renamed in place).
  hier.sort_atoms_in_place()

  his = [e for e in changes if e.rule.startswith("his-")]
  donors_fixed = [e for e in changes if not e.rule.startswith("his-")]
  if donors_fixed:
    print(f"  deprotonated {len(donors_fixed)} metal-coordinating donor(s):",
          file=log)
    for e in donors_fixed:
      ch, rs = e.h_id[0], e.h_id[1]
      print(f"    {e.h_id[4]} {ch}{rs} {e.donor_atom.name.strip()} "
            f"(d={e.distance:.2f} A): {_describe(e)}", file=log)
  if his:
    print(f"  corrected {len(his)} metal-coordinating His:", file=log)
    for e in his:
      print(f"    HIS {e.h_id[0]}{e.h_id[1]}: {_describe(e)}", file=log)
  return len(donors_fixed), len(his)


def prepare_model(model, metals=None, water=None, coord_cutoff=None,
                  log=None):
  """Protonate a raw model for qmi and fix metal-donor protonation, returning a
  fresh protonated model ready for cluster extraction.

  *metals* and *coord_cutoff* scope the deprotonation (see
  :func:`deprotonate_metal_donors`); *water* goes to :func:`add_hydrogens`.
  """
  import mmtbx.model
  if log is None:
    log = null_out()
  model = add_hydrogens(model, water=water, log=log)
  hier = model.get_hierarchy()
  deprotonate_metal_donors(hier, metals=metals, coord_cutoff=coord_cutoff,
                           log=log)
  return mmtbx.model.manager(
    model_input=None, pdb_hierarchy=hier,
    crystal_symmetry=model.crystal_symmetry())
