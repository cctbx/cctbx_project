"""mmtbx.endoexo — QM region builder with BFS expansion and hydrogen capping.

Grows a QM region around each seed site (metal atoms by default, or a
user-supplied selection) by breadth-first traversal of the covalent graph,
optionally caps dangling bonds with hydrogen atoms, estimates the net charge
of the surrounding region, and writes a PDB file, an mmCIF file, and a
sidecar PHIL file per seed.

Usage (via dispatcher):

    mmtbx.endoexo model.pdb
    mmtbx.endoexo model.pdb selection="chain A and resseq 100" radius=5.0
    mmtbx.endoexo model.pdb selection="chain A and resseq 100" \\
                             selection="chain B and resseq 200"

Programmatic use::

    from mmtbx.programs.endoexo import Program
    from iotbx.cli_parser import run_program
    run_program(program_class=Program)
"""

from __future__ import absolute_import, division, print_function

try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate

from libtbx.utils import Sorry

from mmtbx.geometry_restraints.endoexo.builder import QMRegionBuilder


master_phil_str = """
selection = None
  .type = str
  .multiple = True
  .help = "Atom selection string(s) for the initial seed region(s). Each entry produces one QM region output file. \
May be specified multiple times (e.g. selection='chain A and resseq 100' selection='chain B and resseq 200'). \
If no selection is given, all metals in the structure are used as seeds (one output per metal)."
metal_element = None
  .type = str
  .multiple = True
  .help = "Restrict the metal-scan seed search to atoms of these element(s) \
(e.g. metal_element=Fe, or metal_element=Fe metal_element=Cu). Element symbols \
are case-insensitive. Only consulted when no `selection` is given; if any \
selection is provided, the explicit selection wins. Default (no values): \
seed on every metal found by mmtbx.qmi.metals."
altloc = auto
  .type = str
  .help = "Which altloc letter to retain per residue.  'auto' (default) picks the \
highest mean-occupancy non-blank altloc per residue.  A specific letter (e.g. 'A', 'B') \
keeps that letter, falling back to the highest-occupancy altloc with a warning if the \
letter is absent from a residue.  'all' disables altloc filtering."
radius = 5.0
  .type = float(value_min=0)
  .help = "Radius of the buffer region around the selected scatterer."
skip_radius_search = False
  .type = bool
  .help = "If True, the initial radius search is skipped and only the seed atoms themselves seed the QM region (BFS expansion still applies)."
metal_ligand_cutoff = 3.0
  .type = float(value_min=0)
  .help = "Distance cutoff in Angstrom used to add fallback metal-ligand edges when bond proxies do not include coordination bonds."
max_depth = 3
  .type = int(value_min=0)
  .help = "Maximum BFS depth from any atom within the QM region."
use_preferred_cuts = True
  .type = bool
  .help = "Whether to use preferred cut atoms for each residue type when identifying candidate bonds for capping, instead of relying on heuristics alone."
preferred_cut_fallback = False
  .type = bool
  .help = "Only effective when use_preferred_cuts=True. When a preferred-cut \
bond ends up with both endpoints inside the region (e.g. the radius search \
seeded atoms on both sides of it, so the preferred cut can no longer be \
made), fall back to the geometric C-C heuristic to re-cut inward of it and \
trim the resulting backbone overgrowth. Never trims atoms inside 'radius'. \
No effect when use_preferred_cuts=False (the geometric heuristic already \
applies to every bond)."
include_waters_in_convex_hull = True
  .type = bool
  .help = "Whether to check for water molecules inside the convex hull of the selected QM region and add them to the QM region if found."
do_capping = True
  .type = bool
  .help = "Whether to perform capping of boundary atoms based on heuristics. If False, the output QM region will have uncapped dangling bonds."
residues_to_include = None
  .type = str
  .help = "Selection string for residues to almost always include in the output (depending on chain specified in the selection string),
  regardless of the sidechain rules. For example: 'chain A and resseq 50-100'."
include_terminal_charges = False
  .type = bool
  .help = "If True, estimate free peptide termini charges inside the truncated QM region."
write_files = True
  .type = bool
  .help = "If True (default), write a PDB, an mmCIF, and a sidecar PHIL \
file per seed to the current working directory. Set to False when calling \
the program in-memory (e.g. via Program(...).run() + get_results()) and \
the per-seed Model objects are consumed directly without a disk round-trip."
n_terminus_charge = 1
  .type = int(value_min=0, value_max=1)
  .help = "Charge assigned to each detected free N-terminus when include_terminal_charges=True."
c_terminus_charge = -1
  .type = int(value_min=-1, value_max=0)
  .help = "Charge assigned to each detected free C-terminus when include_terminal_charges=True."
"""


class Program(ProgramTemplate):
  """Extract QM regions with BFS expansion and hydrogen capping.

  Seeds the QM region either from all metals in the structure (default) or
  from a user-supplied CCTBX selection string (``selection`` parameter).
  The heavy lifting is delegated to
  :class:`mmtbx.geometry_restraints.endoexo.builder.QMRegionBuilder`.
  """

  description = '''
  Grows a QM region around each seed site by BFS, optionally caps dangling
  bonds with hydrogen atoms, estimates the net charge, and writes a PDB
  file, an mmCIF file, and a sidecar PHIL file per seed.  Seeds are all
  metals in the structure unless a custom selection string is provided
  via the ``selection`` parameter.
  '''

  datatypes = ['model', 'phil']
  master_phil_str = master_phil_str

  def validate(self):
    """Validate user inputs before :meth:`run` is called."""
    if not self.data_manager.has_models():
      raise Sorry('No model provided. Please supply a PDB or mmCIF file.')

    model = self.data_manager.get_model()
    selection_strings = [s for s in (self.params.selection or []) if s]
    for sel_str in selection_strings:
      try:
        model.selection(sel_str)
      except Exception as e:
        raise Sorry(
          f"Invalid selection string '{sel_str}': {e}"
        )

  def run(self):
    """Build the QM region(s) by delegating to :class:`QMRegionBuilder`."""
    model = self.data_manager.get_model()
    self._builder = QMRegionBuilder(self.params, logger=self.logger)
    self._results = self._builder.run(
      model,
      model_name=self.data_manager.get_default_model_name(),
      default_output_filename=self.get_default_output_filename,
    )

  def get_results(self):
    """Return per-seed result dicts produced during :meth:`run`.

    Returns
    -------
    list of dict
        Each dict contains:

        file_name : str
            Output filename stem (no extension); the PDB and mmCIF copies
            are written as ``file_name + '.pdb'`` and ``file_name + '.cif'``.
        n_atoms : int
            Number of atoms in the QM region.
        charge_summary : dict
            As returned by :meth:`ChargeEstimator.calculate`.
        model : mmtbx.model.manager
            Truncated sub-model with caps applied and the parent's
            restraints manager attached. Suitable for direct in-memory
            consumption by downstream tools (e.g., ``mmtbx.qmi``) without
            the disk round-trip through the written PDB/mmCIF.
        seed_iseqs : list of int
            Sorted 0-based positional indices of the seed atoms inside
            ``model``.
        cap_iseqs : list of int
            Sorted 0-based positional indices of the cap atoms inside
            ``model`` (empty when capping is disabled).
        selection_string : str or None
            The CCTBX selection string that produced this seed group, or
            ``None`` for metal-scan groups.
    """
    return self._results
