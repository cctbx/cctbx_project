"""mmtbx.endo_exo: QM region builder with hydrogen capping.

Grows a QM region around each seed site (metal atoms by default, or a
user-supplied selection) by breadth-first traversal of the covalent graph,
optionally caps dangling bonds with hydrogen atoms, and writes a PDB file, an
mmCIF file, and a sidecar PHIL file per seed.

Usage (via dispatcher):

    mmtbx.endo_exo model.pdb
    mmtbx.endo_exo model.pdb selection="chain A and resseq 100" buffer.radius=5.0
    mmtbx.endo_exo model.pdb selection="chain A and resseq 100" \\
                             selection="chain B and resseq 200"

Programmatic use::

    from mmtbx.programs.endo_exo import Program
    from iotbx.cli_parser import run_program
    run_program(program_class=Program)
"""

from __future__ import absolute_import, division, print_function

try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate

from libtbx.utils import Sorry

from mmtbx.geometry_restraints.endo_exo.builder import QMRegionBuilder


master_phil_str = """
selection = None
  .type = str
  .multiple = True
  .help = "Atom selection string for a seed region; each entry produces one \
QM region output. If no selection is given, every metal in the structure is \
used as a seed, one output per metal."
element_filter = None
  .type = str
  .multiple = True
  .help = "Restrict the seed search to atoms of these element(s). Element \
symbols are case-insensitive and need not be metals. Consulted only when no \
`selection` is given."
altloc = auto
  .type = str
  .help = "Which altloc to retain per residue. 'auto' keeps the non-blank \
altloc of highest mean occupancy, a letter keeps that altloc and falls back \
to the highest-occupancy one where the letter is absent, and 'all' retains \
every altloc."
buffer {
  radius = 5.0
    .type = float(value_min=0)
    .help = "Radius of the buffer region around the selected scatterer."
  skip_search = False
    .type = bool
    .help = "Skip the initial radius search, so the region starts from the \
seed atoms alone. Expansion along bonds still applies."
}
# contact_cutoff removed: its seed-contact edges let the BFS drift across the
# lattice forever; the radius search covers what it did.
# contact_cutoff = 3.0
#   .type = float(value_min=0)
#   .help = "Atoms within this distance (Angstrom) of any metal or selected atom are treated as bonded to it, even when the model has no such bond (e.g. metal-ligand coordination)."
max_search_depth = 3
  .type = int(value_min=0)
  .help = "Maximum number of bonds traversed outwards from any atom already \
in the QM region."
capping {
  enable = True
    .type = bool
    .help = "Whether to cap boundary atoms with hydrogens. If False, the QM \
region is written with dangling bonds."
  preferred_cuts = True
    .type = bool
    .help = "Consult the per-residue table of preferred cut bonds before the \
geometric cut heuristic."
}
include_hbond_partners = False
  .type = bool
  .help = "Whether to seed the QM region with atoms hydrogen bonded to it, so \
that a hydrogen bond is not left with only one of its two partners. The added \
seeds also shift where residues are cut, so the region can lose atoms as well \
as gain them. Requires hydrogens."
include_waters_in_convex_hull = True
  .type = bool
  .help = "Whether to add water molecules lying inside the convex hull of the \
selected QM region."
residues_to_include
  .help = "Residues to include in the output whole, exempt from the sidechain \
cut rules. Leave 'selection' unset to disable."
{
  selection = None
    .type = str
    .help = "CCTBX selection string, e.g. 'chain A and resseq 50-100'. \
Matches are expanded to whole residue groups."
  scope = *per_seed global
    .type = choice
    .help = "per_seed: add an included residue only to the region of a seed \
it lies within 'proximity' of. global: add every included residue to every \
seed region."
  proximity = 5.0
    .type = float(value_min=0)
    .help = "Distance (Angstrom) within which an included residue is added to \
a seed's region. Used only when scope=per_seed."
}
write_files = True
  .type = bool
  .help = "Whether to write a PDB, an mmCIF, and a sidecar PHIL file per seed \
to the current working directory."
"""


class Program(ProgramTemplate):
  """Extract QM regions with hydrogen capping.

  Seeds the QM region either from all metals in the structure (default) or
  from a user-supplied CCTBX selection string (``selection`` parameter).
  The work is delegated to
  :class:`mmtbx.geometry_restraints.endo_exo.builder.QMRegionBuilder`.
  """

  description = '''
  Grows a QM region around each seed site, optionally caps dangling bonds
  with hydrogen atoms, and writes a PDB file, an mmCIF file, and a sidecar
  PHIL file per seed.  Seeds are all metals in the structure unless a custom
  selection string is provided via the ``selection`` parameter.
  '''

  datatypes = ['model', 'phil']
  master_phil_str = master_phil_str

  def validate(self):
    """Validate user inputs before :meth:`run` is called."""
    if not self.data_manager.has_models():
      raise Sorry('No model provided. Please supply a PDB or mmCIF file.')
    if self.params.selection is None:
      raise Sorry('No selection="" given.' )

    model = self.data_manager.get_model()
    selection_strings = [s for s in (self.params.selection or []) if s]
    for sel_str in selection_strings:
      try:
        model.selection(sel_str)
      except Exception as e:
        raise Sorry(
          f"Invalid selection string '{sel_str}': {e}"
        )

    if not model.has_hd():
      raise Sorry('Model has no hydrogens. Add them and run again.')

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
        model : mmtbx.model.manager
            Truncated sub-model with caps applied, carrying no restraints
            manager. For in-memory use.
        seed_iseqs : list of int
            Sorted 0-based positional indices of the seed atoms inside
            ``model``.
        cap_iseqs : list of int
            Sorted 0-based positional indices of the boundary atoms inside
            ``model``.  Populated whether or not capping is enabled; with
            ``capping.enable=False`` these atoms keep their original element
            and position, leaving dangling bonds.
        cap_original_elements : list of str
            Element each cap atom carried before capping, parallel to
            ``cap_iseqs``.
        cap_anchor_iseqs : list of int
            Sorted, unique positional indices of the QM-region heavy atoms
            the caps are bonded to.  In-memory only; not written to the
            sidecar.
        sym_image_provenance : dict
            ``{i_seq: ((chain, resseq, resname, name, altloc),
            symmetry_operation_xyz)}`` for each materialized symmetry-image
            atom, so a bond to one can be restrained against its ASU parent.
            In-memory only; not written to the sidecar.
        selection_string : str or None
            The CCTBX selection string that produced this seed group, or
            ``None`` for metal-scan groups.
    """
    return self._results
