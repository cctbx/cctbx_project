"""mmtbx.endo_exo: QM region builder with BFS expansion and hydrogen capping.

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
  .help = "Atom selection string(s) for the initial seed region(s). Each entry produces one QM region output file. \
May be specified multiple times (e.g. selection='chain A and resseq 100' selection='chain B and resseq 200'). \
If no selection is given, all metals in the structure are used as seeds (one output per metal)."
element_filter = None
  .type = str
  .multiple = True
  .help = "Restrict the element-scan seed search to atoms of these element(s) \
(e.g. element_filter=Fe, or element_filter=Fe element_filter=Cu). Element \
symbols are case-insensitive and need not be metals. Only consulted when no \
`selection` is given; if any selection is provided, the explicit selection \
wins. Default (no values): seed on every metal found by \
mmtbx.geometry_restraints.qmi.metals."
altloc = auto
  .type = str
  .help = "Which altloc letter to retain per residue.  'auto' (default) picks the \
highest mean-occupancy non-blank altloc per residue.  A specific letter (e.g. 'A', 'B') \
keeps that letter, falling back to the highest-occupancy altloc with a warning if the \
letter is absent from a residue.  'all' disables altloc filtering."
buffer {
  radius = 5.0
    .type = float(value_min=0)
    .help = "Radius of the buffer region around the selected scatterer."
  skip_search = False
    .type = bool
    .help = "If True, the initial radius search is skipped and only the seed atoms themselves seed the QM region (BFS expansion still applies)."
}
# contact_cutoff removed: its seed-contact edges let the BFS drift across the
# lattice forever; the radius search covers what it did.
# contact_cutoff = 3.0
#   .type = float(value_min=0)
#   .help = "Atoms within this distance (Angstrom) of any metal or selected atom are treated as bonded to it, even when the model has no such bond (e.g. metal-ligand coordination)."
max_search_depth = 3
  .type = int(value_min=0)
  .help = "Maximum BFS depth from any atom within the QM region."
capping {
  enable = True
    .type = bool
    .help = "Whether to perform capping of boundary atoms based on heuristics. If False, the output QM region will have uncapped dangling bonds."
  preferred_cuts = True
    .type = bool
    .help = "Consult the per-residue table of preferred cut bonds before the geometric heuristic, and allow the heuristic only on bonds nearer the backbone than the table entry."
}
include_hbond_partners = False
  .type = bool
  .help = "Whether to seed the QM region with atoms hydrogen-bonded to it, so a bond is not left with only one of its two partners. The added atoms are cut and capped like any other, so a cut can land further out and leave an atom capped that was not before. Requires hydrogens."
include_waters_in_convex_hull = True
  .type = bool
  .help = "Whether to check for water molecules inside the convex hull of the selected QM region and add them to the QM region if found."
residues_to_include
  .help = "Residues to include in the output whole, exempt from the sidechain \
cut rules. Leave 'selection' unset to disable."
{
  selection = None
    .type = str
    .help = "CCTBX selection string, e.g. 'chain A and resseq 50-100'. \
Expanded to whole residue groups, so a partial match still pulls in the \
complete residue."
  scope = *per_seed global
    .type = choice
    .help = "per_seed: add an included residue only to the region of a seed \
it lies within 'proximity' of. global: add every included residue to every \
seed region."
  proximity = 5.0
    .type = float(value_min=0)
    .help = "per_seed only. A residue is included in a seed's region if any \
of its atoms is within this distance (Angstrom) of any seed atom in that \
group. Ignored when scope=global."
}
write_files = True
  .type = bool
  .help = "If True (default), write a PDB, an mmCIF, and a sidecar PHIL \
file per seed to the current working directory. Set to False when calling \
the program in-memory (e.g. via Program(...).run() + get_results()) and \
the per-seed Model objects are consumed directly without a disk round-trip."
"""


class Program(ProgramTemplate):
  """Extract QM regions with BFS expansion and hydrogen capping.

  Seeds the QM region either from all metals in the structure (default) or
  from a user-supplied CCTBX selection string (``selection`` parameter).
  The heavy lifting is delegated to
  :class:`mmtbx.geometry_restraints.endo_exo.builder.QMRegionBuilder`.
  """

  description = '''
  Grows a QM region around each seed site by BFS, optionally caps dangling
  bonds with hydrogen atoms, and writes a PDB file, an mmCIF file, and a
  sidecar PHIL file per seed.  Seeds are all metals in the structure unless
  a custom selection string is provided via the ``selection`` parameter.
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
