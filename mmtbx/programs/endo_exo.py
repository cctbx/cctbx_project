"""mmtbx.endo_exo: QM region builder with hydrogen capping.

Grows a QM region around each seed site (metal atoms by default, or a
user-supplied selection) by breadth-first traversal of the covalent graph,
optionally caps dangling bonds with hydrogen atoms, and writes a model file
and a sidecar PHIL file per seed.

Usage (via dispatcher):

    mmtbx.development.endo_exo model.pdb
    mmtbx.development.endo_exo model.pdb selection="chain A and resseq 100" \\
      buffer.radius=5.0
    mmtbx.development.endo_exo model.pdb selection="chain A and resseq 100" \\
      selection="chain B and resseq 200"

Programmatic use::

    from mmtbx.programs.endo_exo import Program
    from iotbx.cli_parser import run_program
    results = run_program(program_class=Program, args=['model.pdb'])
"""

from libtbx import Auto
from libtbx.program_template import ProgramTemplate
from libtbx.utils import Sorry

from mmtbx.geometry_restraints.endo_exo.builder import (
  QMRegionBuilder, format_sidecar)


master_phil_str = """
selection = None
  .type = atom_selection
  .multiple = True
  .help = "Atoms that seed a QM region; each entry produces one region. If no \
selection is given, every metal in the structure seeds its own region."
element_filter = None
  .type = str
  .multiple = True
  .help = "Restrict the seed search to atoms of these element(s). Element \
symbols are case-insensitive and need not be metals. Consulted only when no \
selection is given."
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
    .type = atom_selection
    .help = "Atoms whose residue groups are included whole in a region."
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
  .help = "Whether to write a model file and a sidecar PHIL file per seed."
verbose = False
  .type = bool
  .help = "Whether to log each bond cut, cap and repair with its atoms."
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
  with hydrogen atoms, and writes a model file and a sidecar PHIL file per
  seed.  Seeds are all metals in the structure unless a custom
  selection string is provided via the selection parameter.
  '''

  datatypes = ['model', 'phil']
  master_phil_str = master_phil_str

  def custom_init(self):
    """Start with no results, so :meth:`get_results` works before
    :meth:`run`."""
    self._results = []

  def validate(self):
    """Validate user inputs before :meth:`run` is called.

    Raises
    ------
    Sorry
        When there is not exactly one model, a selection does not parse,
        altloc is not auto, all or one character, or the model has no
        hydrogens.
    """
    self.data_manager.has_models(
      raise_sorry=True, expected_n=1, exact_count=True)
    model = self.data_manager.get_model()
    n_models = model.get_number_of_models()
    if n_models > 1:
      raise Sorry(f'Multi-model files are not supported ({n_models} models).')
    selection_strings = [('selection', s)
                         for s in (self.params.selection or []) if s]
    include_selection = self.params.residues_to_include.selection
    if include_selection:
      selection_strings.append(
        ('residues_to_include.selection', include_selection))
    for name, sel_str in selection_strings:
      try:
        model.selection(sel_str)
      except Exception as e:
        raise Sorry(
          f"Invalid {name} string '{sel_str}': {e}"
        )

    altloc = self.params.altloc
    if altloc is not Auto and altloc is not None:
      if (altloc.strip().lower() not in ('auto', 'all')
          and len(altloc.strip()) != 1):
        raise Sorry(
          f"altloc must be auto, all or one altloc character, not '{altloc}'.")

    if not model.has_hd():
      raise Sorry('Model has no hydrogens. Add them and run again.')

  def run(self):
    """Build the QM region(s) and write each one with its sidecar.

    Raises
    ------
    Sorry
        When two regions would be written under one file name, or an output
        file exists and overwriting is off.
    """
    model = self.data_manager.get_model()
    self._builder = QMRegionBuilder(self.params, logger=self.logger)
    self._results = self._builder.run(
      model,
      model_name=self.data_manager.get_default_model_name(),
      default_output_filename=self._region_filename,
    )
    if not self.params.write_files:
      return
    file_names = [result['file_name'] for result in self._results]
    if len(set(file_names)) < len(file_names):
      raise Sorry(
        'Several QM regions would be written under one file name; set '
        'output.prefix rather than output.filename.')
    for result in self._results:
      model_file = self.data_manager.write_model_file(
        result['model'], filename=result['file_name'])
      phil_file = self.data_manager.write_phil_file(
        format_sidecar(result), filename=result['file_name'] + '.phil')
      print(f'Wrote QM region to {model_file} ({result["n_atoms"]} atoms) '
            f'and its sidecar to {phil_file}', file=self.logger)

  def _region_filename(self, prefix, suffix):
    """Return the output file name stem of one region.

    Parameters
    ----------
    prefix : str
        The model file's name without its extension; ``output.prefix``
        replaces it when set.
    suffix : str
        The region's own suffix; ``output.suffix`` is appended when set.

    Returns
    -------
    str
    """
    if self.params.output.prefix is not None:
      prefix = self.params.output.prefix
    if self.params.output.suffix is not None:
      suffix += self.params.output.suffix
    return self.get_default_output_filename(prefix=prefix, suffix=suffix)

  def get_results(self):
    """Return per-seed result dicts produced during :meth:`run`.

    Returns
    -------
    list of dict
        Each dict contains:

        file_name : str
            Output filename stem (no extension), unique to the region; the
            model file adds ``.pdb`` or ``.cif`` and the sidecar ``.phil``.
            Set also when ``write_files`` is False.
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
        completed_iseqs : list of int
            Sorted positional indices of the atoms added to repair a valence
            the input model left open; they are not caps.
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
