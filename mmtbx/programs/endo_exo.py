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

import os

from libtbx import Auto
from libtbx.program_template import ProgramTemplate
from libtbx.utils import Sorry

from mmtbx.geometry_restraints.endo_exo.builder import (
  QMRegionBuilder, format_sidecar)


master_phil_str = """
include scope mmtbx.geometry_restraints.endo_exo.builder.master_phil_str
write_files = True
  .type = bool
  .help = "Whether to write a model file and a sidecar PHIL file per seed."
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
    """Build the QM region(s), name them, and write each one with its
    sidecar.

    Raises
    ------
    Sorry
        When two regions would be written under one file name, or an output
        file exists and overwriting is off.
    """
    model = self.data_manager.get_model()
    self._builder = QMRegionBuilder(self.params, logger=self.logger)
    self._results = self._builder.run(model)
    model_stem = os.path.splitext(
      os.path.basename(self.data_manager.get_default_model_name()))[0]
    for region in self._results:
      region.file_name = self._region_filename(model_stem, region.label)
    if not self.params.write_files:
      return
    file_names = [region.file_name for region in self._results]
    if len(set(file_names)) < len(file_names):
      raise Sorry(
        'Several QM regions would be written under one file name; set '
        'output.prefix rather than output.filename.')
    for region in self._results:
      model_file = self.data_manager.write_model_file(
        region.model, filename=region.file_name)
      phil_file = self.data_manager.write_phil_file(
        format_sidecar(region), filename=region.file_name + '.phil')
      print(f'Wrote QM region to {model_file} '
            f'({region.model.get_number_of_atoms()} atoms) '
            f'and its sidecar to {phil_file}', file=self.logger)

  def _region_filename(self, model_stem, label):
    """Return the output file name stem of one region.

    ``{model_stem}_{label}_within{radius}A_depth{depth}``, then the
    template's serial; ``output.prefix`` replaces the model stem and
    ``output.suffix`` follows the depth when they are set, and
    ``output.filename`` replaces the whole name.

    Parameters
    ----------
    model_stem : str
        The model file's name without its extension.
    label : str
        The region's label.

    Returns
    -------
    str
    """
    prefix = model_stem
    if self.params.output.prefix is not None:
      prefix = self.params.output.prefix
    suffix = (f'_{label}_within{self.params.buffer.radius:.2f}A'
              f'_depth{self.params.max_search_depth}')
    if self.params.output.suffix is not None:
      suffix += self.params.output.suffix
    return self.get_default_output_filename(prefix=prefix, suffix=suffix)

  def get_results(self):
    """Return the regions built during :meth:`run`.

    Returns
    -------
    list of libtbx.group_args
        The records :meth:`QMRegionBuilder.run` returns, each with the
        ``file_name`` this program gives it: the stem (no extension) the
        model file, with ``.pdb`` or ``.cif``, and the ``.phil`` sidecar are
        written under, set also when ``write_files`` is False.
    """
    return self._results
