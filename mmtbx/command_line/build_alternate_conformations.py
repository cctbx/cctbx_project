
# XXX has phenix dependency (imports inline)

from __future__ import division
from __future__ import print_function
from mmtbx.building.alternate_conformations import single_residue
from mmtbx.building.alternate_conformations import sliding_window
from mmtbx.building import alternate_conformations
from mmtbx.disorder import analyze_model
from libtbx.str_utils import make_header, make_sub_header
from libtbx.utils import Sorry, create_run_directory
from libtbx import Auto, adopt_init_args
import shutil
import time
import os
import sys

def master_phil():
  from mmtbx.command_line import generate_master_phil_with_inputs
  return generate_master_phil_with_inputs(
    phil_string="""
alt_confs {
  selection = None
    .type = atom_selection
  macro_cycles = 1
    .type = int(value_min=1)
    .optional = False
  include scope libtbx.easy_mp.parallel_phil_str_no_threading
  include scope mmtbx.building.alternate_conformations.single_residue.master_phil_str
#  sliding_window {
#    include scope mmtbx.building.alternate_conformations.sliding_window.master_params_str
#  }
  refinement {
    include scope phenix.automation.refinement.refine_hires_phil_str
    constrain_correlated_occupancies = True
      .type = bool
  }
  merging {
    include scope mmtbx.building.alternate_conformations.rejoin_phil
  }
}
output {
  prefix = alternates
    .type = str
  output_dir = None
    .type = path
  create_dir = True
    .type = bool
  directory_number = None
    .type = int
  debug = 0
    .type = int
  verbose = True
    .type = bool
  remove_hydrogens = False
    .type = bool
}
""")

class build_and_refine(object):
  def __init__(self,
      fmodel,
      pdb_hierarchy,
      params=None,
      processed_pdb_file=None,
      geometry_restraints_manager=None,
      cif_objects=(),
      cif_files=(), # XXX bug
      debug=None,
      verbose=True,
      out=sys.stdout):
    adopt_init_args(self, locals())
    if (self.params is None):
      self.params = master_phil.extract().alt_confs
    self.extract_selection()
    self.refine_cycle = 1
    self.map_file = None
    self.r_work_start = fmodel.r_work()
    self.r_free_start = fmodel.r_free()
    t_start = time.time()
    for i_cycle in range(params.macro_cycles):
      n_alts = self.build_residue_conformers(stop_if_none=(i_cycle==0))
      if (n_alts == 0):
        if (i_cycle == 0):
          raise Sorry("No alternate conformations found.")
      else :
        self.refine(constrain_occupancies=False)
        refine_again = self.params.refinement.constrain_correlated_occupancies
        if (self.rejoin()):
          refine_again = True
        self.refine(title="Refining final model")
    make_header("Finished", out=out)
    from mmtbx.validation import molprobity
    validation = molprobity.molprobity(
      pdb_hierarchy=self.pdb_hierarchy,
      outliers_only=False)
    print("", file=self.out)
    validation.show_summary(out=self.out, prefix="  ")
    make_sub_header("Analyzing final model", out=out)
    analyze_model.process_pdb_hierarchy(
      pdb_hierarchy=self.pdb_hierarchy,
      validation=validation,
      log=self.out).show(out=out, verbose=self.verbose)
    print("", file=self.out)
    print("Start:  r_work=%6.4f  r_free=%6.4f" % \
      (self.r_work_start, self.r_free_start), file=self.out)
    print("Final:  r_work=%6.4f  r_free=%6.4f" % \
      (self.fmodel.r_work(), self.fmodel.r_free()), file=self.out)
    t_end = time.time()
    print("", file=self.out)
    print("Total runtime: %d s" % int(t_end - t_start), file=self.out)
    print("", file=self.out)

  def extract_selection(self):
    self.selection = None
    if (self.params.selection is not None):
      sele_cache = self.pdb_hierarchy.atom_selection_cache()
      self.selection = sele_cache.selection(self.params.selection)
      assert (self.selection.count(True) > 0)

  def build_residue_conformers(self, stop_if_none=False):
    self.extract_selection()
    print("", file=self.out)
    #self.fmodel.info().show_targets(out=self.out, text="starting model")
    make_sub_header("Fitting individual residues", out=self.out)
    t1 = time.time()
    params = self.params
    self.pdb_hierarchy, n_alternates = single_residue.build_cycle(
      pdb_hierarchy = self.pdb_hierarchy,
      fmodel = self.fmodel,
      geometry_restraints_manager = self.geometry_restraints_manager,
      params = params,
      cif_objects=self.cif_objects,
      selection=params.selection,
      nproc=params.nproc,
      verbose=self.verbose,
      debug=self.debug,
      out=self.out)
    if (n_alternates == 0) and (stop_if_none):
      raise Sorry("No new conformations generated.")
    return n_alternates

  def build_window_conformers(self, stop_if_none=False):
    self.extract_selection()
    print("", file=self.out)
    self.fmodel.info().show_targets(out=self.out, text="starting model")
    make_header("Sampling sliding windows", out=self.out)
    t1 = time.time()
    driver = sliding_window.fragment_refinement_driver(
      fmodel=self.fmodel,
      pdb_hierarchy=self.pdb_hierarchy,
      processed_pdb_file=self.processed_pdb_file,
      params=self.params.sliding_window,
      mp_params=self.params,
      selection=self.selection,
      cif_objects=self.cif_objects,
      debug=self.debug,
      verbose=self.verbose,
      out=self.out)
    t2 = time.time()
    print("sampling time: %.3fs" % (t2-t1), file=self.out)
    n_ensembles = driver.n_ensembles()
    if (n_ensembles == 0) and (stop_if_none):
      raise Sorry("No new conformations generated.")
    self.pdb_hierarchy = driver.assemble(out=self.out)
    self.processed_pdb_file = None # needs to be reset
    return n_ensembles

  # XXX should self.selection also apply here and in rejoin()?
  def refine(self, title="Refining multi-conformer model",
      constrain_occupancies=Auto):
    make_sub_header(title, out=self.out)
    t1 = time.time()
    extra_args = []
    if constrain_occupancies :
      if (self.params.refinement.constrain_correlated_occupancies):
        extra_args.append("constrain_correlated_3d_groups=True")
    else :
      print("  Correlated occupancies will *not* be constrained", file=self.out)
    from phenix.automation import refinement
    refined = refinement.refine_hires_simple(
      pdb_hierarchy=self.pdb_hierarchy,
      crystal_symmetry=self.fmodel.xray_structure,
      fmodel=self.fmodel,
      params=self.params.refinement,
      cif_files=self.cif_files,
      cycle=self.refine_cycle,
      extra_args=extra_args,
      out=self.out) # TODO need a verbosity flag
    t2 = time.time()
    print("  refinement time: %.3fs" % (t2-t1), file=self.out)
    print("", file=self.out)
    self.pdb_hierarchy = refined.pdb_hierarchy
    self.fmodel = refined.fmodel
    self.fmodel.info().show_targets(out=self.out, text="refined model")
    self.map_file = refined.map_file
    self.refine_cycle += 1

  def rejoin(self):
    make_sub_header("Re-joining identical conformers", out=self.out)
    pdb_hierarchy = self.pdb_hierarchy.deep_copy()
    n_modified = alternate_conformations.rejoin_split_single_conformers(
      pdb_hierarchy=pdb_hierarchy,
      crystal_symmetry=self.fmodel.xray_structure,
      model_error_ml=self.fmodel.model_error_ml(),
      params=self.params.merging,
      reset_occupancies=self.params.refinement.constrain_correlated_occupancies,
      verbose=self.verbose,
      log=self.out)
    if (n_modified > 0):
      self.pdb_hierarchy = pdb_hierarchy
      xray_structure = self.pdb_hierarchy.extract_xray_structure(
        crystal_symmetry=self.fmodel.xray_structure)
      self.fmodel.update_xray_structure(xray_structure)
      self.map_file = None
    alternate_conformations.finalize_model(
      pdb_hierarchy=self.pdb_hierarchy,
      xray_structure=self.pdb_hierarchy.extract_xray_structure(
        crystal_symmetry=self.fmodel.xray_structure),
      set_b_iso=None,
      convert_to_isotropic=False)
    return (n_modified > 0)

  def write_pdb_file(self, file_name, remove_hd=False):
    if (remove_hd):
      self.pdb_hierarchy.remove_hd()
    self.pdb_hierarchy.write_pdb_file(
      file_name=file_name,
      crystal_symmetry=self.fmodel.xray_structure.crystal_symmetry())
    print("wrote model to %s" % file_name, file=self.out)

  def write_map_file(self, file_name):
    import mmtbx.maps.utils
    import iotbx.map_tools
    if (self.map_file is None):
      two_fofc_coeffs, fofc_coeffs = mmtbx.maps.utils.get_maps_from_fmodel(
        self.fmodel)
      iotbx.map_tools.write_map_coeffs(
        fwt_coeffs=two_fofc_coeffs,
        delfwt_coeffs=fofc_coeffs,
        file_name=file_name)
    else :
      shutil.copyfile(self.map_file, file_name)

def run(args, out=None, driver_class=None,
    require_single_conformer_starting_model=False):
  import mmtbx.utils
  if (out is None) : out = sys.stdout
  usage_string = """phenix.alternator model.pdb data.mtz [selection=...] [options]"""
  cmdline = mmtbx.command_line.load_model_and_data(
    args=args,
    master_phil=master_phil(),
    process_pdb_file=True,
    usage_string=usage_string,
    out=out,
    create_log_buffer=True)
  params = cmdline.params
  if (params.output.output_dir is not None):
    os.chdir(params.output.output_dir)
  dir_name = os.getcwd()
  if (params.output.create_dir):
    dir_name = create_run_directory("alt_confs",
      default_directory_number=params.output.directory_number)
    os.chdir(dir_name)
  log = cmdline.start_log_file("%s.log" % params.output.prefix)
  print("Output will be in %s" % dir_name, file=log)
  #working_phil = master_phil.format(python_object=params)
  #make_sub_header("Final input parameters", out=log)
  #master_phil.fetch_diff(source=working_phil).show(out=log)
  if (driver_class is None):
    driver_class = build_and_refine
  multi_conf_selection = alternate_conformations.multi_conformer_selection(
    pdb_hierarchy=cmdline.pdb_hierarchy)
  if ((len(multi_conf_selection) != 0) and
      require_single_conformer_starting_model):
    atoms = cmdline.pdb_hierarchy.select(multi_conf_selection).atoms()
    print("First %d atoms with alternate conformations:" % min(10,
      len(atoms)), file=log)
    for atom in atoms[0:10] :
      print(atom.format_atom_record(), file=log)
    raise Sorry("Existing alternate conformations detected - this program "+
      "can only be run on a single-conformer model at present.")
  driver = driver_class(
    params=params.alt_confs,
    fmodel=cmdline.fmodel,
    pdb_hierarchy=cmdline.pdb_hierarchy,
    processed_pdb_file=cmdline.processed_pdb_file,
    geometry_restraints_manager=cmdline.geometry,
    cif_objects=[ o for (f,o) in cmdline.cif_objects ],
    cif_files=params.input.monomers.file_name,
    debug=params.output.debug,
    verbose=params.output.verbose,
    out=log)
  output_file_base = os.path.join(os.getcwd(), params.output.prefix)
  driver.write_pdb_file(output_file_base + ".pdb",
    remove_hd=params.output.remove_hydrogens)
  driver.write_map_file(output_file_base + ".mtz")
  # TODO final result object
  return driver.pdb_hierarchy

if (__name__ == "__main__"):
  try :
    import phenix.automation.refinement
  except ImportError :
    print("phenix is required to run this program.", file=sys.stderr)
  else :
    run(sys.argv[1:])
