
"""
This module provides a way to identify a suitably isomorphous initial model
for directly refining against experimental data.  This is commonly done when
solving a series of ligand-bound structures, partly to save time, and partly
to avoid changing the frame of reference when running molecular replacement.
"""

from __future__ import absolute_import, division, print_function
from libtbx.str_utils import make_sub_header
from libtbx import slots_getstate_setstate
from libtbx.utils import null_out, Sorry
from libtbx import easy_mp
import libtbx.phil
import os.path
import sys
from six.moves import zip

master_phil = libtbx.phil.parse("""
d_min = None
  .type = float
rigid_body_refine = True
  .type = bool
optimize_b_factors = False
  .type = bool
max_cell_angle_rmsd=1.0
  .type = float
max_cell_edge_rmsd=1.0
  .type = float
max_r_free = 0.4
  .type = float
""")

class evaluate_model(slots_getstate_setstate):
  """
  Create an fmodel object (including bulk solvent correction) and calculate
  R-factors, with or without optional rigid-body refinement.
  """
  __slots__ = [ "r_work", "r_free", "r_work_start", "r_free_start",
    "xray_structure" ]
  def __init__(self,
      xray_structure,
      pdb_hierarchy,
      f_obs,
      r_free_flags,
      rigid_body_refine=False,
      optimize_b_factors=False,
      skip_twin_detection=False,
      scattering_table="n_gaussian"):
    self.r_work = None
    self.r_free = None
    self.xray_structure = None
    from mmtbx.utils import fmodel_simple
    from cctbx import crystal
    combined_symmetry = crystal.symmetry(
      unit_cell=f_obs.unit_cell(),
      space_group=xray_structure.space_group())
    xray_structure = xray_structure.customized_copy(
      crystal_symmetry=combined_symmetry)
    f_obs = f_obs.customized_copy(
      crystal_symmetry=combined_symmetry).eliminate_sys_absent()
    r_free_flags = r_free_flags.customized_copy(
      crystal_symmetry=combined_symmetry).eliminate_sys_absent()
    fmodel = fmodel_simple(
      f_obs=f_obs,
      r_free_flags=r_free_flags,
      xray_structures=[xray_structure],
      skip_twin_detection=skip_twin_detection,
      scattering_table=scattering_table)
    self.r_work_start = fmodel.r_work()
    self.r_free_start = fmodel.r_free()
    if (not rigid_body_refine):
      self.r_work = self.r_work_start
      self.r_free = self.r_free_start
      self.xray_structure = xray_structure
    else :
      from mmtbx.refinement import rigid_body
      selection_strings = rigid_body.rigid_groups_from_pdb_chains(
        pdb_hierarchy=pdb_hierarchy,
        xray_structure=xray_structure,
        group_all_by_chain=True,
        check_for_atoms_on_special_positions=True,
        log=null_out())
      selections = []
      for sele_str in selection_strings :
        sele = pdb_hierarchy.atom_selection_cache().selection(sele_str)
        selections.append(sele.iselection())
      refined = rigid_body.manager(
        fmodel=fmodel,
        selections=selections,
        params=rigid_body.master_params.extract(),
        log=null_out())
      self.xray_structure = refined.fmodel.xray_structure
      self.r_work = refined.fmodel.r_work()
      self.r_free = refined.fmodel.r_free()

  def show(self, out=sys.stdout, prefix=""):
    print(prefix + "r_work = %6.4f" % self.r_work, file=out)
    print(prefix + "r_free = %6.4f" % self.r_free, file=out)

def ucf(unit_cell):
  return "%g %g %g %g %g %g" % unit_cell.parameters()

class select_model(object):
  def __init__(self,
      model_names,
      model_data,
      f_obs,
      r_free_flags,
      params=None,
      skip_twin_detection=False,
      nproc=1,
      log=sys.stdout):
    if (params is None):
      params = master_phil.extract()
    self.model_names = model_names
    if (model_data is None):
      import iotbx.pdb
      model_data = []
      for file_name in model_names :
        if (not os.path.isfile(file_name)):
          raise RuntimeError("model_data is None, but %s is not a file." %
            file_name)
        model_in = iotbx.pdb.input(file_name)
        pdb_hierarchy = model_in.construct_hierarchy()
        xray_structure = model_in.xray_structure_simple()
        model_data.append((pdb_hierarchy, xray_structure))
    self.model_symmetries = []
    self.models_accepted = []
    self.model_r_frees = []
    self.f_obs = f_obs.resolution_filter(d_min=params.d_min)
    self.r_free_flags = r_free_flags.common_set(other=self.f_obs)
    self.skip_twin_detection = skip_twin_detection
    self.params = params
    self.evaluations = None
    self.best_xray_structure = None
    self.best_pdb_hierarchy = None
    self.best_result = None
    self.best_model_name = None
    from mmtbx.pdb_symmetry import rms_difference
    data_symmetry = f_obs.crystal_symmetry()
    data_space_group = data_symmetry.space_group()
    data_point_group = data_space_group.build_derived_point_group()
    data_unit_cell = data_symmetry.unit_cell()
    data_cell_edges = data_unit_cell.parameters()[0:3]
    data_cell_angles = data_unit_cell.parameters()[3:6]
    make_sub_header("Evaluating models", out=log)
    print("Experimental data:", file=log)
    print("  space group:  %s" % data_space_group.info(), file=log)
    print("  unit cell:    %s" % ucf(data_unit_cell), file=log)
    pdb_hierarchies = []
    xray_structures = []
    for k, file_name in enumerate(model_names):
      pdb_hierarchy, xray_structure = model_data[k]
      pdb_hierarchy.atoms().reset_i_seq()
      pdb_hierarchies.append(pdb_hierarchy)
      model_symmetry = xray_structure.crystal_symmetry()
      self.model_symmetries.append(model_symmetry)
      if (model_symmetry is None):
        print("Model %d is missing symmetry records:" % (k+1), file=log)
        print("  source:  %s" % file_name, file=log)
        xray_structures.append(None)
        continue
      model_unit_cell = model_symmetry.unit_cell()
      model_space_group = model_symmetry.space_group()
      is_compatible_sg = False
      if (model_space_group == data_space_group):
        is_compatible_sg = True
      else :
        model_point_group = model_space_group.build_derived_point_group()
        if (data_point_group == model_point_group):
          is_compatible_sg = True
      if (not is_compatible_sg):
        print("Model %d has incompatible space group:" % (k+1), file=log)
        print("  source:  %s" % file_name, file=log)
        print("  space group: %s" % model_space_group.info(), file=log)
        xray_structures.append(None)
        continue
      is_similar_cell = False
      if (model_unit_cell.is_similar_to(data_unit_cell)):
        is_similar_cell = True
      else :
        model_cell_edges = model_unit_cell.parameters()[0:3]
        model_cell_angles = model_unit_cell.parameters()[3:6]
        cell_edge_rmsd = rms_difference(model_cell_edges, data_cell_edges)
        cell_angle_rmsd = rms_difference(model_cell_angles, data_cell_angles)
        if ((cell_edge_rmsd <= params.max_cell_edge_rmsd) and
            (cell_angle_rmsd <= params.max_cell_angle_rmsd)):
          is_similar_cell = True
      if (not is_similar_cell):
        print("Model %d has incompatible space group:" % (k+1), file=log)
        print("  source: %s" % file_name, file=log)
        print("  model:  %s" % ucf(model_unit_cell), file=log)
        xray_structures.append(None)
        continue
      else :
        xray_structures.append(xray_structure)
    if (xray_structures.count(None) != len(xray_structures)):
      print("", file=log)
      print("Calculating R-factors - will use %s processors." % nproc, file=log)
      evaluations = easy_mp.parallel_map(
        func=self.evaluate_model,
        iterable=list(zip(xray_structures, pdb_hierarchies)),
        processes=nproc)
      passed = []
      for k, result in enumerate(evaluations):
        if (result is not None):
          if (result.r_free <= params.max_r_free):
            passed.append((k, result))
      if (len(passed) > 0):
        passed.sort(key=lambda element: element[1].r_free)
        i_result, result = passed[0]
        self.evaluations = passed
        self.best_xray_structure = result.xray_structure
        self.best_pdb_hierarchy = pdb_hierarchies[i_result]
        self.best_result = result
        self.best_model_name = self.model_names[i_result]
    self.show(out=log, verbose=True)

  def show(self, out=sys.stdout, verbose=False):
    if (self.best_result is None):
      print("No models accepted - will need to run MR.", file=out)
    else :
      print("", file=out)
      print("Best starting model:", file=out)
      print("  source: %s" % self.best_model_name, file=out)
      self.best_result.show(out=out, prefix="    ")
      print("", file=out)
      if (verbose) and (len(self.evaluations) > 1):
        print("Other suitable models:", file=out)
        for i_other, other in self.evaluations[1:] :
          print("  source: %s" % self.model_names[i_other], file=out)
          other.show(out=out, prefix="    ")
        print("", file=out)

  def success(self):
    return self.best_pdb_hierarchy is not None

  def r_free(self):
    return getattr(self.best_result, "r_free", None)

  def r_work(self):
    return getattr(self.best_result, "r_work", None)

  def get_best_model(self, update_structure=True):
    if (self.best_result is None):
      return None
    xray_structure = self.best_xray_structure
    pdb_hierarchy = self.best_pdb_hierarchy
    if (update_structure):
      pdb_hierarchy.adopt_xray_structure(xray_structure)
    return xray_structure, pdb_hierarchy

  def save_best_model(self, file_name="best_model.pdb"):
    assert self.success()
    xray_structure, pdb_hierarchy = self.get_best_model()
    f = open(file_name, "w")
    f.write("REMARK original PDB file:\n")
    f.write("REMARK   %s\n" % self.best_model_name)
    f.write(pdb_hierarchy.as_pdb_string(crystal_symmetry=xray_structure))
    f.close()

  def save_updated_data(self, file_name="best_model_data.mtz"):
    assert self.success()
    xray_structure, pdb_hierarchy = self.get_best_model()
    f_obs = self.f_obs.customized_copy(
      crystal_symmetry=xray_structure).eliminate_sys_absent()
    r_free_flags = self.r_free_flags.customized_copy(
      crystal_symmetry=xray_structure).eliminate_sys_absent()
    mtz_data = f_obs.as_mtz_dataset(column_root_label="F")
    mtz_data.add_miller_array(r_free_flags, column_root_label="FreeR_flag")
    mtz_data.mtz_object().write(file_name)

  def evaluate_model(self, args):
    xray_structure, pdb_hierarchy = args
    if (xray_structure is None):
      return None
    return evaluate_model(
      xray_structure=xray_structure,
      pdb_hierarchy=pdb_hierarchy,
      f_obs=self.f_obs,
      r_free_flags=self.r_free_flags,
      rigid_body_refine=self.params.rigid_body_refine,
      skip_twin_detection=self.skip_twin_detection)

  def space_group_info(self):
    xray_structure, pdb_hierarchy = self.get_best_model()
    return xray_structure.space_group_info()

strip_model_params = """
  remove_waters = True
    .type = bool
    .help = Remove all water molecules (HOH)
  remove_hydrogens = True
    .type = bool
    .help = Remove explicit hydrogen atoms
  remove_alt_confs = True
    .type = bool
    .help = Remove alternate conformations
  convert_semet_to_met = True
    .type = bool
    .help = Change MSE residues to MET
  convert_to_isotropic = True
    .type = bool
    .help = Convert atoms to anisotropic
  reset_occupancies = True
    .type = bool
    .help = Set occupancies to 1.0
  remove_ligands = False
    .type = bool
    .help = Remove all ligands
  reset_hetatm_flag = False
    .type = bool
    .help = Change HETATM records to ATOM
"""

def strip_model(
    pdb_hierarchy=None,
    xray_structure=None,
    file_name=None,
    params=None,
    remove_waters=True,
    remove_hydrogens=True,
    remove_alt_confs=True,
    convert_semet_to_met=True,
    convert_to_isotropic=True,
    reset_occupancies=True,
    remove_ligands=False,
    reset_hetatm_flag=False,
    preserve_remarks=False,
    preserve_symmetry=True,
    add_remarks=None,
    output_file=None,
    log=None):
  """
  Utility for removing extraneous records from a model intended for use in
  molecular replacement, etc., including waters, alternate conformations,
  and other features specific to a particular dataset.
  """
  if (params is not None):
    remove_waters = params.remove_waters
    remove_hydrogens = params.remove_hydrogens
    remove_alt_confs = params.remove_alt_confs
    convert_semet_to_met = params.convert_semet_to_met
    convert_to_isotropic = params.convert_to_isotropic
    reset_occupancies = params.reset_occupancies
    remove_ligands = params.remove_ligands
    reset_hetatm_flag = params.reset_hetatm_flag
  if (log is None):
    log = null_out()
  make_sub_header("Processing input model", out=log)
  remarks = None
  if (file_name is not None):
    print("Reading model from %s" % file_name, file=log)
    assert ([pdb_hierarchy, xray_structure] == [None, None])
    import iotbx.pdb
    pdb_in = iotbx.pdb.input(file_name)
    remarks = pdb_in.remark_section()
    pdb_hierarchy = pdb_in.construct_hierarchy()
    xray_structure = pdb_in.xray_structure_simple()
  else :
    # XXX work with copies, not the original structure
    pdb_hierarchy = pdb_hierarchy.deep_copy()
    xray_structure = xray_structure.deep_copy_scatterers()
  pdb_hierarchy.atoms().reset_i_seq()
  if (len(pdb_hierarchy.models()) > 1):
    raise Sorry("Multiple models not supported.")
  if (remove_hydrogens):
    sele = ~(xray_structure.hd_selection())
    n_hd = sele.count(False)
    if (n_hd > 0):
      pdb_hierarchy = pdb_hierarchy.select(sele)
      xray_structure = xray_structure.select(sele)
      print("  removed %d hydrogens" % n_hd, file=log)
      pdb_hierarchy.atoms().reset_i_seq()
  if (remove_waters):
    sele = pdb_hierarchy.atom_selection_cache().selection("not (resname HOH)")
    n_wat = sele.count(False)
    if (n_wat > 0):
      pdb_hierarchy = pdb_hierarchy.select(sele)
      xray_structure = xray_structure.select(sele)
      print("  removed %d waters" % n_wat, file=log)
      pdb_hierarchy.atoms().reset_i_seq()
  if (remove_alt_confs):
    n_atoms_start = xray_structure.scatterers().size()
    pdb_hierarchy.remove_alt_confs(always_keep_one_conformer=False)
    i_seqs = pdb_hierarchy.atoms().extract_i_seq()
    n_atoms_end = i_seqs.size()
    if (n_atoms_end != n_atoms_start):
      print("  removed %d atoms in alternate conformations" % \
        (n_atoms_end - n_atoms_start), file=log)
    xray_structure = xray_structure.select(i_seqs)
    pdb_hierarchy.atoms().reset_i_seq()
  if (convert_semet_to_met):
    # XXX need to start from a copy here because the atom-parent relationship
    # seems to be messed up otherwise.  this is probably a bug.
    pdb_hierarchy = pdb_hierarchy.deep_copy()
    pdb_hierarchy.convert_semet_to_met()
  if (convert_to_isotropic):
    xray_structure.convert_to_isotropic()
    pdb_hierarchy.adopt_xray_structure(xray_structure)
    print("  converted all atoms to isotropic B-factors", file=log)
  if (reset_occupancies):
    assert (remove_alt_confs)
    xray_structure.adjust_occupancy(occ_max=1.0, occ_min=1.0)
    pdb_hierarchy.adopt_xray_structure(xray_structure)
    print("  reset occupancy to 1.0 for all atoms", file=log)
  if (reset_hetatm_flag):
    for atom in pdb_hierarchy.atoms():
      atom.hetero = False
  if (remove_ligands):
    pdb_hierarchy.atoms().reset_i_seq()
    model = pdb_hierarchy.only_model()
    for chain in model.chains():
      if (not chain.is_protein()) and (not chain.is_na()):
        print("  removing %d ligand atoms in chain '%s'" % \
          (len(chain.atoms()), chain.id), file=log)
        model.remove_chain(chain)
    i_seqs = pdb_hierarchy.atoms().extract_i_seq()
    xray_structure = xray_structure.select(i_seqs)
    pdb_hierarchy.atoms().reset_i_seq()
  assert xray_structure.scatterers().size() == pdb_hierarchy.atoms_size()
  if (output_file is not None):
    f = open(output_file, "w")
    if (add_remarks is not None):
      f.write("\n".join(add_remarks))
      f.write("\n")
    if (preserve_remarks) and (remarks is not None):
      f.write("\n".join(remarks))
      f.write("\n")
    symm = None
    if (preserve_symmetry):
      symm = xray_structure
    f.write(pdb_hierarchy.as_pdb_string(crystal_symmetry=symm))
    f.close()
    print("  wrote model to %s" % output_file, file=log)
  return pdb_hierarchy, xray_structure
