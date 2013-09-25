
"""
This module provides a way to identify a suitably isomorphous initial model
for directly refining against experimental data.  This is commonly done when
solving a series of ligand-bound structures, partly to save time, and partly
to avoid changing the frame of reference when running molecular replacement.
"""

from __future__ import division
from libtbx.str_utils import make_sub_header
from libtbx import slots_getstate_setstate
from libtbx.utils import null_out
from libtbx import easy_mp
import libtbx.phil
import sys

master_phil = libtbx.phil.parse("""
d_min = None
  .type = float
rigid_body_refine = False
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

class evaluate_model (slots_getstate_setstate) :
  """
  Create an fmodel object (including bulk solvent correction) and calculate
  R-factors, with or without optional rigid-body refinement.
  """
  __slots__ = [ "r_work", "r_free", "r_work_start", "r_free_start",
    "xray_structure" ]
  def __init__ (self,
      xray_structure,
      pdb_hierarchy,
      f_obs,
      r_free_flags,
      rigid_body_refine=False,
      optimize_b_factors=False,
      skip_twin_detection=False,
      scattering_table="n_gaussian") :
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
      scattering_table=scattering_table,
      update_f_part1_for="refinement")
    self.r_work_start = fmodel.r_work()
    self.r_free_start = fmodel.r_free()
    if (not rigid_body_refine) :
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

  def show (self, out=sys.stdout, prefix="") :
    print >> out, prefix + "r_work = %6.4f" % self.r_work
    print >> out, prefix + "r_free = %6.4f" % self.r_free

def ucf (unit_cell) :
  return "%g %g %g %g %g %g" % unit_cell.parameters()

class select_model (object) :
  def __init__ (self,
      model_file_names,
      f_obs,
      r_free_flags,
      params=None,
      skip_twin_detection=False,
      nproc=1,
      log=sys.stdout) :
    if (params is None) :
      params = master_phil.extract()
    self.model_file_names = model_file_names
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
    self.best_model_file = None
    from mmtbx.pdb_symmetry import rms_difference
    from iotbx import file_reader
    data_symmetry = f_obs.crystal_symmetry()
    data_space_group = data_symmetry.space_group()
    data_point_group = data_space_group.build_derived_point_group()
    data_unit_cell = data_symmetry.unit_cell()
    data_cell_edges = data_unit_cell.parameters()[0:3]
    data_cell_angles = data_unit_cell.parameters()[3:6]
    make_sub_header("Evaluating models", out=log)
    print >> log, "Experimental data:"
    print >> log, "  space group:  %s" % data_space_group.info()
    print >> log, "  unit cell:    %s" % ucf(data_unit_cell)
    pdb_hierarchies = []
    xray_structures = []
    for k, file_name in enumerate(model_file_names) :
      model_in = file_reader.any_file(file_name,
        force_type="pdb",
        raise_sorry_if_errors=True)
      pdb_hierarchy = model_in.file_object.construct_hierarchy()
      pdb_hierarchy.atoms().reset_i_seq()
      pdb_hierarchies.append(pdb_hierarchy)
      xray_structure = model_in.file_object.xray_structure_simple()
      model_symmetry = xray_structure.crystal_symmetry()
      self.model_symmetries.append(model_symmetry)
      if (model_symmetry is None) :
        print >> log, "Model %d is missing symmetry records:" % (k+1)
        print >> log, "  file:  %s" % file_name
        xray_structures.append(None)
        continue
      model_unit_cell = model_symmetry.unit_cell()
      model_space_group = model_symmetry.space_group()
      is_compatible_sg = False
      if (model_space_group == data_space_group) :
        is_compatible_sg = True
      else :
        model_point_group = model_space_group.build_derived_point_group()
        if (data_point_group == model_point_group) :
          is_compatible_sg = True
      if (not is_compatible_sg) :
        print >> log, "Model %d has incompatible space group:" % (k+1)
        print >> log, "  file:  %s" % file_name
        print >> log, "  space group: %s" % model_space_group.info()
        xray_structures.append(None)
        continue
      is_similar_cell = False
      if (model_unit_cell.is_similar_to(data_unit_cell)) :
        is_similar_cell = True
      else :
        model_cell_edges = model_unit_cell.parameters()[0:3]
        model_cell_angles = model_unit_cell.parameters()[3:6]
        cell_edge_rmsd = rms_difference(model_cell_edges, data_cell_edges)
        cell_angle_rmsd = rms_difference(model_cell_angles, data_cell_angles)
        if ((cell_edge_rmsd <= params.max_cell_edge_rmsd) and
            (cell_angle_rmsd <= params.max_cell_angle_rmsd)) :
          is_similar_cell = True
      if (not is_similar_cell) :
        print >> log, "Model %d has incompatible space group:" % (k+1)
        print >> log, "  file:  %s" % file_name
        print >> log, "  model: %s" % ucf(model_unit_cell.parameters())
        xray_structures.append(None)
        continue
      else :
        xray_structures.append(xray_structure)
    if (xray_structures.count(None) != len(xray_structures)) :
      print >> log, ""
      print >> log, "Calculating R-factors - will use %s processors." % nproc
      evaluations = easy_mp.parallel_map(
        func=self.evaluate_model,
        iterable=zip(xray_structures, pdb_hierarchies),
        processes=nproc)
      passed = []
      for k, result in enumerate(evaluations) :
        if (result is not None) :
          if (result.r_free <= params.max_r_free) :
            passed.append((k, result))
      if (len(passed) > 0) :
        passed.sort(lambda a,b: cmp(a[1].r_free, b[1].r_free))
        i_result, result = passed[0]
        self.evaluations = passed
        self.best_xray_structure = result.xray_structure
        self.best_pdb_hierarchy = pdb_hierarchies[i_result]
        self.best_result = result
        self.best_model_file = self.model_file_names[i_result]
    self.show(out=log, verbose=True)

  def show (self, out=sys.stdout, verbose=False) :
    if (self.best_result is None) :
      print >> out, "No models accepted - will need to run MR."
    else :
      print >> out, ""
      print >> out, "Best starting model:"
      print >> out, "  file: %s" % self.best_model_file
      self.best_result.show(out=out, prefix="    ")
      print >> out, ""
      if (verbose) and (len(self.evaluations) > 1) :
        print >> out, "Other suitable models:"
        for i_other, other in self.evaluations[1:] :
          print >> out, "  file: %s" % self.model_file_names[i_other]
          other.show(out=out, prefix="    ")
        print >> out, ""

  def success (self) :
    return self.best_pdb_hierarchy is not None

  def get_best_model (self, update_structure=True) :
    if (self.best_result is None) :
      return None
    xray_structure = self.best_xray_structure
    pdb_hierarchy = self.best_pdb_hierarchy
    if (update_structure) :
      pdb_hierarchy.adopt_xray_structure(xray_structure)
    return xray_structure, pdb_hierarchy

  def save_best_model (self, file_name="best_model.pdb") :
    assert self.success()
    xray_structure, pdb_hierarchy = self.get_best_model()
    f = open(file_name, "w")
    f.write("REMARK original PDB file:\n")
    f.write("REMARK   %s\n" % self.best_model_file)
    f.write(pdb_hierarchy.as_pdb_string(crystal_symmetry=xray_structure))
    f.close()

  def save_updated_data (self, file_name="best_model_data.mtz") :
    assert self.success()
    xray_structure, pdb_hierarchy = self.get_best_model()
    f_obs = self.f_obs.customized_copy(
      crystal_symmetry=xray_structure).eliminate_sys_absent()
    r_free_flags = self.r_free_flags.customized_copy(
      crystal_symmetry=xray_structure).eliminate_sys_absent()
    mtz_data = f_obs.as_mtz_dataset(column_root_label="F")
    mtz_data.add_miller_array(r_free_flags, column_root_label="FreeR_flag")
    mtz_data.mtz_object().write(file_name)

  def evaluate_model (self, args) :
    xray_structure, pdb_hierarchy = args
    if (xray_structure is None) :
      return None
    return evaluate_model(
      xray_structure=xray_structure,
      pdb_hierarchy=pdb_hierarchy,
      f_obs=self.f_obs,
      r_free_flags=self.r_free_flags,
      rigid_body_refine=self.params.rigid_body_refine,
      skip_twin_detection=self.skip_twin_detection)
