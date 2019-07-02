# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cspad.cbf_metrology
#
from __future__ import absolute_import, division, print_function
from six.moves import range
import os, sys, random
from iotbx.phil import parse
from libtbx import easy_run
from libtbx.utils import Sorry
import six
from six.moves import zip

phil_scope = parse("""
  method = *hierarchical expanding
    .type = choice
  reflections = reindexedstrong *indexed integrated
    .type = choice
    .help = Which subset of reflections
  tag = cspad
    .type = str
    .help = Name of this refinement run. Output filenames will use this tag.
  start_at_hierarchy_level = 0
    .type = int
    .help = Start refinement at this hierarchy level
  refine_to_hierarchy_level = 2
    .type = int
    .help = maximum level to refine cspad to
  refine_distance = True
    .type = bool
    .help = If true, allow root hierarchy level to refine in Z. Otherwise fix this \
            axis. Regardless, higher hierarchy levels will refine in Z.
  refine_energy = False
    .type = bool
    .help = If true, when refining level 0, also refine beam energy. Subsequent hierarchy \
            levels will fix the energy in place.
  flat_refinement = False
    .type = bool
    .help = If True, do not refine tilt (Tau2 and Tau3) when refining panel positions. Further, \
            don't refine distance at levels 1 or higher (respects refine_distance for level 0).
  flat_refinement_with_distance = False
    .type = bool
    .help = If True, and if using flat refinement, then use constraints to allow disance \
            to refine at levels 1 and higher.
  n_subset = None
    .type = int
    .help = Refine a random subset of the provided files
  split_dataset = False
    .type = bool
    .help = After refining the full set of images, if split_dataset is True, the \
            data will be split in two using odd and even file numbers and each half \
            will be refined independently. For each half, _1 or _2 is appended to \
            the tag. If used with n_subset, each half will have n_subset/2 images.
  data_phil = None
    .type = str
    .help = Optional phil file with all experiments and reflections for use during \
            refinement.  If not provided, the program will use whatever directories \
            were specified.
  rmsd_filter {
    enable = True
      .type = bool
      .help = If enabled, between each round of hierarchical refinement, filter \
              the images by positional RMSD
    iqr_multiplier = 1.5
      .type = float
      .help = Interquartile multiplier
  }
  n_subset_method = *random n_refl significance_filter
    .type = choice
    .help = Algorithm to be used for choosing the n_subset images/experiments for \
            refinement.  n_refl chooses the set with the largest numbers of reflections \
            listed in the reflection table files, thus giving maximal coverage of the detector tiles \
            with the fewest refineable parameters. Significance_filter chooses the subset of \
            images with maximum reflections above an I/sigI cutoff
  n_refl_panel_list = None
    .type = ints
    .help = If n_subset_method is n_refl, specify which panels to search on.
  panel_filter = None
    .type = ints
    .help = Specify a list of panels to include during refinement. Default (None) is to use \
            all panels.
  output_lcls_geometry = True
    .type = bool
    .help = If True, convert final refined geometry to LCLS format
""", process_includes=True)

refine_defaults_scope = parse("""
output.include_unused_reflections=False
refinement {
  refinery.engine = SparseLevMar
  parameterisation {
    beam.fix=all
    auto_reduction {
      action=remove
      min_nref_per_parameter=3
    }
  }
  reflections {
    outlier {
      algorithm=sauter_poon
      separate_panels=True
      separate_experiments=False
    }
  }
}
""")

def is_even(filename):
  import re
  return int(re.findall(r'\d+', filename)[-1][-1]) % 2 == 0

refine_scope = parse("""
  include scope dials.command_line.refine.phil_scope
""", process_includes=True)

def run(args):
  print("Parsing input...")
  if "-c" in args or "-h" in args or "--help" in args:
    phil_scope.show(attributes_level=2)
    return
  user_phil = []
  paths = []
  refine_phil_file = None
  for arg in args:
    if os.path.isfile(arg):
      try:
        if os.path.splitext(arg)[1] == ".phil":
          refine_phil_file = arg
          continue
      except Exception as e:
        raise Sorry("Unrecognized file %s"%arg)
    if os.path.isdir(arg):
      paths.append(arg)
    else:
      try:
        user_phil.append(parse(arg))
      except Exception as e:
        raise Sorry("Unrecognized argument: %s"%arg)

  params = phil_scope.fetch(sources=user_phil).extract()

  merged_scope = refine_scope.fetch(refine_defaults_scope)
  if refine_phil_file is not None:
    merged_scope = merged_scope.fetch(parse(file_name = refine_phil_file))

  print("Gathering file names...")
  all_exp = []
  all_ref = []

  if params.data_phil is None:
    for path in paths:
      exp, ref = find_files(path, params.reflections)
      all_exp.extend(exp)
      all_ref.extend(ref)

    if params.split_dataset:
      even_exp = []
      odd_exp = []
      even_ref = []
      odd_ref = []
      for exp, ref in zip(all_exp, all_ref):
        if is_even(exp):
          even_exp.append(exp)
          even_ref.append(ref)
        else:
          odd_exp.append(exp)
          odd_ref.append(ref)

      base_tag = params.tag
      base_n_subset = params.n_subset

      params.n_subset = base_n_subset // 2

      params.tag = base_tag + "_1"
      odd_combine_phil = write_combine_phil(params, odd_exp, odd_ref)
      params.tag = base_tag + "_2"
      even_combine_phil = write_combine_phil(params, even_exp, even_ref)

      params.tag = base_tag
      params.n_subset = base_n_subset
      full_combine_phil = write_combine_phil(params, odd_exp+even_exp, odd_ref+even_ref)

      print("Refining full dataset using tag", params.tag)
      refine(params, merged_scope, full_combine_phil)

      params.tag = base_tag + "_1"
      print("Refining odd numbered data using tag", params.tag)
      refine(params, merged_scope, odd_combine_phil)

      params.tag = base_tag + "_2"
      print("Refining even numbered data using tag", params.tag)
      refine(params, merged_scope, even_combine_phil)
    else:
      combine_phil = write_combine_phil(params, all_exp, all_ref)
      refine(params, merged_scope, combine_phil)
  else:
    assert len(paths) == 0
    assert params.n_subset is None
    print("Refining full dataset using tag", params.tag)
    refine(params, merged_scope, params.data_phil)
    if params.split_dataset:

      input_scope = parse("""
        input {
          experiments = None
            .type = str
            .multiple = True
            .help = "The experiment list file path"
          reflections = None
            .type = str
            .multiple = True
            .help = "The reflection table file path"
        }
        """)
      input_params = input_scope.fetch(parse(file_name = params.data_phil)).extract()
      even_exp = []
      odd_exp = []
      even_ref = []
      odd_ref = []
      for f in input_params.input.experiments:
        if is_even(f):
          even_exp.append(f)
        else:
          odd_exp.append(f)
      for f in input_params.input.reflections:
        if is_even(f):
          even_ref.append(f)
        else:
          odd_ref.append(f)
      base_tag = params.tag
      params.tag = base_tag + "_1"
      odd_combine_phil = write_combine_phil(params, odd_exp, odd_ref)
      params.tag = base_tag + "_2"
      even_combine_phil = write_combine_phil(params, even_exp, even_ref)

      params.tag = base_tag + "_1"
      print("Refining odd numbered data using tag", params.tag)
      refine(params, merged_scope, odd_combine_phil)

      params.tag = base_tag + "_2"
      print("Refining even numbered data using tag", params.tag)
      refine(params, merged_scope, even_combine_phil)

def find_files(path, reflections):
  all_exp = []
  all_ref = []
  for filename in os.listdir(path):
    if reflections in filename:
      extension = os.path.splitext(filename)[1]
      if extension not in ['.pickle', '.mpack']: continue
      exp_path = os.path.join(path, filename.rstrip("_%s%s"%(reflections, extension)) + "_refined_experiments.json")
      if not os.path.exists(exp_path):
        exp_path = os.path.join(path, filename.rstrip("_%s%s"%(reflections, extension)) + "_experiments.json")
        if not os.path.exists(exp_path): continue
      all_exp.append(exp_path)
      all_ref.append(os.path.join(path, filename))
  return all_exp, all_ref

def write_combine_phil(params, all_exp, all_ref):
  combine_phil = "%s_combine.phil"%params.tag
  f = open(combine_phil, 'w')
  for exp_path, ref_path in zip(all_exp, all_ref):
    f.write("input {\n")
    f.write("  experiments = %s\n"%exp_path)
    f.write("  reflections = %s\n"%ref_path)
    f.write("}\n")
  f.close()

  return combine_phil

def refine(params, merged_scope, combine_phil):
  print("Combining experiments...")
  command = "dials.combine_experiments reference_from_experiment.average_detector=True reference_from_experiment.average_hierarchy_level=0 output.experiments_filename=%s_combined_experiments.json output.reflections_filename=%s_combined_reflections.pickle %s"%(params.tag, params.tag, combine_phil)
  if params.n_subset is not None:
    command += " n_subset=%d n_subset_method=%s"%(params.n_subset, params.n_subset_method)
    if params.n_refl_panel_list is not None:
      command += " n_refl_panel_list=%s"%(",".join(["%d"%p for p in params.n_refl_panel_list]))

  if params.refine_energy:
    command += " reference_from_experiment.beam=0"
  print(command)
  result = easy_run.fully_buffered(command=command).raise_if_errors()
  result.show_stdout()

  if params.method == 'hierarchical':
    refine_hierarchical(params, merged_scope, combine_phil)
  elif params.method == 'expanding':
    refine_expanding(params, merged_scope, combine_phil)

def refine_hierarchical(params, merged_scope, combine_phil):
  if params.panel_filter is not None:
    from libtbx import easy_pickle
    print("Filtering out all reflections except those on panels %s"%(", ".join(["%d"%p for p in params.panel_filter])))
    combined_path = "%s_combined_reflections.pickle"%params.tag
    data = easy_pickle.load(combined_path)
    sel = None
    for panel_id in params.panel_filter:
      if sel is None:
        sel = data['panel'] == panel_id
      else:
        sel |= data['panel'] == panel_id
    print("Retaining", len(data.select(sel)), "out of", len(data), "reflections")
    easy_pickle.dump(combined_path, data.select(sel))

  for i in range(params.start_at_hierarchy_level, params.refine_to_hierarchy_level+1):
    if params.rmsd_filter.enable:
      input_name = "filtered"
    else:
      if i == params.start_at_hierarchy_level:
        input_name = "combined"
      else:
        input_name = "refined"

    if params.rmsd_filter.enable:
      command = "cctbx.xfel.filter_experiments_by_rmsd %s %s output.filtered_experiments=%s output.filtered_reflections=%s"
      if i == params.start_at_hierarchy_level:
        command = command%("%s_combined_experiments.json"%params.tag, "%s_combined_reflections.pickle"%params.tag,
                           "%s_filtered_experiments.json"%params.tag, "%s_filtered_reflections.pickle"%params.tag)
      else:
        command = command%("%s_refined_experiments_level%d.json"%(params.tag, i-1), "%s_refined_reflections_level%d.pickle"%(params.tag, i-1),
                           "%s_filtered_experiments_level%d.json"%(params.tag, i-1), "%s_filtered_reflections_level%d.pickle"%(params.tag, i-1))
      command += " iqr_multiplier=%f"%params.rmsd_filter.iqr_multiplier
      print(command)
      result = easy_run.fully_buffered(command=command).raise_if_errors()
      result.show_stdout()

    print("Refining at hierarchy level", i)
    refine_phil_file = "%s_refine_level%d.phil"%(params.tag, i)
    if i == 0:
      fix_list = ['Tau1'] # fix detector rotz
      if not params.refine_distance:
        fix_list.append('Dist')
      if params.flat_refinement:
        fix_list.extend(['Tau2','Tau3'])

      diff_phil = "refinement.parameterisation.detector.fix_list=%s\n"%",".join(fix_list)
      if params.refine_energy:
        diff_phil += " refinement.parameterisation.beam.fix=in_spindle_plane+out_spindle_plane\n" # allow energy to refine
    else:
      # Note, always need to fix something, so pick a panel group and fix its Tau1 (rotation around Z) always
      if params.flat_refinement and params.flat_refinement_with_distance:
        diff_phil = "refinement.parameterisation.detector.fix_list=Group1Tau1,Tau2,Tau3\n" # refine distance, rotz and xy translation
        diff_phil += "refinement.parameterisation.detector.constraints.parameter=Dist\n" # constrain distance to be refined identically for all panels at this hierarchy level
      elif params.flat_refinement:
        diff_phil = "refinement.parameterisation.detector.fix_list=Dist,Group1Tau1,Tau2,Tau3\n" # refine only rotz and xy translation
      else:
        diff_phil = "refinement.parameterisation.detector.fix_list=Group1Tau1\n" # refine almost everything

    if i == params.start_at_hierarchy_level:
      command = "dials.refine %s %s_%s_experiments.json %s_%s_reflections.pickle"%(refine_phil_file, params.tag, input_name, params.tag, input_name)
    else:
      command = "dials.refine %s %s_%s_experiments_level%d.json %s_%s_reflections_level%d.pickle"%(refine_phil_file, params.tag, input_name, i-1, params.tag, input_name, i-1)

    diff_phil += "refinement.parameterisation.detector.hierarchy_level=%d\n"%i

    command += " output.experiments=%s_refined_experiments_level%d.json output.reflections=%s_refined_reflections_level%d.pickle"%( \
      params.tag, i, params.tag, i)

    scope = merged_scope.fetch(parse(diff_phil))
    f = open(refine_phil_file, 'w')
    f.write(refine_scope.fetch_diff(scope).as_str())
    f.close()

    print(command)
    result = easy_run.fully_buffered(command=command).raise_if_errors()
    result.show_stdout()

  output_geometry(params)

def refine_expanding(params, merged_scope, combine_phil):
  assert params.start_at_hierarchy_level == 0
  if params.rmsd_filter.enable:
    input_name = "filtered"
    command = "cctbx.xfel.filter_experiments_by_rmsd %s %s output.filtered_experiments=%s output.filtered_reflections=%s"
    command = command%("%s_combined_experiments.json"%params.tag, "%s_combined_reflections.pickle"%params.tag,
                       "%s_filtered_experiments.json"%params.tag, "%s_filtered_reflections.pickle"%params.tag)
    command += " iqr_multiplier=%f"%params.rmsd_filter.iqr_multiplier
    print(command)
    result = easy_run.fully_buffered(command=command).raise_if_errors()
    result.show_stdout()
  else:
    input_name = "combined"
  # --------------------------
  if params.panel_filter is not None:
    from libtbx import easy_pickle
    print("Filtering out all reflections except those on panels %s"%(", ".join(["%d"%p for p in params.panel_filter])))
    combined_path = "%s_combined_reflections.pickle"%params.tag
    data = easy_pickle.load(combined_path)
    sel = None
    for panel_id in params.panel_filter:
      if sel is None:
        sel = data['panel'] == panel_id
      else:
        sel |= data['panel'] == panel_id
    print("Retaining", len(data.select(sel)), "out of", len(data), "reflections")
    easy_pickle.dump(combined_path, data.select(sel))
  # ----------------------------------
  # this is the order to refine the CSPAD in
  steps = {}
  steps[0] = [2, 3]
  steps[1] = steps[0] + [0, 1]
  steps[2] = steps[1] + [14, 15]
  steps[3] = steps[2] + [6, 7]
  steps[4] = steps[3] + [4, 5]
  steps[5] = steps[4] + [12, 13]
  steps[6] = steps[5] + [8, 9]
  steps[7] = steps[6] + [10, 11]

  for s, panels in six.iteritems(steps):
    rest = []
    for p in panels:
      rest.append(p+16)
      rest.append(p+32)
      rest.append(p+48)
    panels.extend(rest)

  levels = {0: (0,1)} # levels 0 and 1
  for i in range(7):
    levels[i+1] = (2,) # level 2

  previous_step_and_level = None
  for j in range(8):
    from libtbx import easy_pickle
    print("Filtering out all reflections except those on panels %s"%(", ".join(["%d"%p for p in steps[j]])))
    combined_path = "%s_%s_reflections.pickle"%(params.tag, input_name)
    output_path = "%s_reflections_step%d.pickle"%(params.tag, j)
    data = easy_pickle.load(combined_path)
    sel = None
    for panel_id in steps[j]:
      if sel is None:
        sel = data['panel'] == panel_id
      else:
        sel |= data['panel'] == panel_id
    print("Retaining", len(data.select(sel)), "out of", len(data), "reflections")
    easy_pickle.dump(output_path, data.select(sel))

    for i in levels[j]:
      print("Step", j , "refining at hierarchy level", i)
      refine_phil_file = "%s_refine_step%d_level%d.phil"%(params.tag, j, i)
      if i == 0:
        if params.refine_distance:
          diff_phil = "refinement.parameterisation.detector.fix_list=Tau1" # fix detector rotz
        else:
          diff_phil = "refinement.parameterisation.detector.fix_list=Dist,Tau1" # fix detector rotz, distance
        if params.flat_refinement:
          diff_phil += ",Tau2,Tau3" # Also fix x and y rotations
        diff_phil += "\n"
        if params.refine_energy:
          diff_phil += "refinement.parameterisation.beam.fix=in_spindle_plane+out_spindle_plane\n" # allow energy to refine
      else:
        # Note, always need to fix something, so pick a panel group and fix its Tau1 (rotation around Z) always
        if params.flat_refinement and params.flat_refinement_with_distance:
          diff_phil = "refinement.parameterisation.detector.fix_list=Group1Tau1,Tau2,Tau3\n" # refine distance, rotz and xy translation
          diff_phil += "refinement.parameterisation.detector.constraints.parameter=Dist\n" # constrain distance to be refined identically for all panels at this hierarchy level
        elif params.flat_refinement:
          diff_phil = "refinement.parameterisation.detector.fix_list=Dist,Group1Tau1,Tau2,Tau3\n" # refine only rotz and xy translation
        else:
          diff_phil = "refinement.parameterisation.detector.fix_list=Group1Tau1\n" # refine almost everything

      if previous_step_and_level is None:
        command = "dials.refine %s %s_%s_experiments.json %s_reflections_step%d.pickle"%( \
          refine_phil_file, params.tag, input_name, params.tag, j)
      else:
        p_step, p_level = previous_step_and_level
        if p_step == j:
          command = "dials.refine %s %s_refined_experiments_step%d_level%d.json %s_refined_reflections_step%d_level%d.pickle"%( \
            refine_phil_file, params.tag, p_step, p_level, params.tag, p_step, p_level)
        else:
          command = "dials.refine %s %s_refined_experiments_step%d_level%d.json %s_reflections_step%d.pickle"%( \
            refine_phil_file, params.tag, p_step, p_level, params.tag, j)


      diff_phil += "refinement.parameterisation.detector.hierarchy_level=%d\n"%i

      output_experiments = "%s_refined_experiments_step%d_level%d.json"%(params.tag, j, i)
      command += " output.experiments=%s output.reflections=%s_refined_reflections_step%d_level%d.pickle"%( \
        output_experiments, params.tag, j, i)

      scope = merged_scope.fetch(parse(diff_phil))
      f = open(refine_phil_file, 'w')
      f.write(refine_scope.fetch_diff(scope).as_str())
      f.close()

      print(command)
      result = easy_run.fully_buffered(command=command).raise_if_errors()
      result.show_stdout()

      # In expanding mode, if using flat refinement with distance, after having refined this step as a block, unrefined
      # panels will have been left behind.  Read back the new metrology, compute the shift applied to the panels refined
      # in this step,and apply that shift to the unrefined panels in this step
      if params.flat_refinement and params.flat_refinement_with_distance and i > 0:
        from dxtbx.model.experiment_list import ExperimentListFactory, ExperimentListDumper
        from xfel.command_line.cspad_detector_congruence import iterate_detector_at_level, iterate_panels
        from scitbx.array_family import flex
        from scitbx.matrix import col
        from libtbx.test_utils import approx_equal
        experiments = ExperimentListFactory.from_json_file(output_experiments, check_format=False)
        assert len(experiments.detectors()) == 1
        detector = experiments.detectors()[0]
        # Displacements: deltas along the vector normal to the detector
        displacements = flex.double()
        # Iterate through the panel groups at this level
        for panel_group in iterate_detector_at_level(detector.hierarchy(), 0, i):
          # Were there panels refined in this step in this panel group?
          if params.panel_filter:
            test = [list(detector).index(panel) in steps[j] for panel in iterate_panels(panel_group) if list(detector).index(panel) in params.panel_filter]
          else:
            test = [list(detector).index(panel) in steps[j] for panel in iterate_panels(panel_group)]
          if not any(test): continue
          # Compute the translation along the normal of this panel group.  This is defined as distance in dials.refine
          displacements.append(col(panel_group.get_local_fast_axis()).cross(col(panel_group.get_local_slow_axis())).dot(col(panel_group.get_local_origin())))

        # Even though the panels are constrained to move the same amount, there is a bit a variation.
        stats = flex.mean_and_variance(displacements)
        displacement = stats.mean()
        print("Average displacement along normals: %f +/- %f"%(stats.mean(), stats.unweighted_sample_standard_deviation()))

        # Verify the variation isn't significant
        for k in range(1, len(displacements)):
          assert approx_equal(displacements[0], displacements[k])
        # If all of the panel groups in this level moved, no need to do anything.
        if len(displacements) != len(list(iterate_detector_at_level(detector.hierarchy(), 0, i))):
          for panel_group in iterate_detector_at_level(detector.hierarchy(), 0, i):
            if params.panel_filter:
              test = [list(detector).index(panel) in steps[j] and list(detector).index(panel) in params.panel_filter for panel in iterate_panels(panel_group)]
            else:
              test = [list(detector).index(panel) in steps[j] for panel in iterate_panels(panel_group)]
            # If any of the panels in this panel group moved, no need to do anything
            if any(test): continue

            # None of the panels in this panel group moved in this step, so need to apply displacement from other panel
            # groups at this level
            fast = col(panel_group.get_local_fast_axis())
            slow = col(panel_group.get_local_slow_axis())
            ori = col(panel_group.get_local_origin())
            normal = fast.cross(slow)
            panel_group.set_local_frame(fast, slow, (ori.dot(fast)*fast) + (ori.dot(slow)*slow) + (normal*displacement))

        # Check the new displacements. Should be the same across all panels.
        displacements = []
        for panel_group in iterate_detector_at_level(detector.hierarchy(), 0, i):
          displacements.append(col(panel_group.get_local_fast_axis()).cross(col(panel_group.get_local_slow_axis())).dot(col(panel_group.get_local_origin())))

        for k in range(1, len(displacements)):
          assert approx_equal(displacements[0], displacements[k])

        dump = ExperimentListDumper(experiments)
        dump.as_json(output_experiments)

      previous_step_and_level = j,i

  output_geometry(params)

def output_geometry(params):
  print("Creating files to deploy to psana calibration directory...")
  if params.refine_to_hierarchy_level > 2:
    deploy_level = 2
  else:
    deploy_level = params.refine_to_hierarchy_level

  if params.method == 'hierarchical':
    command = "cxi.experiment_json_to_cbf_def %s_refined_experiments_level%d.json output_def_file=%s_refined_detector_level%d.def"%(params.tag, deploy_level, params.tag, deploy_level)
  elif params.method == 'expanding':
    command = "cxi.experiment_json_to_cbf_def %s_refined_experiments_step7_level%d.json output_def_file=%s_refined_detector_level%d.def"%(params.tag, deploy_level, params.tag, deploy_level)
  print(command)
  result = easy_run.fully_buffered(command=command).raise_if_errors()
  result.show_stdout()

  if params.output_lcls_geometry:
    command = "cxi.cbfheader2slaccalib cbf_header=%s_refined_detector_level%d.def out_metrology_file=0-end.data.%s"%(params.tag, deploy_level, params.tag)
    print(command)
    result = easy_run.fully_buffered(command=command)
    errmsg = "\n".join(result.stderr_lines)
    if "ImportError" in errmsg and "PSCalib.GeometryAccess" in errmsg:
      print("Not converting to LCLS geometry as PSDM is not available")
      print("Done.")
    else:
      result.raise_if_errors()
      result.show_stdout()
      print("Done. Soft link 0-end.data.%s to 0-end.data in the geometry folder of your calibration folder for your experiment to deploy this metrology."%params.tag)

if __name__ == "__main__":
  run(sys.argv[1:])
