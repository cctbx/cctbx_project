# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cspad.cbf_metrology
#
from __future__ import division
import os, sys, random
from iotbx.phil import parse
from libtbx import easy_run
from libtbx.utils import Sorry

phil_scope = parse("""
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
  refine_distance = False
    .type = bool
    .help = If true, allow root hierarchy level to refine in Z. Otherwise fix this \
            axis. Regardless, higher hierarchy levels will refine in Z.
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
  n_subset_method = *random n_refl
    .type = choice
    .help = Algorithm to be used for choosing the n_subset images/experiments for \
            refinement.  n_refl chooses the set with the largest numbers of reflections \
            listed in the pickle files, thus giving maximal coverage of the detector tiles \
            with the fewest refineable parameters.
  doit = False
    .type = bool
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

refine_scope = parse("""
  include scope dials.command_line.refine.phil_scope
""", process_includes=True)

def run(args):
  print "Parsing input..."
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
      except Exception, e:
        raise Sorry("Unrecognized file %s"%arg)
    if os.path.isdir(arg):
      paths.append(arg)
    else:
      try:
        user_phil.append(parse(arg))
      except Exception, e:
        raise Sorry("Unrecognized argument: %s"%arg)

  params = phil_scope.fetch(sources=user_phil).extract()

  merged_scope = refine_scope.fetch(refine_defaults_scope)
  if refine_phil_file is not None:
    merged_scope = merged_scope.fetch(parse(file_name = refine_phil_file))

  print "Gathering file names..."
  all_exp = []
  all_ref = []

  if params.data_phil is None:
    for path in paths:
      for filename in os.listdir(path):
        if params.reflections in filename:
          exp_path = os.path.join(path, filename.rstrip("_%s.pickle"%params.reflections) + "_refined_experiments.json")
          if not os.path.exists(exp_path): continue
          all_exp.append(exp_path)
          all_ref.append(os.path.join(path, filename))

    if params.split_dataset:
      import re
      even_exp = []
      odd_exp = []
      even_ref = []
      odd_ref = []
      for exp, ref in zip(all_exp, all_ref):
        if int(re.findall(r'\d+', exp)[-1][-1]) % 2 == 0:
          even_exp.append(exp)
          even_ref.append(ref)
        else:
          odd_exp.append(exp)
          odd_ref.append(ref)

      base_tag = params.tag
      base_n_subset = params.n_subset

      params.n_subset = base_n_subset // 2
      odd_exp, odd_ref = generate_exp_list(params, odd_exp, odd_ref)
      even_exp, even_ref = generate_exp_list(params, even_exp, even_ref)

      params.tag = base_tag + "_1"
      odd_combine_phil = write_combine_phil(params, odd_exp, odd_ref)
      params.tag = base_tag + "_2"
      even_combine_phil = write_combine_phil(params, even_exp, even_ref)

      params.tag = base_tag
      params.n_subset = base_n_subset
      full_combine_phil = write_combine_phil(params, odd_exp+even_exp, odd_ref+even_ref)

      print "Refining full dataset using tag", params.tag
      refine(params, merged_scope, full_combine_phil)

      params.tag = base_tag + "_1"
      print "Refining odd numbered data using tag", params.tag
      refine(params, merged_scope, odd_combine_phil)

      params.tag = base_tag + "_2"
      print "Refining even numbered data using tag", params.tag
      refine(params, merged_scope, even_combine_phil)
    else:
      all_exp, all_ref = generate_exp_list(params, all_exp, all_ref)
      combine_phil = write_combine_phil(params, all_exp, all_ref)
      refine(params, merged_scope, combine_phil)
  else:
    assert len(paths) == 0
    assert params.n_subset is None
    assert not params.split_dataset
    refine(params, merged_scope, params.data_phil)

def generate_exp_list(params, all_exp, all_ref):
  if params.n_subset is not None:
    subset_all_exp = []
    subset_all_ref = []
    n_picked = 0
    if params.n_subset_method=="random":
      while n_picked < params.n_subset:
        idx = random.randint(0, len(all_exp)-1)
        subset_all_exp.append(all_exp.pop(idx))
        subset_all_ref.append(all_ref.pop(idx))
        n_picked += 1
    elif params.n_subset_method=="n_refl":
      from dials.array_family import flex
      import cPickle as pickle
      len_all_ref = flex.size_t(
        [ len(pickle.load(open(A,"rb"))) for A in all_ref ]
      )
      sort_order = flex.sort_permutation(len_all_ref,reverse=True)
      for idx in sort_order[:params.n_subset]:
        subset_all_exp.append(all_exp[idx])
        subset_all_ref.append(all_ref[idx])
      print "Selecting a subset of %d images with highest n_refl out of %d total."%(
        params.n_subset, len(len_all_ref))

    all_exp = subset_all_exp
    all_ref = subset_all_ref
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
  print "Combining experiments..."
  command = "dials.combine_experiments reference_from_experiment.average_detector=True reference_from_experiment.average_hierarchy_level=0 output.experiments_filename=%s_combined_experiments.json output.reflections_filename=%s_combined_reflections.pickle %s"%(params.tag, params.tag, combine_phil)
  print command
  result = easy_run.fully_buffered(command=command).raise_if_errors()
  result.show_stdout()
  for i in xrange(params.start_at_hierarchy_level, params.refine_to_hierarchy_level+1):
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
      print command
      result = easy_run.fully_buffered(command=command).raise_if_errors()
      result.show_stdout()

    print "Refining at hierarchy level", i
    refine_phil_file = "%s_refine_level%d.phil"%(params.tag, i)
    if i == 0:
      if params.refine_distance:
        diff_phil = "refinement.parameterisation.detector.fix_list=None\n" # allow full freedom to refine
      else:
        diff_phil = "refinement.parameterisation.detector.fix_list=Tau1\n" # fix detector rotz
    else:
      diff_phil = "refinement.parameterisation.detector.fix_list=None\n" # allow full freedom to refine

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

    print command
    result = easy_run.fully_buffered(command=command).raise_if_errors()
    result.show_stdout()

  print "Creating files to deploy to psana calibration directory..."
  if params.refine_to_hierarchy_level > 2:
    deploy_level = 2
  else:
    deploy_level = params.refine_to_hierarchy_level

  command = "cxi.experiment_json_to_cbf_def %s_refined_experiments_level%d.json output_def_file=%s_refined_detector_level%d.def"%(params.tag, deploy_level, params.tag, deploy_level)
  print command
  result = easy_run.fully_buffered(command=command).raise_if_errors()
  result.show_stdout()

  command = "cxi.cbfheader2slaccalib cbf_header=%s_refined_detector_level%d.def out_metrology_file=0-end.data.%s"%(params.tag, deploy_level, params.tag)
  print command
  result = easy_run.fully_buffered(command=command).raise_if_errors()
  result.show_stdout()

  print "Done. Soft link 0-end.data.%s to 0-end.data in the geometry folder of your calibration folder for your experiment to deploy this metrology."%params.tag

if __name__ == "__main__":
  run(sys.argv[1:])
