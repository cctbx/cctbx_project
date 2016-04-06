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
  refine_to_hierarchy_level = 3
    .type = int
    .help = maximum level to refine cspad to
  n_subset = None
    .type = int
    .help = Refine a random subset of the provided files
  split_dataset = False
    .type = bool
    .help = Whether to split the data in two using odd and even file numbers. Each \
            half is refined seperately and _1 or _2 is appended to the tag. If used \
            with n_subset, each half will have n_subset images.
  data_phil = None
    .type = str
    .help = Optional phil file with all experiments and reflections for use during \
            refinement.  If not provided, the program will use whatever directories \
            were specified.
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

      params.tag = base_tag + "_1"
      print "Refining odd numbered data using tag", params.tag
      combine_phil = generate_combine_phil(params, odd_exp, odd_ref)
      refine(params, merged_scope, combine_phil)

      params.tag = base_tag + "_2"
      print "Refining even numbered data using tag", params.tag
      combine_phil = generate_combine_phil(params, even_exp, even_ref)
      refine(params, merged_scope, combine_phil)
    else:
      combine_phil = generate_combine_phil(params, all_exp, all_ref)
      refine(params, merged_scope, combine_phil)
  else:
    assert len(paths) == 0
    assert params.n_subset is None
    assert not params.split_dataset
    refine(params, merged_scope, params.data_phil)

def generate_combine_phil(params, all_exp, all_ref):
  if params.n_subset is not None:
    subset_all_exp = []
    subset_all_ref = []
    n_picked = 0

    while n_picked < params.n_subset:
      idx = random.randint(0, len(all_exp)-1)
      subset_all_exp.append(all_exp.pop(idx))
      subset_all_ref.append(all_ref.pop(idx))
      n_picked += 1

    all_exp = subset_all_exp
    all_ref = subset_all_ref

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

  for i in xrange(params.refine_to_hierarchy_level+1):
    print "Refining at hierarchy level", i
    refine_phil_file = "%s_refine_level%d.phil"%(params.tag, i)
    if i == 0:
      diff_phil = "refinement.parameterisation.detector.fix_list=Tau1\n" # fix detector rotz
      command = "dials.refine %s %s_combined_experiments.json %s_combined_reflections.pickle"%(refine_phil_file, params.tag, params.tag)
    else:
      diff_phil = "refinement.parameterisation.detector.fix_list=None\n" # allow full freedom to refine
      command = "dials.refine %s %s_refined_experiments_level%d.json %s_refined_reflections_level%d.pickle"%(refine_phil_file, params.tag, i-1, params.tag, i-1)
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
