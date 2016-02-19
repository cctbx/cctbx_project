from __future__ import division
import os, sys, random
from iotbx.phil import parse
from libtbx import easy_run
from xfel.cftbx.detector.cspad_cbf_tbx import write_cspad_cbf, map_detector_to_basis_dict
from dials.phil import ExperimentListConverters

phil_scope = parse("""
  tag = cspad
    .type = str
    .help = Name of this refinement run. Output filenames will use this tag.
  refine_to_hierarchy_level = 3
    .type = int
    .help = maximum level to refine cspad to
  n_subset = None
    .type = int
    .help = Refine a random subset of the provided files
""", process_includes=True)

def run(args):
  print "Parsing input..."
  user_phil = []
  paths = []
  refine_phil_file = None
  for arg in args:
    try:
      if os.path.splitext(arg)[1] == ".phil":
        refine_phil_file = arg
        continue
    except Exception, e:
      pass
    if os.path.isdir(arg):
      paths.append(arg)
    else:
      try:
        user_phil.append(parse(arg))
      except Exception, e:
        print "Unrecognized argument:", arg
        return

  params = phil_scope.fetch(sources=user_phil).extract()

  print "Gathering file names..."
  all_exp = []
  all_ref = []

  for path in paths:
    for filename in os.listdir(path):
      if "indexed" in filename:
        exp_path = os.path.join(path, filename.rstrip("_indexed.pickle") + "_refined_experiments.json")
        if not os.path.exists(exp_path): continue
        all_exp.append(exp_path)
        all_ref.append(os.path.join(path, filename))

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

  f = open("%s_combine.phil"%params.tag, 'w')
  for exp_path, ref_path in zip(all_exp, all_ref):
    f.write("input {\n")
    f.write("  experiments = %s\n"%exp_path)
    f.write("  reflections = %s\n"%ref_path)
    f.write("}\n")
  f.close()

  print "Combining experiments..."
  command = "dials.combine_experiments reference_from_experiment.average_detector=True reference_from_experiment.average_hierarchy_level=0 output.experiments_filename=%s_combined_experiments.json output.reflections_filename=%s_combined_reflections.pickle %s_combine.phil"%(params.tag, params.tag, params.tag)
  print command
  result = easy_run.fully_buffered(command=command).raise_if_errors()
  result.show_stdout()

  for i in xrange(params.refine_to_hierarchy_level+1):
    print "Refining at hierarchy level", i
    if i == 0:
      command = "dials.refine %s %s_combined_experiments.json %s_combined_reflections.pickle"%(refine_phil_file, params.tag, params.tag)
    else:
      command = "dials.refine %s %s_refined_experiments_level%d.json %s_refined_reflections_level%d.pickle"%(refine_phil_file, params.tag, i-1, params.tag, i-1)
    command += " output.experiments=%s_refined_experiments_level%d.json output.reflections=%s_refined_reflections_level%d.pickle"%( \
      params.tag, i, params.tag, i)
    command += " refinement.parameterisation.detector.hierarchy_level=%d"%i
    print command
    result = easy_run.fully_buffered(command=command).raise_if_errors()
    result.show_stdout()

  print "Creating files to deploy to psana calibration directory..."
  if params.refine_to_hierarchy_level > 2:
    deploy_level = 2
  else:
    deploy_level = params.refine_to_hierarchy_level

  converter = ExperimentListConverters(check_format=False)
  detector = converter.from_string("%s_refined_experiments_level%d.json"%(params.tag, deploy_level)).data[0].detector

  metro = map_detector_to_basis_dict(detector)
  write_cspad_cbf(None, metro, 'cbf', None, '%s_refined_detector_level%d.def'%(params.tag, deploy_level), None, abs(detector.hierarchy().get_distance()), header_only=True)

  command = "cxi.cbfheader2slaccalib cbf_header=%s_refined_detector_level%d.def out_metrology_file=0-end.data.%s"%(params.tag, deploy_level, params.tag)
  result = easy_run.fully_buffered(command=command).raise_if_errors()
  result.show_stdout()

  print "Done. Soft link 0-end.data.%s to 0-end.data in the geometry folder of your calibration folder for your experiment to deploy this metrology."

if __name__ == "__main__":
  run(sys.argv[1:])
