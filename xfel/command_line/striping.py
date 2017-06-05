from __future__ import division
# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.stripe_experiment
#
# Given an LCLS experiment results directory and a trial, group results by
# run group and then distrbute each run group's results into subgroups and run
# dials.combine_experiments (optionally with clustering and selecting clusters).
#
from libtbx.phil import parse
from libtbx.utils import Sorry
from libtbx import easy_run
from xfel.util.dials_file_matcher import match_dials_files
from xfel.util.mp import mp_phil_str as multiprocessing_str
from xfel.util.mp import get_submit_command_chooser

import os, math

multiprocessing_override_str = '''
mp {
  use_mpi = False
}
'''

striping_str = '''
striping {
  results_dir = None
    .type = path
    .help = "LCLS results directory containint runs starting with r."
  trial = None
    .type = int
    .help = "Trial identifier for an XFEL GUI formatted processing trial."
  stripe = False
    .type = bool
    .help = "Enable to select results evenly spaced across each rungroup"
            "(stripes) as opposed to contiguous chunks."
  chunk_size = 1000
    .type = float
    .help = "Maximum number of images per chunk or stripe."
}
'''

combining_str = '''
combine_experiments {
  interactive = False
    .type = bool
    .help = "Overrides any multiprocessing parameters to allow interactive"
    .help = "run. Clustering dendrograms can only be displayed in this mode."
  keep_integrated = False
    .type = bool
    .help = "Combine refined_experiments.json and integrated.pickle files."
    .help = "If False, ignore integrated.pickle files in favor of"
    .help = "indexed.pickle files in preparation for reintegrating."
  include scope dials.command_line.combine_experiments.phil_scope
}
'''

combining_override_str = '''
combine_experiments {
  output {
    experiments_filename = %s_combined_experiments.json
    reflections_filename = %s_combined_reflections.pickle
    delete_shoeboxes = True
  }
  reference_from_experiment {
    detector = 0
  }
}
'''

# future feature: filter experiments by rmsd after combining/clustering
filtering_str = '''
filtering {
  enable = False
}
'''

refinement_str = '''
refinement {
  include scope dials.command_line.refine.phil_scope
  input {
    experiments = None
    reflections = None
  }
}
'''

refinement_override_str = '''
refinement {
  output {
    experiments = %s_refined_experiments_CLUSTER.json
    reflections = %s_refined_reflections_CLUSTER.pickle
    include_unused_reflections = False
    log = %s_refine_CLUSTER.log
    debug_log = %s_refine_CLUSTER.debug.log
  }
  refinement {
    parameterisation {
      auto_reduction {
        action = remove
      }
      beam {
        fix = all
      }
    }
    refinery {
      engine = SparseLevMar
    }
    reflections {
      outlier {
        algorithm = sauter_poon
        minimum_number_of_reflections = 3
        separate_experiments = False
        separate_panels = False
      }
    }
  }
  input {
    experiments = %s_combined_experiments_CLUSTER.json
    reflections = %s_combined_reflections_CLUSTER.pickle
  }
}
'''

# future feature: reintegrate after joint refinement
postprocessing_str = '''
postprocessing {
  enable = False
}
'''

master_defaults_str = multiprocessing_str + striping_str + combining_str + filtering_str + \
                        refinement_str + postprocessing_str

# initialize a master scope from the multiprocessing phil string
master_defaults_scope = parse(master_defaults_str, process_includes=True)
# update master scope with customized and local phil scopes
master_scope = master_defaults_scope.fetch(parse(refinement_override_str, process_includes=True))
master_scope = master_scope.fetch(parse(combining_override_str, process_includes=True))
master_scope = master_scope.fetch(parse(multiprocessing_override_str, process_includes=True))

def allocate_chunks_per_rungroup(results_dir, trial_no, stripe=False, max_size=1000, integrated=False):
  refl_ending = "_integrated.pickle" if integrated else "_indexed.pickle"
  expt_ending = "_refined_experiments.json"
  trial = "%03d" % trial_no
  print "processing trial %s" % trial
  rgs = {} # rungroups and associated runs
  for run in os.listdir(results_dir):
    if not (run.startswith("r") and run.split("r")[1].isdigit()):
      continue
    trgs = [trg for trg in os.listdir(os.path.join(results_dir, run))
            if trg[:6] == trial + "_rg"]
    rungroups = set(map(lambda n: n.split("_")[1], trgs))
    for rg in rungroups:
      if rg not in rgs.keys():
        rgs[rg] = []
      else:
        rgs[rg].append(run)
  rg_ch_nums_sizes = {}
  rg_contents = {}
  for rg, runs in rgs.iteritems():
    n_img = 0
    trg = trial + "_" + rg
    rg_contents[rg] = []
    for run in runs:
      try:
        contents = os.listdir(os.path.join(results_dir, run, trg, "out"))
      except OSError:
        # print "skipping run %s missing out directory" % run
        continue
      abs_contents = [os.path.join(results_dir, run, trg, "out", c)
                      for c in contents]
      rg_contents[rg].extend(abs_contents)
      expts = [c for c in contents if c.endswith(expt_ending)]
      n_img += len(expts)
    if n_img == 0:
      # print "no images found for %s" % rg
      del rg_contents[rg]
      continue
    n_chunks = int(math.ceil(n_img/max_size))
    chunk_size = int(math.ceil(n_img/n_chunks))
    rg_ch_nums_sizes[rg] = (n_chunks, chunk_size)
  if len(rg_contents) == 0:
    raise Sorry, "no DIALS integration results found."
  rg_chunks = {}
  for rg, nst in rg_ch_nums_sizes.iteritems():
    num, size = nst
    rg_chunks[rg] = []
    contents = rg_contents[rg]
    expts = [c for c in contents if c.endswith(expt_ending)]
    refls = [c for c in contents if c.endswith(refl_ending)]
    expts, refls = match_dials_files(expts, refls, expt_ending, refl_ending)
    if stripe:
      for i in xrange(num):
        expts_stripe = expts[i::num]
        refls_stripe = refls[i::num]
        rg_chunks[rg].append((expts_stripe, refls_stripe))
      print "striped %d experiments in %s with %d experiments per stripe and %d stripes" % \
        (len(expts), rg, len(rg_chunks[rg][0][0]), len(rg_chunks[rg]))
    else:
      for i in xrange(num):
        expts_chunk = expts[i*size:(i+1)*size]
        refls_chunk = refls[i*size:(i+1)*size]
        rg_chunks[rg].append((expts_chunk, refls_chunk))
      print "chunked %d experiments in %s with %d experiments per chunk and %d chunks" % \
        (len(expts), rg, len(rg_chunks[rg][0][0]), len(rg_chunks[rg]))
  return rg_chunks

def parse_retaining_scope(args, master_scope=master_scope):
  if "-c" in args:
    master_scope.show(attributes_level=2)
    return
  file_phil = []
  cmdl_phil = []
  for arg in args:
    if os.path.isfile(arg):
      try:
        file_phil.append(parse(file_name=arg))
      except Exception, e:
        raise Sorry("Unrecognized file: %s" % arg)
    else:
      try:
        cmdl_phil.append(parse(arg))
      except Exception, e:
        raise Sorry("Unrecognized argument: %s" % arg)

  run_scope = master_scope.fetch(sources=file_phil)
  run_scope = run_scope.fetch(sources=cmdl_phil)
  return run_scope

def script_to_expand_over_clusters(clustered_json_name,
                                   phil_template_name, command, location):
  """
  Write a bash script to find results of a clustering step and produce customized
  phils and commands to run with each of them. For example, run the command
  dials.refine ...cluster8.json ...cluster8.pickle ...cluster8.phil followed by
  dials.refine ...cluster9.json ...cluster9.pickle ...cluster9.phil.
  clustered_json_name, clustered_refl_name and phil_template_name must each
  contain an asterisk, and substitution in phil_template itself will occur at
  each instance of CLUSTER.
  """
  clj_part_first, clj_part_last = clustered_json_name.split("CLUSTER")
  ph_part_first, ph_part_last = phil_template_name.split("CLUSTER")
  clustered_template_name = "*".join([clj_part_first, clj_part_last])

  bash_str = '''
#! /bin/sh

for file in `ls {clname}`
  do export cluster=`echo $file | sed "s:{cljfirst}::; s:{cljlast}::"`
  export philname="{phfirst}${cluster}{phlast}"
  export outname=`echo $philname | sed "s:.phil:.out:"`
  sed "s:CLUSTER:${cluster}:g" {phtempl} > $philname
  {command} $philname > $outname
done
'''.format(clname=clustered_template_name, phtempl=phil_template_name,
           cljfirst=clj_part_first, cljlast=clj_part_last,
           phfirst=ph_part_first, phlast=ph_part_last,
           command=command, cluster="{cluster}")

  bash_name = "generator".join([ph_part_first, ph_part_last]).split(".phil")[0] + ".sh"
  with open(os.path.join(location, bash_name), "wb") as script:
    script.write(bash_str)
  return bash_name

class Script(object):

  def __init__(self):
    '''Initialise the script.'''

    # The script usage
    self.master_defaults_scope = master_defaults_scope
    self.run_scope = parse_retaining_scope(sys.argv[1:])
    self.params = self.run_scope.extract()

  def run(self):
    '''Execute the script.'''

    rg_chunks = allocate_chunks_per_rungroup(self.params.striping.results_dir,
                                             self.params.striping.trial,
                                             stripe=self.params.striping.stripe,
                                             max_size=self.params.striping.chunk_size,
                                             integrated=self.params.combine_experiments.keep_integrated)
    dirname = "combine_experiments_t%03d" % self.params.striping.trial
    if not os.path.isdir(dirname):
      os.mkdir(dirname)
    cwd = os.getcwd()
    tag = "stripe" if self.params.striping.stripe else "chunk"
    for rg, ch_list in rg_chunks.iteritems():
      for idx in xrange(len(ch_list)):
        chunk = ch_list[idx]

        # get the diff scope
        diff_scope = self.master_defaults_scope.fetch_diff(self.run_scope)

        # set up the file containing input expts and refls
        filename = "t%03d_%s_%s%03d" % (self.params.striping.trial, rg, tag, idx)
        chunk_path = os.path.join(cwd, dirname, filename)
        if os.path.isfile(chunk_path):
          os.remove(chunk_path)
        with open(chunk_path, "wb") as outfile:
          for i in (0, 1): # expts then refls
            outfile.write("\n".join(chunk[i]) + "\n")

        # set up the params for dials.combine_experiments
        combine_diff_str = diff_scope.get("combine_experiments").as_str() % (filename, filename)
        combine_diff_lines = combine_diff_str.split("\n")[:-2]
        combine_diff_lines.pop(0)
        combine_diff_lines.append("  input {")
        for expt_path in chunk[0]:
          combine_diff_lines.append("    experiments = %s" % expt_path)
        for refl_path in chunk[1]:
          combine_diff_lines.append("    reflections = %s" % refl_path)
        combine_diff_lines.append("  }")
        combine_diff_str = "\n".join(combine_diff_lines)
        combine_phil_filename = filename + "_combine.phil"
        combine_phil_path = os.path.join(cwd, dirname, combine_phil_filename)
        if os.path.isfile(combine_phil_path):
          os.remove(combine_phil_path)
        with open(combine_phil_path, "wb") as phil_outfile:
          phil_outfile.write(combine_diff_str + "\n")

        # set up the params for dials.refine (to be customized per cluster at execution time)
        refine_diff_str = diff_scope.get("refinement").as_str() % \
          (filename, filename, filename, filename, filename, filename)
        refine_diff_parts = refine_diff_str.split("\n")[:-2]
        refine_diff_parts.pop(0)
        refine_diff_str = "\n".join(refine_diff_parts)
        if self.params.combine_experiments.clustering.use:
          refine_phil_filename = filename + "_refine_CLUSTER.phil"
          refine_phil_path = os.path.join(cwd, dirname, refine_phil_filename)
          if os.path.isfile(refine_phil_path):
            os.remove(refine_phil_path)
          with open(refine_phil_path, "wb") as phil_outfile:
            phil_outfile.write(refine_diff_str + "\n")
          script = script_to_expand_over_clusters(
            self.params.refinement.input.experiments[0] % filename,
            refine_phil_filename,
            "dials.refine",
            dirname)
          refine_command = ". %s" % os.path.join(cwd, dirname, script)
        else:
          refine_diff_str = refine_diff_str.replace('_CLUSTER', '')
          refine_phil_filename = filename + "_refine.phil"
          refine_phil_path = os.path.join(cwd, dirname, refine_phil_filename)
          if os.path.isfile(refine_phil_path):
            os.remove(refine_phil_path)
          with open(refine_phil_path, "wb") as phil_outfile:
            phil_outfile.write(refine_diff_str + "\n")
          refine_command = "dials.refine %s" % refine_phil_filename

        # submit queued job from appropriate directory
        os.chdir(dirname)
        submit_path = os.path.join(cwd, dirname, "combine_%s.sh" % filename)
        command = "dials.combine_experiments %s && %s" % \
          (combine_phil_filename, refine_command)
        # import pdb; pdb.set_trace()
        if self.params.combine_experiments.interactive:
          easy_run.fully_buffered(command).raise_if_errors().show_stdout()
        else:
          submit_command = get_submit_command_chooser(command, submit_path, dirname, self.params.mp,
            log_name=(submit_path.split(".sh")[0] + ".out"))
          print "executing command: %s" % submit_command
          try:
            easy_run.fully_buffered(submit_command).raise_if_errors().show_stdout()
          except Exception as e:
            if not "Warning: job being submitted without an AFS token." in str(e):
              raise e
        os.chdir(cwd)
        # go back to working directory

if __name__ == "__main__":
  from dials.util import halraiser
  import sys
  if "-c" in sys.argv[1:]:
    if "-e" in sys.argv[1:]:
      expert_level = int(sys.argv[sys.argv.index("-e") + 1])
    else:
      expert_level = 0
    master_defaults_scope.fetch_diff(master_scope).show(attributes_level=expert_level)
    # master_scope.show(attributes_level=expert_level)
    with open("striping_defaults.phil", "wb") as defaults:
      defaults.write(master_scope.as_str())
    exit()
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
