from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME cctbx.stripe_experiment
from libtbx.phil import parse
from libtbx.utils import Sorry
from xfel.util.dials_file_matcher import match_dials_files

import os, math

phil_scope = parse('''
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
  combine_experiments {
    run_serial = False
      .type = bool
      .help = "Run dials.combine_experiments on the resulting files upon"
              "successful completion of striping."
    phil = None
      .type = path
      .help = "Path to a phil file to be used with dials.combine_experiments."
    keep_integrated = False
      .type = bool
      .help = "Combine refined_experiments.json and integrated.pickle files."
      .help = "If False, ignore integrated.pickle files in favor of"
      .help = "indexed.pickle files in preparation for reintegrating."
  }
''', process_includes=True)

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
        expts_stripe = expts[i::size]
        refls_stripe = expts[i::size]
        rg_chunks[rg].append((expts_stripe, refls_stripe))
      print "striped %s with %d experiments per stripe and %d stripes" % \
        (rg, len(rg_chunks[rg][0][0]), len(rg_chunks[rg]))
    else:
      for i in xrange(num):
        expts_chunk = expts[i*size:(i+1)*size]
        refls_chunk = refls[i*size:(i+1)*size]
        rg_chunks[rg].append((expts_chunk, refls_chunk))
      print "chunked %s with %d experiments per chunk and %d chunks" % \
        (rg, len(rg_chunks[rg][0][0]), len(rg_chunks[rg]))
  return rg_chunks


class Script(object):

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    import libtbx.load_env

    # The script usage
    usage  = "usage: %s [options] [param.phil] " \
             "*refined_experiments.json [*integrated.pickle]" \
             "[*indexed.pickle]" \
             % libtbx.env.dispatcher_name

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope)

  def run(self):
    '''Execute the script.'''

    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=False)

    # Set up template phil file for dials.combine_experiments
    import shutil
    template = "template_combine_experiments_phil"
    if params.combine_experiments.phil is not None:
      shutil.copyfile(params.combine_experiments.phil, template)
    else:
      open(template, "wb").close()

    if params.combine_experiments.run_serial:
      from libtbx import easy_run

    rg_chunks = allocate_chunks_per_rungroup(params.results_dir,
      params.trial, stripe=params.stripe, max_size=params.chunk_size,
      integrated=params.combine_experiments.keep_integrated)
    dirname = "combine_experiments_t%03d" % params.trial
    if not os.path.isdir(dirname):
      os.mkdir(dirname)
    tag = "stripe" if params.stripe else "chunk"
    for rg, ch_list in rg_chunks.iteritems():
      for idx in xrange(len(ch_list)):
        chunk = ch_list[idx]
        filename = "t%03d_%s_%s%03d" % (params.trial, rg, tag, idx)
        chunk_path = os.path.join(dirname, filename)
        if os.path.isfile(chunk_path):
          os.remove(chunk_path)
        with open(chunk_path, "wb") as outfile:
          for i in (0, 1): # expts then refls
            outfile.write("\n".join(chunk[i]) + "\n")
        new_phil_filename = filename + ".phil"
        new_phil_path = os.path.join(dirname, new_phil_filename)
        if os.path.isfile(new_phil_path):
          os.remove(new_phil_path)
        shutil.copyfile(template, new_phil_path)
        with open(os.path.join(dirname, new_phil_filename), "ab") as phil_outfile:
          phil_outfile.write("input {\n")
          for expt_path in chunk[0]:
            phil_outfile.write("  experiments = %s\n" % expt_path)
          for refl_path in chunk[1]:
            phil_outfile.write("  reflections = %s\n" % refl_path)
          phil_outfile.write("}\n")
          phil_outfile.write("output {\n")
          phil_outfile.write("  experiments_filename = %s_combined_experiments.json\n" % filename)
          phil_outfile.write("  reflections_filename = %s_combined_reflections.pickle\n" % filename)
          if params.combine_experiments.keep_integrated:
            phil_outfile.write("  delete_shoeboxes = True\n")
          phil_outfile.write("}\n")
        if params.combine_experiments.run_serial:
          working_dir = os.getcwd()
          os.chdir(dirname)
          command = "dials.combine_experiments %s" % new_phil_filename
          print "executing command \"%s\" in location %s" % (command, dirname)
          easy_run.fully_buffered(command).raise_if_errors().show_stdout()
          os.chdir(working_dir)
    os.remove(template)

if __name__ == "__main__":
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
