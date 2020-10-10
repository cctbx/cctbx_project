from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME simtbx.diffBragg.stage_one

from libtbx.mpi4py import MPI
from dxtbx.model import ExperimentList
import os
import time
from simtbx.nanoBragg.utils import H5AttributeGeomWriter
from simtbx.diffBragg.utils import image_data_from_expt

COMM = MPI.COMM_WORLD

if COMM.rank > 0:
  import warnings
  warnings.filterwarnings("ignore")

import sys
import json
import numpy as np
from copy import deepcopy

from simtbx.diffBragg.phil import philz
from simtbx.diffBragg import refine_launcher
from dxtbx.model.experiment_list import ExperimentListFactory
from dials.array_family import flex
from libtbx.phil import parse
from dials.util import show_mail_on_error

help_message = "stage 1 (per-shot) diffBragg refinement"

script_phil = """
exper_id = None
  .type = int
  .help = index of an experiment to process
exper_refls_file = None
  .type = str
  .help = path to two column text file specifying exp json files and refls files as pairs
  .help = to be read and processed together
usempi = False
  .type = bool
  .help = process using mpi
show_timing = True
  .type = bool
  .help = print a refinement duration for each iteration experiment
output {
  directory = .
    .type = str
    .help = path where output files and folders will be written
  save {
    images = None* model model_and_data
      .type = choice
      .help = if model, output an image of the model pixels at the ROIs that
      .help = were used during refinement. If model_and_data, then output
      .help = model, data and model-data
    reflections = False
      .type = bool
      .help = if True, output a refined reflections table with xyz.calc 
      .help = computed by the diffBragg model
    pandas = False       
      .type = bool
      .help = whether to save a pandas output file
    experiments = False
      .type = bool
      .help = whether to save a refined experients output file
  } 
  tag {
    images = after_stage_ine
      .type = str
      .help = output file tag for model images
    reflections = after_stage_one
      .type = str
      .help = output file tag for reflections 
    pandas = info
      .type = str
      .help = output file tag for pandas
    experiments = after_stage_one
      .type = str
      .help = output file tag for experiments
  }
}
"""

philz = script_phil + philz
phil_scope = parse(philz)

class Script:

  def __init__(self):
    from dials.util.options import OptionParser

    self.parser = None
    if COMM.rank==0:
      self.parser = OptionParser(
        usage="", #stage 1 (per-shot) diffBragg refinement",
        sort_options=True,
        phil=phil_scope,
        read_experiments=True,
        read_reflections=True,
        check_format=False,
        epilog=help_message)
    self.parser = COMM.bcast(self.parser)

    self.refls = None
    self.explist = None
    self.has_loaded_expers = False

  def _load_exp_refls_fnames(self):
    self.input_expnames, self.input_reflnames = [], []
    if COMM.rank == 0:
      with open(self.params.exper_refls_file, "r") as exper_refls_file:
        for line in exper_refls_file:
          line = line.strip()
          line_split = line.split()
          if len(line_split) != 2:
            print("Weird input line %s, continuing" % line)
            continue
          exp, ref = line_split
          self.input_expnames.append(exp)
          self.input_reflnames.append(ref)
    self.input_expnames = COMM.bcast(self.input_expnames)
    self.input_reflnames = COMM.bcast(self.input_reflnames)

  def _generate_exp_refl_pairs(self):
    if self.has_loaded_expers:
      for i_exp, exper in enumerate(self.explist):
        if self.params.usempi and i_exp % COMM.size != COMM.rank:
          continue
        elif self.params.exper_id is not None and self.params.exper_id != i_exp:
          continue
        refls_for_exper = self.refls.select(self.refls['id'] == i_exp)
        
        # little hack to check the format now 
        El = ExperimentList()
        El.append(exper)
        El = ExperimentListFactory.from_dict(El.to_dict())
        exp_filename = self.params.input.experiments[i_exp].filename
        yield exp_filename, El[0], refls_for_exper

    elif self.params.exper_refls_file is not None:
      self._load_exp_refls_fnames()
      count = 0
      for exp_f, refls_f in zip(self.input_expnames, self.input_reflnames):
        refls = flex.reflection_table.from_file(refls_f)
        nexper_in_refls = len(set(refls['id']))
        for i_exp in range(nexper_in_refls):
          if self.params.usempi and count % COMM.size != COMM.rank:
            count += 1
            continue
          exper = self._exper_json_single_file(exp_f, i_exp)
          refls_for_exper = refls.select(refls['id'] == i_exp)
          count += 1
          yield exp_f, exper, refls_for_exper

  def _exper_json_single_file(self, exp_file, i_exp=0):
    """
    load a single experiment from an exp_file
    If working with large combined experiment files, we only want to load
    one image at a time on each MPI rank, otherwise at least one rank would need to
    load the entire file into memory.
    :param exp_file:
    :param i_exp:
    :return:
    """
    exper_json = json.load(open(exp_file))
    nexper = len(exper_json["experiment"])
    assert 0 <= i_exp < nexper

    this_exper = exper_json["experiment"][i_exp]

    new_json = {'__id__': "ExperimentList", "experiment": [deepcopy(this_exper)]}

    for model in ['beam', 'detector', 'crystal', 'imageset']:
      if model in this_exper:
        model_index = this_exper[model]
        new_json[model] = [exper_json[model][model_index]]
        new_json["experiment"][0][model] = 0
    new_json["scan"] = []
    new_json["goniometer"] = []
    new_json["profile"] = []
    new_json["scaling_model"] = []
    explist = ExperimentListFactory.from_dict(new_json)
    assert len(explist) == 1
    return explist[0]

  def run(self):
    from dials.util.options import flatten_experiments, flatten_reflections
    self.params = None
    if COMM.rank == 0:
      self.params, _ = self.parser.parse_args(show_diff_phil=True)
    self.params = COMM.bcast(self.params)
    if COMM.size > 1 and not self.params.usempi:
      if COMM.rank == 0:
        print("Using %d MPI ranks, but usempi option is set to False, please try again with usempi=True" % COMM.size)
      sys.exit()

    self.explist = flatten_experiments(self.params.input.experiments)
    reflist = flatten_reflections(self.params.input.reflections)

    if self.explist and self.params.exper_refls_file is not None:
      print("Can't an input-experiments-glob AND exper_refls_file, please to one or the other.")
      sys.exit()

    self.has_loaded_expers = False
    if reflist:
      self.refls = reflist[0]
      for r in reflist[1:]:
        self.refls.extend(r)
      self.has_loaded_expers = True
    if not self.has_loaded_expers and self.params.exper_refls_file is None:
      print("No experiments to process")
      sys.exit()

    i_processed = 0
    for exper_filename, exper, refls_for_exper in self._generate_exp_refl_pairs():
      
      assert len(set(refls_for_exper['id'])) == 1
      
      exp_id = refls_for_exper['id'][0]
      
      print(exper_filename, COMM.rank, set(refls_for_exper['id']))
      
      if self.params.output.save.reflections:
        self.params.refiner.record_xy_calc = True
      
      refine_starttime = time.time()
      refiner = refine_launcher.local_refiner_from_parameters(refls_for_exper, exper, self.params)
      if self.params.show_timing:
        print("Time to refine experiment: %f" % (time.time()- refine_starttime))
      
      basename,_ = os.path.splitext(os.path.basename(exper_filename))


      # Save model image
      if self.params.output.save.images is not None:
        images_outdir = os.path.join(self.params.output.directory, "model_images", "rank%d" % COMM.rank)
        if not os.path.exists(images_outdir):
          os.makedirs(images_outdir)
        img_path = os.path.join(images_outdir, "%s_%s_%d.h5" % (self.params.output.tag.images, basename, i_processed))
        panel_Xdim, panel_Ydim = exper.detector[0].get_image_size()
        img_shape = len(exper.detector), panel_Xdim, panel_Ydim
        writer_args = {"filename": img_path , 
            "image_shape": img_shape, 
            "num_images":1 if self.params.output.save.images=="model" else 3, 
            "detector": exper.detector , "beam": exper.beam}
        model_img = refiner.get_model_image()
        with H5AttributeGeomWriter(**writer_args) as writer:
          if self.params.output.save.images=="model_and_data":
              model_img *= self.params.refiner.adu_per_photon
              data = image_data_from_expt(exper) 
              pids, ys, xs = np.where(model_img==0)
              model_img[pids, ys, xs] = data[pids, ys, xs]
              writer.add_image(model_img)
              writer.add_image(data)
              writer.add_image(model_img-data)
          else:
              writer.add_image(model_img)

      
      # Save reflections 
      if self.params.output.save.reflections:
        refined_refls = refiner.get_refined_reflections(refls_for_exper)
        #NOTE do we really need to reset the id ? 
        refined_refls['id'] = flex.int(len(refined_refls), 0)
        refls_outdir = os.path.join(self.params.output.directory, "reflections_after_stage1", "rank%d" % COMM.rank)
        if not os.path.exists(refls_outdir):
          os.makedirs(refls_outdir)
        refls_path = os.path.join(refls_outdir, "%s_%s_%d.refl" % (self.params.output.tag.reflections, basename, i_processed))
        refined_refls.as_file(refls_path)

      # save pandas
      if self.params.output.save.pandas:
        pandas_outdir = os.path.join(self.params.output.directory, "pandas_pickles", "rank%d" % COMM.rank)
        if not os.path.exists(pandas_outdir):
          os.makedirs(pandas_outdir)
        outpath = os.path.join(pandas_outdir, "%s_%s_%d.pkl" % (self.params.output.tag.pandas,basename, i_processed))
        refiner.save_lbfgs_x_array_as_dataframe(outname=outpath)

      # save experiment
      if self.params.output.save.experiments:
        exp_outdir = os.path.join(self.params.output.directory, "experiments_after_stage1", "rank%d" % COMM.rank)
        if not os.path.exists(exp_outdir):
          os.makedirs(exp_outdir)

        exp_path = os.path.join(exp_outdir, "%s_%s_%d.expt" % (self.params.output.tag.experiments, basename, exp_id))
        exper.crystal = refiner.get_corrected_crystal(i_shot=0)
        exper.detector = refiner.get_optimized_detector()
        new_exp_list = ExperimentList()
        new_exp_list.append(exper)
        new_exp_list.as_file(exp_path)

      i_processed += 1

if __name__ == '__main__':
  with show_mail_on_error():
    script = Script()
    script.run()
