from __future__ import absolute_import, division, print_function
from six.moves import range
#!/usr/bin/env python
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# recompute_mosaicity.py
#
#  Copyright (C) 2017 Lawrence Berkeley National Laboratory (LBNL)
#
#  Author: Aaron Brewster
#
#  This code is distributed under the X license, a copy of which is
#  included in the root directory of this package.
#
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.recompute_mosaicity
#
from dxtbx.model.experiment_list import ExperimentListDumper
from dials.algorithms.indexing.nave_parameters import nave_parameters
from dials.array_family import flex
import libtbx.load_env
from libtbx.phil import parse

help_message = '''
Recompute mosaic parameters for a set of experiments and apply outlier rejection.

Example:

  %s refined_experiments.json refined_reflections.pickle
''' % libtbx.env.dispatcher_name

# Create the phil parameters
phil_scope = parse("""
plot_changes = False
  .type = bool
  .help = If True, plot the change in the mosaic parameters
output {
  experiments = refined_experiments.json
    .type = str
    .help = Name of output  experiments file
  reflections = refined_reflections.pickle
    .type = str
    .help = Name of output reflections file
}
""")

class Script(object):
  ''' Class to parse the command line options. '''

  def __init__(self):
    ''' Set the expected options. '''
    from dials.util.options import OptionParser
    import libtbx.load_env

    # Create the option parser
    usage = "usage: %s refined_experiments.json refined_reflections.pickle" % libtbx.env.dispatcher_name
    self.parser = OptionParser(
      usage=usage,
      sort_options=True,
      phil=phil_scope,
      read_experiments=True,
      read_reflections=True,
      check_format=False,
      epilog=help_message)

  def run(self):
    ''' Parse the options. '''
    from dials.util.options import flatten_experiments, flatten_reflections
    # Parse the command line arguments
    params, options = self.parser.parse_args(show_diff_phil=True)
    experiments = flatten_experiments(params.input.experiments)
    reflections = flatten_reflections(params.input.reflections)
    assert len(reflections) == 1
    reflections = reflections[0]

    domain_size = flex.double()
    mosaic_angle = flex.double()
    filtered_reflections = flex.reflection_table()

    for i in range(len(experiments)):
      refls = reflections.select(reflections['id'] == i)
      try:
        nv = nave_parameters(params = None, experiments=experiments[i:i+1], reflections=refls, refinery=None, graph_verbose=False)
        crystal_model_nv = nv()
      except Exception as e:
        continue
      domain_size.append(experiments[i].crystal.get_domain_size_ang() - crystal_model_nv.get_domain_size_ang())
      mosaic_angle.append(experiments[i].crystal.get_half_mosaicity_deg() - crystal_model_nv.get_half_mosaicity_deg())
      experiments[i].crystal = crystal_model_nv

      refls = refls.select(nv.nv_acceptance_flags)
      filtered_reflections.extend(refls)

    print("Saving new experiments as %s"%params.output.experiments)
    dump = ExperimentListDumper(experiments)
    dump.as_json(params.output.experiments)

    print("Removed %d out of %d reflections as outliers"%(len(reflections) - len(filtered_reflections), len(reflections)))
    print("Saving filtered reflections as %s"%params.output.experiments)
    filtered_reflections.as_pickle(params.output.reflections)

    if params.plot_changes:
      from matplotlib import pyplot as plt
      domain_size = domain_size.select((domain_size >= -10) & (domain_size <= 10))
      mosaic_angle = mosaic_angle.select((mosaic_angle >= -0.1) & (mosaic_angle <= 0.1))

      for d in [domain_size, mosaic_angle]:
        f = plt.figure()
        plt.hist(d, bins=30)
      plt.show()

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
