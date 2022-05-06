#!/usr/bin/env python
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# filter_experiments_by_rmsd.py
#
#  Copyright (C) 2016 Lawrence Berkeley National Laboratory (LBNL)
#
#  Author: Aaron Brewster and David Waterman
#
#  This code is distributed under the X license, a copy of which is
#  included in the root directory of this package.
#
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.filter_experiments_by_rmsd
#
from __future__ import absolute_import, division, print_function
from six.moves import range
from dials.array_family import flex
from dials.util import show_mail_on_error
from scitbx.matrix import col
from libtbx.phil import parse
import libtbx.load_env
import math

try:
  from matplotlib import pyplot as plt
except ImportError: # Can happen on non-GUI nodes
  pass

help_message = '''
Filter a set of experiments and reflections from a multi-experiment job by overall RMSD
using Tukey's rule of thumb.  I.E, for each experiment, determine the RMSD of the
differences between preditions - observations. Then compute the five number summary of
this set of per-image RMSDs. Then, filter outliers more than iqr_multiplier times the
interquartile range from the third quartile. When x=1.5, this is Tukey's rule.

Example:

  %s combined.expt combined.refl
''' % libtbx.env.dispatcher_name

# Create the phil parameters
phil_scope = parse('''
iqr_multiplier = 1.5
  .type = float
  .help = Interquartile multiplier
show_plots = False
  .type = bool
  .help = Show some plots
max_delta = None
  .type = float
  .help = Before filtering, throw out all reflections with obs-pred greater \
          that max_delta mm.
detector = None
  .type = int
  .help = If not None, only filter experiments matching this detector number
output {
  filtered_experiments = filtered.expt
    .type = str
    .help = Name of output filtered experiments file
  filtered_reflections = filtered.refl
    .type = str
    .help = Name of output filtered reflections file
}
delta_psi_filter = None
  .type = float(value_min=0)
  .help = After RMSD filter, filter remaining reflections by delta psi \
          angle (degrees).
''')

def run_with_preparsed(experiments, reflections, params):
  from dxtbx.model import ExperimentList
  from scitbx.math import five_number_summary

  print("Found", len(reflections), "reflections", "and", len(experiments), "experiments")

  filtered_reflections = flex.reflection_table()
  filtered_experiments = ExperimentList()

  skipped_reflections = flex.reflection_table()
  skipped_experiments = ExperimentList()

  if params.detector is not None:
    culled_reflections = flex.reflection_table()
    culled_experiments = ExperimentList()
    detector = experiments.detectors()[params.detector]
    for expt_id, experiment in enumerate(experiments):
      refls = reflections.select(reflections['id']==expt_id)
      if experiment.detector is detector:
        culled_experiments.append(experiment)
        refls['id'] = flex.int(len(refls), len(culled_experiments)-1)
        culled_reflections.extend(refls)
      else:
        skipped_experiments.append(experiment)
        refls['id'] = flex.int(len(refls), len(skipped_experiments)-1)
        skipped_reflections.extend(refls)

    print("RMSD filtering %d experiments using detector %d, out of %d"%(len(culled_experiments), params.detector, len(experiments)))
    reflections = culled_reflections
    experiments = culled_experiments

  difference_vector_norms = (reflections['xyzcal.mm']-reflections['xyzobs.mm.value']).norms()

  if params.max_delta is not None:
    sel = difference_vector_norms <= params.max_delta
    reflections = reflections.select(sel)
    difference_vector_norms = difference_vector_norms.select(sel)

  data = flex.double()
  counts = flex.double()
  for i in range(len(experiments)):
    dvns = difference_vector_norms.select(reflections['id']==i)
    counts.append(len(dvns))
    if len(dvns) == 0:
      data.append(0)
      continue
    rmsd = math.sqrt(flex.sum_sq(dvns)/len(dvns))
    data.append(rmsd)
  data *= 1000
  subset = data.select(counts > 0)
  print(len(subset), "experiments with > 0 reflections")

  if params.show_plots:
    h = flex.histogram(subset, n_slots=40)
    fig = plt.figure()
    ax = fig.add_subplot('111')
    ax.plot(h.slot_centers().as_numpy_array(), h.slots().as_numpy_array(), '-')
    plt.title("Histogram of %d image RMSDs"%len(subset))

    fig = plt.figure()
    plt.boxplot(subset, vert=False)
    plt.title("Boxplot of %d image RMSDs"%len(subset))
    plt.show()

  outliers = counts == 0
  min_x, q1_x, med_x, q3_x, max_x = five_number_summary(subset)
  print("Five number summary of RMSDs (microns): min %.1f, q1 %.1f, med %.1f, q3 %.1f, max %.1f"%(min_x, q1_x, med_x, q3_x, max_x))
  iqr_x = q3_x - q1_x
  cut_x = params.iqr_multiplier * iqr_x
  outliers.set_selected(data > q3_x + cut_x, True)
  #outliers.set_selected(col < q1_x - cut_x, True) # Don't throw away the images that are outliers in the 'good' direction!

  for i in range(len(experiments)):
    if outliers[i]:
      continue
    refls = reflections.select(reflections['id']==i)
    refls['id'] = flex.int(len(refls), len(filtered_experiments))
    filtered_reflections.extend(refls)
    filtered_experiments.append(experiments[i])

  #import IPython;IPython.embed()
  zeroes = counts == 0
  n_zero = len(counts.select(zeroes))
  print("Removed %d bad experiments and %d experiments with zero reflections, out of %d (%%%.1f)"%(
    len(experiments)-len(filtered_experiments)-n_zero,
    n_zero,
    len(experiments),
    100*((len(experiments)-len(filtered_experiments))/len(experiments))))

  if params.detector is not None:
    crystals = filtered_experiments.crystals()
    for expt_id, experiment in enumerate(skipped_experiments):
      if experiment.crystal in crystals:
        filtered_experiments.append(experiment)
        refls = skipped_reflections.select(skipped_reflections['id'] == expt_id)
        refls['id'] = flex.int(len(refls), len(filtered_experiments)-1)
        filtered_reflections.extend(refls)

  if params.delta_psi_filter is not None:
    delta_psi = filtered_reflections['delpsical.rad']*180/math.pi
    sel = (delta_psi <= params.delta_psi_filter) & (delta_psi >= -params.delta_psi_filter)
    l = len(filtered_reflections)
    filtered_reflections = filtered_reflections.select(sel)
    print("Filtering by delta psi, removing %d out of %d reflections"%(l - len(filtered_reflections), l))

  print("Final experiment count", len(filtered_experiments))
  return filtered_experiments, filtered_reflections

class Script(object):
  ''' Class to parse the command line options. '''

  def __init__(self):
    ''' Set the expected options. '''
    from dials.util.options import ArgumentParser
    import libtbx.load_env

    # Create the option parser
    usage = "usage: %s combined.expt combined.refl" % libtbx.env.dispatcher_name
    self.parser = ArgumentParser(
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
    reflections = flatten_reflections(params.input.reflections)[0]

    filtered_experiments, filtered_reflections = run_with_preparsed(experiments, reflections, params)
    filtered_experiments.as_file(params.output.filtered_experiments)
    filtered_reflections.as_pickle(params.output.filtered_reflections)

if __name__ == '__main__':
  with show_mail_on_error():
    script = Script()
    script.run()
