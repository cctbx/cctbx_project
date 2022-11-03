from __future__ import absolute_import, division, print_function
from xfel.ui.db.job import AveragingJob


def AveragingCommand(run, params, trial_num, skip_images, num_images, address, raw):
  job = AveragingJob(run.app)
  job.run = run
  if run.app.params.facility.name == 'lcls':
    job.address = address
  job.trial = job.run.app.create_trial(trial=trial_num)
  job.skip_images = skip_images
  job.num_images = num_images
  job.submit()
