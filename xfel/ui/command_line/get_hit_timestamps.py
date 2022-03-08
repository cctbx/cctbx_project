from __future__ import absolute_import, division, print_function

from libtbx.phil import parse
from libtbx.utils import Sorry
from xfel.ui import db_phil_str
from xfel.ui.db.xfel_db import xfel_db_application
import sys, os

"""
Jiffy script to print out all the hits from the database. If the user was dumping all shots, they can
use this script to print all the file paths of all the hits.
"""

phil_str = """
  trial = None
    .type = int
  hit_cutoff = 16
    .type = int
    .help = Number of reflections to consider an image a hit. Estimate by looking at plot of strong reflections/image.
  output_folder = None
    .type = str
    .help = Folder with results from cctbx.xfel processing
"""
phil_scope = parse(phil_str + db_phil_str)

def run(args):
  user_phil = []
  for arg in args:
    try:
      user_phil.append(parse(arg))
    except Exception as e:
      raise Sorry("Unrecognized argument %s"%arg)
  params = phil_scope.fetch(sources=user_phil).extract()

  app = xfel_db_application(params)
  trial = app.get_trial(trial_number = params.trial)
  jobs = app.get_all_jobs()
  for rungroup in trial.rungroups:
    for run in rungroup.runs:
      if params.output_folder is None:
        job = None
        all_path = None
      else:
        job_found = None
        for job in jobs:
          if job.trial.id == trial.id and job.rungroup.id == rungroup.id and job.run.id == run.id:
            all_path = os.path.join(os.path.sep.join(job.get_log_path().split(os.path.sep)[:-2]), "all")
            assert job_found is None
            job_found = job
        assert job_found is not None
        job = job_found

      where = "WHERE event.n_strong >= %d" % params.hit_cutoff
      events = app.get_all_events(trial = trial, runs = [run], only_indexed = False, where = where)
      for e in events:
        t = e.timestamp
        ts = t[0:4] + t[5:7] + t[8:10] + t[11:13] + t[14:16] + t[17:19] + t[20:23]
        if all_path is None:
          print("%04d"%run.run, ts)
        else:
          print(os.path.join(all_path, "shot-" + ts + ".pickle"))

if __name__ == "__main__":
  run(sys.argv[1:])
