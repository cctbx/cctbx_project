
import libtbx.tracking
import libtbx.phil
from libtbx.utils import Sorry
import sys, os, time


def format_job (job_id=None, app_id=None, title=None, directory=None,
    input_files=()) :
  assert ([job_id, app_id].count(None) == 0)
  files_phil = []
  for (file_name, file_type, param_name) in input_files :
    if file_name is None :
      continue
    files_phil.append("""\
  input_file {
    file_name = %s
    file_type = %s
    param_name = %s
  }""" % (file_name, str(file_type), str(param_name)))
  input_files_formatted = "\n".join(files_phil)
  working_phil = libtbx.phil.parse("""
job {
  job_id = %d
  app_id = %s
  title = %s
  directory = %s
%s
}
""" % (job_id, app_id, str(title), str(directory), input_files_formatted))
  return working_phil

class job_info (libtbx.tracking.container) :
  master_phil = libtbx.phil.read_default(__file__)

  def initialize (self, *args, **kwds) :
    return format_job(*args, **kwds)

  def finish_job (self, statistics=(), output_files=()) :
    info_phil = []
    for (file_name, file_type) in output_files :
      if file_name is None :
        continue
      info_phil.append("""\
  output_file {
    file_name = %s
    file_type = %s
  }""" % (file_name, str(file_type)))
    for (stat_name, stat_value) in statistics :
      if stat_name is None :
        continue
      assert type(stat_value).__name__ in ["str", "int", "float"]
      info_phil.append("""\
  statistic {
    stat_name = %s
    stat_value = %s
  }""" % (stat_name, str(stat_value)))
    output_phil = libtbx.phil.parse("""\
job {
  %s
}""" % ("\n".join(info_phil)))
    self.working_phil = self.master_phil.fetch(
      sources=[self.working_phil, output_phil])
    self.save_file()
    return True

  def add_comments (self, comments) :
    comments_phil = libtbx.phil.parse("""job.user_comments = "%s"\n""" %
      comments)
    self.working_phil = self.master_phil.fetch(
      sources=[self.working_phil,comments_phil])
    self.save_file()
    return True

#-----------------------------------------------------------------------
def exercise (detailed_timings=False, verbose=False) :
  import time, string
  job_file = "exercise_job.phil"
  t1 = time.time()
  new_job = job_info(
    file_name=job_file,
    job_id=15,
    app_id="Test",
    title="exercise of job info I/O",
    directory="/home/nat/Refine_15",
    input_files=[("/home/nat/data.mtz", "X-ray data", "data"),
                 ("/home/nat/model.pdb", "PDB file", "model"),
                 ("/home/nat/XX1.cif", "Restraints", "monomer"),
                 ("/home/nat/XX2.cif", "Restraints", "monomer"),
                 ("/home/nat/phases.mtz", "Experimental phases", "phases")])
  t2 = time.time()
  new_job.finish_job(
    statistics=[("R-work", 0.2456),
                ("R-free", "0.2890"),
                ("RMS(bonds)", 0.016),
                ("RMS(angles)", "1.59"),
                ("Residues placed", 457),
                ("Waters", "100")],
    output_files=[("model_refine_15.pdb", "Output model"),
                  ("model_refine_15_map_coeffs.mtz", "Map coefficients"),
                  ("model_refine_15.log", "Log file"),
                  ("model_refine_15.geo", "Summary of geometry restraints"),
                  ("model_refine_15.eff", "Effective final parameters")])
  t3 = time.time()
  new_job.add_comments(string.uppercase * 10)
  new_job.add_comments(string.lowercase * 10)
  t4 = time.time()
  info = job_info(file_name=job_file)
  job_params = info.extract()
  t5 = time.time()
  assert (job_params.job.user_comments == (string.lowercase * 10))
  if detailed_timings :
    print "start_job:        %5.1fms" % ((t2 - t1) * 1000)
    print "finish_job:       %5.1fms" % ((t3 - t2) * 1000)
    print "read and extract: %5.1fms" % ((t5 - t4) * 1000)
  if verbose :
    info.working_phil.show()

if __name__ == "__main__" :
  exercise(detailed_timings=True, verbose=("--debug" in sys.argv))
  print "OK"
