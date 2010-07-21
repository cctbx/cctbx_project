
import libtbx.tracking
import libtbx.phil
import sys, os, time

class job_history (libtbx.tracking.container) :
  master_phil = libtbx.phil.read_default(__file__)

  def __init__ (self, *args, **kwds) :
    libtbx.tracking.container.__init__(self, *args, **kwds)
    self.params = self.working_phil.extract()

  def initialize (self, *args, **kwds) :
    return self.master_phil.fetch()

  def get_next_job_id (self) :
    return len(self.params.job_summary)

  def start_job (self, app_id, program_name, program_args, config_file,
      title=None, directory=None) :
    job_id = self.get_next_job_id()
    raw_phil = libtbx.phil.parse("""\
      job_summary {
        job_id = %d
        app_id = %s
        program_name = %s
        program_args = %s
        config_file = %s
        title = %s
        time_started = %.1f
        result_directory = %s
      }""" % (job_id, app_id, program_name, " ".join(program_args),
              config_file, str(title), time.time(), str(directory)))
    job_phil = self.master_phil.fetch(source=raw_phil)
    self.working_phil.objects.append(job_phil.objects[1])
    self.params.job_summary.append(job_phil.extract().job_summary[0])
    self.save_file()
    return job_id

  def replace_job_params (self, current_job) :
    job_id = current_job.job_id
    params_stub = self.master_phil.fetch().extract()
    params_stub.job_summary.append(current_job)
    final_phil = self.master_phil.format(python_object=params_stub)
    self.working_phil.objects[job_id + 1] = final_phil.objects[1]
    self.save_file()

  def finish_job (self, job_id, status="complete", r_free=None) :
    assert (status in ["failed", "aborted", "complete", "deleted"])
    assert (r_free is None or isinstance(r_free, float))
    current_job = self.params.job_summary[job_id]
    assert (current_job.job_id == job_id)
    current_job.time_finished = time.time()
    current_job.status = status
    current_job.r_free = r_free
    self.replace_job_params(current_job)

  def set_job_status (self, job_id, status) :
    current_job = self.params.job_summary[job_id]
    assert (current_job.job_id == job_id)
    if current_job.time_finished is None :
      current_job.time_finished = time.time()
    current_job.status = status
    self.replace_job_params(current_job)
    return current_job

  def abort_job (self, job_id) :
    self.set_job_status(job_id, "aborted")

  def job_error (self, job_id) :
    self.set_job_status(job_id, "failed")

  def delete_job (self, job_id, remove_directory=True) :
    current_job = self.set_job_status(job_id, "deleted")
    if remove_directory and (current_job.result_directory is not None) :
      shutil.rmtree(current_job.result_directory)

def exercise (detailed_timings=False) :
  if os.path.isfile("jobs.phil") :
    os.remove("jobs.phil")
  jobs = job_history("jobs.phil")
  n_jobs = 100
  t1 = time.time()
  for i in range(n_jobs) :
    job_id = jobs.get_next_job_id()
    jobs.start_job(
      app_id="Refine",
      program_name="phenix.refine",
      program_args=["/home/nat/refine.eff", "--overwrite"],
      config_file="/home/nat/refine.eff",
      title="Refinement of hen egg lysozyme",
      directory="/home/nat/Refine_%d" % job_id)
    jobs.finish_job(
      job_id=job_id,
      status="complete",
      r_free=0.1234)
  t2 = time.time()
  jobs = job_history("jobs.phil")
  assert (jobs.get_next_job_id() == n_jobs)
  t3 = time.time()
  if detailed_timings :
    print "time to create %d jobs: %5.1fms" % (n_jobs, (t2 - t1) * 1000)
    print "time to reload %d jobs: %5.1fms" % (n_jobs, (t3 - t2) * 1000)

if __name__ == "__main__" :
  exercise(detailed_timings=True)
  print "OK"
