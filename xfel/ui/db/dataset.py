from __future__ import absolute_import, division, print_function
from xfel.ui.db import db_proxy

class Dataset(db_proxy):
  def __init__(self, app, dataset_id = None, **kwargs):
    db_proxy.__init__(self, app, "%s_dataset" % app.params.experiment_tag, id = dataset_id, **kwargs)
    self.dataset_id = self.id

  def __getattr__(self, name):
    # Called only if the property cannot be found
    if name == "tasks":
      return self.app.get_dataset_tasks(self.id)
    elif name == "tags":
      return self.app.get_dataset_tags(self.id)
    elif name == "versions":
      return self.app.get_dataset_versions(self.id)
    elif name == "latest_version":
      versions = self.app.get_dataset_versions(self.id, latest = True)
      if versions: return versions[0]
      else: return None
    else:
      return super(Dataset, self).__getattr__(name)

  def __setattr__(self, name, value):
    assert name not in ["tasks", "tags", "versions", "latest_version"]
    super(Dataset, self).__setattr__(name, value)

  def add_task(self, task):
    query = "INSERT INTO `%s_dataset_task` (dataset_id, task_id) VALUES (%d, %d)" % (
      self.app.params.experiment_tag, self.id, task.id)
    self.app.execute_query(query, commit=True)

  def remove_task(self, task):
    query = "DELETE FROM `%s_dataset_task` WHERE dataset_id = %d AND task_id = %s" % (
      self.app.params.experiment_tag, self.id, task.id)
    self.app.execute_query(query, commit=True)

  def add_tag(self, tag):
    query = "INSERT INTO `%s_dataset_tag` (dataset_id, tag_id) VALUES (%d, %d)" % (
      self.app.params.experiment_tag, self.id, tag.id)
    self.app.execute_query(query, commit=True)

  def remove_tag(self, tag):
    query = "DELETE FROM `%s_dataset_tag` WHERE dataset_id = %d AND tag_id = %s" % (
      self.app.params.experiment_tag, self.id, tag.id)
    self.app.execute_query(query, commit=True)

class DatasetVersion(db_proxy):
  def __init__(self, app, dataset_version_id = None, **kwargs):
    db_proxy.__init__(self, app, "%s_dataset_version" % app.params.experiment_tag, id = dataset_version_id, **kwargs)
    self.dataset_version_id = self.id

  def __getattr__(self, name):
    # Called only if the property cannot be found
    if name == "jobs":
      return self.app.get_dataset_version_jobs(self.id)
    else:
      return super(DatasetVersion, self).__getattr__(name)

  def __setattr__(self, name, value):
    assert name not in ["jobs"]
    super(DatasetVersion, self).__setattr__(name, value)

  def add_job(self, job):
    query = "INSERT INTO `%s_dataset_version_job` (dataset_version_id, job_id) VALUES (%d, %d)" % (
      self.app.params.experiment_tag, self.id, job.id)
    self.app.execute_query(query, commit=True)

  def remove_job(self, job):
    query = "DELETE FROM `%s_dataset_version_job` WHERE dataset_version_id = %d AND job_id = %s" % (
      self.app.params.experiment_tag, self.id, job.id)
    self.app.execute_query(query, commit=True)
