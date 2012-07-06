
from libtbx.program_utils import statistics_info
from libtbx import adopt_init_args
import os

class program_result (object) :
  def __init__ (self,
                program_name,
                job_title,
                directory=None,
                log_file=None,
                input_files=(),
                pdb_files=(),
                map_file=None,
                data_file=None,
                cif_files=(),
                phil_files=(),
                other_files=(),
                statistics={}) :
    adopt_init_args(self, locals())

  def get_output_dir (self) :
    return self.directory

  def get_statistic (self, name) :
    return self.statistics.get(name, None)

  def get_pdb_files (self) :
    return [ self.get_file_path(fn) for fn in self.pdb_files ]

  def get_pdb_file (self) :
    assert (len(self.pdb_files) <= 1)
    if (len(self.pdb_files) == 0) :
      return None
    return self.get_file_path(self.pdb_files[0])

  def show_summary (self, out) :
    pass

  def get_file_path (self, file_name) :
    if (file_name is None) :
      return None
    elif (os.path.isabs(file_name)) :
      return file_name
    else :
      assert (self.directory is not None)
      return os.path.join(self.directory, file_name)

  def get_map_file (self) :
    return self.get_file_path(self.map_file)

  def get_data_file (self) :
    return self.get_file_path(self.data_file)

  def get_phil_files (self) :
    return [ self.get_file_path(fn) for fn in self.phil_files ]

  def r_free (self) :
    return self.get_statistic("r_free")

  def format_statistics (self, stat_keys) :
    formatted = []
    for stat_name in stat_keys :
      stat_value = self.get_statistic(stat_name)
      if (stat_value is not None) :
        stat_label = statistics_info.keys_and_labels.get(stat_name, stat_name)
        format = statistics_info.get_format(stat_label)
        formatted.append((stat_label, format % stat_value))
    return formatted

  def get_pdb_file_caption (self) :
    return "Model"

  def get_map_file_caption (self) :
    return "Map coefficients"

  def finish_job (self) :
    output_files = []
    for file_name in self.pdb_files :
      output_files.append((file_name, self.get_pdb_file_caption()))
    for file_name in self.map_file :
      output_files.append((file_name, self.get_map_file_caption()))
    stats = self.get_final_stats()
    return output_files, stats

  def get_final_stats (self) :
    return []
