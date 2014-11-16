from __future__ import division
import libtbx.load_env
from libtbx import easy_run
import os

class test_cxi_merge(object):
  def __init__(self):

    self.dials_regression = libtbx.env.find_in_repositories(
        relative_path="xfel_regression",
        test=os.path.isdir)

  def test_merge_themolysin(self):
    if self.dials_regression is None:
      print "Skipping test_merge_themolysin: xfel_regression not present"
      return

    merging_dir = os.path.join(self.dials_regression, "merging_test_data")

    if not os.path.exists(merging_dir):
      print "Skipping test_merge_themolysin: merging directory %s not present"%merging_dir
      return

    from libtbx.test_utils import open_tmp_directory
    cwd = os.path.abspath(os.curdir)
    tmp_dir = open_tmp_directory(suffix="test_cxi_merge")
    os.chdir(tmp_dir)

    print "Merging results will be found at", tmp_dir

    command = os.path.join(merging_dir, "merge_thermo.csh")
    result = easy_run.fully_buffered(command=command).raise_if_errors()


    noanom_log = "thermonoanom_2tli_mark0.log"
    assert os.path.exists(noanom_log)

    f = open(noanom_log)
    result = f.readlines()[-6]
    assert float(result.split()[1].strip("%")) == 42.2 # CC1/2


    anom_log = "thermoanom_2tli_mark0.log"
    assert os.path.exists(noanom_log)

    f = open(anom_log)
    result = f.readlines()[-6]
    assert float(result.split()[1].strip("%")) == 33.0 # CC1/2


    os.chdir(cwd)

    print "OK"

  def run_all(self):
    self.test_merge_themolysin()


if __name__ == '__main__':
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()

  tester = test_cxi_merge()
  tester.run_all()

