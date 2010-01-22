from libtbx.utils import Sorry
from libtbx import runtime_utils
from libtbx import easy_pickle, easy_run
from libtbx import adopt_init_args
import sys, os, time

def exercise () :
  run = runtime_utils.simple_run()
  params = runtime_utils.process_master_phil.extract()
  output_dir = os.path.join(os.getcwd(), "simple_run")
  os.makedirs(output_dir)
  params.tmp_dir = output_dir
  params.buffer_stdout = False
  driver = runtime_utils.detached_process_driver(output_dir, run)
  params.run_file = os.path.join(output_dir, "run.pkl")
  easy_pickle.dump(params.run_file, driver)
  easy_run.call("libtbx.start_process %s &" % params.run_file)
  client = runtime_utils.simple_client(params)
  client.run()
  assert (client.out.getvalue() == """\
current is 44444.444444
current is 50000.000000
current is 57142.857143
current is 66666.666667
""")
  assert client.n_cb >= 5 # this is variable!
  assert ([ cb.message for cb in client._accumulated_callbacks ] ==
          ['run 0', 'run 1', 'run 2', 'run 3'])
  print "OK"

if __name__ == "__main__" :
  exercise()
