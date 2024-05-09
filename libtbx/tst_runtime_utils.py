from __future__ import absolute_import, division, print_function
from libtbx import runtime_utils
from libtbx import easy_pickle
from libtbx import easy_run
import time
import os
import sys

def exercise():
  params = runtime_utils.process_master_phil.extract()
  i = 0
  while True :
    output_dir = os.path.join(os.getcwd(), "simple_run%d" % i)
    if os.path.exists(output_dir):
      i += 1
    else :
      os.makedirs(output_dir)
      break
  run = runtime_utils.simple_run(output_dir)
  params.output_dir = output_dir
  params.buffer_stdout = False
  params.tmp_dir = output_dir
#  driver = runtime_utils.detached_process_driver(output_dir, run)
  params.run_file = os.path.join(output_dir, "run.pkl")
  eff_file = os.path.join(output_dir, "run.eff")
  working_phil = runtime_utils.process_master_phil.format(python_object=params)
  with open(eff_file, "w") as f:
    working_phil.show(out=f)
  easy_pickle.dump(params.run_file, run)
  assert not easy_run.call("libtbx.start_process %s &" % eff_file) #params.run_file)
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

def exercise2():
  f = runtime_utils.simple_func(666)
  easy_pickle.dump("myfunc.pkl", f)
  f_out = easy_run.fully_buffered(
    "libtbx.run_pickled_function myfunc.pkl").stdout_lines
  assert (f_out[0] == "666")

# queueing system support
def exercise3():
  i = 0
  while True :
    output_dir = os.path.join(os.getcwd(), "simple_run%d" % i)
    if os.path.exists(output_dir):
      i += 1
    else :
      os.makedirs(output_dir)
      break
  run = runtime_utils.simple_run(output_dir)
  params = runtime_utils.process_master_phil.extract()
  server = runtime_utils.detached_process_server(run, params)
  from libtbx.queuing_system_utils import generic as queuing
  job = queuing.qsub(
    name="tst_runtime_utils",
    target=server)
  job.start()
  client = runtime_utils.simple_client(params)
  client.run()
  assert (client.out.getvalue() == """\
current is 44444.444444
current is 50000.000000
current is 57142.857143
current is 66666.666667
""")
  assert client.n_cb >= 5 # this is variable!
  time.sleep(1)
  assert ([ cb.message for cb in client._accumulated_callbacks ] ==
          ['run 0', 'run 1', 'run 2', 'run 3'])

if __name__ == "__main__" :
  exercise()
  exercise2()
  if ("-q" in sys.argv):
    print("Testing queueing system support...")
    exercise3()
  print("OK")
