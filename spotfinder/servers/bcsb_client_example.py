from __future__ import absolute_import, division, print_function
from six.moves import range
import os,time
import threading
output_lock = threading.RLock()

def do_main_apache(filepath, host, port):
  absfile = os.path.abspath(filepath)
  base_url = "http://%s:%d/spotfinder/distl.signal_strength_bcsb?distl.image=%s"%(host,port,absfile)
  if len(DISTL_OPTIONS) > 0:
    base_url = base_url + "&" + "&".join(DISTL_OPTIONS)
  from six.moves import urllib
  try:
    Response = urllib.request.urlopen(base_url)
    log = Response.read()
    Response.close()
    return log
  except Exception as e:
    return str(e)

def single_thread(idx):
  for x in range(1):
    filepath, host, port = [ABS_DATA_TEMPLATE%idx,HOST,PORT]
    port = int(port)
    if OUTPUT_FILE is not None:
      output_lock.acquire(blocking=True)
      outfile = open(OUTPUT_FILE,"ab")
      print(do_main_apache(filepath, host, port), file=outfile)
      outfile.close()
      output_lock.release()
    else:
      print(do_main_apache(filepath, host, port))

def multi_thread():
  from multiprocessing import Pool
  pool = Pool(processes=N_CLIENT_THREADS)
  results = []
  for x in IMAGE_RANGE:
    if SERVER_TYPE=="Python": time.sleep(TIME_DELAY)
    result = pool.apply_async(single_thread, [x,])
    results.append(result)
  for j,item in enumerate(results):
    item.wait()

if __name__=="__main__":
  ABS_DATA_TEMPLATE = "/net/cci/dials/from_sunbird/sauter/rawdata/pilatus/ssrl_P6/all/I3_1_%04d.cbf"
  HOST = "viper"
  PORT = "8125"
  N_CLIENT_THREADS = 48
  IMAGE_RANGE = range(1,721)
  DISTL_OPTIONS = ["distl.res.outer=3.3","distl.bins.verbose=True"]
  SERVER_TYPE = ["Python","Apache"][1] # choose 0=Python, 1=Apache mod-python
  TIME_DELAY = 0.10 # seconds per-image throughput, depends on server
  OUTPUT_FILE = "/net/viper/raid1/sauter/apache/bcsb_results.txt"
  multi_thread()
