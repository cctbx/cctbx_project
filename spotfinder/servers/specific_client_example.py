import os,time
import threading
output_lock = threading.RLock()

def do_main_apache(filepath, host, port):
  absfile = os.path.abspath(filepath)
  base_url = "http://%s:%d/spotfinder/distl.signal_strength?distl.image=%s"%(host,port,absfile)
  if len(DISTL_OPTIONS) > 0:
    base_url = base_url + "&" + "&".join(DISTL_OPTIONS)
  import urllib2
  try:
    Response = urllib2.urlopen(base_url)
    log = Response.read()
    Response.close()
    return log
  except Exception,e:
    return str(e)

def single_thread(idx):
  for x in xrange(1):
    filepath, host, port = [ABS_DATA_TEMPLATE%idx,HOST,PORT]
    port = int(port)
    if OUTPUT_FILE is not None:
      output_lock.acquire(blocking=True)
      outfile = open(OUTPUT_FILE,"ab")
      print >>outfile, do_main_apache(filepath, host, port)
      outfile.close()
      output_lock.release()
    else:
      print do_main_apache(filepath, host, port)

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
  ABS_DATA_TEMPLATE = "/net/sunbird/raid1/sauter/rawdata/pilatus/ssrl_P6/all/I3_1_%04d.cbf"
  HOST = "viper"
  PORT = "8125"
  N_CLIENT_THREADS = 48
  IMAGE_RANGE = xrange(1,721)
  DISTL_OPTIONS = ["distl.res.outer=3.3","distl.bins.verbose=True"]
  SERVER_TYPE = ["Python","Apache"][1] # choose 0=Python, 1=Apache mod-python
  TIME_DELAY = 0.10 # seconds per-image throughput, depends on server
  OUTPUT_FILE = "/net/viper/raid1/sauter/apache/spotfinder_results.txt"
  multi_thread()
