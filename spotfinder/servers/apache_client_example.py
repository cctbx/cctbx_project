import os

def do_main_apache(filepath, host, port):
  absfile = os.path.abspath(filepath)
  import urllib2
  Response = urllib2.urlopen(
   "http://%s:%d/spotfinder/distl.signal_strength?filename=%s"%(
   host,port,absfile))
  log = Response.read()
  Response.close()
  return log

def single_thread(idx):
  for x in xrange(1):
    filepath, host, port = ["/net/sunbird/raid1/sauter/rawdata/pilatus/ssrl_P6/all/I3_1_%04d.cbf"%idx,
    "apache.lbl.gov","8125"]
    port = int(port)
    print do_main_apache(filepath, host, port)

def multi_thread():
  from multiprocessing import Pool
  pool = Pool(processes=20)
  results = []
  for x in xrange(1,721):
    result = pool.apply_async(single_thread, [x,])
    results.append(result)

  for j,item in enumerate(results):
    item.wait()

if __name__=="__main__":
  multi_thread()
