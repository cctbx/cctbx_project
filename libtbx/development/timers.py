from __future__ import absolute_import, division, print_function
import time,sys

# Restore compatibility to Python >= 3.8, where time.clock() is gone
py_version = sys.version_info
if py_version[0]==3 and py_version[1]>=7:
  #work_clock = time.perf_counter # this is an elapsed time, no good for OpenMP
  work_clock = time.process_time # includes the extra work done with OpenMP
else:
  work_clock = time.clock

class Timer:
  """While the class instance is in scope it accumulates CPU and wall clock
     elapsed time.  A report is printed when the instances goes out of scope"""
  def __init__(self,message):
    self.start = work_clock()
    self.start_el = time.time()
    self.message = message
    print("start timing %s"%self.message)

  def tick(self):
    return work_clock() - self.start

  def __del__(self):
    self.end = work_clock()
    self.end_el = time.time()
    print("time for %s: CPU, %8.3fs; elapsed, %8.3fs"%(
      self.message,self.end-self.start,self.end_el-self.start_el))

class DebuggingTimer:
  def __init__(self,message,filename,filemode = 'a'):
    self.start = work_clock()
    self.start_el = time.time()
    self.message = message
    self.k = open(filename,filemode)
    self.current_out = sys.stdout
    self.current_err = sys.stderr
    sys.stdout = self.k
    sys.stderr = self.k
    print("start timing %s"%self.message)


  def __del__(self):
    self.end = work_clock()
    self.end_el = time.time()
    print("time for %s: CPU, %8.3fs; elapsed, %8.3fs"%(
      self.message,self.end-self.start,self.end_el-self.start_el))
    sys.stdout = self.current_out
    sys.stderr = self.current_err
    self.k.flush()
    self.k.close()

class cumulative:
  def __init__(self,tag):
    self.tag = tag
    self.total = 0.0
  def add(self,amt):
    self.total+=amt
  def __del__(self):
    print("tag",self.tag,"total",self.total)

class singleton_data(dict):
  def __str__(self):
    message = []
    for key in self:
      message.append("time for %30s: CPU, %7.3fs; elapsed, %7.3fs, averaging %7.3fms #calls:%4d"%(key,
        self[key][0],self[key][1],1000.*self[key][1]/self[key][2],self[key][2]))
    return "\n".join(message)

  def __del__(self):
    Nkeys = len(self)
    if Nkeys > 0: print("Exiting profiler")
    total_tm =0.
    total_el =0.
    for key in self:
      print("time for %30s: CPU, %8.3fs; elapsed, %8.3fs, #calls:%4d"%(key,
        self[key][0],self[key][1],self[key][2]))
      total_tm+=self[key][0]
      total_el+=self[key][1]
    if Nkeys > 0:
      print("TOTAL    %30s: CPU, %8.3fs; elapsed, %8.3fs"%("",
        total_tm,total_el))


timing_singleton=singleton_data()

class Profiler:
  """Each time the class is instantiated with the same 'message' it turns on a
     timer to accumulate CPU-elapsed and wall clock-elapsed time, until
     the instance goes out of scope.  When
     the program executable goes out of scope a final report is printed,
     explaining how many instantiations took place and the total time
     accumulated.  The Profiler will separately keep track of timings that
     are instantiated with different 'message's."""
  def __init__(self,message):
    self.start = work_clock()
    self.start_el = time.time()
    self.message = message

  def __enter__(self):pass
  def __exit__(self,exception_type,exception_value,traceback):pass

  def __del__(self):
    self.end = work_clock()
    self.end_el = time.time()
    if self.message not in timing_singleton:
      timing_singleton[self.message]=[0.,0.,0]
    timing_singleton[self.message][0]+=self.end-self.start
    timing_singleton[self.message][1]+=self.end_el-self.start_el
    timing_singleton[self.message][2]+=1

    print("individual call time for %s: CPU, %8.3fs; elapsed, %8.3fs"%(
      self.message,self.end-self.start,self.end_el-self.start_el))

class SlimProfiler(Profiler):
  def __del__(self):
    self.end = work_clock()
    self.end_el = time.time()
    if self.message not in timing_singleton:
      timing_singleton[self.message]=[0.,0.,0]
    timing_singleton[self.message][0]+=self.end-self.start
    timing_singleton[self.message][1]+=self.end_el-self.start_el
    timing_singleton[self.message][2]+=1
