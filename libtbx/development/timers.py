from __future__ import absolute_import, division, print_function
import time,sys

class Timer:
  """While the class instance is in scope it accumulates CPU and wall clock
     elapsed time.  A report is printed when the instances goes out of scope"""
  def __init__(self,message):
    self.start = time.clock()
    self.start_el = time.time()
    self.message = message
    print("start timing %s"%self.message)

  def tick(self):
    return time.clock() - self.start

  def __del__(self):
    self.end = time.clock()
    self.end_el = time.time()
    print("time for %s: CPU, %8.3fs; elapsed, %8.3fs"%(
      self.message,self.end-self.start,self.end_el-self.start_el))

class DebuggingTimer:
  def __init__(self,message,filename,filemode = 'a'):
    self.start = time.clock()
    self.start_el = time.time()
    self.message = message
    self.k = open(filename,filemode)
    self.current_out = sys.stdout
    self.current_err = sys.stderr
    sys.stdout = self.k
    sys.stderr = self.k
    print("start timing %s"%self.message)


  def __del__(self):
    self.end = time.clock()
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
    self.start = time.clock()
    self.start_el = time.time()
    self.message = message

  def __del__(self):
    self.end = time.clock()
    self.end_el = time.time()
    if self.message not in timing_singleton:
      timing_singleton[self.message]=[0.,0.,0]
    timing_singleton[self.message][0]+=self.end-self.start
    timing_singleton[self.message][1]+=self.end_el-self.start_el
    timing_singleton[self.message][2]+=1

    print("individual call time for %s: CPU, %8.3fs; elapsed, %8.3fs"%(
      self.message,self.end-self.start,self.end_el-self.start_el))
